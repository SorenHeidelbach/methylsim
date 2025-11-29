use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};

#[derive(Debug, Clone)]
pub struct ReadRecord {
    pub name: String,
    pub comment: Option<String>,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

impl ReadRecord {
    /// Get sequence as string (for backward compatibility)
    pub fn sequence_str(&self) -> &str {
        std::str::from_utf8(&self.sequence).unwrap_or("")
    }

    /// Get quality as string (for backward compatibility)
    pub fn quality_str(&self) -> &str {
        std::str::from_utf8(&self.quality).unwrap_or("")
    }
}

enum FastxFormat {
    Fasta,
    Fastq,
}

impl ReadRecord {
    pub fn new(name: String, sequence: Vec<u8>, quality: Vec<u8>) -> Self {
        Self {
            name,
            comment: None,
            sequence,
            quality,
        }
    }

    pub fn with_comment(
        name: String,
        comment: Option<String>,
        sequence: Vec<u8>,
        quality: Vec<u8>,
    ) -> Self {
        Self {
            name,
            comment,
            sequence,
            quality,
        }
    }
}

fn detect_format(path: &Path) -> Result<FastxFormat> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open FASTA/FASTQ file '{}'", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    loop {
        line.clear();
        if reader
            .read_line(&mut line)
            .with_context(|| format!("Failed to read from '{}'", path.display()))?
            == 0
        {
            bail!("Input file '{}' is empty", path.display());
        }
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        return match trimmed.chars().next() {
            Some('>') => Ok(FastxFormat::Fasta),
            Some('@') => Ok(FastxFormat::Fastq),
            _ => bail!(
                "Unable to detect FASTA/FASTQ format for '{}'; first non-empty line must start with '>' or '@'",
                path.display()
            ),
        };
    }
}

pub fn read_fastx(path: &Path) -> Result<Vec<ReadRecord>> {
    match detect_format(path)? {
        FastxFormat::Fasta => read_fasta_records(path),
        FastxFormat::Fastq => read_fastq_records(path),
    }
}

pub fn read_fasta_sequences(path: &Path) -> Result<Vec<String>> {
    Ok(read_fasta_records(path)?
        .into_iter()
        .map(|record| String::from_utf8_lossy(&record.sequence).into_owned())
        .collect())
}

pub fn read_fastx_with_progress(path: &Path, show_progress: bool) -> Result<Vec<ReadRecord>> {
    match detect_format(path)? {
        FastxFormat::Fasta => read_fasta_records(path),
        FastxFormat::Fastq => read_fastq_records_with_progress(path, show_progress),
    }
}

fn read_fasta_records(path: &Path) -> Result<Vec<ReadRecord>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open FASTA file '{}'", path.display()))?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line.with_context(|| format!("Failed reading '{}'", path.display()))?;
        if line.starts_with('>') {
            if let Some(name) = current_name.take() {
                if current_seq.is_empty() {
                    bail!(
                        "FASTA record '{}' has no sequence in '{}'",
                        name,
                        path.display()
                    );
                }
                // Make sequence uppercase
                for base in current_seq.iter_mut() {
                    base.make_ascii_uppercase();
                }
                let quality = vec![b'I'; current_seq.len()];
                records.push(ReadRecord::new(name, current_seq, quality));
                current_seq = Vec::new();
            }
            current_name = Some(line[1..].trim().to_string());
        } else {
            current_seq.extend_from_slice(line.trim().as_bytes());
        }
    }

    if let Some(name) = current_name {
        if current_seq.is_empty() {
            bail!(
                "FASTA record '{}' has no sequence in '{}'",
                name,
                path.display()
            );
        }
        for base in current_seq.iter_mut() {
            base.make_ascii_uppercase();
        }
        let quality = vec![b'I'; current_seq.len()];
        records.push(ReadRecord::new(name, current_seq, quality));
    }

    if records.is_empty() {
        bail!("No FASTA records could be read from '{}'", path.display());
    }
    Ok(records)
}

pub fn read_fastq_records(path: &Path) -> Result<Vec<ReadRecord>> {
    read_fastq_records_with_progress(path, false)
}

pub fn read_fastq_records_with_progress(path: &Path, show_progress: bool) -> Result<Vec<ReadRecord>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open FASTQ file '{}'", path.display()))?;
    let mut reader = BufReader::with_capacity(64 * 1024, file);
    let mut records = Vec::new();

    let mut header = String::new();
    let mut sequence = String::new();
    let mut plus = String::new();
    let mut quality = String::new();
    let mut record_count = 0;
    const PROGRESS_INTERVAL: usize = 10_000;

    loop {
        // Read header line
        header.clear();
        let bytes_read = reader.read_line(&mut header)
            .with_context(|| format!("Failed reading FASTQ from '{}'", path.display()))?;

        if bytes_read == 0 {
            break; // EOF
        }

        if !header.starts_with('@') {
            // Skip empty lines at end of file
            if header.trim().is_empty() {
                continue;
            }
            bail!(
                "Invalid FASTQ header (expected '@'): {} in '{}'",
                header.trim_end(),
                path.display()
            );
        }

        // Parse header - remove '@' and trim newline
        let header_content = header[1..].trim_end_matches(|c| c == '\r' || c == '\n');

        // Split off MM/ML tags
        let (raw_identifier, tag_suffix) = split_methylation_suffix(header_content);

        if raw_identifier.is_empty() || raw_identifier.chars().all(|ch| matches!(ch, ' ' | '\t')) {
            bail!("FASTQ header missing read name in '{}'", path.display());
        }

        let name = sanitize_read_identifier(raw_identifier);
        if name.is_empty() {
            bail!("FASTQ header missing read name in '{}'", path.display());
        }

        let comment = if tag_suffix.is_empty() {
            None
        } else {
            Some(tag_suffix.to_string())
        };

        // Read sequence line
        sequence.clear();
        if reader.read_line(&mut sequence)
            .with_context(|| format!("Failed reading sequence for '{}' in '{}'", name, path.display()))? == 0 {
            bail!("Unexpected EOF while reading sequence for '{}' in '{}'", name, path.display());
        }

        // Read '+' separator line
        plus.clear();
        if reader.read_line(&mut plus)
            .with_context(|| format!("Failed reading separator for '{}' in '{}'", name, path.display()))? == 0 {
            bail!("Unexpected EOF while reading separator for '{}' in '{}'", name, path.display());
        }
        if !plus.starts_with('+') {
            bail!("Invalid FASTQ separator (expected '+') for '{}' in '{}'", name, path.display());
        }

        // Read quality line
        quality.clear();
        if reader.read_line(&mut quality)
            .with_context(|| format!("Failed reading quality for '{}' in '{}'", name, path.display()))? == 0 {
            bail!("Unexpected EOF while reading quality for '{}' in '{}'", name, path.display());
        }

        // Convert to bytes and uppercase sequence
        let seq_trimmed = sequence.trim_end_matches(|c| c == '\r' || c == '\n');
        let qual_trimmed = quality.trim_end_matches(|c| c == '\r' || c == '\n');

        if seq_trimmed.len() != qual_trimmed.len() {
            bail!(
                "FASTQ sequence/quality length mismatch for '{}' in '{}' ({} vs {})",
                name,
                path.display(),
                seq_trimmed.len(),
                qual_trimmed.len()
            );
        }

        let mut seq_bytes = seq_trimmed.as_bytes().to_vec();
        let qual_bytes = qual_trimmed.as_bytes().to_vec();

        // Make sequence uppercase
        for base in seq_bytes.iter_mut() {
            base.make_ascii_uppercase();
        }

        records.push(ReadRecord::with_comment(name, comment, seq_bytes, qual_bytes));

        record_count += 1;
        if show_progress && record_count % PROGRESS_INTERVAL == 0 {
            eprint!("\rRead {} records...", record_count);
        }
    }

    if show_progress && record_count > 0 {
        eprintln!("\rRead {} records total", record_count);
    }

    if records.is_empty() {
        bail!("No FASTQ records could be read from '{}'", path.display());
    }
    Ok(records)
}

/// Process FASTQ records one at a time using a callback, minimizing memory usage
pub fn process_fastq_streaming<F>(path: &Path, show_progress: bool, mut callback: F) -> Result<()>
where
    F: FnMut(&ReadRecord) -> Result<()>,
{
    let file = File::open(path)
        .with_context(|| format!("Failed to open FASTQ file '{}'", path.display()))?;
    let mut reader = BufReader::with_capacity(64 * 1024, file);

    let mut header = String::new();
    let mut sequence = String::new();
    let mut plus = String::new();
    let mut quality = String::new();
    let mut record_count = 0;
    const PROGRESS_INTERVAL: usize = 10_000;

    loop {
        // Read header line
        header.clear();
        let bytes_read = reader.read_line(&mut header)
            .with_context(|| format!("Failed reading FASTQ from '{}'", path.display()))?;

        if bytes_read == 0 {
            break; // EOF
        }

        if !header.starts_with('@') {
            // Skip empty lines at end of file
            if header.trim().is_empty() {
                continue;
            }
            bail!(
                "Invalid FASTQ header (expected '@'): {} in '{}'",
                header.trim_end(),
                path.display()
            );
        }

        // Parse header - remove '@' and trim newline
        let header_content = header[1..].trim_end_matches(|c| c == '\r' || c == '\n');

        // Split off MM/ML tags
        let (raw_identifier, tag_suffix) = split_methylation_suffix(header_content);

        if raw_identifier.is_empty() || raw_identifier.chars().all(|ch| matches!(ch, ' ' | '\t')) {
            bail!("FASTQ header missing read name in '{}'", path.display());
        }

        let name = sanitize_read_identifier(raw_identifier);
        if name.is_empty() {
            bail!("FASTQ header missing read name in '{}'", path.display());
        }

        let comment = if tag_suffix.is_empty() {
            None
        } else {
            Some(tag_suffix.to_string())
        };

        // Read sequence line
        sequence.clear();
        if reader.read_line(&mut sequence)
            .with_context(|| format!("Failed reading sequence for '{}' in '{}'", name, path.display()))? == 0 {
            bail!("Unexpected EOF while reading sequence for '{}' in '{}'", name, path.display());
        }

        // Read '+' separator line
        plus.clear();
        if reader.read_line(&mut plus)
            .with_context(|| format!("Failed reading separator for '{}' in '{}'", name, path.display()))? == 0 {
            bail!("Unexpected EOF while reading separator for '{}' in '{}'", name, path.display());
        }
        if !plus.starts_with('+') {
            bail!("Invalid FASTQ separator (expected '+') for '{}' in '{}'", name, path.display());
        }

        // Read quality line
        quality.clear();
        if reader.read_line(&mut quality)
            .with_context(|| format!("Failed reading quality for '{}' in '{}'", name, path.display()))? == 0 {
            bail!("Unexpected EOF while reading quality for '{}' in '{}'", name, path.display());
        }

        // Convert to bytes and uppercase sequence
        let seq_trimmed = sequence.trim_end_matches(|c| c == '\r' || c == '\n');
        let qual_trimmed = quality.trim_end_matches(|c| c == '\r' || c == '\n');

        if seq_trimmed.len() != qual_trimmed.len() {
            bail!(
                "FASTQ sequence/quality length mismatch for '{}' in '{}' ({} vs {})",
                name,
                path.display(),
                seq_trimmed.len(),
                qual_trimmed.len()
            );
        }

        let mut seq_bytes = seq_trimmed.as_bytes().to_vec();
        let qual_bytes = qual_trimmed.as_bytes().to_vec();

        // Make sequence uppercase
        for base in seq_bytes.iter_mut() {
            base.make_ascii_uppercase();
        }

        let record = ReadRecord::with_comment(name, comment, seq_bytes, qual_bytes);

        // Process the record immediately via callback
        callback(&record)?;

        record_count += 1;
        if show_progress && record_count % PROGRESS_INTERVAL == 0 {
            eprint!("\rProcessed {} records...", record_count);
        }
    }

    if show_progress && record_count > 0 {
        eprintln!("\rProcessed {} records total", record_count);
    }

    if record_count == 0 {
        bail!("No FASTQ records could be read from '{}'", path.display());
    }
    Ok(())
}

fn sanitize_read_identifier(identifier: &str) -> String {
    let mut sanitized = String::with_capacity(identifier.len());
    for ch in identifier.chars() {
        match ch {
            ' ' => sanitized.push(','),
            '\t' => sanitized.push(';'),
            _ => sanitized.push(ch),
        }
    }
    sanitized
}

fn split_methylation_suffix(header: &str) -> (&str, &str) {
    let mm_idx = header.find("\tMM");
    let ml_idx = header.find("\tML");
    let split_idx = match (mm_idx, ml_idx) {
        (Some(mm), Some(ml)) => Some(mm.min(ml)),
        (Some(mm), None) => Some(mm),
        (None, Some(ml)) => Some(ml),
        (None, None) => None,
    };
    if let Some(idx) = split_idx {
        header.split_at(idx)
    } else {
        (header, "")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn sanitize_identifier_replaces_spaces_and_tabs() {
        let original = "read id\twith space";
        let sanitized = sanitize_read_identifier(original);
        assert_eq!(sanitized, "read,id;with,space");
    }

    #[test]
    fn fastq_reader_sanitizes_identifier_and_preserves_tags() {
        let mut tmp = NamedTempFile::new().expect("create temp file");
        writeln!(tmp, "@read id extra\tMM:Z:Z\nACGT\n+\n!!!!").unwrap();
        tmp.flush().unwrap();
        let records = read_fastq_records(tmp.path()).expect("parse fastq");
        assert_eq!(records.len(), 1);
        let record = &records[0];
        assert_eq!(record.name, "read,id,extra");
        assert_eq!(record.comment.as_deref(), Some("\tMM:Z:Z"));
    }
}
