use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};

#[derive(Debug, Clone)]
pub struct ReadRecord {
    pub name: String,
    pub comment: Option<String>,
    pub sequence: String,
    pub quality: String,
}

enum FastxFormat {
    Fasta,
    Fastq,
}

impl ReadRecord {
    pub fn new(name: String, sequence: String, quality: String) -> Self {
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
        sequence: String,
        quality: String,
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
        .map(|record| record.sequence)
        .collect())
}

fn read_fasta_records(path: &Path) -> Result<Vec<ReadRecord>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open FASTA file '{}'", path.display()))?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = String::new();

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
                let sequence = current_seq.to_ascii_uppercase();
                let quality = "I".repeat(sequence.len());
                records.push(ReadRecord::new(name, sequence, quality));
                current_seq.clear();
            }
            current_name = Some(line[1..].trim().to_string());
        } else {
            current_seq.push_str(line.trim());
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
        let sequence = current_seq.to_ascii_uppercase();
        let quality = "I".repeat(sequence.len());
        records.push(ReadRecord::new(name, sequence, quality));
    }

    if records.is_empty() {
        bail!("No FASTA records could be read from '{}'", path.display());
    }
    Ok(records)
}

pub fn read_fastq_records(path: &Path) -> Result<Vec<ReadRecord>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open FASTQ file '{}'", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if !line.starts_with('@') {
            bail!(
                "Expected FASTQ header line starting with '@' in '{}'",
                path.display()
            );
        }
        let header = line[1..].trim_end();
        let header = header.trim_start();
        if header.is_empty() {
            bail!("FASTQ header missing read name in '{}'", path.display());
        }
        let (raw_identifier, tag_suffix) = split_methylation_suffix(header);
        if raw_identifier.chars().all(|ch| matches!(ch, ' ' | '\t')) {
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
        let mut seq = String::new();
        reader
            .read_line(&mut seq)
            .map_err(|e| anyhow!("Failed reading sequence for '{}' ({})", name, e))?;
        let sequence = seq.trim_end().to_ascii_uppercase();

        line.clear();
        reader
            .read_line(&mut line)
            .map_err(|e| anyhow!("Failed reading separator for '{}' ({})", name, e))?;
        if !line.starts_with('+') {
            bail!(
                "FASTQ separator '+' missing for read '{}' in '{}'",
                name,
                path.display()
            );
        }

        let mut qual = String::new();
        reader
            .read_line(&mut qual)
            .map_err(|e| anyhow!("Failed reading quality for '{}' ({})", name, e))?;
        let quality = qual.trim_end().to_string();
        if quality.len() != sequence.len() {
            bail!(
                "FASTQ sequence/quality length mismatch for '{}' in '{}'",
                name,
                path.display()
            );
        }
        records.push(ReadRecord::with_comment(name, comment, sequence, quality));
    }

    if records.is_empty() {
        bail!("No FASTQ records could be read from '{}'", path.display());
    }
    Ok(records)
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
