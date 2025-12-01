use std::{collections::HashMap, fs, path::Path};

use anyhow::{anyhow, bail, Context, Result};
use csv::{ReaderBuilder, StringRecord, Trim};

#[derive(Debug, Clone)]
pub struct MotifDefinition {
    pub motif: String,
    pub canonical_base: char,
    pub canonical_offset: usize,
    pub mod_code: String,
    pub strand: char,
    allowed: Vec<Vec<char>>,
}

fn iupac_mapping() -> HashMap<char, Vec<char>> {
    [
        ('A', vec!['A']),
        ('C', vec!['C']),
        ('G', vec!['G']),
        ('T', vec!['T']),
        ('R', vec!['A', 'G']),
        ('Y', vec!['C', 'T']),
        ('S', vec!['G', 'C']),
        ('W', vec!['A', 'T']),
        ('K', vec!['G', 'T']),
        ('M', vec!['A', 'C']),
        ('B', vec!['C', 'G', 'T']),
        ('D', vec!['A', 'G', 'T']),
        ('H', vec!['A', 'C', 'T']),
        ('V', vec!['A', 'C', 'G']),
        ('N', vec!['A', 'C', 'G', 'T']),
    ]
    .into_iter()
    .collect()
}

fn mod_code_to_base(code: &str) -> Option<char> {
    let lower = code.to_ascii_lowercase();
    match lower.as_str() {
        "m" | "5mc" => Some('C'),
        "a" | "6ma" => Some('A'),
        "21839" | "4mc" => Some('C'),
        _ => lower
            .chars()
            .rev()
            .find(|c| matches!(c, 'a' | 'c' | 'g' | 't'))
            .map(|c| c.to_ascii_uppercase()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn load_motif_file_reads_basic_rows() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "motif\tmod_type\tmod_position\nGATC\t6mA\t1\nCCWGG\tm\t0"
        )
        .unwrap();
        let specs = load_motif_file(file.path()).unwrap();
        assert_eq!(
            specs,
            vec!["GATC_6mA_1".to_string(), "CCWGG_5mC_0".to_string()]
        );
    }

    #[test]
    fn load_motif_file_handles_complements() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "motif,mod_type,mod_position,motif_complement,mod_position_complement\n\
             GATC,6mA,1,GATC,2\n\
             CCWGG,m,0,CCWGG,4"
        )
        .unwrap();
        let specs = load_motif_file(file.path()).unwrap();
        assert_eq!(
            specs,
            vec![
                "GATC_6mA_1".to_string(),
                "GATC_6mA_2".to_string(),
                "CCWGG_5mC_0".to_string(),
                "CCWGG_5mC_4".to_string()
            ]
        );
    }

    #[test]
    fn load_motif_file_errors_on_empty_input() {
        let file = NamedTempFile::new().unwrap();
        let err = load_motif_file(file.path()).unwrap_err();
        assert!(err.to_string().contains("Motif file") && err.to_string().contains("is empty"));
    }

    #[test]
    fn load_motif_file_handles_extra_columns_and_reference_prefix() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "reference\tmotif\tmod_position\tmod_type\tn_mod\tn_nomod\tmotif_type\tmotif_complement\tmod_position_complement\tn_mod_complement\tn_nomod_complement\n\
             bin_1\tCCWGG\t1\tm\t51\t0\tpalindrome\tCCWGG\t1\t51\t0\n\
             bin_1\tGCWGC\t1\tm\t256\t0\tpalindrome\tGCWGC\t1\t256\t0"
        )
        .unwrap();
        let specs = load_motif_file(file.path()).unwrap();
        assert_eq!(
            specs,
            vec!["CCWGG_5mC_1".to_string(), "GCWGC_5mC_1".to_string(),]
        );
    }
}

pub fn load_motif_file(path: &Path) -> Result<Vec<String>> {
    let data = fs::read(path)
        .with_context(|| format!("Failed to read motif file '{}'", path.display()))?;
    if data.is_empty() {
        bail!("Motif file '{}' is empty", path.display());
    }
    let delimiter = detect_delimiter(&data);
    let mut reader = ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(true)
        .trim(Trim::All)
        .flexible(true)
        .from_reader(&data[..]);
    let headers = reader
        .headers()
        .with_context(|| format!("Failed to read headers from '{}'", path.display()))?
        .clone();
    let motif_idx = find_column(&headers, "motif")
        .ok_or_else(|| anyhow!("motif column missing in '{}'", path.display()))?;
    let mod_type_idx = find_column(&headers, "mod_type")
        .ok_or_else(|| anyhow!("mod_type column missing in '{}'", path.display()))?;
    let mod_position_idx = find_column(&headers, "mod_position")
        .ok_or_else(|| anyhow!("mod_position column missing in '{}'", path.display()))?;
    let motif_comp_idx = find_column(&headers, "motif_complement");
    let mod_pos_comp_idx = find_column(&headers, "mod_position_complement");
    if motif_comp_idx.is_some() ^ mod_pos_comp_idx.is_some() {
        bail!(
            "motif_complement and mod_position_complement columns must both be present in '{}'",
            path.display()
        );
    }
    let mut specs = Vec::new();
    for record in reader.records() {
        let record = record
            .with_context(|| format!("Failed to read motif file record in '{}'", path.display()))?;
        let motif = record
            .get(motif_idx)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| anyhow!("Row missing motif value in '{}'", path.display()))?;
        let motif_len = motif.chars().count();
        let mod_label = record
            .get(mod_type_idx)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| anyhow!("Row missing mod_type value in '{}'", path.display()))
            .and_then(normalize_mod_label)?;
        let position_raw = record
            .get(mod_position_idx)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| anyhow!("Row missing mod_position value in '{}'", path.display()))?;
        let mod_position: usize = position_raw.parse().with_context(|| {
            format!(
                "Invalid mod_position value '{}' in '{}'",
                position_raw,
                path.display()
            )
        })?;
        if mod_position >= motif_len {
            bail!(
                "mod_position {} exceeds motif length {} in '{}'",
                mod_position,
                motif_len,
                path.display()
            );
        }
        let spec = format!("{}_{}_{}", motif, mod_label, mod_position);
        specs.push(spec.clone());
        if let (Some(motif_idx), Some(pos_idx)) = (motif_comp_idx, mod_pos_comp_idx) {
            let comp_motif = record
                .get(motif_idx)
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string());
            let comp_position = record.get(pos_idx).map(str::trim).filter(|s| !s.is_empty());
            if let (Some(comp_motif), Some(comp_position)) = (comp_motif, comp_position) {
                let position: usize = comp_position.parse().with_context(|| {
                    format!(
                        "Invalid mod_position_complement value '{}' in '{}'",
                        comp_position,
                        path.display()
                    )
                })?;
                let comp_len = comp_motif.chars().count();
                if position >= comp_len {
                    bail!(
                        "mod_position_complement {} exceeds motif length {} in '{}'",
                        position,
                        comp_len,
                        path.display()
                    );
                }
                let comp_spec = format!("{}_{}_{}", comp_motif, mod_label, position);
                if comp_spec != spec {
                    specs.push(comp_spec);
                }
            }
        }
    }
    if specs.is_empty() {
        bail!(
            "Motif file '{}' did not contain any usable entries",
            path.display()
        );
    }
    Ok(specs)
}

fn detect_delimiter(data: &[u8]) -> u8 {
    let first_line = data.split(|b| *b == b'\n').next().unwrap_or(data);
    if first_line.contains(&b'\t') {
        b'\t'
    } else {
        b','
    }
}

fn find_column(headers: &StringRecord, name: &str) -> Option<usize> {
    headers.iter().position(|h| h.eq_ignore_ascii_case(name))
}

fn normalize_mod_label(raw: &str) -> Result<String> {
    let token = raw.trim();
    if token.is_empty() {
        bail!("mod_type value is empty");
    }
    if token.len() == 1 {
        if let Some(label) = code_to_label(token.chars().next().unwrap()) {
            return Ok(label.to_string());
        }
    }
    let lowered = token.to_ascii_lowercase();
    if let Some(label) = name_to_label(&lowered) {
        return Ok(label.to_string());
    }
    Ok(token.to_string())
}

fn code_to_label(code: char) -> Option<&'static str> {
    match code.to_ascii_lowercase() {
        'a' => Some("6mA"),
        'm' => Some("5mC"),
        'h' => Some("5hmC"),
        'f' => Some("5fC"),
        'g' => Some("5gmC"),
        'c' => Some("4mC"),
        _ => None,
    }
}

fn name_to_label(name: &str) -> Option<&'static str> {
    match name {
        "6ma" => Some("6mA"),
        "5mc" => Some("5mC"),
        "5hmc" => Some("5hmC"),
        "5fc" => Some("5fC"),
        "5gmc" => Some("5gmC"),
        "4mc" => Some("4mC"),
        _ => None,
    }
}

fn canonical_mod_code(code: &str) -> String {
    match code.to_ascii_lowercase().as_str() {
        "5mc" => "m".to_string(),
        "6ma" => "a".to_string(),
        other => {
            if other == "m" || other == "a" || other == "21839" || other == "4mc" {
                other.to_string()
            } else {
                code.to_string()
            }
        }
    }
}

impl MotifDefinition {
    pub fn parse(spec: &str) -> Result<Self> {
        let trimmed = spec.trim();
        if trimmed.is_empty() {
            bail!("Motif specification cannot be empty");
        }
        let mapping = iupac_mapping();
        if trimmed.contains('_') && !trimmed.contains(':') {
            Self::parse_underscore(trimmed, &mapping)
        } else {
            Self::parse_colon(trimmed, &mapping)
        }
    }

    fn parse_colon(spec: &str, mapping: &HashMap<char, Vec<char>>) -> Result<Self> {
        let parts: Vec<_> = spec.split(':').collect();
        let motif_raw = parts.get(0).map(|s| s.trim()).unwrap_or("");
        if motif_raw.is_empty() {
            bail!("Motif string missing in specification '{}'", spec);
        }
        let motif = motif_raw.to_ascii_uppercase();
        let canonical_base = parts
            .get(1)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| s.to_ascii_uppercase().chars().next().unwrap())
            .or_else(|| {
                parts
                    .get(3)
                    .map(|s| s.trim())
                    .filter(|s| !s.is_empty())
                    .and_then(|code| mod_code_to_base(code))
            })
            .or_else(|| motif.chars().find(|c| matches!(c, 'A' | 'C' | 'G' | 'T')))
            .ok_or_else(|| anyhow!("Unable to determine canonical base for '{}'", spec))?;
        if !matches!(canonical_base, 'A' | 'C' | 'G' | 'T') {
            bail!(
                "Canonical base must be one of A/C/G/T (received '{}')",
                canonical_base
            );
        }

        let canonical_offset = parts
            .get(2)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|offset| {
                offset.parse::<usize>().map_err(|_| {
                    anyhow!(
                        "Unable to parse canonical offset '{}' in '{}'",
                        offset,
                        spec
                    )
                })
            })
            .transpose()?
            .unwrap_or_else(|| {
                motif
                    .chars()
                    .position(|c| {
                        mapping
                            .get(&c)
                            .map_or(false, |allowed| allowed.contains(&canonical_base))
                    })
                    .unwrap_or(0)
            });

        let raw_mod_code = parts
            .get(3)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .unwrap_or("m")
            .to_string();
        let mod_code = canonical_mod_code(&raw_mod_code);
        let strand = parts
            .get(4)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| s.chars().next().unwrap())
            .unwrap_or('+');
        Self::build(
            motif,
            canonical_base,
            canonical_offset,
            mod_code,
            strand,
            mapping,
            spec,
        )
    }

    fn parse_underscore(spec: &str, mapping: &HashMap<char, Vec<char>>) -> Result<Self> {
        let parts: Vec<_> = spec.split('_').map(|s| s.trim()).collect();
        if parts.len() < 3 {
            bail!(
                "Motif specification '{}' must follow sequence_modtype_offset format",
                spec
            );
        }
        let motif = parts[0].to_ascii_uppercase();
        if motif.is_empty() {
            bail!("Motif string missing in specification '{}'", spec);
        }
        let mod_code_raw = parts[1];
        if mod_code_raw.is_empty() {
            bail!("Modification type missing in specification '{}'", spec);
        }
        let canonical_offset = parts[2].parse::<usize>().map_err(|_| {
            anyhow!(
                "Unable to parse canonical offset '{}' in '{}'",
                parts[2],
                spec
            )
        })?;
        let strand = if parts.len() > 3 && !parts[3].is_empty() {
            let strand_char = parts[3].chars().next().unwrap();
            if strand_char != '+' && strand_char != '-' {
                bail!("Strand must be '+' or '-' (received '{}')", strand_char);
            }
            strand_char
        } else {
            '+'
        };
        let canonical_base = mod_code_to_base(mod_code_raw)
            .or_else(|| motif.chars().find(|c| matches!(c, 'A' | 'C' | 'G' | 'T')))
            .ok_or_else(|| anyhow!("Unable to determine canonical base for '{}'", spec))?;
        let mod_code = canonical_mod_code(mod_code_raw);
        Self::build(
            motif,
            canonical_base,
            canonical_offset,
            mod_code,
            strand,
            mapping,
            spec,
        )
    }

    fn build(
        motif: String,
        canonical_base: char,
        canonical_offset: usize,
        mod_code: String,
        strand: char,
        mapping: &HashMap<char, Vec<char>>,
        spec: &str,
    ) -> Result<Self> {
        if canonical_offset >= motif.len() {
            bail!(
                "Canonical offset {} out of bounds for motif '{}'",
                canonical_offset,
                motif
            );
        }
        if strand != '+' && strand != '-' {
            bail!("Strand must be '+' or '-' (received '{}')", strand);
        }
        if !matches!(canonical_base, 'A' | 'C' | 'G' | 'T') {
            bail!(
                "Canonical base must be one of A/C/G/T (received '{}')",
                canonical_base
            );
        }
        let allowed: Vec<Vec<char>> = motif
            .chars()
            .map(|c| {
                mapping
                    .get(&c)
                    .cloned()
                    .ok_or_else(|| anyhow!("Unsupported motif character '{}' in '{}'", c, spec))
            })
            .collect::<Result<_>>()?;
        Ok(MotifDefinition {
            motif,
            canonical_base,
            canonical_offset,
            mod_code,
            strand,
            allowed,
        })
    }

    pub fn find_matches(&self, sequence: &str) -> Vec<usize> {
        let seq_bytes: Vec<u8> = sequence
            .as_bytes()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        if seq_bytes.len() < self.allowed.len() {
            return Vec::new();
        }
        let mut hits = Vec::new();
        for start in 0..=seq_bytes.len() - self.allowed.len() {
            let mut valid = true;
            for (offset, allowed) in self.allowed.iter().enumerate() {
                let base = seq_bytes[start + offset] as char;
                if !allowed.contains(&base) {
                    valid = false;
                    break;
                }
            }
            if valid {
                let canonical_index = start + self.canonical_offset;
                if seq_bytes[canonical_index] as char == self.canonical_base {
                    hits.push(canonical_index);
                }
            }
        }
        hits
    }

    pub fn model_key(&self) -> String {
        format!("{}_{}_{}", self.motif, self.mod_code, self.canonical_offset)
    }
}
