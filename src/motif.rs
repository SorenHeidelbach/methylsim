use std::collections::HashMap;

use anyhow::{Result, anyhow, bail};

#[derive(Debug, Clone)]
pub struct MotifDefinition {
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
}
