use std::collections::HashMap;
use std::fs;
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};
use rand::Rng;
use rand::rngs::StdRng;
use rand_distr::{Beta, Distribution};
use serde::{Deserialize, Serialize};

use crate::fastx::ReadRecord;
use crate::motif::MotifDefinition;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BetaParams {
    pub alpha: f64,
    pub beta: f64,
    pub n: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EmpiricalDistribution {
    pub bin_width: f64,
    pub bins: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MotifModel {
    pub key: String,
    pub motif: String,
    pub canonical_base: char,
    pub canonical_offset: usize,
    pub mod_code: String,
    pub strand: char,
    pub total_sites: usize,
    pub covered_sites: usize,
    pub methylated_fraction: f64,
    pub empirical: EmpiricalDistribution,
    pub methylated_beta: Option<BetaParams>,
    pub unmethylated_beta: Option<BetaParams>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LearnedModel {
    pub version: u32,
    pub threshold: f64,
    pub motifs: Vec<MotifModel>,
}

#[derive(Debug, Clone)]
pub struct MotifSampler {
    methylated_fraction: f64,
    methylated_beta: Beta<f64>,
    unmethylated_beta: Beta<f64>,
}

struct TaggedEvent {
    position: usize,
    canonical_base: char,
    mod_code: String,
    strand: char,
    probability: f64,
}

type EventMap = HashMap<(char, String, char), HashMap<usize, f64>>;

pub fn learn_model_from_reads(
    reads: &[ReadRecord],
    motifs: &[MotifDefinition],
    threshold: f64,
) -> Result<LearnedModel> {
    if !(0.0..=1.0).contains(&threshold) {
        bail!("Model threshold must be between 0 and 1 (received {threshold})");
    }
    if motifs.is_empty() {
        bail!("At least one motif is required to learn a model");
    }
    let mut values: HashMap<String, Vec<f64>> = HashMap::new();
    let mut totals: HashMap<String, usize> = HashMap::new();
    for motif in motifs {
        values.insert(motif.model_key(), Vec::new());
        totals.insert(motif.model_key(), 0);
    }

    for read in reads {
        let parsed = parse_tagged_events(read)?;
        for motif in motifs {
            let key = motif.model_key();
            let hits = motif.find_matches(&read.sequence);
            *totals.entry(key.clone()).or_default() += hits.len();
            let Some(event_map) = parsed.as_ref() else {
                continue;
            };
            let Some(probs) =
                event_map.get(&(motif.canonical_base, motif.mod_code.clone(), motif.strand))
            else {
                continue;
            };
            for pos in hits {
                if let Some(prob) = probs.get(&pos) {
                    values.entry(key.clone()).or_default().push(*prob);
                }
            }
        }
    }

    let mut motif_models = Vec::new();
    for motif in motifs {
        let key = motif.model_key();
        let covered: Vec<f64> = values.remove(&key).unwrap_or_default();
        let total_sites = totals.remove(&key).unwrap_or_default();
        let methylated: Vec<f64> = covered
            .iter()
            .copied()
            .filter(|p| *p >= threshold)
            .collect();
        let unmethylated: Vec<f64> = covered.iter().copied().filter(|p| *p < threshold).collect();
        let meth_beta = fit_beta(&methylated);
        let unmeth_beta = fit_beta(&unmethylated);
        let methylated_fraction = if covered.is_empty() {
            0.0
        } else {
            methylated.len() as f64 / covered.len() as f64
        };
        let empirical = EmpiricalDistribution::from_values(&covered, 0.05);
        motif_models.push(MotifModel {
            key,
            motif: motif.motif.clone(),
            canonical_base: motif.canonical_base,
            canonical_offset: motif.canonical_offset,
            mod_code: motif.mod_code.clone(),
            strand: motif.strand,
            total_sites,
            covered_sites: covered.len(),
            methylated_fraction,
            empirical,
            methylated_beta: meth_beta,
            unmethylated_beta: unmeth_beta,
        });
    }

    Ok(LearnedModel {
        version: 1,
        threshold,
        motifs: motif_models,
    })
}

impl LearnedModel {
    pub fn save(&self, path: &Path) -> Result<()> {
        let data = serde_json::to_vec_pretty(self)?;
        fs::write(path, data)
            .with_context(|| format!("Failed to write model to '{}'", path.display()))
    }

    pub fn load(path: &Path) -> Result<Self> {
        let data = fs::read(path)
            .with_context(|| format!("Failed to read model file '{}'", path.display()))?;
        let model: LearnedModel = serde_json::from_slice(&data)
            .with_context(|| format!("Failed to deserialize model from '{}'", path.display()))?;
        Ok(model)
    }

    pub fn build_samplers(&self) -> Result<HashMap<String, MotifSampler>> {
        let mut map = HashMap::new();
        for motif in &self.motifs {
            let meth_params = motif.methylated_beta.clone().unwrap_or(BetaParams {
                alpha: 1.0,
                beta: 1.0,
                n: 0,
            });
            let unmeth_params = motif.unmethylated_beta.clone().unwrap_or(BetaParams {
                alpha: 1.0,
                beta: 1.0,
                n: 0,
            });
            let methylated_beta =
                Beta::new(meth_params.alpha.max(1e-3), meth_params.beta.max(1e-3)).map_err(
                    |e| anyhow!("Invalid methylated beta params for {}: {}", motif.key, e),
                )?;
            let unmethylated_beta =
                Beta::new(unmeth_params.alpha.max(1e-3), unmeth_params.beta.max(1e-3)).map_err(
                    |e| anyhow!("Invalid unmethylated beta params for {}: {}", motif.key, e),
                )?;
            let fraction = motif.methylated_fraction.clamp(0.0, 1.0);
            map.insert(
                motif.key.clone(),
                MotifSampler {
                    methylated_fraction: fraction,
                    methylated_beta,
                    unmethylated_beta,
                },
            );
        }
        Ok(map)
    }
}

impl EmpiricalDistribution {
    pub fn from_values(values: &[f64], bin_width: f64) -> Self {
        let bins = if bin_width > 0.0 {
            (1.0 / bin_width).ceil() as usize
        } else {
            0
        };
        if bins == 0 {
            return EmpiricalDistribution {
                bin_width: 0.0,
                bins: Vec::new(),
            };
        }
        let mut counts = vec![0usize; bins];
        for v in values {
            let mut idx = (v.clamp(0.0, 1.0) / bin_width).floor() as usize;
            if idx >= bins {
                idx = bins - 1;
            }
            counts[idx] += 1;
        }
        EmpiricalDistribution {
            bin_width,
            bins: counts,
        }
    }
}

impl MotifSampler {
    pub fn sample_probability(&self, rng: &mut StdRng) -> f64 {
        let pick_methylated = rng.gen_bool(self.methylated_fraction.clamp(0.0, 1.0));
        let prob = if pick_methylated {
            self.methylated_beta.sample(rng)
        } else {
            self.unmethylated_beta.sample(rng)
        };
        prob.clamp(0.0, 1.0)
    }
}

fn fit_beta(values: &[f64]) -> Option<BetaParams> {
    if values.is_empty() {
        return None;
    }
    let n = values.len();
    let mean = values.iter().copied().sum::<f64>() / n as f64;
    let variance = if n > 1 {
        let mean_sq = values.iter().map(|v| v * v).sum::<f64>() / n as f64;
        (mean_sq - mean * mean).abs()
    } else {
        0.0
    };
    let mean = mean.clamp(1e-6, 1.0 - 1e-6);
    let variance = variance.max(1e-6);
    let temp = mean * (1.0 - mean) / variance - 1.0;
    if !temp.is_finite() || temp <= 0.0 {
        let concentration = 10.0;
        let alpha = (mean * concentration).max(1e-3);
        let beta = ((1.0 - mean) * concentration).max(1e-3);
        return Some(BetaParams { alpha, beta, n });
    }
    let alpha = (mean * temp).max(1e-3);
    let beta = ((1.0 - mean) * temp).max(1e-3);
    Some(BetaParams { alpha, beta, n })
}

fn parse_tagged_events(read: &ReadRecord) -> Result<Option<EventMap>> {
    let comment = match read.comment.as_deref() {
        Some(c) if !c.is_empty() => c,
        _ => return Ok(None),
    };
    let (mm_tag, ml_tag) = extract_mm_ml(comment);
    let Some(mm_tag) = mm_tag else {
        return Ok(None);
    };
    let ml_values = parse_ml_values(ml_tag)?;
    let mut ml_iter = ml_values.into_iter();
    let segments = normalize_mm_segments(mm_tag);
    if segments.is_empty() {
        return Ok(None);
    }
    let mut events: Vec<TaggedEvent> = Vec::new();
    for segment in segments {
        let Some((base, strand, mod_code, deltas)) = parse_mm_segment(&segment) else {
            continue;
        };
        let positions = decode_positions(&read.sequence, base, &deltas);
        for pos in positions {
            if let Some(ml) = ml_iter.next() {
                events.push(TaggedEvent {
                    position: pos,
                    canonical_base: base,
                    mod_code: mod_code.clone(),
                    strand,
                    probability: (ml as f64) / 255.0,
                });
            }
        }
    }
    if events.is_empty() {
        return Ok(None);
    }
    let mut grouped: EventMap = HashMap::new();
    for event in events {
        grouped
            .entry((event.canonical_base, event.mod_code.clone(), event.strand))
            .or_default()
            .insert(event.position, event.probability);
    }
    Ok(Some(grouped))
}

fn extract_mm_ml(comment: &str) -> (Option<&str>, Option<&str>) {
    let mut mm = None;
    let mut ml = None;
    for token in comment.split_whitespace() {
        if token.starts_with("MM:Z:") {
            mm = Some(token);
        } else if token.starts_with("ML:B:") {
            ml = Some(token);
        }
    }
    (mm, ml)
}

fn parse_ml_values(tag: Option<&str>) -> Result<Vec<u8>> {
    let Some(tag) = tag else {
        return Ok(Vec::new());
    };
    let (_, values) = tag
        .split_once(',')
        .ok_or_else(|| anyhow!("Malformed ML tag '{}'", tag))?;
    let mut result = Vec::new();
    for entry in values.split(',') {
        if entry.is_empty() {
            continue;
        }
        let value: u8 = entry
            .parse()
            .with_context(|| format!("Invalid ML value '{}'", entry))?;
        result.push(value);
    }
    Ok(result)
}

fn normalize_mm_segments(mm_tag: &str) -> Vec<String> {
    let trimmed = mm_tag.trim_start_matches("MM:Z:");
    trimmed
        .split(';')
        .filter_map(|seg| {
            let s = seg.trim();
            if s.is_empty() {
                None
            } else {
                Some(s.to_string())
            }
        })
        .collect()
}

fn parse_mm_segment(segment: &str) -> Option<(char, char, String, Vec<usize>)> {
    let marker_idx = segment.find(".,")?;
    let prefix = &segment[..marker_idx];
    let mut chars = prefix.chars();
    let base = chars.next()?;
    let strand = chars.next()?;
    let mod_code: String = chars.collect();
    let mut deltas = Vec::new();
    for part in segment[marker_idx + 2..].split(',') {
        if part.is_empty() {
            continue;
        }
        if let Ok(delta) = part.parse::<usize>() {
            deltas.push(delta);
        }
    }
    Some((base, strand, mod_code, deltas))
}

fn decode_positions(sequence: &str, base: char, deltas: &[usize]) -> Vec<usize> {
    let seq_bytes: Vec<u8> = sequence
        .as_bytes()
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let target = base.to_ascii_uppercase() as u8;
    let positions: Vec<usize> = seq_bytes
        .iter()
        .enumerate()
        .filter_map(|(idx, b)| if *b == target { Some(idx) } else { None })
        .collect();
    if positions.is_empty() {
        return Vec::new();
    }
    let mut decoded = Vec::new();
    let mut cursor = 0usize;
    for delta in deltas {
        cursor = cursor.saturating_add(*delta);
        if cursor >= positions.len() {
            break;
        }
        decoded.push(positions[cursor]);
        cursor = cursor.saturating_add(1);
    }
    decoded
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode_positions_matches_encoder_logic() {
        let seq = "ACGATCACGATC";
        // Deltas count intervening canonical bases between modification events.
        let deltas = vec![0, 1, 0];
        let positions = decode_positions(seq, 'A', &deltas);
        assert_eq!(positions, vec![0, 6, 9]);
    }

    #[test]
    fn empirical_bins_cover_full_range() {
        let values = vec![0.0, 0.1, 0.2, 0.99];
        let emp = EmpiricalDistribution::from_values(&values, 0.2);
        assert_eq!(emp.bins.len(), 5);
        assert_eq!(emp.bins.iter().sum::<usize>(), values.len());
    }
}
