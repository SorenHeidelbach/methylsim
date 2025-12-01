use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};
use log::{info, debug};
use rand::Rng;
use rand::rngs::StdRng;
use rand::Rng;
use rand_distr::{Beta, Distribution};
use serde::{Deserialize, Serialize};

use crate::fastx::{process_fastq_streaming, ReadRecord};
use crate::motif::MotifDefinition;

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct ModelKey(String);

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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strand: Option<char>,
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

type EventMap = HashMap<(char, String, char), HashMap<usize, f64>>;

const MAX_PARSE_WARNINGS: usize = 20;

impl ModelKey {
    pub fn normalize(raw: &str) -> Self {
        if let Some(stripped) = raw.strip_suffix("_+") {
            return ModelKey(stripped.to_string());
        }
        if let Some(stripped) = raw.strip_suffix("_-") {
            return ModelKey(stripped.to_string());
        }
        ModelKey(raw.to_string())
    }

    #[allow(dead_code)]
    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub fn into_string(self) -> String {
        self.0
    }
}

struct ParseWarningTracker {
    count: usize,
    max_warnings: usize,
}

impl ParseWarningTracker {
    fn new() -> Self {
        Self {
            count: 0,
            max_warnings: MAX_PARSE_WARNINGS,
        }
    }

    fn record_failure(&mut self, read_id: &str, error: &anyhow::Error) {
        if self.count < self.max_warnings {
            eprintln!(
                "warning: failed to parse modifications for read {}: {}",
                read_id, error
            );
        } else if self.count == self.max_warnings {
            eprintln!("warning: additional modification parse failures suppressed");
        }
        self.count += 1;
    }

    fn report(&self) {
        if self.count > 0 {
            eprintln!(
                "warning: skipped {} reads with incompatible MM/ML tags (counts suppressed after first {})",
                self.count, self.max_warnings
            );
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SkipMode {
    /// Explicit mode (?): only report modifications, not unmodified canonical bases
    Explicit,
    /// Implicit mode (.): report both modifications and unmodified canonical bases
    ImplicitUnmodified,
    /// Default: behave like ImplicitUnmodified
    DefaultImplicitUnmodified,
}

impl SkipMode {
    fn is_implicit(self) -> bool {
        matches!(
            self,
            SkipMode::ImplicitUnmodified | SkipMode::DefaultImplicitUnmodified
        )
    }
}

pub fn learn_model_from_reads(
    reads_path: &Path,
    motifs: &[MotifDefinition],
    threshold: f64,
    max_reads: Option<usize>,
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

    let mut reads_with_tags = 0;
    let mut reads_without_tags = 0;
    let mut total_reads = 0;
    let mut warning_tracker = ParseWarningTracker::new();

    // Process reads in streaming fashion to minimize memory usage
    process_fastq_streaming(reads_path, true, max_reads, |read| {
        total_reads += 1;

        // Find all motif positions in this read first
        let mut motif_positions: HashMap<(char, String, char), HashSet<usize>> = HashMap::new();
        for motif in motifs {
            let key = motif.model_key();
            let hits = motif.find_matches(read.sequence_str());
            *totals.entry(key.clone()).or_default() += hits.len();

            // Store positions for this motif's base/mod/strand combination
            let group_key = (motif.canonical_base, motif.mod_code.clone(), motif.strand);
            motif_positions.entry(group_key).or_default().extend(hits);
        }

        // Parse tags only for positions that match motifs - gracefully handle errors
        let parsed = match parse_tagged_events_filtered(read, &motif_positions) {
            Ok(Some(events)) => {
                reads_with_tags += 1;
                Some(events)
            }
            Ok(None) => {
                reads_without_tags += 1;
                None
            }
            Err(e) => {
                warning_tracker.record_failure(&read.name, &e);
                reads_without_tags += 1;
                None
            }
        };

        // Collect probabilities for each motif
        for motif in motifs {
            let key = motif.model_key();
            let Some(event_map) = parsed.as_ref() else {
                continue;
            };
            let Some(probs) =
                event_map.get(&(motif.canonical_base, motif.mod_code.clone(), motif.strand))
            else {
                continue;
            };
            for prob in probs.values() {
                values.entry(key.clone()).or_default().push(*prob);
            }
        }

        Ok(())
    })?;

    // Report parse failures summary
    warning_tracker.report();

    eprintln!(
        "  - {} reads with MM/ML tags ({:.1}%)",
        reads_with_tags,
        if total_reads > 0 {
            (reads_with_tags as f64 / total_reads as f64) * 100.0
        } else {
            0.0
        }
    );
    eprintln!(
        "  - {} reads without MM/ML tags ({:.1}%)",
        reads_without_tags,
        if total_reads > 0 {
            (reads_without_tags as f64 / total_reads as f64) * 100.0
        } else {
            0.0
        }
    );

    let mut motif_models = Vec::new();
    eprintln!("\nFitting model for {} motif(s)", motifs.len());
    for motif in motifs {
        eprintln!("Processing motif: {}", motif.motif);
        let key = motif.model_key();
        let covered: Vec<f64> = values.remove(&key).unwrap_or_default();
        let total_sites = totals.remove(&key).unwrap_or_default();

        eprintln!(
            "  - Found {} total motif sites across all reads",
            total_sites
        );
        eprintln!(
            "  - {} sites have methylation probability data ({:.1}%)",
            covered.len(),
            if total_sites > 0 {
                (covered.len() as f64 / total_sites as f64) * 100.0
            } else {
                0.0
            }
        );

        let methylated: Vec<f64> = covered
            .iter()
            .copied()
            .filter(|p| *p >= threshold)
            .collect();

        eprintln!(
            "  - {} sites considered methylated (probability >= {:.2})",
            methylated.len(),
            threshold
        );
        let unmethylated: Vec<f64> = covered.iter().copied().filter(|p| *p < threshold).collect();
        let meth_beta = fit_beta(&methylated);
        let unmeth_beta = fit_beta(&unmethylated);
        let methylated_fraction = if covered.is_empty() {
            0.0
        } else {
            methylated.len() as f64 / covered.len() as f64
        };
        let empirical = EmpiricalDistribution::from_values(&covered, 0.05);

        // Log per-motif statistics
        let coverage_pct = if total_sites > 0 {
            covered.len() as f64 / total_sites as f64 * 100.0
        } else {
            0.0
        };
        info!(
            "  {} ({}): {} total sites, {} covered ({:.1}%), {:.1}% methylated",
            motif.motif,
            key,
            total_sites,
            covered.len(),
            coverage_pct,
            methylated_fraction * 100.0
        );
        debug!(
            "    Methylated: {} sites, Unmethylated: {} sites",
            methylated.len(),
            unmethylated.len()
        );

        motif_models.push(MotifModel {
            key,
            motif: motif.motif.clone(),
            canonical_base: motif.canonical_base,
            canonical_offset: motif.canonical_offset,
            mod_code: motif.mod_code.clone(),
            strand: None,
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
    /// Extract motif specifications from the model
    pub fn extract_motif_specs(&self) -> Vec<String> {
        self.motifs
            .iter()
            .map(|m| {
                // Format: motif:canonical:offset:mod:strand
                format!(
                    "{}:{}:{}:{}",
                    m.motif, m.canonical_base, m.canonical_offset, m.mod_code
                )
            })
            .collect()
    }

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
            let key = ModelKey::normalize(&motif.key);
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
                key.into_string(),
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

/// Parse tagged events only for positions that match motifs (memory-optimized)
fn parse_tagged_events_filtered(
    read: &ReadRecord,
    motif_positions: &HashMap<(char, String, char), HashSet<usize>>,
) -> Result<Option<EventMap>> {
    let comment = match read.comment.as_deref() {
        Some(c) if !c.is_empty() => c,
        _ => return Ok(None),
    };
    let (mm_tag, ml_tag) = extract_mm_ml(comment);
    let Some(mm_tag) = mm_tag else {
        return Ok(None);
    };
    // Parse ML values - return error for warning tracker
    let ml_values = parse_ml_values(ml_tag).context("Invalid ML tag format")?;

    let mut ml_iter = ml_values.into_iter();
    let segments = normalize_mm_segments(mm_tag);
    if segments.is_empty() {
        return Ok(None);
    }

    let sequence = read.sequence_str();
    let mut grouped: EventMap = HashMap::new();
    let mut total_modifications = 0;

    for segment in segments {
        let (base, strand, mod_code, deltas, skip_mode) = match parse_mm_segment(&segment) {
            Some(parsed) => parsed,
            None => {
                bail!("Invalid MM segment format: {}", segment);
            }
        };

        // Check if we care about this base/strand/mod combination
        let group_key = (base, mod_code.clone(), strand);
        let is_target_group = motif_positions.contains_key(&group_key);

        // Get all canonical positions for tracking used positions
        let canonical_positions = get_canonical_positions(sequence, base);

        // Decode modification positions with validation
        let mod_positions = decode_positions_validated(sequence, base, &deltas, &segment)?;

        total_modifications += mod_positions.len();

        // Track which canonical positions are explicitly modified
        let mut used_canonical_indices = vec![false; canonical_positions.len()];
        for &mod_pos in &mod_positions {
            if let Some(canonical_idx) = canonical_positions.iter().position(|&p| p == mod_pos) {
                used_canonical_indices[canonical_idx] = true;
            }
        }

        if !is_target_group {
            // Skip ML values for this segment since we don't care about it
            for _ in &mod_positions {
                ml_iter.next();
            }
            continue;
        }

        let target_positions = motif_positions.get(&group_key).unwrap();

        // Process explicit modifications
        for pos in mod_positions {
            match ml_iter.next() {
                Some(ml) => {
                    // Only store if this position matches a motif
                    if target_positions.contains(&pos) {
                        grouped
                            .entry(group_key.clone())
                            .or_default()
                            .insert(pos, (ml as f64) / 255.0);
                    }
                }
                None => {
                    bail!(
                        "ML tag has fewer values than MM modifications (expected at least {} more)",
                        total_modifications - (ml_iter.len() + 1)
                    );
                }
            }
        }

        // Handle implicit unmodified bases if skip mode is implicit
        if skip_mode.is_implicit() {
            for (canonical_idx, &pos) in canonical_positions.iter().enumerate() {
                if !used_canonical_indices[canonical_idx] && target_positions.contains(&pos) {
                    // This canonical base is unmodified (implicitly prob=0)
                    grouped
                        .entry(group_key.clone())
                        .or_default()
                        .insert(pos, 0.0);
                }
            }
        }
    }

    // Validate ML tag length
    let remaining_ml = ml_iter.count();
    if remaining_ml > 0 {
        bail!(
            "ML tag contains {} extra values beyond {} MM modifications",
            remaining_ml,
            total_modifications
        );
    }
    if grouped.is_empty() {
        Ok(None)
    } else {
        Ok(Some(grouped))
    }
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

fn parse_mm_segment(segment: &str) -> Option<(char, char, String, Vec<usize>, SkipMode)> {
    // Split by comma to separate header from deltas
    let mut parts = segment.split(',');
    let header = parts.next()?;

    // Parse header: BASE+STRAND+CODE[SKIPMODE]
    let mut chars = header.chars();
    let base = chars.next()?;
    let strand = chars.next()?;
    let mut mod_code: String = chars.collect();

    // Parse skip mode from modification code suffix
    let mut skip_mode = SkipMode::DefaultImplicitUnmodified;
    if let Some(last) = mod_code.chars().last() {
        match last {
            '?' => {
                skip_mode = SkipMode::Explicit;
                mod_code.pop();
            }
            '.' => {
                skip_mode = SkipMode::ImplicitUnmodified;
                mod_code.pop();
            }
            _ => {}
        }
    }

    if mod_code.is_empty() {
        return None;
    }

    // Parse deltas (may be empty)
    let mut deltas = Vec::new();
    for part in parts {
        let trimmed = part.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Ok(delta) = trimmed.parse::<usize>() {
            deltas.push(delta);
        }
    }

    Some((base, strand, mod_code, deltas, skip_mode))
}

#[allow(dead_code)]
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

/// Get canonical base positions in sequence
fn get_canonical_positions(sequence: &str, base: char) -> Vec<usize> {
    let target = base.to_ascii_uppercase() as u8;
    sequence
        .as_bytes()
        .iter()
        .enumerate()
        .filter_map(|(idx, b)| {
            if b.to_ascii_uppercase() == target {
                Some(idx)
            } else {
                None
            }
        })
        .collect()
}

/// Decode positions with validation to provide better error messages
fn decode_positions_validated(
    sequence: &str,
    base: char,
    deltas: &[usize],
    segment: &str,
) -> Result<Vec<usize>> {
    let positions = get_canonical_positions(sequence, base);

    if positions.is_empty() {
        bail!(
            "MM segment for base {} but no occurrences found in sequence",
            base
        );
    }

    let mut decoded = Vec::new();
    let mut cursor: isize = -1;

    for delta in deltas {
        cursor += *delta as isize + 1;

        if cursor < 0 {
            bail!(
                "Negative modification offset encountered in MM segment: {}",
                segment
            );
        }

        let canonical_idx = cursor as usize;
        if canonical_idx >= positions.len() {
            bail!(
                "MM segment for base {} references canonical index {} but only {} occurrences found",
                base,
                canonical_idx,
                positions.len()
            );
        }

        decoded.push(positions[canonical_idx]);
    }

    Ok(decoded)
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

    #[test]
    fn validation_detects_out_of_range_positions() {
        let seq = "ACGT";
        let deltas = vec![0, 10]; // Second position is way beyond what exists
        let result = decode_positions_validated(seq, 'A', &deltas, "A+m,0,10");
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("only 1 occurrences"));
    }

    #[test]
    fn validation_detects_missing_canonical_base() {
        let seq = "TGCT";
        let deltas = vec![0];
        let result = decode_positions_validated(seq, 'A', &deltas, "A+m,0");
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("no occurrences found"));
    }

    #[test]
    fn validation_succeeds_for_valid_positions() {
        let seq = "ACGATCG";
        // Delta 0 means modify first A (pos 0), delta 0 means skip 0 and modify next A (pos 3)
        let deltas = vec![0, 0];
        let result = decode_positions_validated(seq, 'A', &deltas, "A+a,0,0");
        assert!(result.is_ok(), "Expected Ok, got: {:?}", result.err());
        let positions = result.unwrap();
        // Should get both A positions: 0 and 3
        assert_eq!(positions, vec![0, 3]);
    }

    #[test]
    fn parse_segment_detects_explicit_skip_mode() {
        let result = parse_mm_segment("C+m?,0,1");
        assert!(result.is_some());
        let (base, strand, mod_code, deltas, skip_mode) = result.unwrap();
        assert_eq!(base, 'C');
        assert_eq!(strand, '+');
        assert_eq!(mod_code, "m");
        assert_eq!(deltas, vec![0, 1]);
        assert_eq!(skip_mode, SkipMode::Explicit);
    }

    #[test]
    fn parse_segment_detects_implicit_skip_mode() {
        let result = parse_mm_segment("C+m.,0,1");
        assert!(result.is_some());
        let (_, _, _, _, skip_mode) = result.unwrap();
        assert_eq!(skip_mode, SkipMode::ImplicitUnmodified);
    }

    #[test]
    fn parse_segment_defaults_to_implicit_skip_mode() {
        let result = parse_mm_segment("C+m,0,1");
        assert!(result.is_some());
        let (_, _, _, _, skip_mode) = result.unwrap();
        assert_eq!(skip_mode, SkipMode::DefaultImplicitUnmodified);
        assert!(skip_mode.is_implicit());
    }

    #[test]
    fn parse_segment_handles_no_deltas_with_skip_mode() {
        // Real-world case: C+21839. means 4mC implicit mode with no modifications
        let result = parse_mm_segment("C+21839.");
        assert!(result.is_some());
        let (base, strand, mod_code, deltas, skip_mode) = result.unwrap();
        assert_eq!(base, 'C');
        assert_eq!(strand, '+');
        assert_eq!(mod_code, "21839");
        assert_eq!(deltas, Vec::<usize>::new());
        assert_eq!(skip_mode, SkipMode::ImplicitUnmodified);
    }

    #[test]
    fn parse_segment_handles_no_deltas_default() {
        let result = parse_mm_segment("C+m");
        assert!(result.is_some());
        let (_, _, mod_code, deltas, skip_mode) = result.unwrap();
        assert_eq!(mod_code, "m");
        assert_eq!(deltas, Vec::<usize>::new());
        assert_eq!(skip_mode, SkipMode::DefaultImplicitUnmodified);
    }
}
