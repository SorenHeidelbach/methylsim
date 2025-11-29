use std::collections::{BTreeMap, VecDeque};

use anyhow::{Result, bail};
use rand::{Rng, SeedableRng, rngs::StdRng};
use rand_distr::{Distribution, Normal};

use crate::motif::MotifDefinition;

#[derive(Debug, Clone, Default)]
pub struct TagResult {
    pub mm_tag: Option<String>,
    pub ml_tag: Option<String>,
    pub motif_hit_count: usize,
    pub motif_high_count: usize,
}

impl TagResult {}

struct ModificationEvent {
    position: usize,
    ml_value: u8,
    canonical_base: char,
    mod_code: String,
    strand: char,
}

struct MotifGroup {
    canonical_base: char,
    mod_code: String,
    strand: char,
    motifs: Vec<MotifDefinition>,
}

pub struct MethylationTagger {
    groups: Vec<MotifGroup>,
    motif_high_prob: f64,
    background_high_prob: f64,
    rng: StdRng,
    high_ml_mean: f64,
    high_ml_normal: Option<Normal<f64>>,
    low_ml_mean: f64,
    low_ml_normal: Option<Normal<f64>>,
}

impl MethylationTagger {
    pub fn new(
        motifs: Vec<MotifDefinition>,
        motif_high_prob: f64,
        background_high_prob: f64,
        high_ml_mean: f64,
        high_ml_std: f64,
        low_ml_mean: f64,
        low_ml_std: f64,
        seed: u64,
    ) -> Result<Self> {
        if motifs.is_empty() {
            bail!("At least one motif definition must be provided");
        }
        if !(0.0..=1.0).contains(&motif_high_prob) {
            bail!("Motif high probability must be between 0 and 1 (received {motif_high_prob})");
        }
        if !(0.0..=1.0).contains(&background_high_prob) {
            bail!(
                "Background high probability must be between 0 and 1 (received {background_high_prob})"
            );
        }
        if !(0.0..=255.0).contains(&high_ml_mean) {
            bail!("High ML mean must be between 0 and 255 (received {high_ml_mean})");
        }
        if high_ml_std.is_sign_negative() {
            bail!("High ML standard deviation must be non-negative");
        }
        if !(0.0..=255.0).contains(&low_ml_mean) {
            bail!("Low ML mean must be between 0 and 255 (received {low_ml_mean})");
        }
        if low_ml_std.is_sign_negative() {
            bail!("Low ML standard deviation must be non-negative");
        }
        let high_ml_normal = if high_ml_std > 0.0 {
            Some(Normal::new(high_ml_mean, high_ml_std)?)
        } else {
            None
        };
        let low_ml_normal = if low_ml_std > 0.0 {
            Some(Normal::new(low_ml_mean, low_ml_std)?)
        } else {
            None
        };
        let mut grouped: BTreeMap<(char, String, char), Vec<MotifDefinition>> = BTreeMap::new();
        for motif in motifs {
            let key = (motif.canonical_base, motif.mod_code.clone(), motif.strand);
            grouped.entry(key).or_default().push(motif);
        }
        let groups = grouped
            .into_iter()
            .map(|((base, mod_code, strand), motifs)| MotifGroup {
                canonical_base: base,
                mod_code,
                strand,
                motifs,
            })
            .collect();
        Ok(Self {
            groups,
            motif_high_prob,
            background_high_prob,
            rng: StdRng::seed_from_u64(seed),
            high_ml_mean,
            high_ml_normal,
            low_ml_mean,
            low_ml_normal,
        })
    }

    pub fn annotate(&mut self, sequence: &str) -> TagResult {
        let seq_bytes: Vec<u8> = sequence
            .as_bytes()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        let mut events = Vec::new();
        let mut motif_hit_count = 0usize;
        let mut motif_high_count = 0usize;
        for idx in 0..self.groups.len() {
            let match_flags =
                Self::collect_match_flags(&self.groups[idx], sequence, seq_bytes.len());
            self.collect_events(
                idx,
                &seq_bytes,
                &match_flags,
                &mut events,
                &mut motif_hit_count,
                &mut motif_high_count,
            );
        }
        if events.is_empty() {
            return TagResult::default();
        }
        let mut result = self.build_tags(sequence, events);
        result.motif_hit_count = motif_hit_count;
        result.motif_high_count = motif_high_count;
        result
    }

    fn collect_match_flags(group: &MotifGroup, sequence: &str, seq_len: usize) -> Vec<bool> {
        let mut flags = vec![false; seq_len];
        for motif in &group.motifs {
            for position in motif.find_matches(sequence) {
                if position < seq_len {
                    flags[position] = true;
                }
            }
        }
        flags
    }

    fn collect_events(
        &mut self,
        group_index: usize,
        seq_bytes: &[u8],
        match_flags: &[bool],
        events: &mut Vec<ModificationEvent>,
        motif_hit_count: &mut usize,
        motif_high_count: &mut usize,
    ) {
        let (canonical_base, strand, mod_code) = {
            let group = &self.groups[group_index];
            (group.canonical_base, group.strand, group.mod_code.clone())
        };
        let target_base = canonical_base.to_ascii_uppercase() as u8;
        for (idx, base) in seq_bytes.iter().enumerate() {
            if *base != target_base {
                continue;
            }
            let is_motif = match_flags.get(idx).copied().unwrap_or(false);
            let probability = if is_motif {
                self.motif_high_prob
            } else {
                self.background_high_prob
            };
            let is_high = self.sample_high(probability);
            let ml_value = self.sample_ml_value(is_high);
            if is_motif {
                *motif_hit_count += 1;
                if is_high {
                    *motif_high_count += 1;
                }
            }
            events.push(ModificationEvent {
                position: idx,
                ml_value,
                canonical_base,
                mod_code: mod_code.clone(),
                strand,
            });
        }
    }

    fn sample_high(&mut self, probability: f64) -> bool {
        if probability <= 0.0 {
            false
        } else if probability >= 1.0 {
            true
        } else {
            self.rng.gen_bool(probability)
        }
    }

    fn sample_ml_value(&mut self, is_high: bool) -> u8 {
        let (mean, normal) = if is_high {
            (self.high_ml_mean, self.high_ml_normal.as_ref())
        } else {
            (self.low_ml_mean, self.low_ml_normal.as_ref())
        };
        let value = if let Some(dist) = normal {
            dist.sample(&mut self.rng)
        } else {
            mean
        };
        value.clamp(0.0, 255.0).round() as u8
    }

    fn build_tags(&mut self, sequence: &str, events: Vec<ModificationEvent>) -> TagResult {
        let mut groups: BTreeMap<(char, String, char), Vec<(usize, u8)>> = BTreeMap::new();
        for event in events {
            let key = (event.canonical_base, event.mod_code.clone(), event.strand);
            groups
                .entry(key)
                .or_default()
                .push((event.position, event.ml_value));
        }
        let mut grouped: Vec<_> = groups.into_iter().collect();
        grouped.sort_by_key(|(_, evts)| evts.iter().map(|(pos, _)| *pos).min().unwrap_or(0));
        let seq_bytes: Vec<u8> = sequence
            .as_bytes()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();

        let mut mm_segments = Vec::new();
        let mut ml_values: Vec<u8> = Vec::new();

        for ((base, mod_code, strand), mut evts) in grouped {
            evts.sort_by_key(|(pos, _)| *pos);
            let mut queue: VecDeque<(usize, u8)> = VecDeque::from(evts);
            let target_base = base.to_ascii_uppercase() as u8;
            let mut deltas: Vec<usize> = Vec::new();
            let mut skip = 0usize;
            let total_events = queue.len();
            let mut processed = 0usize;

            for (idx, seq_base) in seq_bytes.iter().enumerate() {
                if *seq_base != target_base {
                    continue;
                }
                if processed >= total_events {
                    break;
                }
                if let Some(&(pos, _)) = queue.front() {
                    if pos == idx {
                        while let Some((pos_inner, prob)) = queue.pop_front() {
                            if pos_inner != idx {
                                queue.push_front((pos_inner, prob));
                                break;
                            }
                            deltas.push(skip);
                            ml_values.push(prob);
                            skip = 0;
                            processed += 1;
                        }
                        continue;
                    }
                }
                skip += 1;
            }

            if !deltas.is_empty() {
                let delta_str = deltas
                    .iter()
                    .map(|d| d.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                mm_segments.push(format!("{base}{strand}{mod_code}.,{delta_str};"));
            }
        }

        if mm_segments.is_empty() {
            TagResult::default()
        } else {
            let mm_tag = format!("MM:Z:{}", mm_segments.join(""));
            let ml_tag = if ml_values.is_empty() {
                None
            } else {
                Some(format!(
                    "ML:B:C,{}",
                    ml_values
                        .iter()
                        .map(|v| v.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                ))
            };
            TagResult {
                mm_tag: Some(mm_tag),
                ml_tag,
                motif_hit_count: 0,
                motif_high_count: 0,
            }
        }
    }
}
