use std::collections::{BTreeMap, HashMap, VecDeque};

use anyhow::{bail, Result};
use rand::{rngs::StdRng, Rng, SeedableRng};
use rand_distr::{Beta, Distribution};

use crate::model::MotifSampler;
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

#[derive(Clone)]
struct EventCandidate {
    position: usize,
    probability: f64,
    canonical_base: char,
    mod_code: String,
    strand: char,
    is_motif: bool,
}

struct MotifGroup {
    canonical_base: char,
    mod_code: String,
    strand: char,
    motifs: Vec<MotifDefinition>,
    sampler: Option<MotifSampler>,
}

fn normalize_mod_code_for_mm(code: &str) -> String {
    let lower = code.to_ascii_lowercase();
    match lower.as_str() {
        "5mc" | "m" => "m".to_string(),
        "6ma" | "a" => "a".to_string(),
        "4mc" | "21839" => "21839".to_string(),
        _other => code.to_string(),
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum SiteKind {
    Motif,
    Background,
}

pub trait ProbabilitySampler {
    fn sample_probability(
        &mut self,
        site: SiteKind,
        model_sampler: Option<&MotifSampler>,
        mod_code: &str,
        rng: &mut StdRng,
    ) -> f64;
}

struct MixtureProbabilitySampler {
    high_beta: Beta<f64>,
    low_betas: HashMap<String, Beta<f64>>,
    default_low_beta: Beta<f64>,
    motif_high_prob: f64,
    background_high_prob: f64,
}

impl MixtureProbabilitySampler {
    fn new(
        high_beta: Beta<f64>,
        low_betas: HashMap<String, Beta<f64>>,
        default_low_beta: Beta<f64>,
        motif_high_prob: f64,
        background_high_prob: f64,
    ) -> Self {
        Self {
            high_beta,
            low_betas,
            default_low_beta,
            motif_high_prob,
            background_high_prob,
        }
    }
}

pub struct MethylationTagger {
    groups: Vec<MotifGroup>,
    probability_sampler: Box<dyn ProbabilitySampler + Send>,
    rng: StdRng,
}

#[derive(Clone, Debug)]
struct DefaultBetas {
    high: (f64, f64),
    fallback_low: (f64, f64),
    mod_low: HashMap<String, (f64, f64)>,
}

impl Default for DefaultBetas {
    fn default() -> Self {
        let high = (2.8685813093772063, 0.055071333309938346);
        let fallback_low = (0.3451168730448373, 5.040038112000441); // 6mA defaults

        let mut mod_low = HashMap::new();
        // Adenine methylation (6mA)
        let sixma = (0.3451168730448373, 5.040038112000441);
        mod_low.insert("6mA".to_ascii_lowercase(), sixma);
        mod_low.insert("a".to_string(), sixma);
        mod_low.insert("6ma".to_string(), sixma);
        // Cytosine methylation (5mC / 4mC)
        let fivemc = (1.3799717363106534, 8.927216036532014);
        mod_low.insert("5mC".to_ascii_lowercase(), fivemc);
        mod_low.insert("5mc".to_string(), fivemc);
        mod_low.insert("m".to_string(), fivemc);

        let forumc = (0.8124532876545941, 9.107388497092732);
        mod_low.insert("4mC".to_ascii_lowercase(), forumc);
        mod_low.insert("4mc".to_string(), forumc);
        mod_low.insert("c".to_string(), forumc);
        mod_low.insert("21839".to_string(), forumc);

        Self {
            high,
            fallback_low,
            mod_low,
        }
    }
}

impl DefaultBetas {
    fn high(&self) -> (f64, f64) {
        self.high
    }

    fn low_for(&self, mod_code: &str, fallback: (f64, f64)) -> (f64, f64) {
        let key = mod_code.to_ascii_lowercase();
        self.mod_low.get(&key).copied().unwrap_or(fallback)
    }

    fn fallback_low(&self) -> (f64, f64) {
        self.fallback_low
    }
}

pub struct MethylationTaggerBuilder {
    motifs: Vec<MotifDefinition>,
    motif_high_prob: f64,
    background_high_prob: f64,
    sampler_map: Option<HashMap<String, MotifSampler>>,
    seed: u64,
    high_beta_alpha: f64,
    high_beta_beta: f64,
    low_beta_alpha: f64,
    low_beta_beta: f64,
}

impl MethylationTaggerBuilder {
    pub fn new(motifs: Vec<MotifDefinition>) -> Self {
        Self {
            motifs,
            motif_high_prob: 0.95,
            background_high_prob: 0.01,
            sampler_map: None,
            seed: 1,
            high_beta_alpha: DefaultBetas::default().high().0,
            high_beta_beta: DefaultBetas::default().high().1,
            low_beta_alpha: DefaultBetas::default().fallback_low().0,
            low_beta_beta: DefaultBetas::default().fallback_low().1,
        }
    }

    pub fn with_probs(mut self, motif_high: f64, background_high: f64) -> Self {
        self.motif_high_prob = motif_high;
        self.background_high_prob = background_high;
        self
    }

    pub fn with_sampler_map(mut self, sampler_map: Option<HashMap<String, MotifSampler>>) -> Self {
        self.sampler_map = sampler_map;
        self
    }

    pub fn with_betas(
        mut self,
        high_alpha: f64,
        high_beta: f64,
        low_alpha: f64,
        low_beta: f64,
    ) -> Self {
        self.high_beta_alpha = high_alpha;
        self.high_beta_beta = high_beta;
        self.low_beta_alpha = low_alpha;
        self.low_beta_beta = low_beta;
        self
    }

    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }

    pub fn build(self) -> Result<MethylationTagger> {
        if self.motifs.is_empty() {
            bail!("At least one motif definition must be provided");
        }
        if !(0.0..=1.0).contains(&self.motif_high_prob) {
            bail!(
                "Motif high probability must be between 0 and 1 (received {})",
                self.motif_high_prob
            );
        }
        if !(0.0..=1.0).contains(&self.background_high_prob) {
            bail!(
                "Background high probability must be between 0 and 1 (received {})",
                self.background_high_prob
            );
        }

        let defaults = DefaultBetas::default();
        let high_beta = build_beta(self.high_beta_alpha, self.high_beta_beta)?;
        let mut low_betas: HashMap<String, Beta<f64>> = HashMap::new();
        for motif in &self.motifs {
            let params =
                defaults.low_for(&motif.mod_code, (self.low_beta_alpha, self.low_beta_beta));
            low_betas
                .entry(motif.mod_code.clone())
                .or_insert(build_beta(params.0, params.1)?);
        }
        let default_low_beta = build_beta(self.low_beta_alpha, self.low_beta_beta)?;

        let probability_sampler: Box<dyn ProbabilitySampler + Send> =
            Box::new(MixtureProbabilitySampler::new(
                high_beta,
                low_betas,
                default_low_beta,
                self.motif_high_prob,
                self.background_high_prob,
            ));

        let model_map = self.sampler_map.unwrap_or_default();
        let mut grouped: BTreeMap<(char, String, char), Vec<MotifDefinition>> = BTreeMap::new();
        for motif in self.motifs {
            let key = (motif.canonical_base, motif.mod_code.clone(), motif.strand);
            grouped.entry(key).or_default().push(motif);
        }
        let groups = grouped
            .into_iter()
            .map(|((base, mod_code, strand), motifs)| {
                let sampler = motifs
                    .get(0)
                    .and_then(|motif| model_map.get(&motif.model_key()).cloned());
                MotifGroup {
                    canonical_base: base,
                    mod_code,
                    strand,
                    motifs,
                    sampler,
                }
            })
            .collect();

        Ok(MethylationTagger {
            groups,
            probability_sampler,
            rng: StdRng::seed_from_u64(self.seed),
        })
    }
}

impl MethylationTagger {
    #[allow(dead_code)]
    pub fn new(
        motifs: Vec<MotifDefinition>,
        motif_high_prob: f64,
        background_high_prob: f64,
        high_beta_alpha: f64,
        high_beta_beta: f64,
        low_beta_alpha: f64,
        low_beta_beta: f64,
        seed: u64,
        motif_models: Option<HashMap<String, MotifSampler>>,
    ) -> Result<Self> {
        MethylationTaggerBuilder::new(motifs)
            .with_probs(motif_high_prob, background_high_prob)
            .with_betas(
                high_beta_alpha,
                high_beta_beta,
                low_beta_alpha,
                low_beta_beta,
            )
            .with_seed(seed)
            .with_sampler_map(motif_models)
            .build()
    }

    pub fn annotate(&mut self, sequence: &str) -> TagResult {
        let seq_bytes: Vec<u8> = sequence
            .as_bytes()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        let mut candidates = Vec::new();
        for idx in 0..self.groups.len() {
            let match_flags =
                Self::collect_match_flags(&self.groups[idx], sequence, seq_bytes.len());
            self.collect_events(idx, &seq_bytes, &match_flags, &mut candidates);
        }
        if candidates.is_empty() {
            return TagResult::default();
        }

        // Deduplicate per canonical position: keep highest-probability event
        use std::collections::HashMap;
        let mut best: HashMap<usize, EventCandidate> = HashMap::new();
        for cand in candidates {
            let entry = best.entry(cand.position).or_insert_with(|| cand.clone());
            if cand.probability > entry.probability {
                *entry = cand;
            }
        }
        let mut final_events = Vec::with_capacity(best.len());
        let mut motif_hit_count = 0usize;
        let mut motif_high_count = 0usize;
        for cand in best.into_values() {
            if cand.is_motif {
                motif_hit_count += 1;
                motif_high_count += 1; // all kept candidates are modified
            }
            final_events.push(ModificationEvent {
                position: cand.position,
                ml_value: (cand.probability * 255.0).round().clamp(0.0, 255.0) as u8,
                canonical_base: cand.canonical_base,
                mod_code: cand.mod_code,
                strand: cand.strand,
            });
        }

        let mut result = self.build_tags(sequence, final_events);
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
        candidates: &mut Vec<EventCandidate>,
    ) {
        let (canonical_base, strand, mod_code, sampler) = {
            let group = &self.groups[group_index];
            (
                group.canonical_base,
                group.strand,
                normalize_mod_code_for_mm(&group.mod_code),
                group.sampler.clone(),
            )
        };
        let target_base = canonical_base.to_ascii_uppercase() as u8;
        for (idx, base) in seq_bytes.iter().enumerate() {
            if *base != target_base {
                continue;
            }
            let is_motif = match_flags.get(idx).copied().unwrap_or(false);
            let probability = self.probability_sampler.sample_probability(
                if is_motif {
                    SiteKind::Motif
                } else {
                    SiteKind::Background
                },
                sampler.as_ref(),
                &mod_code,
                &mut self.rng,
            );
            let is_modified = self.sample_high(probability);
            if is_motif {
                // stats computed after deduplication
            }
            if !is_modified {
                continue;
            }

            candidates.push(EventCandidate {
                position: idx,
                probability,
                canonical_base,
                mod_code: mod_code.clone(),
                strand,
                is_motif,
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
        grouped.sort_by(|(key_a, ev_a), (key_b, ev_b)| {
            let (base_a, mod_a, _) = key_a;
            let (base_b, mod_b, _) = key_b;
            base_a.cmp(base_b).then(mod_a.cmp(mod_b)).then(
                ev_a.iter()
                    .map(|(pos, _)| *pos)
                    .min()
                    .unwrap_or(0)
                    .cmp(&ev_b.iter().map(|(pos, _)| *pos).min().unwrap_or(0)),
            )
        });

        // If multiple mod codes share the same base/strand, use explicit skip mode ('?') to avoid inference conflicts
        let mut base_counts: std::collections::HashMap<(char, char), usize> =
            std::collections::HashMap::new();
        for ((base, _, strand), evts) in &grouped {
            if !evts.is_empty() {
                *base_counts.entry((*base, *strand)).or_insert(0) += 1;
            }
        }
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
                let skip_char = if base_counts.get(&(base, strand)).copied().unwrap_or(0) > 1 {
                    '?'
                } else {
                    '.'
                };
                let delta_str = deltas
                    .iter()
                    .map(|d| d.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                mm_segments.push(format!("{base}{strand}{mod_code}{skip_char},{delta_str};"));
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

fn build_beta(alpha: f64, beta: f64) -> Result<Beta<f64>> {
    Ok(Beta::new(alpha.max(1e-6), beta.max(1e-6))?)
}

impl ProbabilitySampler for MixtureProbabilitySampler {
    fn sample_probability(
        &mut self,
        site: SiteKind,
        model_sampler: Option<&MotifSampler>,
        mod_code: &str,
        rng: &mut StdRng,
    ) -> f64 {
        if let Some(model) = model_sampler {
            return model.sample_probability(rng);
        }
        let pick_high = match site {
            SiteKind::Motif => rng.gen_bool(self.motif_high_prob.clamp(0.0, 1.0)),
            SiteKind::Background => rng.gen_bool(self.background_high_prob.clamp(0.0, 1.0)),
        };
        let prob = if pick_high {
            self.high_beta.sample(rng)
        } else {
            self.low_betas
                .get(mod_code)
                .unwrap_or(&self.default_low_beta)
                .sample(rng)
        };
        prob.clamp(0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::motif::MotifDefinition;

    #[test]
    fn annotate_skips_unmodified_sites() {
        struct ZeroSampler;
        impl ProbabilitySampler for ZeroSampler {
            fn sample_probability(
                &mut self,
                _site: SiteKind,
                _model_sampler: Option<&MotifSampler>,
                _mod_code: &str,
                _rng: &mut StdRng,
            ) -> f64 {
                0.0
            }
        }

        let motif = MotifDefinition::parse("A:A:0:m:+").expect("parse motif");
        let mut tagger = MethylationTaggerBuilder::new(vec![motif])
            .with_seed(1)
            .build()
            .expect("build tagger");
        tagger.probability_sampler = Box::new(ZeroSampler);

        let result = tagger.annotate("AAAA");
        assert!(
            result.mm_tag.is_none() && result.ml_tag.is_none(),
            "Tags should be absent when no sites are modified"
        );
    }

    #[test]
    fn annotate_uses_numeric_code_for_4mc() {
        struct OneSampler;
        impl ProbabilitySampler for OneSampler {
            fn sample_probability(
                &mut self,
                _site: SiteKind,
                _model_sampler: Option<&MotifSampler>,
                _mod_code: &str,
                _rng: &mut StdRng,
            ) -> f64 {
                1.0
            }
        }

        let motif = MotifDefinition::parse("CCGG:C:1:4mC:+").expect("parse motif");
        let mut tagger = MethylationTaggerBuilder::new(vec![motif])
            .with_seed(1)
            .build()
            .expect("build tagger");
        tagger.probability_sampler = Box::new(OneSampler);

        let result = tagger.annotate("CCGGCCGG");
        let mm = result.mm_tag.expect("mm tag");
        assert!(
            mm.contains("21839"),
            "Expected numeric 21839 code, got {mm}"
        );
    }
}
