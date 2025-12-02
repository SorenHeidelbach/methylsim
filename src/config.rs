use std::collections::HashSet;
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use log::{debug, info, warn};

use crate::cli::{SimulateArgs, SimulatorKind};
use crate::model::{LearnedModel, MotifSampler};
use crate::motif::MotifDefinition;
use crate::presets::ModelPreset;
use crate::simulator::{
    load_references, BadreadsSimulationConfig, BadreadsStrategy, BuiltinSimulationConfig,
    BuiltinStrategy, ErrorProfile, ReadLengthSpec, SimulatorKindStrategy,
};
use crate::tagging::MethylationTaggerBuilder;

pub fn collect_motif_specs(
    motifs: &[String],
    motifs_file: Option<&PathBuf>,
    model: Option<&LearnedModel>,
) -> Result<Vec<String>> {
    debug!("Collecting motif specifications from command line, files, and model");
    let mut specs = motifs.to_vec();

    if let Some(path) = motifs_file {
        info!("Loading motifs from file: {}", path.display());
        let mut file_specs = crate::motif::load_motif_file(path)?;
        debug!("Loaded {} motif specs from file", file_specs.len());
        specs.append(&mut file_specs);
    }

    if specs.is_empty() {
        if let Some(model) = model {
            info!("No motifs specified, extracting from model");
            specs = model.extract_motif_specs();
            info!("Extracted {} motif(s) from model", specs.len());
        }
    }

    let mut filtered = Vec::with_capacity(specs.len());
    for spec in specs {
        if looks_like_motif_header(&spec) {
            warn!("Ignoring motif entry that looks like a header: {}", spec);
            continue;
        }
        filtered.push(spec);
    }
    if filtered.is_empty() {
        bail!("At least one motif specification must be supplied (via --motif, --motifs-file, or --model-in/--model-preset)");
    }
    info!("Using {} motif(s): {}", filtered.len(), filtered.join(", "));
    Ok(filtered)
}

pub fn parse_quantity(
    quantity: &str,
    reference_size: usize,
    mean_read_length: usize,
) -> Result<usize> {
    let quantity = quantity.trim().to_uppercase();
    debug!("Parsing quantity specification: {}", quantity);

    if quantity.ends_with('X') {
        let coverage_str = &quantity[..quantity.len() - 1];
        let coverage: f64 = coverage_str.parse().with_context(|| {
            format!(
                "Invalid coverage value in '{}'. Expected format like '25x'",
                quantity
            )
        })?;

        if coverage <= 0.0 {
            bail!("Coverage must be positive, got: {}", coverage);
        }

        let num_reads =
            ((coverage * reference_size as f64) / mean_read_length as f64).ceil() as usize;

        info!(
            "Coverage-based quantity: {}x coverage → {} reads (ref_size: {}, read_length: {})",
            coverage, num_reads, reference_size, mean_read_length
        );

        Ok(num_reads)
    } else {
        let has_suffix =
            quantity.ends_with('K') || quantity.ends_with('M') || quantity.ends_with('G');
        if has_suffix {
            let (num_str, multiplier) = if quantity.ends_with('K') {
                (&quantity[..quantity.len() - 1], 1_000)
            } else if quantity.ends_with('M') {
                (&quantity[..quantity.len() - 1], 1_000_000)
            } else {
                (&quantity[..quantity.len() - 1], 1_000_000_000)
            };

            let total_bases: f64 = num_str.parse().with_context(|| {
                format!(
                    "Invalid quantity '{}'. Expected format like '250M', '25x', or a plain number",
                    quantity
                )
            })?;

            if total_bases <= 0.0 {
                bail!("Quantity must be positive, got: {}", total_bases);
            }

            let total_bases = (total_bases * multiplier as f64) as usize;
            let num_reads = (total_bases as f64 / mean_read_length as f64).ceil() as usize;

            info!(
                "Absolute quantity: {} bases → {} reads (read_length: {})",
                total_bases, num_reads, mean_read_length
            );

            Ok(num_reads)
        } else {
            let num_reads: usize = quantity.parse().with_context(|| {
                format!(
                    "Invalid quantity '{}'. Expected integer reads, coverage (25x), or bases (250M)",
                    quantity
                )
            })?;
            if num_reads == 0 {
                bail!("Quantity must be greater than zero");
            }
            info!("Quantity specified as reads: {}", num_reads);
            Ok(num_reads)
        }
    }
}

fn looks_like_motif_header(spec: &str) -> bool {
    let lower = spec.to_ascii_lowercase();
    let mut contains = 0usize;
    for token in lower.split_whitespace() {
        if token == "motif" || token == "motif_complement" {
            contains |= 0b001;
        } else if token == "mod_type" {
            contains |= 0b010;
        } else if token.starts_with("mod_position") {
            contains |= 0b100;
        }
    }
    contains == 0b111
}

pub struct LoadedModel {
    pub model: LearnedModel,
}

#[derive(Debug)]
pub enum ModelOrigin {
    File(PathBuf),
    Preset(ModelPreset),
}

impl ModelOrigin {
    fn describe(&self) -> String {
        match self {
            ModelOrigin::File(path) => path.display().to_string(),
            ModelOrigin::Preset(preset) => {
                format!("preset '{}' ({})", preset.cli_token(), preset.label())
            }
        }
    }
}

pub fn load_model_source(
    model_path: Option<&PathBuf>,
    preset: Option<ModelPreset>,
) -> Result<Option<LoadedModel>> {
    if let Some(path) = model_path {
        let origin = ModelOrigin::File(path.clone());
        info!("Loading model from {}", origin.describe());
        let model = LearnedModel::load(path)?;
        return Ok(Some(LoadedModel { model }));
    }

    if let Some(preset) = preset {
        let origin = ModelOrigin::Preset(preset);
        info!("Loading built-in model {}", origin.describe());
        let model = preset
            .load_model()
            .with_context(|| format!("Failed to load preset '{}'", preset.cli_token()))?;
        return Ok(Some(LoadedModel { model }));
    }

    Ok(None)
}

pub fn build_tagger_from_args(
    motifs: Vec<MotifDefinition>,
    args: &SimulateArgs,
    sampler_map: Option<std::collections::HashMap<String, MotifSampler>>,
    model: Option<&LearnedModel>,
) -> Result<crate::tagging::MethylationTagger> {
    if let Some(model) = model {
        let model_keys: HashSet<_> = model.motifs.iter().map(|m| m.key.clone()).collect();
        for motif in &motifs {
            let key = motif.model_key();
            if !model_keys.contains(&key) {
                bail!("Motif '{}' is not present in the supplied model file", key);
            }
        }
    }
    MethylationTaggerBuilder::new(motifs)
        .with_probs(args.motif_high_prob, args.non_motif_high_prob)
        .with_betas(
            args.high_beta_alpha,
            args.high_beta_beta,
            args.low_beta_alpha,
            args.low_beta_beta,
        )
        .with_seed(args.seed)
        .with_sampler_map(sampler_map)
        .build()
}

pub fn build_simulator_strategy(args: &SimulateArgs) -> Result<Option<SimulatorKindStrategy>> {
    let Some(reference_path) = args.reference.as_ref() else {
        return Ok(None);
    };

    let references = load_references(reference_path)?;
    let reference_size: usize = references.iter().map(|r| r.len()).sum();
    info!(
        "Loaded reference: {} sequences, {} total bases",
        references.len(),
        reference_size
    );

    let length_spec = ReadLengthSpec::Mode(args.read_length);
    let num_reads = parse_quantity(&args.quantity, reference_size, args.read_length)?;

    match args.simulator {
        SimulatorKind::Builtin => {
            let profile = ErrorProfile {
                substitution_rate: args.substitution_rate,
                insertion_rate: args.insertion_rate,
                deletion_rate: args.deletion_rate,
            };
            info!("Starting builtin simulator");
            debug!(
                "Error profile: sub={:.3}, ins={:.3}, del={:.3}",
                profile.substitution_rate, profile.insertion_rate, profile.deletion_rate
            );
            let config = BuiltinSimulationConfig {
                references,
                profile,
                length_spec,
                num_reads,
                name_prefix: args.name_prefix.clone(),
                seed: args.seed,
            };
            let strategy = BuiltinStrategy::new(config)?;
            Ok(Some(SimulatorKindStrategy::Builtin(strategy)))
        }
        SimulatorKind::Badreads => {
            info!("Starting badread simulation");
            let config = BadreadsSimulationConfig {
                reference: reference_path.clone(),
                num_reads,
                read_length: args.read_length,
                executable: args.badreads_exec.clone(),
                extra_args: args.badreads_extra.clone(),
                name_prefix: args.name_prefix.clone(),
                seed: args.seed,
            };
            Ok(Some(SimulatorKindStrategy::Badreads(
                BadreadsStrategy::new(config),
            )))
        }
    }
}
