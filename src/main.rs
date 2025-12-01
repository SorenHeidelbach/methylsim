mod fastx;
mod model;
mod motif;
mod presets;
mod simulator;
mod tagging;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use anyhow::{anyhow, bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use log::{debug, info, warn};

use fastx::{read_fastx, ReadRecord};
use model::{learn_model_from_reads, LearnedModel, MotifSampler};
use motif::{load_motif_file, MotifDefinition};
use presets::ModelPreset;
use simulator::{
    load_references, BadreadsSimulationConfig, BadreadsStrategy, BuiltinSimulationConfig,
    BuiltinStrategy, ErrorProfile, ReadLengthSpec, SimulatorKindStrategy, SimulatorStrategy,
};
use tagging::{MethylationTagger, MethylationTaggerBuilder};

#[derive(Copy, Clone, Debug, ValueEnum, PartialEq, Eq)]
enum SimulatorKind {
    Builtin,
    Badreads,
}

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Simulate/tag reads with MM/ML tags or fit methylation models",
    long_about = None
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Simulate reads or tag existing FASTQ/FASTA/BAM with methylation marks
    ///
    /// This command can:
    /// - Simulate new reads from a reference FASTA (--reference)
    /// - Tag existing FASTQ/FASTA files with methylation (--reads)
    /// - Tag BAM files using precomputed tags (--bam-input + --tags-tsv-input)
    Simulate(SimulateArgs),

    /// Learn a methylation model from reads with MM/ML tags
    ///
    /// Takes tagged FASTQ (with MM/ML in headers) and learns empirical distributions
    /// and Beta parameters for each motif. Output model can be used with simulate --model-in.
    FitModel(FitModelArgs),
}

#[derive(Parser, Debug)]
struct SimulateArgs {
    // ========== Input Mode Selection (pick one) ==========
    /// Tag existing FASTQ/FASTA file with methylation (untagged inputs allowed)
    #[arg(long, help_heading = "Input Mode (choose one)")]
    reads: Option<PathBuf>,

    /// Reference FASTA for simulating new reads (use with --num-reads)
    #[arg(long, help_heading = "Input Mode (choose one)")]
    reference: Option<PathBuf>,

    /// Input BAM file to tag directly with methylation (reads SEQ field, adds MM/ML tags)
    #[arg(long, help_heading = "Input Mode (choose one)")]
    bam_input: Option<PathBuf>,

    /// Output BAM path (required when using --bam-input)
    #[arg(long, help_heading = "BAM Tagging Mode")]
    bam_output: Option<PathBuf>,

    // ========== Required: Motif Definitions ==========
    /// Motif specification(s): 'motif[:canonical[:offset[:mod[:strand]]]]' or 'sequence_modtype_offset[_strand]'
    ///
    /// Examples: 'GATC:A:1:a:+' or 'CG_m_0' or 'CCWGG_5mC_1'
    /// Can specify multiple times or use comma-separated values.
    /// If --model-in is provided and no motifs specified, uses motifs from model file.
    #[arg(
        long = "motif",
        num_args = 1..,
        value_delimiter = ',',
        required_unless_present_any = ["motifs_file", "model_in", "model_preset"],
        help_heading = "Motif Definitions"
    )]
    motifs: Vec<String>,

    /// File with motif specifications (one per line, '#' for comments)
    #[arg(long, help_heading = "Motif Definitions")]
    motifs_file: Option<PathBuf>,

    // ========== Simulation Options (when using --reference) ==========
    /// Choose simulator: 'builtin' (fast) or 'badreads' (realistic errors, requires badread installed)
    #[arg(long, default_value_t = SimulatorKind::Builtin, value_enum, help_heading = "Simulation Options")]
    simulator: SimulatorKind,

    /// Amount of sequence to generate: absolute like '250M' (bases) or coverage like '25x'
    ///
    /// Examples: '100M' = 100 million bases, '50x' = 50x coverage of reference
    /// If not specified, uses --num-reads instead
    #[arg(
        long,
        help_heading = "Simulation Options",
        conflicts_with = "num_reads"
    )]
    quantity: Option<String>,

    /// Number of reads to generate (deprecated: use --quantity instead)
    #[arg(long, default_value_t = 100, help_heading = "Simulation Options")]
    num_reads: usize,

    /// Target read length in bp (mode of distribution)
    #[arg(long, default_value_t = 5000, help_heading = "Simulation Options")]
    read_length: usize,

    /// Alternative: specify mean read length (overrides --read-length)
    #[arg(long, help_heading = "Simulation Options")]
    read_length_mean: Option<usize>,

    /// Alternative: specify read length N50 (overrides --read-length)
    #[arg(long, help_heading = "Simulation Options")]
    read_length_n50: Option<usize>,

    /// Prefix for simulated read names (e.g., 'methylsim_000001')
    #[arg(long, default_value = "methylsim", help_heading = "Simulation Options")]
    name_prefix: String,

    /// Random seed for reproducibility
    #[arg(long, default_value_t = 1, help_heading = "Simulation Options")]
    seed: u64,

    /// Path to badreads executable (default: searches $PATH)
    #[arg(long, help_heading = "Badreads Simulator Options")]
    badreads_exec: Option<PathBuf>,

    /// Extra arguments passed to badreads (e.g., '--quantity 50x --error_model nanopore2023')
    #[arg(
        long,
        help_heading = "Badreads Simulator Options",
        allow_hyphen_values = true
    )]
    badreads_extra: Option<String>,

    /// Substitution error rate (0.0-1.0) for builtin simulator
    #[arg(
        long,
        default_value_t = 0.03,
        help_heading = "Builtin Simulator Error Model"
    )]
    substitution_rate: f64,

    /// Insertion error rate (0.0-1.0) for builtin simulator
    #[arg(
        long,
        default_value_t = 0.01,
        help_heading = "Builtin Simulator Error Model"
    )]
    insertion_rate: f64,

    /// Deletion error rate (0.0-1.0) for builtin simulator
    #[arg(
        long,
        default_value_t = 0.01,
        help_heading = "Builtin Simulator Error Model"
    )]
    deletion_rate: f64,

    // ========== Methylation Model Parameters ==========
    /// Load a pre-trained model JSON (from 'fit-model' command). Overrides simple probability parameters below
    #[arg(long, help_heading = "Methylation Model")]
    model_in: Option<PathBuf>,

    /// Use a bundled preset model (e.g., 'ecoli' for Dam/Dcm methylation)
    #[arg(
        long,
        value_enum,
        conflicts_with = "model_in",
        help_heading = "Methylation Model"
    )]
    model_preset: Option<ModelPreset>,

    /// Probability a motif site is methylated (0.0-1.0) [used when --model-in not provided]
    #[arg(
        long,
        default_value_t = 0.95,
        help_heading = "Methylation Model (simple mode)"
    )]
    motif_high_prob: f64,

    /// Background probability for non-motif sites (0.0-1.0)
    #[arg(
        long,
        default_value_t = 0.01,
        help_heading = "Methylation Model (simple mode)"
    )]
    non_motif_high_prob: f64,

    /// Alpha parameter for the high (methylated) Beta distribution when not using a learned model
    #[arg(
        long,
        default_value_t = 2.8685813093772063,
        help_heading = "Methylation Model (simple mode)"
    )]
    high_beta_alpha: f64,

    /// Beta parameter for the high (methylated) Beta distribution when not using a learned model
    #[arg(
        long,
        default_value_t = 0.055071333309938346,
        help_heading = "Methylation Model (simple mode)"
    )]
    high_beta_beta: f64,

    /// Alpha parameter for the low (unmethylated) Beta distribution when not using a learned model
    #[arg(
        long,
        default_value_t = 2.460887754108385,
        help_heading = "Methylation Model (simple mode)"
    )]
    low_beta_alpha: f64,

    /// Beta parameter for the low (unmethylated) Beta distribution when not using a learned model
    #[arg(
        long,
        default_value_t = 6.6833122937987115,
        help_heading = "Methylation Model (simple mode)"
    )]
    low_beta_beta: f64,

    // ========== Output Options ==========
    /// Output FASTQ file path
    #[arg(
        long,
        default_value = "methylsim.fastq",
        help_heading = "Output Options"
    )]
    output_fastq: PathBuf,

    /// Optional TSV output with read_id, MM, ML columns
    #[arg(long, help_heading = "Output Options")]
    tags_tsv: Option<PathBuf>,
}

#[derive(Parser, Debug)]
struct FitModelArgs {
    // ========== Required Inputs ==========
    /// FASTQ file with MM/ML tags in read headers (from basecaller or prior tagging)
    #[arg(long, help_heading = "Required")]
    reads: PathBuf,

    /// Output path for learned model JSON (use with 'simulate --model-in')
    #[arg(long, help_heading = "Required")]
    model_out: PathBuf,

    // ========== Motif Definitions (at least one required) ==========
    /// Motif specification(s): 'motif[:canonical[:offset[:mod[:strand]]]]' or 'sequence_modtype_offset[_strand]'
    ///
    /// Examples: 'GATC:A:1:a:+' or 'CG_m_0' or 'CCWGG_5mC_1'
    /// Must match the modification types present in input MM/ML tags
    #[arg(
        long = "motif",
        num_args = 1..,
        value_delimiter = ',',
        required_unless_present = "motifs_file",
        help_heading = "Motif Definitions (at least one required)"
    )]
    motifs: Vec<String>,

    /// File with motif specifications (one per line, '#' for comments)
    #[arg(long, help_heading = "Motif Definitions (at least one required)")]
    motifs_file: Option<PathBuf>,

    // ========== Model Fitting Parameters ==========
    /// Probability threshold (0.0-1.0) to classify sites as methylated vs unmethylated
    ///
    /// ML values are converted to probabilities (ML/255), then compared to this threshold
    /// to fit separate Beta distributions for methylated and unmethylated populations
    #[arg(long, default_value_t = 0.5, help_heading = "Model Fitting Parameters")]
    model_threshold: f64,

    /// Limit the number of reads used for fitting (processes the first N reads)
    #[arg(long, help_heading = "Model Fitting Parameters")]
    n_reads: Option<usize>,
}

#[derive(Default)]
struct WriteStats {
    total_reads: usize,
    tagged_reads: usize,
    reads_with_motif_hits: usize,
    total_motif_hits: usize,
    total_motif_high: usize,
}

struct LoadedModel {
    model: LearnedModel,
}

#[derive(Debug)]
enum ModelOrigin {
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

fn main() -> Result<()> {
    // Initialize logger - control verbosity with RUST_LOG environment variable
    // Examples: RUST_LOG=info, RUST_LOG=debug, RUST_LOG=methylsim=debug
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp(None)
        .format_module_path(false)
        .init();

    let cli = Cli::parse();
    match cli.command {
        Commands::Simulate(args) => run_simulate(args),
        Commands::FitModel(args) => run_fit_model(args),
    }
}

fn collect_motif_specs(
    motifs: &[String],
    motifs_file: Option<&PathBuf>,
    model: Option<&LearnedModel>,
) -> Result<Vec<String>> {
    debug!("Collecting motif specifications from command line, files, and model");
    let mut specs = motifs.to_vec();

    if let Some(path) = motifs_file {
        info!("Loading motifs from file: {}", path.display());
        let mut file_specs = load_motif_file(path)?;
        debug!("Loaded {} motif specs from file", file_specs.len());
        specs.append(&mut file_specs);
    }

    // If no motifs provided but model is available, extract from model
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

/// Parse quantity string (e.g., "250M" or "25x") and calculate number of reads
fn parse_quantity(quantity: &str, reference_size: usize, mean_read_length: usize) -> Result<usize> {
    let quantity = quantity.trim().to_uppercase();
    debug!("Parsing quantity specification: {}", quantity);

    // Check for coverage specification (ends with 'X')
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

        // Calculate number of reads needed for desired coverage
        // coverage = (num_reads * mean_read_length) / reference_size
        // num_reads = (coverage * reference_size) / mean_read_length
        let num_reads =
            ((coverage * reference_size as f64) / mean_read_length as f64).ceil() as usize;

        info!(
            "Coverage-based quantity: {}x coverage → {} reads (ref_size: {}, read_length: {})",
            coverage, num_reads, reference_size, mean_read_length
        );

        Ok(num_reads)
    } else {
        // Parse absolute base count with suffix (K, M, G, or just a number)
        let (num_str, multiplier) = if quantity.ends_with('K') {
            (&quantity[..quantity.len() - 1], 1_000)
        } else if quantity.ends_with('M') {
            (&quantity[..quantity.len() - 1], 1_000_000)
        } else if quantity.ends_with('G') {
            (&quantity[..quantity.len() - 1], 1_000_000_000)
        } else {
            (quantity.as_str(), 1)
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

fn load_model_source(
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

fn build_tagger_from_args(
    motifs: Vec<MotifDefinition>,
    args: &SimulateArgs,
    sampler_map: Option<HashMap<String, MotifSampler>>,
    model: Option<&LearnedModel>,
) -> Result<MethylationTagger> {
    if let Some(model) = model {
        // Ensure every motif requested exists in the model when using a learned model
        let model_keys: std::collections::HashSet<_> =
            model.motifs.iter().map(|m| m.key.clone()).collect();
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

fn build_simulator_strategy(args: &SimulateArgs) -> Result<Option<SimulatorKindStrategy>> {
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

    let length_spec = if let Some(mean) = args.read_length_mean {
        ReadLengthSpec::Mean(mean)
    } else if let Some(n50) = args.read_length_n50 {
        ReadLengthSpec::N50(n50)
    } else {
        ReadLengthSpec::Mode(args.read_length)
    };

    let num_reads = if let Some(quantity) = &args.quantity {
        let mean_length = args
            .read_length_mean
            .unwrap_or_else(|| args.read_length_n50.unwrap_or(args.read_length));
        parse_quantity(quantity, reference_size, mean_length)?
    } else {
        args.num_reads
    };

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
                seed: args.seed,
            };
            Ok(Some(SimulatorKindStrategy::Badreads(
                BadreadsStrategy::new(config),
            )))
        }
    }
}

fn run_simulate(args: SimulateArgs) -> Result<()> {
    validate_simulate_inputs(&args)?;
    let loaded_model = load_model_source(args.model_in.as_ref(), args.model_preset)?;
    let model_ref = loaded_model.as_ref().map(|m| &m.model);

    // Handle direct BAM tagging mode
    if args.reads.is_none() && args.reference.is_none() && args.bam_input.is_some() {
        let bam_in = args.bam_input.as_ref().unwrap();
        let bam_out = args
            .bam_output
            .as_ref()
            .expect("--bam-output must be supplied when using --bam-input");

        let motif_specs = collect_motif_specs(&args.motifs, args.motifs_file.as_ref(), model_ref)?;
        let motifs = motif_specs
            .iter()
            .map(|spec| MotifDefinition::parse(spec))
            .collect::<Result<Vec<_>>>()?;
        let sampler_map = model_ref.map(|model| model.build_samplers()).transpose()?;

        let mut tagger = build_tagger_from_args(motifs, &args, sampler_map, model_ref)?;

        tag_bam_directly(bam_in, bam_out, &mut tagger)?;
        return Ok(());
    }

    let motif_specs = collect_motif_specs(&args.motifs, args.motifs_file.as_ref(), model_ref)?;
    let motifs = motif_specs
        .iter()
        .map(|spec| MotifDefinition::parse(spec))
        .collect::<Result<Vec<_>>>()?;
    let sampler_map = model_ref.map(|model| model.build_samplers()).transpose()?;

    let mut tagger = build_tagger_from_args(motifs, &args, sampler_map, model_ref)?;

    let mut fastq_writer = BufWriter::new(
        File::create(&args.output_fastq)
            .with_context(|| format!("Failed to create '{}'", args.output_fastq.display()))?,
    );
    let mut tags_writer = if let Some(path) = &args.tags_tsv {
        let mut writer = BufWriter::new(
            File::create(path).with_context(|| format!("Failed to create '{}'", path.display()))?,
        );
        writer.write_all(b"read_id\tMM\tML\n")?;
        Some(writer)
    } else {
        None
    };
    let mut stats = WriteStats::default();

    if let Some(reads_path) = &args.reads {
        eprintln!(
            "Annotating existing reads from '{}'...",
            reads_path.display()
        );
        let reads = read_fastx(reads_path)?;
        for read in reads {
            write_read(
                &read,
                &mut tagger,
                &mut fastq_writer,
                &mut tags_writer,
                &mut stats,
            )?;
        }
    } else {
        let mut strategy = build_simulator_strategy(&args)?
            .ok_or_else(|| anyhow!("Simulation requires a reference sequence"))?;
        let reads = strategy.simulate()?;
        for read in reads {
            write_read(
                &read,
                &mut tagger,
                &mut fastq_writer,
                &mut tags_writer,
                &mut stats,
            )?;
        }
    }

    fastq_writer.flush()?;
    if let Some(writer) = tags_writer.as_mut() {
        writer.flush()?;
    }

    println!(
        "Wrote {} reads to {} ({} reads carried MM/ML tags)",
        stats.total_reads,
        args.output_fastq.display(),
        stats.tagged_reads
    );
    if stats.total_reads > 0 {
        let avg_motif_hits = stats.total_motif_hits as f64 / stats.total_reads as f64;
        let avg_high_hits = stats.total_motif_high as f64 / stats.total_reads as f64;
        println!(
            "Motif summary: {} reads contained motif hits | avg motifs/read = {:.2} | avg high-methylated motifs/read = {:.2}",
            stats.reads_with_motif_hits, avg_motif_hits, avg_high_hits
        );
    }
    Ok(())
}

fn run_fit_model(args: FitModelArgs) -> Result<()> {
    info!("Starting model fitting");
    validate_fit_inputs(&args)?;

    let motif_specs = collect_motif_specs(&args.motifs, args.motifs_file.as_ref(), None)?;
    debug!("Parsing {} motif specifications", motif_specs.len());
    let motifs = motif_specs
        .iter()
        .map(|spec| MotifDefinition::parse(spec))
        .collect::<Result<Vec<_>>>()?;

    eprintln!(
        "Learning methylation patterns from reads in: {}",
        args.reads.display()
    );
    debug!("Model threshold: {}", args.model_threshold);
    debug!(
        "Read sampling limit: {}",
        args.n_reads
            .map(|n| n.to_string())
            .unwrap_or_else(|| "all".to_string())
    );
    let model = learn_model_from_reads(&args.reads, &motifs, args.model_threshold, args.n_reads)?;

    info!("Saving model to: {}", args.model_out.display());
    model.save(&args.model_out)?;
    info!(
        "Successfully saved model with {} motif(s)",
        model.motifs.len()
    );

    Ok(())
}

fn validate_simulate_inputs(args: &SimulateArgs) -> Result<()> {
    if args.model_in.is_some() && args.model_preset.is_some() {
        bail!("Only one of --model-in or --model-preset can be supplied");
    }
    if args.reads.is_none() && args.reference.is_none() && args.bam_input.is_none() {
        bail!("Simulate requires either --reads, --reference, or --bam-input");
    }
    if args.bam_input.is_some() && args.bam_output.is_none() {
        bail!("--bam-input requires --bam-output");
    }
    let doing_reads = args.reads.is_some() || args.reference.is_some();
    if doing_reads {
        if args.reads.is_none() && args.reference.is_none() {
            bail!("Either --reads or --reference must be supplied");
        }
        if args.simulator == SimulatorKind::Badreads && args.reads.is_some() {
            eprintln!("Warning: --reads provided; ignoring --simulator badreads");
        }
        if args.reads.is_none()
            && args.simulator == SimulatorKind::Badreads
            && args.reference.is_none()
        {
            bail!("--simulator badreads requires --reference");
        }
        if args.read_length == 0 {
            bail!("--read-length must be greater than zero");
        }
        if let Some(mean) = args.read_length_mean {
            if mean == 0 {
                bail!("--read-length-mean must be greater than zero");
            }
        }
        if let Some(n50) = args.read_length_n50 {
            if n50 == 0 {
                bail!("--read-length-n50 must be greater than zero");
            }
        }
        if args.read_length_mean.is_some() && args.read_length_n50.is_some() {
            bail!("Only one of --read-length-mean or --read-length-n50 can be supplied");
        }
        if !(0.0..=1.0).contains(&args.motif_high_prob) {
            bail!("--motif-high-prob must be between 0 and 1");
        }
        if !(0.0..=1.0).contains(&args.non_motif_high_prob) {
            bail!("--non-motif-high-prob must be between 0 and 1");
        }
        for (name, val) in [
            ("--high-beta-alpha", args.high_beta_alpha),
            ("--high-beta-beta", args.high_beta_beta),
            ("--low-beta-alpha", args.low_beta_alpha),
            ("--low-beta-beta", args.low_beta_beta),
        ] {
            if !val.is_finite() || val <= 0.0 {
                bail!("{name} must be positive and finite");
            }
        }
    }
    Ok(())
}

fn validate_fit_inputs(args: &FitModelArgs) -> Result<()> {
    if !(0.0..=1.0).contains(&args.model_threshold) {
        bail!("--model-threshold must be between 0 and 1");
    }
    if args.motifs.is_empty() && args.motifs_file.is_none() {
        bail!("At least one motif must be supplied via --motif or --motifs-file");
    }
    if let Some(n) = args.n_reads {
        if n == 0 {
            bail!("--n-reads must be greater than zero");
        }
    }
    Ok(())
}

fn write_read(
    read: &ReadRecord,
    tagger: &mut MethylationTagger,
    fastq_writer: &mut BufWriter<File>,
    tags_writer: &mut Option<BufWriter<File>>,
    stats: &mut WriteStats,
) -> Result<()> {
    let tags = tagger.annotate(read.sequence_str());
    let mut header = format!("@{}", read.name);
    if let Some(comment) = &read.comment {
        header.push_str(comment);
    }
    if tags.mm_tag.is_some() || tags.ml_tag.is_some() {
        stats.tagged_reads += 1;
        if let Some(mm) = &tags.mm_tag {
            header.push('\t');
            header.push_str(mm);
        }
        if let Some(ml) = &tags.ml_tag {
            header.push('\t');
            header.push_str(ml);
        }
    }
    writeln!(fastq_writer, "{header}")?;
    writeln!(fastq_writer, "{}", read.sequence_str())?;
    writeln!(fastq_writer, "+")?;
    writeln!(fastq_writer, "{}", read.quality_str())?;

    if let Some(writer) = tags_writer.as_mut() {
        writeln!(
            writer,
            "{}\t{}\t{}",
            read.name,
            tags.mm_tag.unwrap_or_default(),
            tags.ml_tag.unwrap_or_default()
        )?;
    }
    stats.total_reads += 1;
    stats.total_motif_hits += tags.motif_hit_count;
    stats.total_motif_high += tags.motif_high_count;
    if tags.motif_hit_count > 0 {
        stats.reads_with_motif_hits += 1;
    }
    Ok(())
}

/// Tag BAM file directly by reading SEQ field and adding MM/ML tags
fn tag_bam_directly(bam_in: &Path, bam_out: &Path, tagger: &mut MethylationTagger) -> Result<()> {
    info!("Starting direct BAM tagging");
    debug!("Input BAM: {}", bam_in.display());
    debug!("Output BAM: {}", bam_out.display());

    let mut view_proc = Command::new("samtools")
        .arg("view")
        .arg("-h")
        .arg(bam_in)
        .stdout(Stdio::piped())
        .spawn()
        .with_context(|| "Failed to spawn 'samtools view -h'")?;

    let mut write_proc = Command::new("samtools")
        .arg("view")
        .arg("-b")
        .arg("-o")
        .arg(bam_out)
        .arg("-")
        .stdin(Stdio::piped())
        .spawn()
        .with_context(|| "Failed to spawn 'samtools view -b'")?;

    let reader_stdout = view_proc
        .stdout
        .take()
        .ok_or_else(|| anyhow!("Failed to capture samtools stdout"))?;
    let writer_stdin = write_proc
        .stdin
        .take()
        .ok_or_else(|| anyhow!("Failed to capture samtools stdin"))?;

    let mut reader = BufReader::new(reader_stdout);
    let mut writer = BufWriter::new(writer_stdin);
    let mut line = String::new();
    let mut tagged = 0usize;
    let mut total = 0usize;
    let mut headers = 0usize;

    info!("Processing BAM records...");
    loop {
        line.clear();
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }

        // Pass through header lines
        if line.starts_with('@') {
            headers += 1;
            writer.write_all(line.as_bytes())?;
            continue;
        }

        // Log progress periodically
        if total > 0 && total % 10000 == 0 {
            info!(
                "Processed {} reads, tagged {} ({:.1}%)",
                total,
                tagged,
                (tagged as f64 / total as f64) * 100.0
            );
        }

        // Trim newlines
        while line.ends_with('\n') || line.ends_with('\r') {
            line.pop();
        }

        if line.is_empty() {
            writer.write_all(b"\n")?;
            continue;
        }

        let mut fields: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();

        // BAM format: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS...]
        if fields.len() >= 11 {
            let seq = &fields[9];

            // Tag the sequence
            let tags = tagger.annotate(seq);

            // Remove existing MM/ML tags
            fields.retain(|field| !(field.starts_with("MM:Z:") || field.starts_with("ML:B:")));

            // Add new tags if present
            if tags.mm_tag.is_some() || tags.ml_tag.is_some() {
                if let Some(mm) = &tags.mm_tag {
                    fields.push(mm.clone());
                }
                if let Some(ml) = &tags.ml_tag {
                    fields.push(ml.clone());
                }
                tagged += 1;
            }
        }

        total += 1;
        let mut rebuilt = fields.join("\t");
        rebuilt.push('\n');
        writer.write_all(rebuilt.as_bytes())?;
    }

    writer.flush()?;
    drop(writer);

    let write_status = write_proc
        .wait()
        .with_context(|| "samtools view -b terminated unexpectedly")?;
    if !write_status.success() {
        bail!("samtools view -b exited with {}", write_status);
    }

    let view_status = view_proc
        .wait()
        .with_context(|| "samtools view -h terminated unexpectedly")?;
    if !view_status.success() {
        bail!("samtools view -h exited with {}", view_status);
    }

    debug!("Processed {} header lines", headers);
    info!(
        "BAM tagging complete: tagged {} of {} reads ({:.1}%) → {}",
        tagged,
        total,
        if total > 0 {
            (tagged as f64 / total as f64) * 100.0
        } else {
            0.0
        },
        bam_out.display()
    );

    Ok(())
}

#[cfg(test)]
mod cli_tests {
    use super::*;
    use clap::Parser;

    #[test]
    fn badreads_extra_accepts_leading_hyphenated_values() {
        let cli = Cli::try_parse_from([
            "methylsim",
            "simulate",
            "--motif",
            "CG",
            "--reference",
            "ref.fa",
            "--simulator",
            "badreads",
            "--badreads-extra",
            "--quantity 50x",
        ])
        .expect("Cli parsing failed");

        match cli.command {
            Commands::Simulate(args) => {
                assert_eq!(args.badreads_extra.as_deref(), Some("--quantity 50x"));
            }
            _ => panic!("Expected Simulate command"),
        }
    }
}
