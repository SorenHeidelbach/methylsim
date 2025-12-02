mod fastx;
mod model;
mod motif;
mod presets;
mod simulator;
mod tagging;

use anyhow::{anyhow, bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use log::{debug, info, warn};
use noodles::bam::io::{
    reader::Builder as BamReaderBuilder, writer::Builder as BamWriterBuilder, Writer as BamWriter,
};
use noodles::bgzf;
use noodles::sam::{
    self,
    alignment::io::Write as AlignmentWrite,
    alignment::{self, record::data::field::Tag, record::QualityScores, record_buf},
    io::reader::Builder as SamReaderBuilder,
    io::Writer as SamWriter,
};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;

use fastx::{read_fastx, ReadRecord};
use model::{learn_model_from_reads, LearnedModel, MotifSampler};
use motif::{load_motif_file, MotifDefinition};
use presets::ModelPreset;
use simulator::{
    load_references, BadreadsSimulationConfig, BadreadsStrategy, BuiltinSimulationConfig,
    BuiltinStrategy, ErrorProfile, ReadLengthSpec, SimulatorKindStrategy, SimulatorStrategy,
};
use tagging::{MethylationTagger, MethylationTaggerBuilder, TagResult};

#[derive(Copy, Clone, Debug, ValueEnum, PartialEq, Eq)]
enum SimulatorKind {
    Builtin,
    Badreads,
}

#[derive(Copy, Clone, Debug, ValueEnum, PartialEq, Eq)]
enum OutputFormat {
    Fastq,
    Sam,
    Bam,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum InputFormat {
    Fastx,
    Sam,
    Bam,
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
    /// Simulate reads or tag existing reads with methylation marks
    ///
    /// This command can:
    /// - Simulate new reads from a reference FASTA (--reference)
    /// - Tag existing FASTQ/SAM/BAM files with methylation (--reads)
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
    /// Tag existing FASTQ/SAM/BAM file with methylation (untagged inputs allowed)
    #[arg(long, help_heading = "Input Mode (choose one)")]
    reads: Option<PathBuf>,

    /// Reference FASTA for simulating new reads (use with --num-reads)
    #[arg(long, help_heading = "Input Mode (choose one)")]
    reference: Option<PathBuf>,

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

    /// Number of reads to generate, or sequence amount: absolute like '250M' (bases) or coverage like '25x'
    ///
    /// Examples: '100M' = 100 million bases, '50x' = 50x coverage of reference, '10000' = 10k reads
    #[arg(long, default_value = "100", help_heading = "Simulation Options")]
    quantity: String,

    /// Target read length in bp (mode of distribution)
    #[arg(long, default_value_t = 5000, help_heading = "Simulation Options")]
    read_length: usize,

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
    /// Output file path
    #[arg(
        short = 'o',
        long = "out",
        default_value = "methylsim.fastq",
        help_heading = "Output Options"
    )]
    out: PathBuf,

    /// Output format (fastq, sam, bam)
    #[arg(
        long = "out-format",
        default_value_t = OutputFormat::Fastq,
        value_enum,
        help_heading = "Output Options"
    )]
    out_format: OutputFormat,

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

impl OutputFormat {
    fn label(&self) -> &'static str {
        match self {
            OutputFormat::Fastq => "FASTQ",
            OutputFormat::Sam => "SAM",
            OutputFormat::Bam => "BAM",
        }
    }

    fn extensions(&self) -> &'static [&'static str] {
        match self {
            OutputFormat::Fastq => &["fastq", "fq"],
            OutputFormat::Sam => &["sam"],
            OutputFormat::Bam => &["bam"],
        }
    }
}

impl InputFormat {
    fn from_path(path: &Path) -> Result<Self> {
        if let Some(ext) = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|s| s.to_ascii_lowercase())
        {
            match ext.as_str() {
                "bam" => return Ok(InputFormat::Bam),
                "sam" => return Ok(InputFormat::Sam),
                _ => {}
            }
        }
        Ok(InputFormat::Fastx)
    }
}

struct OutputWriter {
    target: OutputTarget,
}

enum OutputTarget {
    Fastq(BufWriter<File>),
    Sam(SamOutput),
    Bam(BamOutput),
}

impl OutputWriter {
    fn new(path: &Path, format: OutputFormat) -> Result<Self> {
        let target = match format {
            OutputFormat::Fastq => {
                let writer = BufWriter::new(
                    File::create(path)
                        .with_context(|| format!("Failed to create '{}'", path.display()))?,
                );
                OutputTarget::Fastq(writer)
            }
            OutputFormat::Sam => OutputTarget::Sam(SamOutput::new(path)?),
            OutputFormat::Bam => OutputTarget::Bam(BamOutput::new(path)?),
        };
        Ok(Self { target })
    }

    fn write(&mut self, read: &ReadRecord, tags: &TagResult) -> Result<()> {
        match &mut self.target {
            OutputTarget::Fastq(writer) => write_fastq_record(writer, read, tags),
            OutputTarget::Sam(writer) => writer.write_record(read, tags),
            OutputTarget::Bam(writer) => writer.write_record(read, tags),
        }
    }

    fn finish(self) -> Result<()> {
        match self.target {
            OutputTarget::Fastq(mut writer) => {
                writer.flush()?;
            }
            OutputTarget::Sam(mut writer) => {
                writer.finish()?;
            }
            OutputTarget::Bam(mut writer) => {
                writer.finish()?;
            }
        }
        Ok(())
    }
}

struct SamOutput {
    writer: SamWriter<BufWriter<File>>,
    header: sam::Header,
}

impl SamOutput {
    fn new(path: &Path) -> Result<Self> {
        let header = default_output_header()?;
        let writer = BufWriter::new(
            File::create(path).with_context(|| format!("Failed to create '{}'", path.display()))?,
        );
        let mut writer = SamWriter::new(writer);
        writer.write_header(&header)?;
        Ok(Self { writer, header })
    }

    fn write_record(&mut self, read: &ReadRecord, tags: &TagResult) -> Result<()> {
        let record = build_alignment_record(read, tags)?;
        self.writer
            .write_alignment_record(&self.header, &record)
            .with_context(|| "Failed to write SAM record")?;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        self.writer.finish(&self.header)?;
        Ok(())
    }
}

struct BamOutput {
    writer: BamWriter<bgzf::Writer<File>>,
    header: sam::Header,
}

impl BamOutput {
    fn new(path: &Path) -> Result<Self> {
        let header = default_output_header()?;
        let mut writer = BamWriterBuilder::default()
            .build_from_path(path)
            .with_context(|| format!("Failed to create '{}'", path.display()))?;
        writer.write_header(&header)?;
        Ok(Self { writer, header })
    }

    fn write_record(&mut self, read: &ReadRecord, tags: &TagResult) -> Result<()> {
        let record = build_alignment_record(read, tags)?;
        self.writer
            .write_alignment_record(&self.header, &record)
            .with_context(|| "Failed to write BAM record")?;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        self.writer.finish(&self.header)?;
        Ok(())
    }
}

fn default_output_header() -> Result<sam::Header> {
    const HEADER: &str = "@HD\tVN:1.6\tSO:unknown\n@PG\tID:methylsim\tPN:methylsim\n";
    sam::Header::from_str(HEADER).context("Failed to build SAM/BAM header")
}

fn ensure_output_extension(path: &Path, format: OutputFormat) -> Result<()> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase())
        .ok_or_else(|| anyhow!("Output file must include an extension"))?;

    if !format.extensions().iter().any(|allowed| *allowed == ext) {
        bail!(
            "Output path '{}' does not match output format {} (expected extension: {})",
            path.display(),
            format.label(),
            format.extensions().join(", ")
        );
    }
    Ok(())
}

fn build_alignment_record(read: &ReadRecord, tags: &TagResult) -> Result<alignment::RecordBuf> {
    let name = record_buf::Name::from(read.name.as_bytes());
    let sequence = record_buf::Sequence::from(read.sequence.clone());
    let quality_scores = record_buf::QualityScores::from(read.quality.clone());
    let data = build_tag_data(tags)?;

    let record = record_buf::RecordBuf::builder()
        .set_name(name)
        .set_sequence(sequence)
        .set_quality_scores(quality_scores)
        .set_data(data)
        .build();

    Ok(record)
}

fn build_tag_data(tags: &TagResult) -> Result<record_buf::Data> {
    use noodles::sam::alignment::record_buf::data::field::{self, Value};

    let mut fields: Vec<(Tag, Value)> = Vec::new();

    if let Some(mm) = &tags.mm_tag {
        let value = mm
            .strip_prefix("MM:Z:")
            .or_else(|| mm.strip_prefix("MM:"))
            .unwrap_or(mm.as_str())
            .to_string();
        fields.push((Tag::new(b'M', b'M'), Value::String(value.into())));
    }

    if let Some(ml) = &tags.ml_tag {
        let values = parse_ml_values(ml)?;
        if let Some(values) = values {
            fields.push((
                Tag::new(b'M', b'L'),
                Value::Array(field::value::Array::UInt8(values)),
            ));
        }
    }

    Ok(fields.into_iter().collect())
}

fn parse_ml_values(raw: &str) -> Result<Option<Vec<u8>>> {
    let trimmed = raw
        .strip_prefix("ML:B:")
        .map(|s| s.to_string())
        .unwrap_or_else(|| raw.to_string());

    let mut parts = trimmed.splitn(2, ',');
    let type_code = parts.next().unwrap_or("C");
    let list = parts.next().unwrap_or("");

    if type_code != "C" {
        warn!(
            "Unsupported ML type '{}' in '{}'; skipping ML tag for this read",
            type_code, raw
        );
        return Ok(None);
    }
    if list.is_empty() {
        return Ok(Some(Vec::new()));
    }

    let mut values = Vec::new();
    for value in list.split(',').filter(|s| !s.is_empty()) {
        let parsed = value.parse::<u8>().with_context(|| {
            format!(
                "Invalid ML value '{}' in '{}'; expected unsigned byte",
                value, raw
            )
        })?;
        values.push(parsed);
    }
    Ok(Some(values))
}

struct LoadedModel {
    model: LearnedModel,
}

fn read_alignment_record_buf(record: &record_buf::RecordBuf, idx: usize) -> ReadRecord {
    let name = record
        .name()
        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
        .unwrap_or_else(|| format!("read_{idx:06}"));

    let mut sequence = record.sequence().as_ref().to_vec();
    for base in sequence.iter_mut() {
        base.make_ascii_uppercase();
    }

    let mut quality: Vec<u8> = record
        .quality_scores()
        .iter()
        .map(|q| q.saturating_add(33))
        .collect();
    if quality.is_empty() {
        quality = vec![b'I'; sequence.len()];
    }

    ReadRecord::new(name, sequence, quality)
}

fn read_sam_reads(path: &Path) -> Result<Vec<ReadRecord>> {
    let mut reader = SamReaderBuilder::default()
        .build_from_path(path)
        .with_context(|| format!("Failed to open SAM file '{}'", path.display()))?;
    let header = reader
        .read_header()
        .with_context(|| format!("Failed to read SAM header from '{}'", path.display()))?;
    let mut reads = Vec::new();
    for (idx, result) in reader.record_bufs(&header).enumerate() {
        let record = result.with_context(|| format!("Failed to read SAM record {}", idx + 1))?;
        reads.push(read_alignment_record_buf(&record, idx));
    }
    Ok(reads)
}

fn read_bam_reads(path: &Path) -> Result<Vec<ReadRecord>> {
    let mut reader = BamReaderBuilder::default()
        .build_from_path(path)
        .with_context(|| format!("Failed to open BAM file '{}'", path.display()))?;
    let header = reader
        .read_header()
        .with_context(|| format!("Failed to read BAM header from '{}'", path.display()))?;
    let mut reads = Vec::new();
    for (idx, result) in reader.record_bufs(&header).enumerate() {
        let record = result.with_context(|| format!("Failed to read BAM record {}", idx + 1))?;
        reads.push(read_alignment_record_buf(&record, idx));
    }
    Ok(reads)
}

fn read_input_reads(path: &Path) -> Result<Vec<ReadRecord>> {
    match InputFormat::from_path(path)? {
        InputFormat::Fastx => read_fastx(path),
        InputFormat::Sam => read_sam_reads(path),
        InputFormat::Bam => read_bam_reads(path),
    }
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
        // Parse absolute base count with suffix (K, M, G) or plain number = reads
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

fn run_simulate(args: SimulateArgs) -> Result<()> {
    validate_simulate_inputs(&args)?;
    ensure_output_extension(&args.out, args.out_format)?;
    let loaded_model = load_model_source(args.model_in.as_ref(), args.model_preset)?;
    let model_ref = loaded_model.as_ref().map(|m| &m.model);

    let motif_specs = collect_motif_specs(&args.motifs, args.motifs_file.as_ref(), model_ref)?;
    let motifs = motif_specs
        .iter()
        .map(|spec| MotifDefinition::parse(spec))
        .collect::<Result<Vec<_>>>()?;
    let sampler_map = model_ref.map(|model| model.build_samplers()).transpose()?;

    let mut tagger = build_tagger_from_args(motifs, &args, sampler_map, model_ref)?;

    let mut output_writer = OutputWriter::new(&args.out, args.out_format)?;
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
        let reads = read_input_reads(reads_path)?;
        for read in reads {
            write_read(
                &read,
                &mut tagger,
                &mut output_writer,
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
                &mut output_writer,
                &mut tags_writer,
                &mut stats,
            )?;
        }
    }

    output_writer.finish()?;
    if let Some(writer) = tags_writer.as_mut() {
        writer.flush()?;
    }

    println!(
        "Wrote {} reads to {} as {} ({} reads carried MM/ML tags)",
        stats.total_reads,
        args.out.display(),
        args.out_format.label(),
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
    if args.reads.is_none() && args.reference.is_none() {
        bail!("Simulate requires either --reads or --reference");
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
    output_writer: &mut OutputWriter,
    tags_writer: &mut Option<BufWriter<File>>,
    stats: &mut WriteStats,
) -> Result<()> {
    let tags = tagger.annotate(read.sequence_str());
    if tags.mm_tag.is_some() || tags.ml_tag.is_some() {
        stats.tagged_reads += 1;
    }
    output_writer.write(read, &tags)?;

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

fn write_fastq_record(
    fastq_writer: &mut BufWriter<File>,
    read: &ReadRecord,
    tags: &TagResult,
) -> Result<()> {
    let mut header = format!("@{}", read.name);
    if let Some(comment) = &read.comment {
        header.push_str(comment);
    }
    if tags.mm_tag.is_some() || tags.ml_tag.is_some() {
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

    #[test]
    fn simulate_defaults_output_flags() {
        let cli = Cli::try_parse_from([
            "methylsim",
            "simulate",
            "--motif",
            "CG",
            "--reads",
            "reads.fastq",
        ])
        .expect("Cli parsing failed");

        match cli.command {
            Commands::Simulate(args) => {
                assert_eq!(args.out, PathBuf::from("methylsim.fastq"));
                assert_eq!(args.out_format, OutputFormat::Fastq);
            }
            _ => panic!("Expected Simulate command"),
        }
    }

    #[test]
    fn simulate_accepts_bam_output_format() {
        let cli = Cli::try_parse_from([
            "methylsim",
            "simulate",
            "--motif",
            "CG",
            "--reads",
            "reads.fastq",
            "--out",
            "custom.bam",
            "--out-format",
            "bam",
        ])
        .expect("Cli parsing failed");

        match cli.command {
            Commands::Simulate(args) => {
                assert_eq!(args.out, PathBuf::from("custom.bam"));
                assert_eq!(args.out_format, OutputFormat::Bam);
            }
            _ => panic!("Expected Simulate command"),
        }
    }
}
