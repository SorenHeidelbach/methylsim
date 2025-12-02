use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

use crate::presets::ModelPreset;

#[derive(Copy, Clone, Debug, ValueEnum, PartialEq, Eq)]
pub enum SimulatorKind {
    Builtin,
    Badreads,
}

#[derive(Copy, Clone, Debug, ValueEnum, PartialEq, Eq)]
pub enum OutputFormat {
    Fastq,
    Sam,
    Bam,
}

impl OutputFormat {
    pub fn label(&self) -> &'static str {
        match self {
            OutputFormat::Fastq => "FASTQ",
            OutputFormat::Sam => "SAM",
            OutputFormat::Bam => "BAM",
        }
    }

    pub fn extensions(&self) -> &'static [&'static str] {
        match self {
            OutputFormat::Fastq => &["fastq", "fq"],
            OutputFormat::Sam => &["sam"],
            OutputFormat::Bam => &["bam"],
        }
    }
}

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Simulate/tag reads with MM/ML tags or fit methylation models",
    long_about = None
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
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
pub struct SimulateArgs {
    // ========== Input Mode Selection (pick one) ==========
    /// Tag existing FASTQ/SAM/BAM file with methylation (untagged inputs allowed)
    #[arg(long, help_heading = "Input Mode (choose one)")]
    pub reads: Option<PathBuf>,

    /// Reference FASTA for simulating new reads (use with --num-reads)
    #[arg(long, help_heading = "Input Mode (choose one)")]
    pub reference: Option<PathBuf>,

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
    pub motifs: Vec<String>,

    /// File with motif specifications (one per line, '#' for comments)
    #[arg(long, help_heading = "Motif Definitions")]
    pub motifs_file: Option<PathBuf>,

    // ========== Simulation Options (when using --reference) ==========
    /// Choose simulator: 'builtin' (fast) or 'badreads' (realistic errors, requires badread installed)
    #[arg(long, default_value_t = SimulatorKind::Builtin, value_enum, help_heading = "Simulation Options")]
    pub simulator: SimulatorKind,

    /// Number of reads to generate, or sequence amount: absolute like '250M' (bases) or coverage like '25x'
    ///
    /// Examples: '100M' = 100 million bases, '50x' = 50x coverage of reference, '10000' = 10k reads
    #[arg(long, default_value = "100", help_heading = "Simulation Options")]
    pub quantity: String,

    /// Target read length in bp (mode of distribution)
    #[arg(long, default_value_t = 5000, help_heading = "Simulation Options")]
    pub read_length: usize,

    /// Prefix for simulated read names (e.g., 'methylsim_000001')
    #[arg(long, default_value = "methylsim", help_heading = "Simulation Options")]
    pub name_prefix: String,

    /// Random seed for reproducibility
    #[arg(long, default_value_t = 1, help_heading = "Simulation Options")]
    pub seed: u64,

    /// Path to badreads executable (default: searches $PATH)
    #[arg(long, help_heading = "Badreads Simulator Options")]
    pub badreads_exec: Option<PathBuf>,

    /// Extra arguments passed to badreads (e.g., '--quantity 50x --error_model nanopore2023')
    #[arg(
        long,
        help_heading = "Badreads Simulator Options",
        allow_hyphen_values = true
    )]
    pub badreads_extra: Option<String>,

    /// Substitution error rate (0.0-1.0) for builtin simulator
    #[arg(
        long,
        default_value_t = 0.03,
        help_heading = "Builtin Simulator Error Model"
    )]
    pub substitution_rate: f64,

    /// Insertion error rate (0.0-1.0) for builtin simulator
    #[arg(
        long,
        default_value_t = 0.01,
        help_heading = "Builtin Simulator Error Model"
    )]
    pub insertion_rate: f64,

    /// Deletion error rate (0.0-1.0) for builtin simulator
    #[arg(
        long,
        default_value_t = 0.01,
        help_heading = "Builtin Simulator Error Model"
    )]
    pub deletion_rate: f64,

    // ========== Methylation Model Parameters ==========
    /// Load a pre-trained model JSON (from 'fit-model' command). Overrides simple probability parameters below
    #[arg(long, help_heading = "Methylation Model")]
    pub model_in: Option<PathBuf>,

    /// Use a bundled preset model (e.g., 'ecoli' for Dam/Dcm methylation)
    #[arg(
        long,
        value_enum,
        conflicts_with = "model_in",
        help_heading = "Methylation Model"
    )]
    pub model_preset: Option<ModelPreset>,

    /// Probability a motif site is methylated (0.0-1.0) [used when --model-in not provided]
    #[arg(
        long,
        default_value_t = 0.95,
        help_heading = "Methylation Model (simple mode)"
    )]
    pub motif_high_prob: f64,

    /// Background probability for non-motif sites (0.0-1.0)
    #[arg(
        long,
        default_value_t = 0.01,
        help_heading = "Methylation Model (simple mode)"
    )]
    pub non_motif_high_prob: f64,

    /// Alpha parameter for the high (methylated) Beta distribution when not using a learned model
    #[arg(
        long,
        default_value_t = 2.8685813093772063,
        help_heading = "Methylation Model (simple mode)"
    )]
    pub high_beta_alpha: f64,

    /// Beta parameter for the high (methylated) Beta distribution when not using a learned model
    #[arg(
        long,
        default_value_t = 0.055071333309938346,
        help_heading = "Methylation Model (simple mode)"
    )]
    pub high_beta_beta: f64,

    /// Alpha parameter for the low (unmethylated) Beta distribution when not using a learned model
    #[arg(
        long,
        default_value_t = 2.460887754108385,
        help_heading = "Methylation Model (simple mode)"
    )]
    pub low_beta_alpha: f64,

    /// Beta parameter for the low (unmethylated) Beta distribution when not using a learned model
    #[arg(
        long,
        default_value_t = 6.6833122937987115,
        help_heading = "Methylation Model (simple mode)"
    )]
    pub low_beta_beta: f64,

    // ========== Output Options ==========
    /// Output file path
    #[arg(
        short = 'o',
        long = "out",
        default_value = "methylsim.fastq",
        help_heading = "Output Options"
    )]
    pub out: PathBuf,

    /// Output format (fastq, sam, bam)
    #[arg(
        long = "out-format",
        default_value_t = OutputFormat::Fastq,
        value_enum,
        help_heading = "Output Options"
    )]
    pub out_format: OutputFormat,

    /// Optional TSV output with read_id, MM, ML columns
    #[arg(long, help_heading = "Output Options")]
    pub tags_tsv: Option<PathBuf>,
}

#[derive(Parser, Debug)]
pub struct FitModelArgs {
    // ========== Required Inputs ==========
    /// FASTQ file with MM/ML tags in read headers (from basecaller or prior tagging)
    #[arg(long, help_heading = "Required")]
    pub reads: PathBuf,

    /// Output path for learned model JSON (use with 'simulate --model-in')
    #[arg(long, help_heading = "Required")]
    pub model_out: PathBuf,

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
    pub motifs: Vec<String>,

    /// File with motif specifications (one per line, '#' for comments)
    #[arg(long, help_heading = "Motif Definitions (at least one required)")]
    pub motifs_file: Option<PathBuf>,

    // ========== Model Fitting Parameters ==========
    /// Probability threshold (0.0-1.0) to classify sites as methylated vs unmethylated
    ///
    /// ML values are converted to probabilities (ML/255), then compared to this threshold
    /// to fit separate Beta distributions for methylated and unmethylated populations
    #[arg(long, default_value_t = 0.5, help_heading = "Model Fitting Parameters")]
    pub model_threshold: f64,

    /// Limit the number of reads used for fitting (processes the first N reads)
    #[arg(long, help_heading = "Model Fitting Parameters")]
    pub n_reads: Option<usize>,
}

pub fn validate_simulate_inputs(args: &SimulateArgs) -> anyhow::Result<()> {
    use anyhow::bail;
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

pub fn validate_fit_inputs(args: &FitModelArgs) -> anyhow::Result<()> {
    use anyhow::bail;
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

#[cfg(test)]
mod cli_tests {
    use super::*;

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
