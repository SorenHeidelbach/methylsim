mod fastx;
mod motif;
mod simulator;
mod tagging;

use std::collections::HashMap;
use std::fs::{File, read_to_string};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use anyhow::{Context, Result, anyhow, bail};
use clap::{Parser, ValueEnum};

use fastx::{ReadRecord, read_fastx};
use motif::MotifDefinition;
use simulator::{BuiltinSimulator, ErrorProfile, ReadLengthSpec, load_references, run_badreads};
use tagging::MethylationTagger;

#[derive(Copy, Clone, Debug, ValueEnum, PartialEq, Eq)]
enum SimulatorKind {
    Builtin,
    Badreads,
}

#[derive(Parser, Debug)]
#[command(author, version, about = "Simulate nanopore reads with MM/ML tags", long_about = None)]
struct Cli {
    /// Existing reads (FASTQ/FASTA) to annotate instead of simulating
    #[arg(long, help_heading = "Input/Annotation")]
    reads: Option<PathBuf>,

    /// Reference FASTA used when simulating reads
    #[arg(long, help_heading = "Input/Annotation")]
    reference: Option<PathBuf>,

    /// One or more motif specifications (motif[:canonical[:offset[:mod[:strand]]]] or sequence_modtype_offset[_strand])
    #[arg(
        long = "motif",
        num_args = 1..,
        value_delimiter = ',',
        required_unless_present = "motifs_file",
        help_heading = "Motif Definitions"
    )]
    motifs: Vec<String>,

    /// File containing motif specifications (one per line, '#' for comments)
    #[arg(long, help_heading = "Motif Definitions")]
    motifs_file: Option<PathBuf>,

    /// Simulation mode for generating reads when --reads is not supplied
    #[arg(long, default_value_t = SimulatorKind::Builtin, value_enum, help_heading = "Simulation")]
    simulator: SimulatorKind,

    /// Number of reads to simulate
    #[arg(long, default_value_t = 100, help_heading = "Simulation")]
    num_reads: usize,

    /// Target read length (treated as the peak of the distribution)
    #[arg(long, default_value_t = 5000, help_heading = "Simulation")]
    read_length: usize,

    /// Mean read length for simulated reads (overrides --read-length)
    #[arg(long, help_heading = "Simulation")]
    read_length_mean: Option<usize>,

    /// Read length N50 for simulated reads (overrides --read-length)
    #[arg(long, help_heading = "Simulation")]
    read_length_n50: Option<usize>,

    /// Prefix used for naming simulated reads
    #[arg(long, default_value = "methylsim", help_heading = "Simulation")]
    name_prefix: String,

    /// Path to badreads executable
    #[arg(long, help_heading = "Badreads Options")]
    badreads_exec: Option<PathBuf>,

    /// Additional arguments forwarded to badreads (quoted string)
    #[arg(long, help_heading = "Badreads Options", allow_hyphen_values = true)]
    badreads_extra: Option<String>,

    /// Probability that a motif hit carries high methylation (true positive rate)
    #[arg(long, default_value_t = 0.95, help_heading = "Methylation Model")]
    motif_high_prob: f64,

    /// Probability that a non-motif canonical base carries high methylation (false positive rate)
    #[arg(long, default_value_t = 0.01, help_heading = "Methylation Model")]
    non_motif_high_prob: f64,

    /// Mean ML tag value (0-255) for high methylation events
    #[arg(long, default_value_t = 230.0, help_heading = "Methylation Model")]
    high_ml_mean: f64,

    /// Standard deviation for high methylation ML values
    #[arg(long, default_value_t = 10.0, help_heading = "Methylation Model")]
    high_ml_std: f64,

    /// Mean ML tag value for low methylation events
    #[arg(long, default_value_t = 20.0, help_heading = "Methylation Model")]
    low_ml_mean: f64,

    /// Standard deviation for low methylation ML values
    #[arg(long, default_value_t = 5.0, help_heading = "Methylation Model")]
    low_ml_std: f64,

    /// Output FASTQ file
    #[arg(long, default_value = "methylsim.fastq", help_heading = "Output")]
    output_fastq: PathBuf,

    /// Optional TSV output listing MM/ML tags per read
    #[arg(long, help_heading = "Output")]
    tags_tsv: Option<PathBuf>,

    /// Existing TSV file of MM/ML tags used when tagging BAMs
    #[arg(long, help_heading = "BAM Tagging")]
    tags_tsv_input: Option<PathBuf>,

    /// RNG seed shared by simulator and methylation annotator
    #[arg(long, default_value_t = 1, help_heading = "Simulation")]
    seed: u64,

    /// Substitution rate used by builtin simulator
    #[arg(long, default_value_t = 0.03, help_heading = "Simulation")]
    substitution_rate: f64,

    /// Insertion rate used by builtin simulator
    #[arg(long, default_value_t = 0.01, help_heading = "Simulation")]
    insertion_rate: f64,

    /// Deletion rate used by builtin simulator
    #[arg(long, default_value_t = 0.01, help_heading = "Simulation")]
    deletion_rate: f64,

    /// Input BAM whose reads should be tagged using --tags-tsv-input
    #[arg(long, help_heading = "BAM Tagging")]
    bam_input: Option<PathBuf>,

    /// Output BAM path when tagging BAM files
    #[arg(long, help_heading = "BAM Tagging")]
    bam_output: Option<PathBuf>,
}

#[derive(Default)]
struct WriteStats {
    total_reads: usize,
    tagged_reads: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    run(cli)
}

fn run(cli: Cli) -> Result<()> {
    validate_inputs(&cli)?;
    if cli.reads.is_none() && cli.reference.is_none() {
        let bam_in = cli
            .bam_input
            .as_ref()
            .expect("--bam-input must be supplied in tag-only mode");
        let bam_out = cli
            .bam_output
            .as_ref()
            .expect("--bam-output must be supplied in tag-only mode");
        let tags_path = cli
            .tags_tsv_input
            .as_ref()
            .or(cli.tags_tsv.as_ref())
            .expect("--bam-input requires --tags-tsv-input or --tags-tsv");
        apply_tags_to_bam(bam_in, bam_out, tags_path)?;
        return Ok(());
    }
    let motif_specs = collect_motif_specs(&cli)?;
    let motifs = motif_specs
        .iter()
        .map(|spec| MotifDefinition::parse(spec))
        .collect::<Result<Vec<_>>>()?;
    let mut tagger = MethylationTagger::new(
        motifs,
        cli.motif_high_prob,
        cli.non_motif_high_prob,
        cli.high_ml_mean,
        cli.high_ml_std,
        cli.low_ml_mean,
        cli.low_ml_std,
        cli.seed,
    )?;

    let mut fastq_writer = BufWriter::new(
        File::create(&cli.output_fastq)
            .with_context(|| format!("Failed to create '{}'", cli.output_fastq.display()))?,
    );
    let mut tags_writer = if let Some(path) = &cli.tags_tsv {
        let mut writer = BufWriter::new(
            File::create(path).with_context(|| format!("Failed to create '{}'", path.display()))?,
        );
        writer.write_all(b"read_id\tMM\tML\n")?;
        Some(writer)
    } else {
        None
    };
    let mut stats = WriteStats::default();

    if let Some(reads_path) = &cli.reads {
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
        match cli.simulator {
            SimulatorKind::Builtin => simulate_builtin(
                &cli,
                &mut tagger,
                &mut fastq_writer,
                &mut tags_writer,
                &mut stats,
            )?,
            SimulatorKind::Badreads => simulate_badreads(
                &cli,
                &mut tagger,
                &mut fastq_writer,
                &mut tags_writer,
                &mut stats,
            )?,
        }
    }

    fastq_writer.flush()?;
    if let Some(writer) = tags_writer.as_mut() {
        writer.flush()?;
    }

    println!(
        "Wrote {} reads to {} ({} reads carried MM/ML tags)",
        stats.total_reads,
        cli.output_fastq.display(),
        stats.tagged_reads
    );
    if let Some(bam_in) = cli.bam_input.as_ref() {
        let bam_out = cli
            .bam_output
            .as_ref()
            .expect("--bam-output must be supplied when --bam-input is used");
        let tags_path = cli
            .tags_tsv_input
            .as_ref()
            .or(cli.tags_tsv.as_ref())
            .expect("--bam-input requires --tags-tsv-input or --tags-tsv");
        apply_tags_to_bam(bam_in, bam_out, tags_path)?;
    }
    Ok(())
}

fn collect_motif_specs(cli: &Cli) -> Result<Vec<String>> {
    let mut specs = cli.motifs.clone();
    if let Some(path) = &cli.motifs_file {
        let contents = read_to_string(path)
            .with_context(|| format!("Failed to read motif file '{}'", path.display()))?;
        for line in contents.lines() {
            let trimmed = line.split('#').next().unwrap_or("").trim();
            if trimmed.is_empty() {
                continue;
            }
            specs.push(trimmed.to_string());
        }
        if specs.is_empty() {
            bail!(
                "Motif file '{}' did not contain any usable entries",
                path.display()
            );
        }
    }
    if specs.is_empty() {
        bail!("At least one motif specification must be supplied");
    }
    Ok(specs)
}

fn validate_inputs(cli: &Cli) -> Result<()> {
    if cli.reads.is_none() && cli.reference.is_none() && cli.bam_input.is_none() {
        bail!("Either --reads, --reference, or --bam-input must be supplied");
    }
    if cli.bam_input.is_some() {
        if cli.bam_output.is_none() {
            bail!("--bam-input requires --bam-output");
        }
        if cli.tags_tsv_input.is_none() && cli.tags_tsv.is_none() {
            bail!("--bam-input requires --tags-tsv-input or --tags-tsv");
        }
    }
    let performing_read_work = cli.reads.is_some() || cli.reference.is_some();
    if !performing_read_work {
        return Ok(());
    }
    if cli.reads.is_none() && cli.reference.is_none() {
        bail!("Either --reads or --reference must be supplied");
    }
    if cli.simulator == SimulatorKind::Badreads && cli.reads.is_some() {
        eprintln!("Warning: --reads provided; ignoring --simulator badreads");
    }
    if cli.reads.is_none() && cli.simulator == SimulatorKind::Badreads && cli.reference.is_none() {
        bail!("--simulator badreads requires --reference");
    }
    if cli.read_length == 0 {
        bail!("--read-length must be greater than zero");
    }
    if let Some(mean) = cli.read_length_mean {
        if mean == 0 {
            bail!("--read-length-mean must be greater than zero");
        }
    }
    if let Some(n50) = cli.read_length_n50 {
        if n50 == 0 {
            bail!("--read-length-n50 must be greater than zero");
        }
    }
    if cli.read_length_mean.is_some() && cli.read_length_n50.is_some() {
        bail!("Only one of --read-length-mean or --read-length-n50 can be supplied");
    }
    if !(0.0..=1.0).contains(&cli.motif_high_prob) {
        bail!("--motif-high-prob must be between 0 and 1");
    }
    if !(0.0..=1.0).contains(&cli.non_motif_high_prob) {
        bail!("--non-motif-high-prob must be between 0 and 1");
    }
    if !(0.0..=255.0).contains(&cli.high_ml_mean) {
        bail!("--high-ml-mean must lie between 0 and 255");
    }
    if cli.high_ml_std.is_sign_negative() {
        bail!("--high-ml-std must be non-negative");
    }
    if !(0.0..=255.0).contains(&cli.low_ml_mean) {
        bail!("--low-ml-mean must lie between 0 and 255");
    }
    if cli.low_ml_std.is_sign_negative() {
        bail!("--low-ml-std must be non-negative");
    }
    Ok(())
}

fn simulate_builtin(
    cli: &Cli,
    tagger: &mut MethylationTagger,
    fastq_writer: &mut BufWriter<File>,
    tags_writer: &mut Option<BufWriter<File>>,
    stats: &mut WriteStats,
) -> Result<()> {
    let reference_path = cli
        .reference
        .as_ref()
        .ok_or_else(|| anyhow!("Builtin simulator requires --reference"))?;
    eprintln!(
        "Running builtin simulator using reference '{}'...",
        reference_path.display()
    );
    let references = load_references(reference_path)?;
    let profile = ErrorProfile {
        substitution_rate: cli.substitution_rate,
        insertion_rate: cli.insertion_rate,
        deletion_rate: cli.deletion_rate,
    };
    let length_spec = if let Some(mean) = cli.read_length_mean {
        ReadLengthSpec::Mean(mean)
    } else if let Some(n50) = cli.read_length_n50 {
        ReadLengthSpec::N50(n50)
    } else {
        ReadLengthSpec::Mode(cli.read_length)
    };
    let mut simulator = BuiltinSimulator::new(references, profile, length_spec, cli.seed)?;
    for idx in 0..cli.num_reads {
        let name = format!("{}_{}", cli.name_prefix, format!("{:06}", idx + 1));
        let read = simulator.generate_read(name);
        write_read(&read, tagger, fastq_writer, tags_writer, stats)?;
    }
    Ok(())
}

fn simulate_badreads(
    cli: &Cli,
    tagger: &mut MethylationTagger,
    fastq_writer: &mut BufWriter<File>,
    tags_writer: &mut Option<BufWriter<File>>,
    stats: &mut WriteStats,
) -> Result<()> {
    let reference_path = cli
        .reference
        .as_ref()
        .ok_or_else(|| anyhow!("badread simulation requires --reference"))?;
    eprintln!(
        "Invoking badread simulator with reference '{}'...",
        reference_path.display()
    );
    let reads = run_badreads(
        reference_path,
        cli.num_reads,
        Some(cli.read_length),
        cli.badreads_exec.as_ref(),
        cli.badreads_extra.as_deref(),
        Some(cli.seed),
    )?;
    for read in reads {
        write_read(&read, tagger, fastq_writer, tags_writer, stats)?;
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
    let tags = tagger.annotate(&read.sequence);
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
    writeln!(fastq_writer, "{}", read.sequence)?;
    writeln!(fastq_writer, "+")?;
    writeln!(fastq_writer, "{}", read.quality)?;

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
    Ok(())
}

#[derive(Debug)]
struct TagStrings {
    mm: Option<String>,
    ml: Option<String>,
}

fn load_tags_table(path: &Path) -> Result<HashMap<String, TagStrings>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open tags TSV '{}'", path.display()))?;
    let reader = BufReader::new(file);
    let mut tags: HashMap<String, TagStrings> = HashMap::new();
    for (idx, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Failed reading '{}'", path.display()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if idx == 0 && trimmed.starts_with("read_id") {
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let read_id = parts[0].to_string();
        let mm = parts[1].trim();
        let ml = parts[2].trim();
        tags.insert(
            read_id,
            TagStrings {
                mm: if mm.is_empty() {
                    None
                } else {
                    Some(mm.to_string())
                },
                ml: if ml.is_empty() {
                    None
                } else {
                    Some(ml.to_string())
                },
            },
        );
    }
    Ok(tags)
}

fn apply_tags_to_bam(bam_in: &Path, bam_out: &Path, tags_path: &Path) -> Result<()> {
    let tags = load_tags_table(tags_path)?;
    if tags.is_empty() {
        bail!(
            "Tags TSV '{}' did not contain any entries",
            tags_path.display()
        );
    }
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
    let mut patched = 0usize;
    let mut missing = 0usize;
    loop {
        line.clear();
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        if line.starts_with('@') {
            writer.write_all(line.as_bytes())?;
            continue;
        }
        while line.ends_with('\n') || line.ends_with('\r') {
            line.pop();
        }
        if line.is_empty() {
            writer.write_all(b"\n")?;
            continue;
        }
        let mut fields: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if let Some(entry) = tags.get(fields[0].as_str()) {
            fields.retain(|field| !(field.starts_with("MM:Z:") || field.starts_with("ML:B:")));
            if let Some(mm) = &entry.mm {
                if !mm.is_empty() {
                    fields.push(mm.clone());
                }
            }
            if let Some(ml) = &entry.ml {
                if !ml.is_empty() {
                    fields.push(ml.clone());
                }
            }
            patched += 1;
        } else {
            missing += 1;
        }
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
    eprintln!(
        "Tagged BAM: updated {} reads ({} had no matching tags) -> '{}'",
        patched,
        missing,
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
        assert_eq!(cli.badreads_extra.as_deref(), Some("--quantity 50x"));
    }
}
