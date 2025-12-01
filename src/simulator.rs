use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use anyhow::{anyhow, bail, Context, Result};
use rand::seq::SliceRandom;
use rand::{rngs::StdRng, Rng, SeedableRng};
use rand_distr::{Distribution, Exp};
use shell_words::split as shell_split;
use tempfile::tempdir;

use crate::fastx::{read_fasta_sequences, read_fastq_records, ReadRecord};

#[derive(Clone, Debug)]
pub struct BuiltinSimulationConfig {
    pub references: Vec<String>,
    pub profile: ErrorProfile,
    pub length_spec: ReadLengthSpec,
    pub num_reads: usize,
    pub name_prefix: String,
    pub seed: u64,
}

#[derive(Clone, Debug)]
pub struct BadreadsSimulationConfig {
    pub reference: PathBuf,
    pub num_reads: usize,
    pub read_length: usize,
    pub executable: Option<PathBuf>,
    pub extra_args: Option<String>,
    pub seed: u64,
}

pub trait SimulatorStrategy {
    fn simulate(&mut self) -> Result<Vec<ReadRecord>>;
}

#[derive(Copy, Clone, Debug)]
pub struct ErrorProfile {
    pub substitution_rate: f64,
    pub insertion_rate: f64,
    pub deletion_rate: f64,
}

#[derive(Copy, Clone, Debug)]
pub enum ReadLengthSpec {
    Mode(usize),
    Mean(usize),
    N50(usize),
}

pub fn load_references(path: &Path) -> Result<Vec<String>> {
    read_fasta_sequences(path)
}

pub fn run_badreads(
    reference: &Path,
    num_reads: usize,
    read_length: Option<usize>,
    executable: Option<&PathBuf>,
    extra_args: Option<&str>,
    seed: Option<u64>,
) -> Result<Vec<ReadRecord>> {
    let badreads_exec = executable
        .cloned()
        .unwrap_or_else(|| PathBuf::from("badread"));
    let tmp_dir = tempdir()?;
    let fastq_path = tmp_dir.path().join("badreads.fastq");
    let quantity_bases = if let Some(length) = read_length {
        (num_reads as u128).saturating_mul(length as u128)
    } else {
        num_reads as u128
    };
    let quantity_arg = quantity_bases.to_string();

    let mut cmd = Command::new(&badreads_exec);
    cmd.arg("simulate")
        .arg("--reference")
        .arg(reference)
        .arg("--quantity")
        .arg(&quantity_arg);
    if let Some(length) = read_length {
        let std_dev = ((length as f64 * 0.3).round() as usize).max(1);
        // badread expects "mean,stdev" for --length; use the requested length as the mean
        // with a modest spread so reads stay near the desired size.
        cmd.arg("--length").arg(format!("{length},{std_dev}"));
    }
    if let Some(seed_val) = seed {
        cmd.arg("--seed").arg(seed_val.to_string());
    }
    if let Some(extra) = extra_args {
        for token in
            shell_split(extra).map_err(|e| anyhow!("Failed to parse --badreads-extra: {e}"))?
        {
            cmd.arg(token);
        }
    }

    let output_file = File::create(&fastq_path)
        .with_context(|| format!("Failed to create temporary file '{}'", fastq_path.display()))?;
    cmd.stdout(Stdio::from(output_file));

    let status = cmd.status().with_context(|| {
        format!(
            "Failed to invoke badread executable '{}'",
            badreads_exec.display()
        )
    })?;
    if !status.success() {
        bail!("badread exited with status {}", status);
    }
    let reads = read_fastq_records(&fastq_path)?;
    Ok(reads)
}

pub enum SimulatorKindStrategy {
    Builtin(BuiltinStrategy),
    Badreads(BadreadsStrategy),
}

impl SimulatorStrategy for SimulatorKindStrategy {
    fn simulate(&mut self) -> Result<Vec<ReadRecord>> {
        match self {
            SimulatorKindStrategy::Builtin(inner) => inner.simulate(),
            SimulatorKindStrategy::Badreads(inner) => inner.simulate(),
        }
    }
}

pub struct BuiltinSimulator {
    references: Vec<String>,
    profile: ErrorProfile,
    rng: StdRng,
    length_sampler: ReadLengthSampler,
}

impl BuiltinSimulator {
    pub fn new(
        references: Vec<String>,
        profile: ErrorProfile,
        length_spec: ReadLengthSpec,
        seed: u64,
    ) -> Result<Self> {
        if references.is_empty() {
            bail!("Reference FASTA did not contain any sequences");
        }
        let length_sampler = ReadLengthSampler::from_spec(length_spec)?;
        Ok(Self {
            references,
            profile,
            rng: StdRng::seed_from_u64(seed),
            length_sampler,
        })
    }

    pub fn generate_read(&mut self, name: String) -> ReadRecord {
        let length = self.length_sampler.sample(&mut self.rng);
        let template = self.sample_template(length);
        let (sequence, quality) = self.introduce_errors(&template);
        ReadRecord::new(name, sequence.into_bytes(), quality.into_bytes())
    }

    fn sample_template(&mut self, length: usize) -> String {
        let template = self
            .references
            .choose(&mut self.rng)
            .expect("references cannot be empty")
            .clone();
        if template.len() <= length {
            template
        } else {
            let start = self.rng.gen_range(0..=template.len() - length);
            template[start..start + length].to_string()
        }
    }

    fn introduce_errors(&mut self, template: &str) -> (String, String) {
        let mut sequence = String::with_capacity(template.len());
        let mut quality = String::with_capacity(template.len());
        for ch in template.chars() {
            let base = ch.to_ascii_uppercase();
            if self
                .rng
                .gen_bool(self.profile.deletion_rate.clamp(0.0, 1.0))
            {
                continue;
            }
            let mut current_base = base;
            let mut mutated = false;
            if self
                .rng
                .gen_bool(self.profile.substitution_rate.clamp(0.0, 1.0))
            {
                current_base = self.random_base(Some(base));
                mutated = true;
            }
            sequence.push(current_base);
            quality.push(Self::quality_char(mutated));

            while self
                .rng
                .gen_bool(self.profile.insertion_rate.clamp(0.0, 1.0))
            {
                sequence.push(self.random_base(None));
                quality.push(Self::quality_char(true));
            }
        }
        if sequence.is_empty() {
            sequence.push(self.random_base(None));
            quality.push(Self::quality_char(true));
        }
        (sequence, quality)
    }

    fn random_base(&mut self, exclude: Option<char>) -> char {
        let mut bases = vec!['A', 'C', 'G', 'T'];
        if let Some(ex) = exclude {
            bases.retain(|b| *b != ex);
        }
        *bases
            .choose(&mut self.rng)
            .expect("List of bases cannot be empty")
    }

    fn quality_char(mutated: bool) -> char {
        let phred = if mutated { 18 } else { 38 };
        (phred + 33) as u8 as char
    }
}

pub struct BuiltinStrategy {
    simulator: BuiltinSimulator,
    num_reads: usize,
    name_prefix: String,
}

impl BuiltinStrategy {
    pub fn new(config: BuiltinSimulationConfig) -> Result<Self> {
        let simulator = BuiltinSimulator::new(
            config.references,
            config.profile,
            config.length_spec,
            config.seed,
        )?;
        Ok(Self {
            simulator,
            num_reads: config.num_reads,
            name_prefix: config.name_prefix,
        })
    }
}

impl SimulatorStrategy for BuiltinStrategy {
    fn simulate(&mut self) -> Result<Vec<ReadRecord>> {
        let mut reads = Vec::with_capacity(self.num_reads);
        let progress_interval = (self.num_reads / 10).max(1).min(1000);
        for idx in 0..self.num_reads {
            let name = format!("{}_{}", self.name_prefix, format!("{:06}", idx + 1));
            let read = self.simulator.generate_read(name);
            if (idx + 1) % progress_interval == 0 || idx + 1 == self.num_reads {
                let pct = (idx + 1) as f64 / self.num_reads as f64 * 100.0;
                log::info!(
                    "Progress: {}/{} reads ({:.1}%)",
                    idx + 1,
                    self.num_reads,
                    pct
                );
            }
            reads.push(read);
        }
        Ok(reads)
    }
}

pub struct BadreadsStrategy {
    config: BadreadsSimulationConfig,
}

impl BadreadsStrategy {
    pub fn new(config: BadreadsSimulationConfig) -> Self {
        Self { config }
    }
}

impl SimulatorStrategy for BadreadsStrategy {
    fn simulate(&mut self) -> Result<Vec<ReadRecord>> {
        run_badreads(
            &self.config.reference,
            self.config.num_reads,
            Some(self.config.read_length),
            self.config.executable.as_ref(),
            self.config.extra_args.as_deref(),
            Some(self.config.seed),
        )
    }
}

struct ReadLengthSampler {
    peak: f64,
    min_length: usize,
    exp: Exp<f64>,
}

impl ReadLengthSampler {
    fn from_spec(spec: ReadLengthSpec) -> Result<Self> {
        let (peak, tail_mean) = match spec {
            ReadLengthSpec::Mode(mode) => {
                if mode == 0 {
                    bail!("Peak read length must be greater than zero");
                }
                let tail_mean = (mode as f64 * 0.3).max(50.0);
                (mode as f64, tail_mean)
            }
            ReadLengthSpec::Mean(mean) => {
                if mean == 0 {
                    bail!("Mean read length must be greater than zero");
                }
                let tail_mean = (mean as f64 * 0.25).max(50.0);
                (mean as f64 + tail_mean, tail_mean)
            }
            ReadLengthSpec::N50(n50) => {
                if n50 == 0 {
                    bail!("Read length N50 must be greater than zero");
                }
                let tail_mean = (n50 as f64 * 0.2).max(50.0);
                let peak = n50 as f64 + tail_mean * std::f64::consts::LN_2;
                (peak, tail_mean)
            }
        };
        let lambda = (1.0 / tail_mean).max(1e-6);
        let exp = Exp::new(lambda)
            .map_err(|_| anyhow!("Failed to configure read length distribution"))?;
        let mut min_length = ((peak * 0.1) as usize).max(100);
        if min_length >= peak as usize {
            min_length = peak.floor().max(1.0) as usize;
        }
        Ok(Self {
            peak,
            min_length: min_length.max(1),
            exp,
        })
    }

    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> usize {
        let delta = self.exp.sample(rng);
        let candidate = (self.peak - delta).round();
        let clamped = candidate.max(self.min_length as f64).max(1.0);
        clamped as usize
    }
}
