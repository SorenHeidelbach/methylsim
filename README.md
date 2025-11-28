# methylsim

`methylsim` creates Oxford Nanopore–style FASTQ files with MM/ML methylation tags or adds tags to reads that already exist. It can:

- Generate synthetic reads with an internal simulator.
- Call [badread](https://github.com/rrwick/Badread) for more detailed error models.
- Read FASTQ/FASTA files and annotate them in place.

## Requirements
- Rust toolchain (1.74+ recommended) for building.
- `badread` executable when using `--simulator badreads`.

## Installation
### Build locally
```bash
git clone git@github.com:SorenHeidelbach/methylsim.git
cd methylsim
cargo build --release
# binary at target/release/methylsim
```

### Install into `$HOME/.cargo/bin`
```bash
git clone git@github.com:SorenHeidelbach/methylsim.git
cd methylsim
cargo install --path .
```

## Basic usage
Simulate 1,000 reads with a G6mATC motif methylated:

```bash
methylsim \
  --reference resources/reference.fasta \
  --motif GATC_6mA_1 \
  --num-reads 1000 \
  --read-length 5000 \
  --output-fastq results/GATC_a_1.fastq \
  --tags-tsv results/GATC_a_1_tags.tsv
```

Important global flags:
- `--reference`: FASTA used for simulation (required unless you pass `--reads`).
- `--reads`: FASTA/FASTQ to annotate instead of simulating.
- `--motif` / `--motifs-file`: motif definition(s); at least one is required.
- `--seed`: RNG seed shared by the simulator and methylation tagger (default `1`).
- `--output-fastq`: FASTQ output path (default `methylsim.fastq`).
- `--tags-tsv`: optional TSV output containing `read_id`, `MM`, `ML`.

## Use cases

### 1. Builtin simulator (fast and dependency-free)
Use when you need synthetic reads without installing extra tools. Provide a reference genome and describe how long the reads should be.

```bash
methylsim \
  --reference references/bacteria.fa \
  --motif GATC:C:3:m \
  --num-reads 5000 \
  --read-length 4000 \
  --substitution-rate 0.03 \
  --insertion-rate 0.015 \
  --deletion-rate 0.015 \
  --output-fastq results/builtin.fastq \
  --tags-tsv results/builtin_tags.tsv
```

Key simulator flags:
- `--num-reads`: number of reads to emit.
- `--read-length`, `--read-length-mean`, `--read-length-n50`: different ways to control the length distribution.
- `--substitution-rate`, `--insertion-rate`, `--deletion-rate`: simple ONT-style error rates.
- `--name-prefix`: prefix for read identifiers (`methylsim_000001`, … by default).

### 2. badread-backed simulation
Choose this when you want badread’s built-in models or coverage presets. `methylsim` still handles motif tagging and output files; it just sources sequences from badread.

```bash
methylsim \
  --reference resources/reference.fasta \
  --motif CCWGG_5mC_1 \
  --simulator badreads \
  --read-length 10000 \
  --badreads-extra "--quantity 20x --error_model nanopore2023" \
  --output-fastq results/badread.fastq \
  --tags-tsv results/badread_tags.tsv
```

Flags realted to this mode:
- `--simulator badreads`: switch from the builtin simulator.
- `--badreads-exec`: path to the `badread` binary (defaults to `badread` on `$PATH`).
- `--badreads-extra`: quoted string that is split with shell rules and sent directly to badread (`--quantity`, `--error_model`, etc.).

If you also pass `--reads`, no simulation happens and the simulator flag is ignored.

### 3. Annotate an existing FASTQ or FASTA
Use this path to tag reads produced by real devices or other simulators. Only motifs and `--reads` are required.

```bash
methylsim \
  --reads data/real_run.fastq.gz \
  --motifs-file configs/cpg_motifs.txt \
  --motif-high-prob 0.8 \
  --non-motif-high-prob 0.02 \
  --output-fastq results/real_run_tagged.fastq \
  --tags-tsv results/real_run_tags.tsv
```

Notes:
- FASTA inputs get dummy high-quality strings so they remain valid FASTQ.
- Existing MM/ML annotations inside FASTQ headers are stripped before new ones are written.


## Motif definitions and methylation model
- Colon syntax: `motif[:canonical[:offset[:mod[:strand]]]]`
  - Example: `GATC:C:2:m:+` targets the cytosine two bases into `GATC`.
- Underscore syntax: `sequence_modtype_offset[_strand]`
  - Example: `CGCG_m_1_-` targets the cytosine at offset 1 on the reverse strand.

Repeat `--motif` or use comma-separated entries (e.g. `--motif CG,GC_m_0`). Use `--motifs-file` for one-per-line specifications (lines starting with `#` are ignored).

Methylation-related flags:

| Flag | Description |
| --- | --- |
| `--motif-high-prob` | Probability that a motif hit is “high” methylated (default 0.95). |
| `--non-motif-high-prob` | Background probability for other canonical bases (default 0.01). |
| `--high-ml-mean`, `--high-ml-std` | Distribution parameters for ML values of methylated bases. |
| `--low-ml-mean`, `--low-ml-std` | Distribution parameters for ML values of unmethylated bases. |

## Outputs
- **FASTQ** (`--output-fastq`, default `methylsim.fastq`): sequences with MM/ML tags appended to the header.
- **Tags TSV** (`--tags-tsv`): columns `read_id`, `MM`, `ML` for reuse or auditing.
