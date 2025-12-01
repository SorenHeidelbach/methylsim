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

`methylsim` has two main subcommands:

### 1. `simulate` - Generate or tag reads with methylation
Simulate reads with a G6mATC motif to achieve 50x coverage:

```bash
methylsim simulate \
  --reference resources/reference.fasta \
  --motif GATC_6mA_1 \
  --quantity 50x \
  --output-fastq results/GATC_a_1.fastq
```

Important global flags:
- `--reference`: FASTA used for simulation (required unless you pass `--reads`).
- `--reads`: FASTA/FASTQ to annotate instead of simulating.
- `--motif` / `--motifs-file`: motif definition(s); the file form expects `motif`, `mod_type`, `mod_position` columns (CSV/TSV) matching the methylation_phasing format (bin/contig columns are ignored).
- `--motif` / `--motifs-file`: motif definition(s); the file form expects `motif`, `mod_type`, `mod_position` columns (CSV/TSV) matching the methylation_phasing format (bin/contig columns are ignored).
- `--seed`: RNG seed shared by the simulator and methylation tagger (default `1`).
- `--output-fastq`: FASTQ output path (default `methylsim.fastq`).

Or specify absolute amount of sequence (e.g., 250 million bases):

```bash
methylsim simulate \
  --reference resources/reference.fasta \
  --motif GATC_6mA_1 \
  --quantity 250M \
  --read-length 5000 \
  --output-fastq results/GATC_a_1.fastq
```

### 2. `fit-model` - Learn methylation model from tagged reads
Extract methylation patterns from real basecalled data:

```bash
methylsim fit-model \
  --reads basecalled_reads.fastq \
  --motif GATC_6mA_1 \
  --model-out models/ecoli_gatc.json
```

To speed up exploratory runs, you can limit fitting to the first N reads:

```bash
methylsim fit-model \
  --reads basecalled_reads.fastq \
  --motif GATC_6mA_1 \
  --model-out models/ecoli_gatc.json \
  --n-reads 1000
```

Then use the learned model to simulate realistic reads:

```bash
methylsim simulate \
  --reference resources/reference.fasta \
  --model-in models/ecoli_gatc.json \
  --motif GATC_6mA_1 \
  --output-fastq results/realistic_reads.fastq
```

### Key flags:
- **Input modes** (choose one):
  - `--reference`: FASTA for simulation
  - `--reads`: Existing FASTQ/FASTA to tag
  - `--bam-input`: BAM to tag directly (reads SEQ field, adds MM/ML tags)
- **Motifs** (required): `--motif` or `--motifs-file`
- **Quantity** (simulation):
  - `--quantity 50x`: Coverage-based (50x coverage)
  - `--quantity 250M`: Absolute bases (250 million bases)
  - `--num-reads`: Legacy option (number of reads)
- **Model**: `--model-in` to use learned model (optional)
- **Output**: `--output-fastq` (FASTQ) and `--tags-tsv` (TSV table)

## Use cases

### 1. Builtin simulator (fast and dependency-free)
Use when you need synthetic reads without installing extra tools. Provide a reference genome and describe how long the reads should be.

```bash
methylsim simulate \
  --reference references/bacteria.fa \
  --motif GATC:C:3:m \
  --quantity 25x \
  --read-length 4000 \
  --substitution-rate 0.03 \
  --insertion-rate 0.015 \
  --deletion-rate 0.015 \
  --output-fastq results/builtin.fastq \
  --tags-tsv results/builtin_tags.tsv
```

Key simulator flags:
- `--quantity`: Coverage (e.g., `25x`) or absolute bases (e.g., `100M`). Calculates number of reads automatically.
- `--num-reads`: (Legacy) Explicit number of reads to generate.
- `--read-length`, `--read-length-mean`, `--read-length-n50`: Different ways to control the length distribution.
- `--substitution-rate`, `--insertion-rate`, `--deletion-rate`: Simple ONT-style error rates.
- `--name-prefix`: Prefix for read identifiers (`methylsim_000001`, … by default).

### 2. badread-backed simulation
Choose this when you want badread's built-in models or coverage presets. `methylsim` still handles motif tagging and output files; it just sources sequences from badread.

```bash
methylsim simulate \
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
methylsim simulate \
  --reads real_run.fastq \
  --motifs-file motifs.txt \
  --output-fastq real_run_tagged.fastq
```

Notes:
- FASTA inputs get dummy high-quality strings so they remain valid FASTQ.

### 4. Learn and reuse methylation models
If you already have reads with MM/ML tags (e.g., real E. coli data with `GATC_6mA_1`), you can learn the per-motif methylation probability distribution and reuse it when simulating or re-tagging reads.

**Input requirements:** FASTQ reads must already carry MM/ML tags in the header (e.g., produced by ONT basecalling and `samtools fastq`). Only FASTQ is supported for learning currently.

#### Single Motif Workflow:

**Step 1:** Learn from tagged reads (here: E. coli with `GATC_6mA_1`)
```bash
methylsim fit-model \
  --reads data/ecoli_tagged.fastq \
  --motif GATC_6mA_1 \
  --model-out models/ecoli_gatc.json
```
- Splits ML values into methylated vs unmethylated using `--model-threshold` (probability = ML/255; default 0.5)
- Stores an empirical histogram plus two Beta fits and the observed methylated fraction

**Step 2:** Use the learned model to simulate realistic reads
```bash
methylsim simulate \
  --reference references/ecoli.fa \
  --motif GATC_6mA_1 \
  --model-in models/ecoli_gatc.json \
  --num-reads 500 \
  --output-fastq results/ecoli_sim.fastq
```
- The learned model drives motif probabilities instead of `--motif-high-prob`
- ML value distributions match the empirical data

#### Multi-Motif Models:

**Models support multiple motifs!** Learn methylation patterns for multiple motifs simultaneously:

```bash
# Learn model for both Dam (GATC) and Dcm (CCWGG) methylation
methylsim fit-model \
  --reads data/ecoli_tagged.fastq \
  --motif GATC_6mA_1,CCWGG_5mC_1 \
  --model-out models/ecoli_multi.json
```

The tool will print per-motif statistics:
```
Fitting model for 2 motif(s):
  GATC (GATC_a_1): 1523 total sites, 1489 covered (97.8%), 94.2% methylated
  CCWGG (CCWGG_m_1): 832 total sites, 798 covered (95.9%), 89.5% methylated
```

Then use the multi-motif model for simulation:
```bash
methylsim simulate \
  --reference references/ecoli.fa \
  --motif GATC_6mA_1,CCWGG_5mC_1 \
  --model-in models/ecoli_multi.json \
  --quantity 50x \
  --output-fastq results/ecoli_sim.fastq
```

**Note:** If your FASTQ lacks MM/ML tags, you must tag it first (e.g., with a methylation caller or by running `methylsim simulate` in default high/low mode) before `fit-model` will learn anything.

### Bundled model presets

The crate ships with ready-to-use models in `models/`. The first preset is an E. coli Dam/Dcm model (`models/e_coli_G6mATC_C5mCWGG.toml`).
Use it with `--model-preset` (motifs are pulled from the preset automatically, so `--motif` is optional):

```bash
methylsim simulate \
  --reference references/ecoli.fa \
  --model-preset ecoli \
  --quantity 25x \
  --output-fastq results/ecoli_preset.fastq
```

### ML tag sampling

ML values come directly from methylation probabilities: a Beta draw yields the probability, then `ML = prob * 255`. Simple mode uses two default Beta templates scaled to `--motif-high-prob` and `--non-motif-high-prob`; learned models keep their own Beta parameters.
You can override the simple-mode Betas with `--high-beta-alpha/--high-beta-beta` and `--low-beta-alpha/--low-beta-beta`; the probability of sampling the high vs low Beta is controlled by `--motif-high-prob` (motif sites) and `--non-motif-high-prob` (background).


## Motif definitions and methylation model
- Colon syntax: `motif[:canonical[:offset[:mod[:strand]]]]`
  - Example: `GATC:C:2:m:+` targets the cytosine two bases into `GATC`.
- Underscore syntax: `sequence_modtype_offset[_strand]`
  - Example: `CGCG_m_1_-` targets the cytosine at offset 1 on the reverse strand.

Repeat `--motif` or use comma-separated entries (e.g. `--motif CG,GC_m_0`). Use `--motifs-file` to load the same TSV/CSV schema used by `methylation_phasing`: required columns `motif`, `mod_type`, `mod_position` plus optional `motif_complement`/`mod_position_complement`. Any `bin`, `id`, or `reference` columns present in shared motif panels are simply ignored.

Example `motifs.tsv`:

```
reference	motif	mod_position	mod_type
bin_1	CCWGG	1	m
bin_1	GCWGC	1	m
```
