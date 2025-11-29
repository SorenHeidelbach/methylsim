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

## Testing
Run the test suite to verify functionality:
```bash
cargo test
```

The test suite includes:
- **9 unit tests** covering motif parsing, model fitting, and FASTQ handling
- **16 integration tests** covering:
  - CLI help messages and argument validation
  - Read simulation with the builtin simulator
  - Coverage-based and absolute quantity specifications
  - Model fitting from tagged reads
  - Using learned models for simulation
  - Various parameter combinations (error rates, methylation probabilities, read lengths)
  - BAM tagging workflow validation

Test data is located in `tests/data/` including:
- `test_reference.fasta` - Small reference genome for testing
- `test_motifs.txt` - Example motif specification file

## Basic usage

`methylsim` has two main subcommands:

### 1. `simulate` - Generate or tag reads with methylation
Simulate reads with a G6mATC motif to achieve 50x coverage:

```bash
methylsim simulate \
  --reference resources/reference.fasta \
  --motif GATC_6mA_1 \
  --quantity 50x \
  --read-length 5000 \
  --output-fastq results/GATC_a_1.fastq \
  --tags-tsv results/GATC_a_1_tags.tsv
```

<<<<<<< HEAD
Important global flags:
- `--reference`: FASTA used for simulation (required unless you pass `--reads`).
- `--reads`: FASTA/FASTQ to annotate instead of simulating.
- `--motif` / `--motifs-file`: motif definition(s); the file form expects `motif`, `mod_type`, `mod_position` columns (CSV/TSV) matching the methylation_phasing format (bin/contig columns are ignored).
- `--motif` / `--motifs-file`: motif definition(s); the file form expects `motif`, `mod_type`, `mod_position` columns (CSV/TSV) matching the methylation_phasing format (bin/contig columns are ignored).
- `--seed`: RNG seed shared by the simulator and methylation tagger (default `1`).
- `--output-fastq`: FASTQ output path (default `methylsim.fastq`).
- `--tags-tsv`: optional TSV output containing `read_id`, `MM`, `ML`.
=======
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

Then use the learned model to simulate realistic reads:

```bash
methylsim simulate \
  --reference resources/reference.fasta \
  --model-in models/ecoli_gatc.json \
  --motif GATC_6mA_1 \
  --num-reads 1000 \
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
>>>>>>> 125d715 (Refactor code structure for improved readability and maintainability)

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
  --reads data/real_run.fastq.gz \
  --motifs-file configs/cpg_motifs.txt \
  --motif-high-prob 0.8 \
  --non-motif-high-prob 0.02 \
  --output-fastq results/real_run_tagged.fastq \
  --tags-tsv results/real_run_tags.tsv
```

Notes:
- FASTA inputs get dummy high-quality strings so they remain valid FASTQ.
- Existing MM/ML annotations inside FASTQ headers are preserved in the header comment; new tags are appended.

### 4. Learn and reuse methylation models
If you already have reads with MM/ML tags (e.g., real E. coli data with `GATC_6mA_1`), you can learn the per-motif methylation probability distribution and reuse it when simulating or re-tagging reads.

**Input requirements:** FASTQ reads must already carry MM/ML tags in the header (e.g., produced by ONT basecalling, `samtools fastq`, or another tool). Only FASTQ is supported for learning today; BAM/CRAM learning is not wired yet. The motifs you pass (`--motif`/`--motifs-file`) must match the MM tags, otherwise the learner will see zero sites.

#### Workflow:

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

**Note:** If your FASTQ lacks MM/ML tags, you must tag it first (e.g., with a methylation caller or by running `methylsim simulate` in default high/low mode) before `fit-model` will learn anything.


## Motif definitions and methylation model
- Colon syntax: `motif[:canonical[:offset[:mod[:strand]]]]`
  - Example: `GATC:C:2:m:+` targets the cytosine two bases into `GATC`.
- Underscore syntax: `sequence_modtype_offset[_strand]`
  - Example: `CGCG_m_1_-` targets the cytosine at offset 1 on the reverse strand.

Repeat `--motif` or use comma-separated entries (e.g. `--motif CG,GC_m_0`). Use `--motifs-file` to load the same TSV/CSV schema used by `methylation_phasing`: required columns `motif`, `mod_type`, `mod_position` plus optional `motif_complement`/`mod_position_complement`. Any `bin`, `id`, or `reference` columns present in shared motif panels are simply ignored.

Example `motifs.tsv`:

```
reference	motif	mod_position	mod_type	n_mod	n_nomod	motif_type	motif_complement	mod_position_complement	n_mod_complement	n_nomod_complement
bin_1	CCWGG	1	m	51	0	palindrome	CCWGG	1	51	0
bin_1	GCWGC	1	m	256	0	palindrome	GCWGC	1	256	0
```


<<<<<<< HEAD
Example `motifs.tsv`:

```
reference	motif	mod_position	mod_type	n_mod	n_nomod	motif_type	motif_complement	mod_position_complement	n_mod_complement	n_nomod_complement
bin_1	CCWGG	1	m	51	0	palindrome	CCWGG	1	51	0
bin_1	GCWGC	1	m	256	0	palindrome	GCWGC	1	256	0
```

# Full help output

```
Simulate nanopore reads with MM/ML tags


Usage: methylsim [OPTIONS]

Options:
  -h, --help     Print help
  -V, --version  Print version

Input/Annotation:
      --reads <READS>          Existing reads (FASTQ/FASTA) to annotate instead of simulating
      --reference <REFERENCE>  Reference FASTA used when simulating reads

Motif Definitions:
      --motif <MOTIFS>...          One or more motif specifications (motif[:canonical[:offset[:mod[:strand]]]] or sequence_modtype_offset[_strand])
      --motifs-file <MOTIFS_FILE>  File containing motif specifications (one per line, '#' for comments)

Simulation:
      --simulator <SIMULATOR>
          Simulation mode for generating reads when --reads is not supplied [default: builtin] [possible values: builtin, badreads]
      --num-reads <NUM_READS>
          Number of reads to simulate [default: 100]
      --read-length <READ_LENGTH>
          Target read length (treated as the peak of the distribution) [default: 5000]
      --read-length-mean <READ_LENGTH_MEAN>
          Mean read length for simulated reads (overrides --read-length)
      --read-length-n50 <READ_LENGTH_N50>
          Read length N50 for simulated reads (overrides --read-length)
      --name-prefix <NAME_PREFIX>
          Prefix used for naming simulated reads [default: methylsim]
      --seed <SEED>
          RNG seed shared by simulator and methylation annotator [default: 1]
      --substitution-rate <SUBSTITUTION_RATE>
          Substitution rate used by builtin simulator [default: 0.03]
      --insertion-rate <INSERTION_RATE>
          Insertion rate used by builtin simulator [default: 0.01]
      --deletion-rate <DELETION_RATE>
          Deletion rate used by builtin simulator [default: 0.01]

Badreads Options:
      --badreads-exec <BADREADS_EXEC>    Path to badreads executable
      --badreads-extra <BADREADS_EXTRA>  Additional arguments forwarded to badreads (quoted string)

Methylation Model:
      --motif-high-prob <MOTIF_HIGH_PROB>
          Probability that a motif hit carries high methylation (true positive rate) [default: 0.95]
      --non-motif-high-prob <NON_MOTIF_HIGH_PROB>
          Probability that a non-motif canonical base carries high methylation (false positive rate) [default: 0.01]
      --high-ml-mean <HIGH_ML_MEAN>
          Mean ML tag value (0-255) for high methylation events [default: 230]
      --high-ml-std <HIGH_ML_STD>
          Standard deviation for high methylation ML values [default: 10]
      --low-ml-mean <LOW_ML_MEAN>
          Mean ML tag value for low methylation events [default: 20]
      --low-ml-std <LOW_ML_STD>
          Standard deviation for low methylation ML values [default: 5]

Output:
      --output-fastq <OUTPUT_FASTQ>  Output FASTQ file [default: methylsim.fastq]
      --tags-tsv <TAGS_TSV>          Optional TSV output listing MM/ML tags per read```
```
=======
**For `simulate` command:**
| Flag | Description |
| --- | --- |
| `--model-in` | Load a previously learned model JSON (from `fit-model`). Overrides simple probability parameters. |
| `--motif-high-prob` | Probability that a motif hit is "high" methylated (default 0.95). Used when no model is loaded. |
| `--non-motif-high-prob` | Background probability for other canonical bases (default 0.01). |
| `--high-ml-mean`, `--high-ml-std` | Distribution parameters for ML values of methylated bases (defaults: 230, 10). |
| `--low-ml-mean`, `--low-ml-std` | Distribution parameters for ML values of unmethylated bases (defaults: 20, 5). |

**For `fit-model` command:**
| Flag | Description |
| --- | --- |
| `--model-out` | Output path for learned model JSON (required). |
| `--model-threshold` | Probability cutoff (0–1) to split methylated vs. unmethylated when fitting Betas (default 0.5). |

## BAM Tagging

Tag BAM files directly without intermediate files. The tool reads the SEQ field from each BAM record, applies methylation tagging, and writes MM/ML tags directly:

```bash
methylsim simulate \
  --bam-input input.bam \
  --bam-output output.bam \
  --motif GATC_6mA_1
```

This workflow:
- Reads SEQ field from each BAM record
- Tags sequences with MM/ML based on motifs
- Writes tagged BAM directly (no intermediate FASTQ/TSV files needed)
- Preserves all other BAM fields and tags
>>>>>>> 125d715 (Refactor code structure for improved readability and maintainability)

## Outputs
- **FASTQ** (`--output-fastq`, default `methylsim.fastq`): Sequences with MM/ML tags appended to the header.
- **Tags TSV** (`--tags-tsv`): Columns `read_id`, `MM`, `ML` for reuse or auditing.
- **BAM** (`--bam-output`): Tagged BAM when using `--bam-input` (direct tagging from SEQ field).
- **Logging stats**: Every run prints how many reads contained motif hits plus the average number of motif hits and high-methylation events per read, making it easy to catch mismatched motif definitions.

## Workflow
- Snakemake and Pixi project files now live under `workflow/` so the Rust crate and workflow dependencies stay separate.
- Update `workflow/config.yml` with your samples, SNP rates, and simulation settings. Paths inside the config are resolved relative to the `workflow/` directory (so `resources/...` points at `workflow/resources/...`).
- Run the workflow from that directory, e.g.

```bash
cd workflow
XDG_CACHE_HOME=$PWD/.cache pixi run snakemake --cores 8
```

- The DAG mutates each reference at every SNP rate, simulates reads, discovers motifs with nanomotif, and finishes by running `methylation_phasing split-reads` using the discovered motifs. Outputs still land under the top-level `results/` directory next to the Rust crate.
