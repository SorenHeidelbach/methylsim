# methylsim

`methylsim` creates Oxford Nanopore–style read files with MM/ML methylation tags or adds tags to reads that already exist.

## Requirements
- Rust toolchain for building.
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

### 0. Sequence source

Methylsim can use seqeunces from three sources:
* Simulated from reference using basic built-in simulator. Use this if length and error validity of sequences are not important.
* Simulate reads from reference using badread. Methylsim assumes badread is available on $PATH.
* Existing reads. Usefull if you already have a sequence set, but just want to simulate the methylation. Both fastq and BAM files are accepted as input.


`methylsim` has two main subcommands:

### 1. `simulate` - Generate or tag reads with methylation
Simulate reads with methylated at G6mATC motif with 50x reference coverage:

```bash
# Built-in
methylsim simulate \
  --reference ref.fasta \
  --motif GATC_6mA_1 \
  --quantity 50x \
  --out sim_GATC_a_1.fastq
```

```bash
# Badread
methylsim simulate \
  --reference ref.fasta \
  --motif GATC_6mA_1 \
  --out sim_GATC_a_1.fastq \
  --simulator badread \
  --badreads-extra "--quantity 50x"
```



### 2. `fit-model` - Learn methylation model from tagged reads
Extract methylation patterns from real basecalled data:

```bash
methylsim fit-model \
  --reads basecalled_reads.fastq \
  --motif GATC_6mA_1 \
  --model-out models/ecoli_gatc.yaml
```

Then you can use the learned model to simulate realistic reads:

```bash
methylsim simulate \
  --reference ref.fasta \
  --model-in models/ecoli_gatc.yaml \
  --motif GATC_6mA_1 \
  --out realistic_reads.fastq
```


#### Multi-Motif Models:

Learn methylation patterns for multiple motifs simultaneously:

```bash
# Learn model for both Dam (GATC) and Dcm (CCWGG) methylation
methylsim fit-model \
  --reads data/ecoli_tagged.fastq \
  --motif GATC_6mA_1,CCWGG_5mC_1 \
  --model-out models/ecoli_multi.json
```

Then use the multi-motif model for simulation:

```bash
methylsim simulate \
  --reference references/ecoli.fa \
  --model-in models/ecoli_multi.json \
  --quantity 50x \
  --output-fastq results/ecoli_sim.fastq
```

### Bundled model presets

Premade models are located in `models/`. Currently only E. coli Dam/Dcm model is available (`models/e_coli_G6mATC_C5mCWGG.toml`).
Use it with `--model-preset` (motifs are pulled from the preset automatically, so `--motif` is optional):

```bash
methylsim simulate \
  --reference references/ecoli.fa \
  --model-preset ecoli \
  --quantity 25x \
  --output-fastq results/ecoli_preset.fastq
```

### ML tag sampling

Each position is partioned into *High methylation probability* or *Low methylation probability*. Position within a specified motifs are assigned to high methylation probability with a frequency specified in `--motif-high-prob`, whereas positions outside a motif are assigned to high methylation probability with a frequency specified in `--non-motif-high-prob`. For high methylation probability positions, the methylaiton probability is sampled from the corresponding beta. The low methylation probability positions are sampled likevise from a beta with parameters.  The parameters can adjusted with the `--high-beta-alpha`, `--high-beta-beta`, `--low-beta-alpha` and `--low-beta-beta`.
Learned models keep their own Beta parameters specified in the model file.


## Motif definitions and methylation model
- Colon syntax: `motif[:canonical[:offset[:mod]]]`
  - Example: `GATC:A:1:a` = G**6mA**TC.
- Underscore syntax: `sequence_modtype_offset`
  - Example: `GCCG_m_1` G**5mC**CG.

Repeat `--motif` or use comma-separated entries (e.g. `--motif CG,GC_m_1`). Use `--motifs-file` to load motifs: required columns `motif`, `mod_type`, `mod_position` plus optional `motif_complement`/`mod_position_complement`.

Example `motifs.tsv`:

```
motif	mod_position	mod_type
CCWGG	1	m
GCWGC	1	m
```
