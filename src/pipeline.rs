use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{anyhow, Context, Result};
use log::info;

use crate::cli::{FitModelArgs, SimulateArgs};
use crate::config::{
    build_simulator_strategy, build_tagger_from_args, collect_motif_specs, load_model_source,
};
use crate::fastx::ReadRecord;
use crate::io_utils::{ensure_output_extension, read_input_reads, OutputWriter};
use crate::model::learn_model_from_reads;
use crate::motif::MotifDefinition;
use crate::simulator::SimulatorStrategy;
use crate::tagging::{MethylationTagger, WriteStats};

pub fn run_simulate(args: SimulateArgs) -> Result<()> {
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

pub fn run_fit_model(args: FitModelArgs) -> Result<()> {
    info!("Starting model fitting");

    let motif_specs = collect_motif_specs(&args.motifs, args.motifs_file.as_ref(), None)?;
    let motifs = motif_specs
        .iter()
        .map(|spec| MotifDefinition::parse(spec))
        .collect::<Result<Vec<_>>>()?;

    eprintln!(
        "Learning methylation patterns from reads in: {}",
        args.reads.display()
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
