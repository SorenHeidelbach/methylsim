use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::str::FromStr;

use anyhow::{anyhow, bail, Context, Result};
use noodles::bam::io::writer::Builder as BamWriterBuilder;
use noodles::bam::io::Writer as BamWriter;
use noodles::bgzf;
use noodles::sam::{
    self,
    alignment::io::Write as AlignmentWrite,
    alignment::{self, record::data::field::Tag, record_buf},
    io::reader::Builder as SamReaderBuilder,
    io::Writer as SamWriter,
};

use crate::cli::OutputFormat;
use crate::fastx::read_fastx;
use crate::fastx::ReadRecord;
use crate::tagging::TagResult;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum InputFormat {
    Fastx,
    Sam,
    Bam,
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

pub struct OutputWriter {
    target: OutputTarget,
}

enum OutputTarget {
    Fastq(BufWriter<File>),
    Sam(SamOutput),
    Bam(BamOutput),
}

impl OutputWriter {
    pub fn new(path: &Path, format: OutputFormat) -> Result<Self> {
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

    pub fn write(&mut self, read: &ReadRecord, tags: &TagResult) -> Result<()> {
        match &mut self.target {
            OutputTarget::Fastq(writer) => write_fastq_record(writer, read, tags),
            OutputTarget::Sam(writer) => writer.write_record(read, tags),
            OutputTarget::Bam(writer) => writer.write_record(read, tags),
        }
    }

    pub fn finish(self) -> Result<()> {
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
    let Some(list) = raw.strip_prefix("ML:B:C,").or_else(|| raw.strip_prefix("ML:B:")) else {
        return Err(anyhow!("Invalid ML tag: {}", raw));
    };

    if list.trim().is_empty() {
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

pub fn write_fastq_record(
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

pub fn ensure_output_extension(path: &Path, format: OutputFormat) -> Result<()> {
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

pub fn read_input_reads(path: &Path) -> Result<Vec<ReadRecord>> {
    match InputFormat::from_path(path)? {
        InputFormat::Fastx => read_fastx(path),
        InputFormat::Sam => read_sam_reads(path),
        InputFormat::Bam => read_bam_reads(path),
    }
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
        .as_ref()
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
    let mut reader = noodles::bam::io::reader::Builder::default()
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
