use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

/// Helper to get the path to the compiled binary
fn get_binary_path() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push("debug");
    path.push("methylsim");
    path
}

/// Helper to get test data path
fn get_test_data_path(filename: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests");
    path.push("data");
    path.push(filename);
    path
}

#[test]
fn test_cli_help() {
    let output = Command::new(get_binary_path())
        .arg("--help")
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("simulate"));
    assert!(stdout.contains("fit-model"));
}

#[test]
fn test_simulate_help() {
    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--help")
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Input Mode"));
    assert!(stdout.contains("Motif Definitions"));
    assert!(stdout.contains("Simulation Options"));
}

#[test]
fn test_fit_model_help() {
    let output = Command::new(get_binary_path())
        .arg("fit-model")
        .arg("--help")
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Required"));
    assert!(stdout.contains("Model Fitting Parameters"));
    assert!(stdout.contains("n-reads"));
}

#[test]
fn test_simulate_basic() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_path = temp_dir.path().join("output.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("5")
        .arg("--read-length")
        .arg("100")
        .arg("--out")
        .arg(&output_path)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute command");

    assert!(
        output.status.success(),
        "Command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output_path.exists(), "Output file was not created");

    // Verify output contains FASTQ records
    let content = fs::read_to_string(&output_path).expect("Failed to read output");
    let lines: Vec<&str> = content.lines().collect();
    assert!(
        lines.len() >= 20,
        "Expected at least 20 lines (5 reads × 4 lines)"
    );

    // Check FASTQ format
    assert!(lines[0].starts_with('@'), "First line should be header");
    assert!(lines[2].starts_with('+'), "Third line should be +");
}

#[test]
fn test_simulate_with_tags_tsv() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_fastq = temp_dir.path().join("output.fastq");
    let output_tsv = temp_dir.path().join("tags.tsv");
    let reference_path = get_test_data_path("test_reference.fasta");

    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("3")
        .arg("--read-length")
        .arg("80")
        .arg("--out")
        .arg(&output_fastq)
        .arg("--tags-tsv")
        .arg(&output_tsv)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute command");

    assert!(
        output.status.success(),
        "Command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output_tsv.exists(), "Tags TSV was not created");

    // Verify TSV format
    let content = fs::read_to_string(&output_tsv).expect("Failed to read TSV");
    let lines: Vec<&str> = content.lines().collect();
    assert!(lines.len() >= 4, "Expected header + at least 3 data lines");
    assert!(
        lines[0].contains("read_id"),
        "Header should contain read_id"
    );
    assert!(lines[0].contains("MM"), "Header should contain MM");
    assert!(lines[0].contains("ML"), "Header should contain ML");
}

#[test]
fn test_simulate_tags_tsv_contents() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_fastq = temp_dir.path().join("output.fastq");
    let output_tsv = temp_dir.path().join("tags.tsv");
    let reference_path = get_test_data_path("test_reference.fasta");

    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("4")
        .arg("--read-length")
        .arg("60")
        .arg("--out")
        .arg(&output_fastq)
        .arg("--tags-tsv")
        .arg(&output_tsv)
        .arg("--seed")
        .arg("5")
        .output()
        .expect("Failed to execute command");

    assert!(
        output.status.success(),
        "Command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output_tsv.exists(), "Tags TSV was not created");

    let tsv = fs::read_to_string(&output_tsv).expect("Failed to read TSV");
    let lines: Vec<&str> = tsv.lines().collect();
    assert_eq!(lines.len(), 5, "TSV lines should match read count + header");
    assert!(
        lines.iter().skip(1).all(|line| line.contains('\t')),
        "Data lines should be tab-delimited"
    );
}

#[test]
fn test_simulate_with_motif_file() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_path = temp_dir.path().join("output.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");
    let motifs_path = get_test_data_path("test_motifs.txt");

    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motifs-file")
        .arg(&motifs_path)
        .arg("--quantity")
        .arg("2")
        .arg("--read-length")
        .arg("100")
        .arg("--out")
        .arg(&output_path)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute command");

    assert!(
        output.status.success(),
        "Command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output_path.exists(), "Output file was not created");

    // Check that output mentions both motifs
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("2 motif(s)"),
        "Should load 2 motifs from file"
    );
}

#[test]
fn test_simulate_with_existing_reads_fastq() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_path = temp_dir.path().join("tagged.fastq");
    let reads_path = get_test_data_path("test_reads.fastq");

    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reads")
        .arg(&reads_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--out")
        .arg(&output_path)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute command");

    assert!(
        output.status.success(),
        "Command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output_path.exists(), "Output file was not created");

    let content = fs::read_to_string(&output_path).expect("Failed to read output");
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len(), 8, "Expected two FASTQ records");
    assert!(
        lines[0].contains("MM:Z:"),
        "Tagged FASTQ headers should include MM tags"
    );
}

#[test]
fn test_simulate_with_existing_reads_sam_input() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_path = temp_dir.path().join("tagged_from_sam.fastq");
    let reads_path = get_test_data_path("test_reads.sam");

    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reads")
        .arg(&reads_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--out")
        .arg(&output_path)
        .arg("--seed")
        .arg("7")
        .output()
        .expect("Failed to execute command");

    assert!(
        output.status.success(),
        "Command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output_path.exists(), "Output file was not created");

    let content = fs::read_to_string(&output_path).expect("Failed to read output");
    let headers: Vec<&str> = content.lines().take(2).collect();
    assert!(
        headers.iter().any(|h| h.contains("MM:Z:")),
        "Tagged output should include MM tags from SAM input"
    );
}

#[test]
fn test_simulate_with_existing_reads_bam_input() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let initial_bam = temp_dir.path().join("initial.bam");
    let final_fastq = temp_dir.path().join("from_bam.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    // First generate BAM output from simulation
    let sim_to_bam = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("3")
        .arg("--read-length")
        .arg("80")
        .arg("--out")
        .arg(&initial_bam)
        .arg("--out-format")
        .arg("bam")
        .arg("--seed")
        .arg("11")
        .output()
        .expect("Failed to execute simulate to BAM");

    assert!(
        sim_to_bam.status.success(),
        "Simulate to BAM failed: {}",
        String::from_utf8_lossy(&sim_to_bam.stderr)
    );
    assert!(initial_bam.exists(), "Initial BAM output missing");

    // Now read BAM as input and tag again (should still succeed)
    let sim_from_bam = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reads")
        .arg(&initial_bam)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--out")
        .arg(&final_fastq)
        .arg("--seed")
        .arg("22")
        .output()
        .expect("Failed to execute simulate from BAM");

    assert!(
        sim_from_bam.status.success(),
        "Simulate from BAM failed: {}",
        String::from_utf8_lossy(&sim_from_bam.stderr)
    );
    assert!(final_fastq.exists(), "Final FASTQ not created");

    let content = fs::read_to_string(&final_fastq).expect("Failed to read final FASTQ");
    assert!(
        content.contains("MM:Z:"),
        "BAM input path should still produce MM tags"
    );
}

#[test]
fn test_simulate_missing_required_args() {
    // Missing motif
    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(get_test_data_path("test_reference.fasta"))
        .arg("--quantity")
        .arg("5")
        .output()
        .expect("Failed to execute command");

    assert!(!output.status.success(), "Should fail without motif");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("motif") || stderr.contains("required"),
        "Error should mention motif requirement"
    );
}

#[test]
fn test_fit_model_basic() {
    // First generate some tagged reads
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let tagged_reads = temp_dir.path().join("tagged.fastq");
    let model_output = temp_dir.path().join("model.json");
    let reference_path = get_test_data_path("test_reference.fasta");

    // Generate tagged reads
    let sim_output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("10")
        .arg("--read-length")
        .arg("100")
        .arg("--out")
        .arg(&tagged_reads)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute simulate command");

    assert!(
        sim_output.status.success(),
        "Simulate failed: {}",
        String::from_utf8_lossy(&sim_output.stderr)
    );

    // Fit model from tagged reads
    let fit_output = Command::new(get_binary_path())
        .arg("fit-model")
        .arg("--reads")
        .arg(&tagged_reads)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--model-out")
        .arg(&model_output)
        .output()
        .expect("Failed to execute fit-model command");

    assert!(
        fit_output.status.success(),
        "Fit-model failed: {}",
        String::from_utf8_lossy(&fit_output.stderr)
    );
    assert!(model_output.exists(), "Model JSON was not created");

    // Verify model JSON structure
    let model_content = fs::read_to_string(&model_output).expect("Failed to read model");
    assert!(
        model_content.contains("version"),
        "Model should contain version"
    );
    assert!(
        model_content.contains("threshold"),
        "Model should contain threshold"
    );
    assert!(
        model_content.contains("motifs"),
        "Model should contain motifs"
    );
    assert!(
        model_content.contains("GATC"),
        "Model should contain GATC motif"
    );
}

#[test]
fn test_fit_model_with_n_reads_limit() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let tagged_reads = temp_dir.path().join("tagged.fastq");
    let model_output = temp_dir.path().join("model_limited.json");
    let reference_path = get_test_data_path("test_reference.fasta");

    // Generate a handful of tagged reads
    let sim_output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("5")
        .arg("--read-length")
        .arg("100")
        .arg("--out")
        .arg(&tagged_reads)
        .arg("--seed")
        .arg("7")
        .output()
        .expect("Failed to execute simulate command");

    assert!(
        sim_output.status.success(),
        "Simulate failed: {}",
        String::from_utf8_lossy(&sim_output.stderr)
    );

    // Fit using only a subset of reads
    let fit_output = Command::new(get_binary_path())
        .arg("fit-model")
        .arg("--reads")
        .arg(&tagged_reads)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--model-out")
        .arg(&model_output)
        .arg("--n-reads")
        .arg("2")
        .output()
        .expect("Failed to execute fit-model command");

    assert!(
        fit_output.status.success(),
        "Fit-model failed: {}",
        String::from_utf8_lossy(&fit_output.stderr)
    );
    assert!(model_output.exists(), "Model JSON was not created");
}

#[test]
fn test_fit_model_missing_required_args() {
    // Missing reads
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let model_output = temp_dir.path().join("model.json");

    let output = Command::new(get_binary_path())
        .arg("fit-model")
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--model-out")
        .arg(&model_output)
        .output()
        .expect("Failed to execute command");

    assert!(!output.status.success(), "Should fail without reads");
}

#[test]
fn test_fit_model_with_motifs_file() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let tagged_reads = temp_dir.path().join("tagged.fastq");
    let model_output = temp_dir.path().join("model_file.json");
    let reference_path = get_test_data_path("test_reference.fasta");
    let motifs_path = get_test_data_path("test_motifs.txt");

    let sim_output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motifs-file")
        .arg(&motifs_path)
        .arg("--quantity")
        .arg("8")
        .arg("--read-length")
        .arg("100")
        .arg("--out")
        .arg(&tagged_reads)
        .arg("--seed")
        .arg("404")
        .output()
        .expect("Failed to simulate with motifs file");
    assert!(sim_output.status.success());

    let fit_output = Command::new(get_binary_path())
        .arg("fit-model")
        .arg("--reads")
        .arg(&tagged_reads)
        .arg("--motifs-file")
        .arg(&motifs_path)
        .arg("--model-out")
        .arg(&model_output)
        .output()
        .expect("Failed to fit model with motifs file");

    assert!(
        fit_output.status.success(),
        "Fit-model with motifs file failed: {}",
        String::from_utf8_lossy(&fit_output.stderr)
    );
    assert!(model_output.exists(), "Model JSON was not created");
}

#[test]
fn test_fit_model_invalid_n_reads() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let model_output = temp_dir.path().join("model_invalid.json");

    let output = Command::new(get_binary_path())
        .arg("fit-model")
        .arg("--reads")
        .arg(get_test_data_path("test_reads.fastq"))
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--model-out")
        .arg(&model_output)
        .arg("--n-reads")
        .arg("0")
        .output()
        .expect("Failed to execute fit-model with invalid n-reads");

    assert!(
        !output.status.success(),
        "fit-model should fail with n-reads=0"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("n-reads") || stderr.contains("greater than zero"),
        "Error should mention invalid n-reads"
    );
}

#[test]
fn test_fit_model_with_custom_threshold() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let tagged_reads = temp_dir.path().join("tagged.fastq");
    let model_output = temp_dir.path().join("model_threshold.json");
    let reference_path = get_test_data_path("test_reference.fasta");

    // Generate tagged reads
    let sim_output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("6")
        .arg("--read-length")
        .arg("90")
        .arg("--out")
        .arg(&tagged_reads)
        .arg("--seed")
        .arg("808")
        .output()
        .expect("Failed to simulate for threshold test");
    assert!(sim_output.status.success());

    // Fit model with non-default threshold
    let fit_output = Command::new(get_binary_path())
        .arg("fit-model")
        .arg("--reads")
        .arg(&tagged_reads)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--model-out")
        .arg(&model_output)
        .arg("--model-threshold")
        .arg("0.2")
        .output()
        .expect("Failed to fit model with custom threshold");

    assert!(
        fit_output.status.success(),
        "Fit-model with custom threshold failed: {}",
        String::from_utf8_lossy(&fit_output.stderr)
    );
    let model_json = fs::read_to_string(&model_output).expect("Failed to read model");
    assert!(
        model_json.contains("\"threshold\": 0.2"),
        "Model JSON should persist custom threshold"
    );
}

#[test]
fn test_simulate_with_learned_model() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let tagged_reads = temp_dir.path().join("tagged.fastq");
    let model_file = temp_dir.path().join("model.json");
    let output_with_model = temp_dir.path().join("output_with_model.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    // Step 1: Generate tagged reads
    let sim1 = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("10")
        .arg("--read-length")
        .arg("100")
        .arg("--out")
        .arg(&tagged_reads)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute simulate");

    assert!(sim1.status.success(), "First simulate failed");

    // Step 2: Fit model
    let fit = Command::new(get_binary_path())
        .arg("fit-model")
        .arg("--reads")
        .arg(&tagged_reads)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--model-out")
        .arg(&model_file)
        .output()
        .expect("Failed to execute fit-model");

    assert!(fit.status.success(), "Fit-model failed");

    // Step 3: Simulate with learned model
    let sim2 = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--model-in")
        .arg(&model_file)
        .arg("--quantity")
        .arg("5")
        .arg("--read-length")
        .arg("100")
        .arg("--out")
        .arg(&output_with_model)
        .arg("--seed")
        .arg("123")
        .output()
        .expect("Failed to execute simulate with model");

    assert!(
        sim2.status.success(),
        "Simulate with model failed: {}",
        String::from_utf8_lossy(&sim2.stderr)
    );
    assert!(
        output_with_model.exists(),
        "Output with model was not created"
    );

    // Verify output
    let content = fs::read_to_string(&output_with_model).expect("Failed to read output");
    assert!(content.contains("MM:Z:"), "Output should contain MM tags");
}

#[test]
fn test_simulate_with_model_preset() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_path = temp_dir.path().join("preset.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--model-preset")
        .arg("ecoli")
        .arg("--quantity")
        .arg("3")
        .arg("--read-length")
        .arg("80")
        .arg("--out")
        .arg(&output_path)
        .arg("--seed")
        .arg("7")
        .output()
        .expect("Failed to execute simulate with preset");

    assert!(
        output.status.success(),
        "Simulate with preset failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output_path.exists(), "Output with preset was not created");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Extracted 2 motif(s) from model"),
        "Should pull motifs from preset"
    );

    let content = fs::read_to_string(&output_path).expect("Failed to read preset output");
    assert!(
        content.contains("MM:Z:"),
        "Preset simulation should contain MM tags"
    );
}

#[test]
fn test_simulate_with_model_only_motifs() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let tagged_reads = temp_dir.path().join("tagged.fastq");
    let model_file = temp_dir.path().join("model.json");
    let final_output = temp_dir.path().join("output.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    // Generate tagged reads
    let sim_output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("6")
        .arg("--read-length")
        .arg("80")
        .arg("--out")
        .arg(&tagged_reads)
        .arg("--seed")
        .arg("123")
        .output()
        .expect("Failed to execute simulate");

    assert!(
        sim_output.status.success(),
        "Simulate failed: {}",
        String::from_utf8_lossy(&sim_output.stderr)
    );

    // Fit model from tagged reads
    let fit_output = Command::new(get_binary_path())
        .arg("fit-model")
        .arg("--reads")
        .arg(&tagged_reads)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--model-out")
        .arg(&model_file)
        .output()
        .expect("Failed to execute fit-model");

    assert!(
        fit_output.status.success(),
        "Fit-model failed: {}",
        String::from_utf8_lossy(&fit_output.stderr)
    );

    // Simulate using only the model (motifs should be inferred from model file)
    let sim_with_model = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--model-in")
        .arg(&model_file)
        .arg("--quantity")
        .arg("3")
        .arg("--read-length")
        .arg("80")
        .arg("--out")
        .arg(&final_output)
        .arg("--seed")
        .arg("999")
        .output()
        .expect("Failed to execute simulate with model-only motifs");

    assert!(
        sim_with_model.status.success(),
        "Simulate with model-only motifs failed: {}",
        String::from_utf8_lossy(&sim_with_model.stderr)
    );
    assert!(final_output.exists(), "Final output was not created");

    let content = fs::read_to_string(&final_output).expect("Failed to read output");
    assert!(
        content.contains("MM:Z:"),
        "Simulation driven by model should emit MM tags"
    );
}

#[test]
fn test_simulate_with_model_and_missing_motif_fails() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let tagged_reads = temp_dir.path().join("tagged.fastq");
    let model_file = temp_dir.path().join("model.json");
    let reference_path = get_test_data_path("test_reference.fasta");

    // Build model for a single motif
    let sim_out = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("5")
        .arg("--read-length")
        .arg("80")
        .arg("--out")
        .arg(&tagged_reads)
        .arg("--seed")
        .arg("101")
        .output()
        .expect("Failed to simulate");
    assert!(sim_out.status.success());

    let fit_out = Command::new(get_binary_path())
        .arg("fit-model")
        .arg("--reads")
        .arg(&tagged_reads)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--model-out")
        .arg(&model_file)
        .output()
        .expect("Failed to fit model");
    assert!(fit_out.status.success());

    // Attempt to simulate with a motif not present in the model should fail
    let sim_fail = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("CCWGG_5mC_1")
        .arg("--model-in")
        .arg(&model_file)
        .arg("--quantity")
        .arg("2")
        .arg("--read-length")
        .arg("50")
        .arg("--out")
        .arg(temp_dir.path().join("should_not_exist.fastq"))
        .output()
        .expect("Failed to execute simulate with mismatched model");

    assert!(
        !sim_fail.status.success(),
        "Simulate should fail when model is missing requested motif"
    );
    let stderr = String::from_utf8_lossy(&sim_fail.stderr);
    assert!(
        stderr.contains("not present") || stderr.contains("missing"),
        "Error should mention missing motif in model"
    );
}

#[test]
fn test_simulate_invalid_quantity_and_read_length() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_path = temp_dir.path().join("invalid.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    // Invalid quantity string
    let bad_quantity = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("abc")
        .arg("--read-length")
        .arg("50")
        .arg("--out")
        .arg(&output_path)
        .output()
        .expect("Failed to execute simulate with bad quantity");
    assert!(
        !bad_quantity.status.success(),
        "Simulate should fail with invalid quantity"
    );
    let stderr = String::from_utf8_lossy(&bad_quantity.stderr);
    assert!(
        stderr.contains("quantity") || stderr.contains("Invalid"),
        "Error should mention quantity parsing"
    );

    // Invalid read length (zero)
    let bad_read_length = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("1")
        .arg("--read-length")
        .arg("0")
        .arg("--out")
        .arg(&output_path)
        .output()
        .expect("Failed to execute simulate with bad read length");
    assert!(
        !bad_read_length.status.success(),
        "Simulate should fail with zero read length"
    );
    let stderr = String::from_utf8_lossy(&bad_read_length.stderr);
    assert!(
        stderr.contains("read-length") || stderr.contains("greater than zero"),
        "Error should mention read-length validation"
    );
}

#[test]
fn test_different_simulators() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_builtin = temp_dir.path().join("builtin.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    // Test builtin simulator
    let builtin = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--simulator")
        .arg("builtin")
        .arg("--quantity")
        .arg("3")
        .arg("--read-length")
        .arg("80")
        .arg("--out")
        .arg(&output_builtin)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute builtin simulator");

    assert!(builtin.status.success(), "Builtin simulator failed");
    assert!(output_builtin.exists(), "Builtin output not created");

    // Note: We skip badreads execution because it requires external dependency
}

#[test]
fn test_simulate_badreads_missing_exec() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_path = temp_dir.path().join("badreads.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--simulator")
        .arg("badreads")
        .arg("--badreads-exec")
        .arg("nonexistent_badread_exec")
        .arg("--quantity")
        .arg("1")
        .arg("--read-length")
        .arg("50")
        .arg("--out")
        .arg(&output_path)
        .arg("--seed")
        .arg("3")
        .output()
        .expect("Failed to execute badreads simulate");

    assert!(
        !output.status.success(),
        "Simulate should fail when badread executable is missing"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("badread") || stderr.contains("executable"),
        "Error should mention badread invocation"
    );
}

#[test]
fn test_error_rates() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_path = temp_dir.path().join("output.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("5")
        .arg("--read-length")
        .arg("100")
        .arg("--substitution-rate")
        .arg("0.05")
        .arg("--insertion-rate")
        .arg("0.02")
        .arg("--deletion-rate")
        .arg("0.02")
        .arg("--out")
        .arg(&output_path)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Command with error rates failed");
    assert!(output_path.exists(), "Output not created");
}

#[test]
fn test_methylation_probabilities() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_path = temp_dir.path().join("output.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("5")
        .arg("--read-length")
        .arg("100")
        .arg("--out")
        .arg(&output_path)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute command");

    assert!(
        output.status.success(),
        "Command with methylation params failed"
    );
    assert!(output_path.exists(), "Output not created");
}

#[test]
fn test_read_length_specifications() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let reference_path = get_test_data_path("test_reference.fasta");

    let output = temp_dir.path().join("output_len.fastq");
    let run = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("3")
        .arg("--read-length")
        .arg("150")
        .arg("--out")
        .arg(&output)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute read length test");

    assert!(run.status.success(), "Read length run failed");
    assert!(output.exists());
}

#[test]
fn test_quantity_coverage() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_path = temp_dir.path().join("output.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("5x")
        .arg("--read-length")
        .arg("100")
        .arg("--out")
        .arg(&output_path)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute with coverage quantity");

    assert!(
        output.status.success(),
        "Coverage quantity failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output_path.exists());

    // Verify stderr mentions coverage calculation
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Coverage-based quantity") && stderr.contains("coverage"),
        "Should mention coverage calculation"
    );
}

#[test]
fn test_quantity_absolute() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let output_path = temp_dir.path().join("output.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--quantity")
        .arg("10K")
        .arg("--read-length")
        .arg("100")
        .arg("--out")
        .arg(&output_path)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute with absolute quantity");

    assert!(
        output.status.success(),
        "Absolute quantity failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output_path.exists());

    // Should generate 10,000 / 100 = 100 reads
    let content = fs::read_to_string(&output_path).expect("Failed to read output");
    let line_count = content.lines().count();
    // FASTQ has 4 lines per read
    assert!(
        line_count >= 400 && line_count <= 400,
        "Should have ~100 reads (400 lines)"
    );
}

#[test]
fn test_multi_motif_model_fitting() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let tagged_reads = temp_dir.path().join("tagged.fastq");
    let model_output = temp_dir.path().join("model.json");
    let reference_path = get_test_data_path("test_reference.fasta");

    // Generate reads with multiple motifs
    let sim_output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1,CCWGG_5mC_1")
        .arg("--quantity")
        .arg("20")
        .arg("--read-length")
        .arg("150")
        .arg("--out")
        .arg(&tagged_reads)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to simulate with multiple motifs");

    assert!(
        sim_output.status.success(),
        "Simulate with multiple motifs failed"
    );

    // Fit model from reads with multiple motifs
    let fit_output = Command::new(get_binary_path())
        .arg("fit-model")
        .arg("--reads")
        .arg(&tagged_reads)
        .arg("--motif")
        .arg("GATC_6mA_1,CCWGG_5mC_1")
        .arg("--model-out")
        .arg(&model_output)
        .output()
        .expect("Failed to fit model with multiple motifs");

    assert!(
        fit_output.status.success(),
        "Fit model with multiple motifs failed: {}",
        String::from_utf8_lossy(&fit_output.stderr)
    );
    assert!(model_output.exists(), "Model JSON was not created");

    // Verify model contains both motifs
    let model_content = fs::read_to_string(&model_output).expect("Failed to read model");
    assert!(
        model_content.contains("GATC"),
        "Model should contain GATC motif"
    );
    assert!(
        model_content.contains("CCWGG"),
        "Model should contain CCWGG motif"
    );

    // Parse JSON and verify structure
    let model: serde_json::Value =
        serde_json::from_str(&model_content).expect("Failed to parse model JSON");

    let motifs = model["motifs"]
        .as_array()
        .expect("Model should have motifs array");
    assert_eq!(motifs.len(), 2, "Model should contain exactly 2 motifs");

    // Verify stderr shows per-motif statistics
    let stderr = String::from_utf8_lossy(&fit_output.stderr);
    assert!(
        stderr.contains("Fitting model for 2 motif(s)"),
        "Should mention 2 motifs"
    );
    assert!(
        stderr.contains("GATC") && stderr.contains("CCWGG"),
        "Should show statistics for both motifs"
    );
}

#[test]
fn test_simulate_with_multi_motif_model() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let tagged_reads = temp_dir.path().join("tagged.fastq");
    let model_file = temp_dir.path().join("model.json");
    let output_with_model = temp_dir.path().join("output.fastq");
    let reference_path = get_test_data_path("test_reference.fasta");

    // Generate reads with multiple motifs
    let sim1 = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1,CCWGG_5mC_1")
        .arg("--quantity")
        .arg("15")
        .arg("--read-length")
        .arg("120")
        .arg("--out")
        .arg(&tagged_reads)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to simulate");

    assert!(sim1.status.success());

    // Fit multi-motif model
    let fit = Command::new(get_binary_path())
        .arg("fit-model")
        .arg("--reads")
        .arg(&tagged_reads)
        .arg("--motif")
        .arg("GATC_6mA_1,CCWGG_5mC_1")
        .arg("--model-out")
        .arg(&model_file)
        .output()
        .expect("Failed to fit model");

    assert!(fit.status.success());

    // Simulate with multi-motif model
    let sim2 = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1,CCWGG_5mC_1")
        .arg("--model-in")
        .arg(&model_file)
        .arg("--quantity")
        .arg("5")
        .arg("--read-length")
        .arg("100")
        .arg("--out")
        .arg(&output_with_model)
        .arg("--seed")
        .arg("123")
        .output()
        .expect("Failed to simulate with multi-motif model");

    assert!(
        sim2.status.success(),
        "Simulate with multi-motif model failed: {}",
        String::from_utf8_lossy(&sim2.stderr)
    );
    assert!(output_with_model.exists());

    // Verify output contains tags for both motif types
    let content = fs::read_to_string(&output_with_model).expect("Failed to read output");
    assert!(content.contains("MM:Z:"), "Output should contain MM tags");
}
