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
        .arg("--num-reads")
        .arg("5")
        .arg("--read-length")
        .arg("100")
        .arg("--output-fastq")
        .arg(&output_path)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Command failed: {}", String::from_utf8_lossy(&output.stderr));
    assert!(output_path.exists(), "Output file was not created");

    // Verify output contains FASTQ records
    let content = fs::read_to_string(&output_path).expect("Failed to read output");
    let lines: Vec<&str> = content.lines().collect();
    assert!(lines.len() >= 20, "Expected at least 20 lines (5 reads × 4 lines)");

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
        .arg("--num-reads")
        .arg("3")
        .arg("--read-length")
        .arg("80")
        .arg("--output-fastq")
        .arg(&output_fastq)
        .arg("--tags-tsv")
        .arg(&output_tsv)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Command failed: {}", String::from_utf8_lossy(&output.stderr));
    assert!(output_tsv.exists(), "Tags TSV was not created");

    // Verify TSV format
    let content = fs::read_to_string(&output_tsv).expect("Failed to read TSV");
    let lines: Vec<&str> = content.lines().collect();
    assert!(lines.len() >= 4, "Expected header + at least 3 data lines");
    assert!(lines[0].contains("read_id"), "Header should contain read_id");
    assert!(lines[0].contains("MM"), "Header should contain MM");
    assert!(lines[0].contains("ML"), "Header should contain ML");
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
        .arg("--num-reads")
        .arg("2")
        .arg("--read-length")
        .arg("100")
        .arg("--output-fastq")
        .arg(&output_path)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Command failed: {}", String::from_utf8_lossy(&output.stderr));
    assert!(output_path.exists(), "Output file was not created");

    // Check that output mentions both motifs
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("2 motif(s)"), "Should load 2 motifs from file");
}

#[test]
fn test_simulate_missing_required_args() {
    // Missing motif
    let output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(get_test_data_path("test_reference.fasta"))
        .arg("--num-reads")
        .arg("5")
        .output()
        .expect("Failed to execute command");

    assert!(!output.status.success(), "Should fail without motif");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("motif") || stderr.contains("required"), "Error should mention motif requirement");
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
        .arg("--num-reads")
        .arg("10")
        .arg("--read-length")
        .arg("100")
        .arg("--output-fastq")
        .arg(&tagged_reads)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute simulate command");

    assert!(sim_output.status.success(), "Simulate failed: {}", String::from_utf8_lossy(&sim_output.stderr));

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

    assert!(fit_output.status.success(), "Fit-model failed: {}", String::from_utf8_lossy(&fit_output.stderr));
    assert!(model_output.exists(), "Model JSON was not created");

    // Verify model JSON structure
    let model_content = fs::read_to_string(&model_output).expect("Failed to read model");
    assert!(model_content.contains("version"), "Model should contain version");
    assert!(model_content.contains("threshold"), "Model should contain threshold");
    assert!(model_content.contains("motifs"), "Model should contain motifs");
    assert!(model_content.contains("GATC"), "Model should contain GATC motif");
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
        .arg("--num-reads")
        .arg("10")
        .arg("--read-length")
        .arg("100")
        .arg("--output-fastq")
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
        .arg("--num-reads")
        .arg("5")
        .arg("--read-length")
        .arg("100")
        .arg("--output-fastq")
        .arg(&output_with_model)
        .arg("--seed")
        .arg("123")
        .output()
        .expect("Failed to execute simulate with model");

    assert!(sim2.status.success(), "Simulate with model failed: {}", String::from_utf8_lossy(&sim2.stderr));
    assert!(output_with_model.exists(), "Output with model was not created");

    // Verify output
    let content = fs::read_to_string(&output_with_model).expect("Failed to read output");
    assert!(content.contains("MM:Z:"), "Output should contain MM tags");
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
        .arg("--num-reads")
        .arg("3")
        .arg("--read-length")
        .arg("80")
        .arg("--output-fastq")
        .arg(&output_builtin)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute builtin simulator");

    assert!(builtin.status.success(), "Builtin simulator failed");
    assert!(output_builtin.exists(), "Builtin output not created");

    // Note: We skip badreads test as it requires external dependency
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
        .arg("--num-reads")
        .arg("5")
        .arg("--read-length")
        .arg("100")
        .arg("--substitution-rate")
        .arg("0.05")
        .arg("--insertion-rate")
        .arg("0.02")
        .arg("--deletion-rate")
        .arg("0.02")
        .arg("--output-fastq")
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
        .arg("--num-reads")
        .arg("5")
        .arg("--read-length")
        .arg("100")
        .arg("--motif-high-prob")
        .arg("0.8")
        .arg("--non-motif-high-prob")
        .arg("0.05")
        .arg("--high-ml-mean")
        .arg("220")
        .arg("--high-ml-std")
        .arg("15")
        .arg("--low-ml-mean")
        .arg("30")
        .arg("--low-ml-std")
        .arg("8")
        .arg("--output-fastq")
        .arg(&output_path)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Command with methylation params failed");
    assert!(output_path.exists(), "Output not created");
}

#[test]
fn test_read_length_specifications() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let reference_path = get_test_data_path("test_reference.fasta");

    // Test with mean
    let output_mean = temp_dir.path().join("output_mean.fastq");
    let mean_output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--num-reads")
        .arg("3")
        .arg("--read-length-mean")
        .arg("150")
        .arg("--output-fastq")
        .arg(&output_mean)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute with mean");

    assert!(mean_output.status.success(), "Mean read length failed");
    assert!(output_mean.exists());

    // Test with N50
    let output_n50 = temp_dir.path().join("output_n50.fastq");
    let n50_output = Command::new(get_binary_path())
        .arg("simulate")
        .arg("--reference")
        .arg(&reference_path)
        .arg("--motif")
        .arg("GATC_6mA_1")
        .arg("--num-reads")
        .arg("3")
        .arg("--read-length-n50")
        .arg("120")
        .arg("--output-fastq")
        .arg(&output_n50)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute with N50");

    assert!(n50_output.status.success(), "N50 read length failed");
    assert!(output_n50.exists());
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
        .arg("--output-fastq")
        .arg(&output_path)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute with coverage quantity");

    assert!(output.status.success(), "Coverage quantity failed: {}", String::from_utf8_lossy(&output.stderr));
    assert!(output_path.exists());

    // Verify stderr mentions coverage calculation
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Calculating reads for") && stderr.contains("coverage"),
            "Should mention coverage calculation");
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
        .arg("--output-fastq")
        .arg(&output_path)
        .arg("--seed")
        .arg("42")
        .output()
        .expect("Failed to execute with absolute quantity");

    assert!(output.status.success(), "Absolute quantity failed: {}", String::from_utf8_lossy(&output.stderr));
    assert!(output_path.exists());

    // Should generate 10,000 / 100 = 100 reads
    let content = fs::read_to_string(&output_path).expect("Failed to read output");
    let line_count = content.lines().count();
    // FASTQ has 4 lines per read
    assert!(line_count >= 400 && line_count <= 400, "Should have ~100 reads (400 lines)");
}
