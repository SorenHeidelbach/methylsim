mod cli;
mod config;
mod fastx;
mod io_utils;
mod model;
mod motif;
mod pipeline;
mod presets;
mod simulator;
mod tagging;

use anyhow::Result;
use clap::Parser;

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp(None)
        .format_module_path(false)
        .init();

    let cli = cli::Cli::parse();
    match cli.command {
        cli::Commands::Simulate(args) => {
            cli::validate_simulate_inputs(&args)?;
            pipeline::run_simulate(args)
        }
        cli::Commands::FitModel(args) => {
            cli::validate_fit_inputs(&args)?;
            pipeline::run_fit_model(args)
        }
    }
}
