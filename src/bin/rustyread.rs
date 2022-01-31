/* crate use */
use anyhow::Result;

use clap::Parser as _;

/* local use */
use rustyread::*;

#[cfg(not(tarpaulin_include))]
fn main() -> Result<()> {
    let params = cli::Command::parse();

    if let Some(level) = cli::i82level(params.verbosity) {
        env_logger::builder()
            .format_timestamp(Some(env_logger::fmt::TimestampPrecision::Millis))
            .filter_level(level.to_level_filter())
            .init();
    } else {
        env_logger::Builder::from_env("BADREAD_LOG")
            .format_timestamp(Some(env_logger::fmt::TimestampPrecision::Millis))
            .init();
    }

    if let Some(threads) = params.threads {
        log::info!("Set number of threads to {}", threads);

        cli::set_nb_threads(threads);
    }

    match params.subcmd {
        cli::SubCommand::Simulate(sub) => simulate::simulate(sub),
    }
}
