//! All stuff relate to error

/* crate use */
use thiserror::Error;

/* module declaration */
pub mod cli;

/* reexport for easiest use */
pub use cli::Cli;

/// Enum to manage error polymorphism
#[derive(Debug, Error)]
pub enum Error {
    /// Error related to command line interface
    #[error(transparent)]
    Cli(#[from] Cli),
}
