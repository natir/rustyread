//! All stuff relate to error

/* crate use */
use thiserror::Error;

/* module declaration */
pub mod cli;
pub mod model;

/* reexport for easiest use */
pub use cli::Cli;
pub use model::Model;

/// Enum to manage error polymorphism
#[derive(Debug, Error)]
pub enum Error {
    /// Error related to command line interface
    #[error(transparent)]
    Cli(#[from] Cli),

    /// Error related to model
    #[error(transparent)]
    Model(#[from] Model),
}
