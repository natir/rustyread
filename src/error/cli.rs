//! Command line interface error

/* crate use */
use thiserror::Error;

/// Enum to manage error polymorphism
#[derive(Debug, Error)]
pub enum Cli {
    /// quantity didn't match to pattern \d+\[KMGx\]?
    #[error("We aren't able to parse quantity, quantity must match with this regex '\\d+[KMGx]?'")]
    CantParseQuantity,

    /// Cant parse a duo of value
    #[error("We aren't able to parse a value you provide as argument for a parameter")]
    CantParseDuo,

    /// Cant parse a trio of value
    #[error("We aren't able to parse a value you provide as argument for a parameter")]
    CantParseTrio,

    /// Cant found model path
    #[error("Can't found model path use qscore_model and error_model with file")]
    CantFoundModelPath,
}
