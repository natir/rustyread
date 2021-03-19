//! Command line interface error

/* crate use */
use thiserror::Error;

/// Enum to manage error polymorphism
#[derive(Debug, Error)]
pub enum Cli {
    /// quantity didn't match to pattern \d+[KMGx]?
    #[error("We aren't able to parse quantity, quantity must match with this regex '\\d+[KMGx]?'")]
    CantParssQuantity,

    /// Cant parse a duo of value
    #[error("We aren't able to parse a value you provide as argument for a parameter")]
    CantParseDuo,

    /// Cant parse a trio of value
    #[error("We aren't able to parse a value you provide as argument for a parameter")]
    CantParseTrio,
}
