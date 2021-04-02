//! Model error

/* crate use */
use thiserror::Error;

/// Enum to manage error polymorphism
#[derive(Debug, Error)]
pub enum Model {
    /// Error durring error model parsing
    #[error("We aren't able to parse error model")]
    ErrorModelParsing,

    /// Error durring quality model parsing
    #[error("We aren't able to parse error model")]
    QualityModelParsing,

    /// Quality model not contains minimal cigar string
    #[error("Quality model not contains minimal cigar string")]
    QualityModelNotMinimalCigarString,

    /// Cigar string length must be odd
    #[error("Cigar string length must be odd")]
    QualityModelCigarLenNotOdd,
}
