//! Model error

/* crate use */
use thiserror::Error;

/// Enum to manage error polymorphism
#[derive(Debug, Error)]
pub enum Model {
    /// Error durring error model parsing
    #[error("We aren't able to parse error model")]
    ErrorParsing,

    /// Error durring quality model parsing
    #[error("We aren't able to parse error model")]
    QualityParsing,

    /// Quality model not contains minimal cigar string
    #[error("Quality model not contains minimal cigar string")]
    QualityNotMinimalCigarString,

    /// Cigar string length must be odd
    #[error("Cigar string length must be odd")]
    QualityCigarLenNotOdd,

    /// Length model parameter must be upper than 0.0
    #[error("Length model parameter must be upper than 0.0")]
    LengthParamMustBeUpperThan0,
}
