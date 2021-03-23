//! Model error

/* crate use */
use thiserror::Error;

/// Enum to manage error polymorphism
#[derive(Debug, Error)]
pub enum Model {
    /// Error durring error model parsing
    #[error("We aren't able to parse error model")]
    ErrorModelParsing,
}
