//! Manage model

/* module declaration */
pub mod error;
pub mod identity;
pub mod length;
pub mod quality;

/* reexport for easiest use */
pub use error::Error;
pub use identity::Identity;
pub use length::Length;
pub use quality::Quality;
