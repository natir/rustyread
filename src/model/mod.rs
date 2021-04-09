//! Manage model

/* module declaration */
pub mod error;
pub mod identity;
pub mod length;
pub mod quality;
pub mod adapter;

/* reexport for easiest use */
pub use error::Error;
pub use identity::Identity;
pub use length::Length;
pub use quality::Quality;
pub use adapter::Adapter;
