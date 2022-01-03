#![warn(clippy::all)]
#![allow(clippy::too_many_arguments)]

mod dft;
mod eos;
pub mod micelles;
mod record;
pub use dft::{GcPcSaftFunctional, GcPcSaftFunctionalParameters};
pub use eos::{GcPcSaft, GcPcSaftEosParameters, GcPcSaftOptions};
pub use record::GcPcSaftRecord;

#[cfg(feature = "python")]
pub mod python;
