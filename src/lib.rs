#![warn(clippy::all)]
#![allow(clippy::too_many_arguments)]

mod dft;
mod eos;
mod parameter;
pub use dft::{GcPcSaftFunctional, GcPcSaftFunctionalParameters};
pub use eos::{GcPcSaft, GcPcSaftEosParameters, GcPcSaftOptions};
pub use parameter::GcPcSaftRecord;

#[cfg(feature = "python")]
pub mod python;
