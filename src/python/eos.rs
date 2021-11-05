use super::parameters::PyGcPcSaftParameters;
use crate::eos::{GcPcSaft, GcPcSaftOptions};
use feos_core::python::{PyContributions, PyVerbosity};
use feos_core::EquationOfState;
use feos_core::*;
use numpy::convert::ToPyArray;
use numpy::{PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use quantity::python::*;
use quantity::si::*;
use std::collections::HashMap;
use std::rc::Rc;

#[pyclass(name = "GcPcSaft", unsendable)]
#[pyo3(text_signature = "(parameters, max_eta, max_iter_cross_assoc, tol_cross_assoc)")]
#[derive(Clone)]
pub struct PyGcPcSaft(pub Rc<GcPcSaft>);

#[pymethods]
impl PyGcPcSaft {
    #[new]
    #[args(
        max_eta = "0.5",
        max_iter_cross_assoc = "50",
        tol_cross_assoc = "1e-10"
    )]
    fn new(
        parameters: PyGcPcSaftParameters,
        max_eta: f64,
        max_iter_cross_assoc: usize,
        tol_cross_assoc: f64,
    ) -> Self {
        let options = GcPcSaftOptions {
            max_eta,
            max_iter_cross_assoc,
            tol_cross_assoc,
        };
        Self(Rc::new(GcPcSaft::with_options(
            parameters.0.clone(),
            options,
        )))
    }
}

impl_state!(GcPcSaft, PyGcPcSaft);
impl_state_molarweight!(GcPcSaft, PyGcPcSaft);
impl_vle_state!(GcPcSaft, PyGcPcSaft);

pub fn gc_pcsaft_eos(_: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyGcPcSaft>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagramPure>()?;
    m.add_class::<PyPhaseDiagramBinary>()?;
    m.add_class::<PyPhaseDiagramHetero>()?;
    m.add_class::<PyPhaseEquilibrium>()?;

    Ok(())
}
