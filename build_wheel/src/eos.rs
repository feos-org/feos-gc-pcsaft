use feos_core::*;
use feos_gc_pcsaft::python::PyGcPcSaftEosParameters;
use feos_gc_pcsaft::{GcPcSaft, GcPcSaftOptions};
use numpy::convert::ToPyArray;
use numpy::{PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use quantity::python::*;
use quantity::si::*;
use std::collections::HashMap;
use std::rc::Rc;

/// Initialize PC-SAFT equation of state.
///
/// Parameters
/// ----------
/// parameters : GcPcSaftParameters
///     The parameters of the gc-PC-Saft equation of state to use.
/// max_eta : float, optional
///     Maximum packing fraction. Defaults to 0.5.
/// max_iter_cross_assoc : unsigned integer, optional
///     Maximum number of iterations for cross association. Defaults to 50.
/// tol_cross_assoc : float
///     Tolerance for convergence of cross association. Defaults to 1e-10.
///
/// Returns
/// -------
/// GcPcSaft
///     The gc-PC-SAFT equation of state that can be used to compute thermodynamic
///     states.
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
        parameters: PyGcPcSaftEosParameters,
        max_eta: f64,
        max_iter_cross_assoc: usize,
        tol_cross_assoc: f64,
    ) -> Self {
        let options = GcPcSaftOptions {
            max_eta,
            max_iter_cross_assoc,
            tol_cross_assoc,
        };
        Self(Rc::new(GcPcSaft::with_options(parameters.0, options)))
    }
}

impl_equation_of_state!(PyGcPcSaft);
impl_virial_coefficients!(PyGcPcSaft);

impl_state!(GcPcSaft, PyGcPcSaft);
impl_state_molarweight!(GcPcSaft, PyGcPcSaft);
impl_phase_equilibrium!(GcPcSaft, PyGcPcSaft);

#[pymodule]
pub fn eos(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyGcPcSaft>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagram>()?;
    m.add_class::<PyPhaseEquilibrium>()?;
    Ok(())
}
