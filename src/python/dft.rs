use super::parameter::*;
use crate::dft::GcPcSaftFunctional;
use feos_core::python::{PyContributions, PyVerbosity};
use feos_core::utils::{
    DataSet, EquilibriumLiquidDensity, Estimator, LiquidDensity, VaporPressure,
};
use feos_core::*;
use feos_dft::adsorption::*;
use feos_dft::interface::*;
use feos_dft::python::*;
use feos_dft::solvation::*;
use feos_dft::*;
use numpy::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use quantity::python::*;
use quantity::si::*;
use std::collections::HashMap;
use std::rc::Rc;

/// PC-SAFT Helmholtz energy functional.
///
/// Parameters
/// ----------
/// parameters: GcPcSaftFunctionalParameters
///     The set of gc-PC-SAFT parameters.
///
/// Returns
/// -------
/// GcPcSaftFunctional
#[pyclass(name = "GcPcSaftFunctional", unsendable)]
#[pyo3(text_signature = "(parameters)")]
#[derive(Clone)]
pub struct PyGcPcSaftFunctional(pub Rc<DFT<GcPcSaftFunctional>>);

#[pymethods]
impl PyGcPcSaftFunctional {
    #[new]
    fn new(parameters: PyGcPcSaftFunctionalParameters) -> Self {
        Self(Rc::new(GcPcSaftFunctional::new(parameters.0.clone())))
    }
}

impl_equation_of_state!(PyGcPcSaftFunctional);

impl_state!(DFT<GcPcSaftFunctional>, PyGcPcSaftFunctional);
impl_state_molarweight!(DFT<GcPcSaftFunctional>, PyGcPcSaftFunctional);
impl_vle_state!(DFT<GcPcSaftFunctional>, PyGcPcSaftFunctional);

impl_estimator!(DFT<GcPcSaftFunctional>, PyGcPcSaftFunctional);

impl_planar_interface!(GcPcSaftFunctional);
impl_surface_tension_diagram!(GcPcSaftFunctional);

impl_pore!(GcPcSaftFunctional, PyGcPcSaftFunctional);
impl_adsorption!(GcPcSaftFunctional, PyGcPcSaftFunctional);

impl_solvation_profile!(GcPcSaftFunctional);

#[pymodule]
pub fn gc_pcsaft_dft(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyGcPcSaftFunctional>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagramPure>()?;
    m.add_class::<PyPhaseDiagramBinary>()?;
    m.add_class::<PyPhaseDiagramHetero>()?;
    m.add_class::<PyPhaseEquilibrium>()?;
    m.add_class::<PyPlanarInterface>()?;
    m.add_class::<PyGeometry>()?;
    m.add_class::<PyPore1D>()?;
    m.add_class::<PyPore3D>()?;
    m.add_class::<PyExternalPotential>()?;
    m.add_class::<PyAdsorption1D>()?;
    m.add_class::<PyAdsorption3D>()?;
    m.add_class::<PySurfaceTensionDiagram>()?;
    m.add_class::<PyDFTSolver>()?;
    m.add_class::<PySolvationProfile>()?;
    m.add_class::<PyFMTVersion>()?;

    let utils = PyModule::new(py, "utils")?;
    utils.add_class::<PyDataSet>()?;
    utils.add_class::<PyEstimator>()?;
    m.add_submodule(utils)?;
    Ok(())
}
