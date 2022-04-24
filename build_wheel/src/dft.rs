use feos_core::*;
use feos_dft::adsorption::*;
use feos_dft::fundamental_measure_theory::FMTVersion;
use feos_dft::interface::*;
use feos_dft::python::*;
use feos_dft::solvation::*;
use feos_dft::*;
use feos_gc_pcsaft::micelles::*;
use feos_gc_pcsaft::python::*;
use feos_gc_pcsaft::*;
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
/// parameters: GcPcSaftParameters
///     The set of PC-SAFT parameters.
/// fmt_version: FMTVersion, optional
///     The specific variant of the FMT term. Defaults to FMTVersion.WhiteBear
/// max_eta : float, optional
///     Maximum packing fraction. Defaults to 0.5.
/// max_iter_cross_assoc : unsigned integer, optional
///     Maximum number of iterations for cross association. Defaults to 50.
/// tol_cross_assoc : float
///     Tolerance for convergence of cross association. Defaults to 1e-10.
///
/// Returns
/// -------
/// GcPcSaftFunctional
#[pyclass(name = "GcPcSaftFunctional", unsendable)]
#[pyo3(
    text_signature = "(parameters, fmt_version, max_eta, max_iter_cross_assoc, tol_cross_assoc)"
)]
#[derive(Clone)]
pub struct PyGcPcSaftFunctional(pub Rc<DFT<GcPcSaftFunctional>>);

#[pymethods]
impl PyGcPcSaftFunctional {
    #[new]
    #[args(
        fmt_version = "FMTVersion::WhiteBear",
        max_eta = "0.5",
        max_iter_cross_assoc = "50",
        tol_cross_assoc = "1e-10"
    )]
    fn new(
        parameters: PyGcPcSaftFunctionalParameters,
        fmt_version: FMTVersion,
        max_eta: f64,
        max_iter_cross_assoc: usize,
        tol_cross_assoc: f64,
    ) -> Self {
        let options = GcPcSaftOptions {
            max_eta,
            max_iter_cross_assoc,
            tol_cross_assoc,
        };
        Self(Rc::new(GcPcSaftFunctional::with_options(
            parameters.0,
            fmt_version,
            options,
        )))
    }
}

impl_equation_of_state!(PyGcPcSaftFunctional);

impl_state!(DFT<GcPcSaftFunctional>, PyGcPcSaftFunctional);
impl_state_molarweight!(DFT<GcPcSaftFunctional>, PyGcPcSaftFunctional);
impl_phase_equilibrium!(DFT<GcPcSaftFunctional>, PyGcPcSaftFunctional);

impl_planar_interface!(GcPcSaftFunctional);
impl_surface_tension_diagram!(GcPcSaftFunctional);

impl_pore!(GcPcSaftFunctional, PyGcPcSaftFunctional);
impl_adsorption!(GcPcSaftFunctional, PyGcPcSaftFunctional);

impl_solvation_profile!(GcPcSaftFunctional);

impl_micelle_profile!(GcPcSaftFunctional);

#[pymodule]
pub fn dft(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyGcPcSaftFunctional>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagram>()?;
    m.add_class::<PyPhaseEquilibrium>()?;
    m.add_class::<PyPlanarInterface>()?;
    m.add_class::<Geometry>()?;
    m.add_class::<PyPore1D>()?;
    m.add_class::<PyPore3D>()?;
    m.add_class::<PyExternalPotential>()?;
    m.add_class::<PyAdsorption1D>()?;
    m.add_class::<PyAdsorption3D>()?;
    m.add_class::<PySurfaceTensionDiagram>()?;
    m.add_class::<PyDFTSolver>()?;
    m.add_class::<PySolvationProfile>()?;
    m.add_class::<FMTVersion>()?;
    m.add_class::<PyMicelleProfile>()?;
    Ok(())
}
