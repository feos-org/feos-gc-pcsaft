use crate::dft::GcPcSaftFunctionalParameters;
use crate::eos::GcPcSaftEosParameters;
use crate::record::GcPcSaftRecord;
use feos_core::joback::JobackRecord;
use feos_core::parameter::{BinaryRecord, IdentifierOption, ParameterError, SegmentRecord};
use feos_core::python::joback::PyJobackRecord;
use feos_core::python::parameter::{PyBinarySegmentRecord, PyChemicalRecord};
use feos_core::{impl_json_handling, impl_parameter_from_segments, impl_segment_record};
use numpy::{PyArray2, ToPyArray};
use pyo3::prelude::*;
use std::convert::TryFrom;
use std::rc::Rc;

#[pyclass(name = "GcPcSaftRecord", unsendable)]
#[pyo3(
    text_signature = "(m, sigma, epsilon_k, mu=None, q=None, kappa_ab=None, epsilon_k_ab=None, na=None, nb=None)"
)]
#[derive(Clone)]
pub struct PyGcPcSaftRecord(GcPcSaftRecord);

#[pymethods]
impl PyGcPcSaftRecord {
    #[new]
    fn new(
        m: f64,
        sigma: f64,
        epsilon_k: f64,
        mu: Option<f64>,
        kappa_ab: Option<f64>,
        epsilon_k_ab: Option<f64>,
        na: Option<f64>,
        nb: Option<f64>,
        psi_dft: Option<f64>,
    ) -> Self {
        Self(GcPcSaftRecord::new(
            m,
            sigma,
            epsilon_k,
            mu,
            kappa_ab,
            epsilon_k_ab,
            na,
            nb,
            psi_dft,
        ))
    }
}

#[pyproto]
impl pyo3::class::basic::PyObjectProtocol for PyGcPcSaftRecord {
    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}

impl_json_handling!(PyGcPcSaftRecord);

impl_segment_record!(
    GcPcSaftRecord,
    PyGcPcSaftRecord,
    JobackRecord,
    PyJobackRecord
);

#[pyclass(name = "GcPcSaftEosParameters", unsendable)]
#[pyo3(
    text_signature = "(pure_records, segmentbinary_records=None, substances=None, search_option='Name')"
)]
#[derive(Clone)]
pub struct PyGcPcSaftEosParameters(pub Rc<GcPcSaftEosParameters>);

impl_parameter_from_segments!(GcPcSaftEosParameters, PyGcPcSaftEosParameters);

#[pymethods]
impl PyGcPcSaftEosParameters {
    fn _repr_markdown_(&self) -> String {
        self.0.to_markdown()
    }
}

#[pyproto]
impl pyo3::class::basic::PyObjectProtocol for PyGcPcSaftEosParameters {
    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}

#[pyclass(name = "GcPcSaftFunctionalParameters", unsendable)]
#[pyo3(
    text_signature = "(pure_records, segmentbinary_records=None, substances=None, search_option='Name')"
)]
#[derive(Clone)]
pub struct PyGcPcSaftFunctionalParameters(pub Rc<GcPcSaftFunctionalParameters>);

impl_parameter_from_segments!(GcPcSaftFunctionalParameters, PyGcPcSaftFunctionalParameters);

#[pymethods]
impl PyGcPcSaftFunctionalParameters {
    fn _repr_markdown_(&self) -> String {
        self.0.to_markdown()
    }

    #[getter]
    fn get_graph(&self, py: Python) -> PyResult<PyObject> {
        let fun: Py<PyAny> = PyModule::from_code(
            py,
            "def f(s): 
                import graphviz
                return graphviz.Source(s.replace('\\\\\"', ''))",
            "",
            "",
        )?
        .getattr("f")?
        .into();
        fun.call1(py, (self.0.graph(),))
    }

    #[getter]
    fn get_k_ij<'py>(&self, py: Python<'py>) -> &'py PyArray2<f64> {
        self.0.k_ij.view().to_pyarray(py)
    }
}

#[pyproto]
impl pyo3::class::basic::PyObjectProtocol for PyGcPcSaftFunctionalParameters {
    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}
