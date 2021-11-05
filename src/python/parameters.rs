use crate::parameters::{GcPcSaftParameters, GcPcSaftRecord};
use feos_core::joback::JobackRecord;
use feos_core::parameter::{IdentifierOption, NoRecord, ParameterError, PureRecord, SegmentRecord};
use feos_core::python::joback::PyJobackRecord;
use feos_core::python::parameter::{
    PyBinarySegmentRecord, PyChemicalRecord, PyIdentifier, PyNoRecord,
};
use feos_core::{
    impl_json_handling, impl_parameter_from_segments, impl_pure_record, impl_segment_record,
};
use pyo3::prelude::*;
use std::convert::TryFrom;

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
        q: Option<f64>,
        kappa_ab: Option<f64>,
        epsilon_k_ab: Option<f64>,
        na: Option<f64>,
        nb: Option<f64>,
    ) -> Self {
        Self(GcPcSaftRecord::new(
            m,
            sigma,
            epsilon_k,
            mu,
            q,
            kappa_ab,
            epsilon_k_ab,
            na,
            nb,
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

impl_pure_record!(NoRecord, PyNoRecord, NoRecord, PyNoRecord);
impl_segment_record!(
    GcPcSaftRecord,
    PyGcPcSaftRecord,
    JobackRecord,
    PyJobackRecord
);

#[pyclass(name = "GcPcSaftParameters", unsendable)]
#[pyo3(
    text_signature = "(pure_records, segmentbinary_records=None, substances=None, search_option='Name')"
)]
#[derive(Clone)]
pub struct PyGcPcSaftParameters(pub GcPcSaftParameters);

impl_parameter_from_segments!(GcPcSaftParameters, PyGcPcSaftParameters);

#[pymethods]
impl PyGcPcSaftParameters {
    fn _repr_markdown_(&self) -> String {
        self.0.to_markdown()
    }
}

#[pyproto]
impl pyo3::class::basic::PyObjectProtocol for PyGcPcSaftParameters {
    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}
