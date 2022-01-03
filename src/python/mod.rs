use feos_core::python::joback::PyJobackRecord;
use feos_core::python::parameter::*;
use feos_core::python::*;
use pyo3::prelude::*;
use pyo3::wrap_pymodule;
use quantity::python::PyInit_quantity;

mod dft;
mod eos;
mod micelles;
mod parameter;
use dft::*;
use eos::*;
use parameter::*;

#[pymodule]
pub fn feos_gc_pcsaft(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyIdentifier>()?;
    m.add_class::<PyVerbosity>()?;
    m.add_class::<PyContributions>()?;
    m.add_class::<PyChemicalRecord>()?;
    m.add_class::<PyChemicalRecord>()?;
    m.add_class::<PyJobackRecord>()?;

    m.add_class::<PyGcPcSaftRecord>()?;
    m.add_class::<PyPureRecord>()?;
    m.add_class::<PySegmentRecord>()?;
    m.add_class::<PyBinaryRecord>()?;
    m.add_class::<PyBinarySegmentRecord>()?;
    m.add_class::<PyGcPcSaftEosParameters>()?;
    m.add_class::<PyGcPcSaftFunctionalParameters>()?;

    let eos = PyModule::new(py, "eos")?;
    gc_pcsaft_eos(py, eos)?;
    m.add_submodule(eos)?;

    let dft = PyModule::new(py, "dft")?;
    gc_pcsaft_dft(py, dft)?;
    m.add_submodule(dft)?;

    m.add_wrapped(wrap_pymodule!(quantity))?;

    py.run(
        "\
import sys
sys.modules['feos_gc_pcsaft.eos'] = eos
sys.modules['feos_gc_pcsaft.dft'] = dft
sys.modules['feos_gc_pcsaft.si'] = quantity
    ",
        None,
        Some(m.dict()),
    )?;
    Ok(())
}
