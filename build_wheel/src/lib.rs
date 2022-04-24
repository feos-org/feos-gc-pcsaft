use feos_core::python::joback::PyJobackRecord;
use feos_core::python::parameter::*;
use feos_core::{Contributions, Verbosity};
use feos_gc_pcsaft::python::*;
use pyo3::prelude::*;
use pyo3::wrap_pymodule;
use quantity::python::__PYO3_PYMODULE_DEF_QUANTITY;

mod dft;
mod eos;
use dft::__PYO3_PYMODULE_DEF_DFT;
use eos::__PYO3_PYMODULE_DEF_EOS;

#[pymodule]
pub fn feos_pcsaft(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyIdentifier>()?;
    m.add_class::<Verbosity>()?;
    m.add_class::<Contributions>()?;
    m.add_class::<PyChemicalRecord>()?;
    m.add_class::<PyJobackRecord>()?;

    m.add_class::<PyGcPcSaftRecord>()?;
    m.add_class::<PySegmentRecord>()?;
    m.add_class::<PyBinaryRecord>()?;
    m.add_class::<PyBinarySegmentRecord>()?;
    m.add_class::<PyGcPcSaftEosParameters>()?;
    m.add_class::<PyGcPcSaftFunctionalParameters>()?;

    m.add_wrapped(wrap_pymodule!(eos))?;
    m.add_wrapped(wrap_pymodule!(dft))?;
    m.add_wrapped(wrap_pymodule!(quantity))?;

    py.run(
        "\
import sys
sys.modules['feos_gc_pcsaft.eos'] = eos
sys.modules['feos_gc_pcsaft.dft'] = dft
quantity.SINumber.__module__ = 'feos_gc_pcsaft.si'
quantity.SIArray1.__module__ = 'feos_gc_pcsaft.si'
quantity.SIArray2.__module__ = 'feos_gc_pcsaft.si'
quantity.SIArray3.__module__ = 'feos_gc_pcsaft.si'
quantity.SIArray4.__module__ = 'feos_gc_pcsaft.si'
sys.modules['feos_gc_pcsaft.si'] = quantity
    ",
        None,
        Some(m.dict()),
    )?;
    Ok(())
}
