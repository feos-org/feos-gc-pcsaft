use approx::assert_relative_eq;
use feos_core::parameter::IdentifierOption;
use feos_core::{EosResult, State};
use feos_gc_pcsaft::{GcPcSaft, GcPcSaftParameters};
use ndarray::arr1;
use quantity::si::{KELVIN, MOL};
use std::rc::Rc;

#[test]
fn test_binary() -> EosResult<()> {
    let parameters_feos = GcPcSaftParameters::from_json_segments(
        &["ethanol", "methanol"],
        "parameters/associating.json",
        "parameters/parameters_hetero_segments.json",
        None,
        IdentifierOption::Name,
    )
    .unwrap();
    let feos = Rc::new(GcPcSaft::new(parameters_feos));
    let cp_feos = State::critical_point(
        &feos,
        Some(&(arr1(&[0.5, 0.5]) * MOL)),
        None,
        Default::default(),
    )?;
    assert_relative_eq!(
        cp_feos.temperature,
        536.4129479522177 * KELVIN,
        max_relative = 1e-14
    );
    println!("{}", cp_feos.temperature);
    Ok(())
}
