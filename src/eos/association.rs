use super::hard_sphere::zeta;
use crate::parameters::GcPcSaftParameters;
use feos_core::{EosError, HelmholtzEnergyDual, StateHD};
use ndarray::*;
use ndarray_linalg::Norm;
use num_dual::*;
use std::fmt;
use std::rc::Rc;

#[derive(Clone)]
pub struct Association {
    pub parameters: Rc<GcPcSaftParameters>,
}

#[derive(Clone)]
pub struct CrossAssociation {
    pub parameters: Rc<GcPcSaftParameters>,
    pub max_iter: usize,
    pub tol: f64,
}

fn association_strength<D: DualNum<f64>>(
    p: &GcPcSaftParameters,
    temperature: D,
    diameter: &Array1<D>,
    n2: D,
    n3i: D,
    i: usize,
    j: usize,
) -> D {
    let ai = p.assoc_segment[i];
    let aj = p.assoc_segment[j];
    let k = diameter[ai] * diameter[aj] / (diameter[ai] + diameter[aj]) * (n2 * n3i);
    n3i * (k * (k / 18.0 + 0.5) + 1.0)
        * p.sigma3_kappa_aibj[(i, j)]
        * (temperature.recip() * p.epsilon_k_aibj[(i, j)]).exp_m1()
}

impl<D: DualNum<f64>> HelmholtzEnergyDual<D> for Association {
    fn helmholtz_energy(&self, state: &StateHD<D>) -> D {
        let p = &self.parameters;
        let c = p.component_index[p.assoc_segment[0]];

        // temperature dependent segment diameter
        let diameter = p.hs_diameter(state.temperature);

        // Packing fractions
        let n2 = zeta(p, &diameter, &state.partial_density, [2])[0] * 6.0;
        let n3 = zeta(p, &diameter, &state.partial_density, [3])[0];
        let n3i = (-n3 + 1.0).recip();

        // association strength
        let deltarho = association_strength(p, state.temperature, &diameter, n2, n3i, 0, 0)
            * state.partial_density[c];

        let na = p.na[0];
        let nb = p.nb[0];
        if nb > 0.0 {
            // no cross association, two association sites
            let xa = assoc_site_frac_ab(deltarho, na, nb);
            let xb = (xa - 1.0) * (na / nb) + 1.0;

            state.moles[c] * ((xa.ln() - xa * 0.5 + 0.5) * na + (xb.ln() - xb * 0.5 + 0.5) * nb)
        } else {
            // no cross association, one association site
            let xa = assoc_site_frac_a(deltarho, na);

            state.moles[c] * (xa.ln() - xa * 0.5 + 0.5) * na
        }
    }
}

impl fmt::Display for Association {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Association")
    }
}

pub(crate) fn assoc_site_frac_ab<D: DualNum<f64>>(deltarho: D, na: f64, nb: f64) -> D {
    if deltarho.re() > f64::EPSILON.sqrt() {
        (((deltarho * (na - nb) + 1.0).powi(2) + deltarho * nb * 4.0).sqrt()
            - (deltarho * (nb - na) + 1.0))
            / (deltarho * na * 2.0)
    } else {
        D::one() + deltarho * nb * (deltarho * (nb + na) - 1.0)
    }
}

pub(crate) fn assoc_site_frac_a<D: DualNum<f64>>(deltarho: D, na: f64) -> D {
    if deltarho.re() > f64::EPSILON.sqrt() {
        ((deltarho * na * 4.0 + 1.0).powi(2) - 1.0).sqrt() / (deltarho * na * 2.0)
    } else {
        D::one() + deltarho * na * (deltarho * na * 2.0 - 1.0)
    }
}

impl<D: DualNum<f64> + ScalarOperand> HelmholtzEnergyDual<D> for CrossAssociation
where
    Array2<D>: SolveDual<D>,
{
    fn helmholtz_energy(&self, state: &StateHD<D>) -> D {
        let p = &self.parameters;

        // temperature dependent segment diameter
        let diameter = p.hs_diameter(state.temperature);

        // Packing fractions
        let n2 = zeta(p, &diameter, &state.partial_density, [2])[0] * 6.0;
        let n3 = zeta(p, &diameter, &state.partial_density, [3])[0];
        let n3i = (-n3 + 1.0).recip();

        // association strength
        let nassoc = p.assoc_segment.len();
        let delta = Array::from_shape_fn([nassoc; 2], |(i, j)| {
            association_strength(p, state.temperature, &diameter, n2, n3i, i, j)
        });

        // extract densities of associating segments
        let rho = Array::from_shape_fn(nassoc, |i| {
            state.partial_density[p.component_index[p.assoc_segment[i]]]
        });

        // cross-association according to Michelsen2006
        // initialize monomer fraction
        let mut x = Array::from_elem(2 * nassoc, 0.2);
        for k in 0..self.max_iter {
            if newton_step_cross_association::<f64>(
                &mut x,
                nassoc,
                &delta.map(D::re),
                &p.na,
                &p.nb,
                &rho.map(D::re),
                self.tol,
            )
            .unwrap()
            {
                break;
            }
            if k == self.max_iter - 1 {
                Err(EosError::NotConverged("Cross association".into())).unwrap()
            }
        }
        // calculate derivatives
        let mut x_dual = x.mapv(D::from);
        for _ in 0..D::NDERIV {
            newton_step_cross_association(
                &mut x_dual,
                nassoc,
                &delta,
                &p.na,
                &p.nb,
                &rho,
                self.tol,
            )
            .unwrap();
        }
        // Helmholtz energy density
        let xa = x_dual.slice(s![..nassoc]);
        let xb = x_dual.slice(s![nassoc..]);
        let f = |x: D| x.ln() - x * 0.5 + 0.5;
        (rho * (xa.mapv(f) * &p.na + xb.mapv(f) * &p.nb)).sum() * state.volume
    }
}

impl fmt::Display for CrossAssociation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Cross-association")
    }
}

fn newton_step_cross_association<D: DualNum<f64> + ScalarOperand>(
    x: &mut Array1<D>,
    nassoc: usize,
    delta: &Array2<D>,
    na: &Array1<f64>,
    nb: &Array1<f64>,
    rho: &Array1<D>,
    tol: f64,
) -> Result<bool, EosError>
where
    Array2<D>: SolveDual<D>,
{
    // gradient
    let mut g: Array1<D> = Array::zeros(2 * nassoc);
    // Hessian
    let mut h: Array2<D> = Array::zeros((2 * nassoc, 2 * nassoc));

    // slice arrays
    let (xa, xb) = x.multi_slice_mut((s![..nassoc], s![nassoc..]));
    let (mut ga, mut gb) = g.multi_slice_mut((s![..nassoc], s![nassoc..]));
    let (mut haa, mut hab, mut hba, mut hbb) = h.multi_slice_mut((
        s![..nassoc, ..nassoc],
        s![..nassoc, nassoc..],
        s![nassoc.., ..nassoc],
        s![nassoc.., nassoc..],
    ));

    // calculate gradients and approximate Hessian
    for i in 0..nassoc {
        let d = &delta.index_axis(Axis(0), i) * rho;

        let dnx = (&xb * nb * &d).sum() + 1.0;
        ga[i] = xa[i].recip() - dnx;
        hab.index_axis_mut(Axis(0), i).assign(&(&d * &(-nb)));
        haa[(i, i)] = -dnx / xa[i];

        let dnx = (&xa * na * &d).sum() + 1.0;
        gb[i] = xb[i].recip() - dnx;
        hba.index_axis_mut(Axis(0), i).assign(&(&d * &(-na)));
        hbb[(i, i)] = -dnx / xb[i];
    }

    // Newton step
    x.assign(&(&*x - &h.solve(&g)?));

    // check convergence
    Ok(g.map(D::re).norm() < tol)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::parameters::test::*;
    use approx::assert_relative_eq;
    use feos_core::EosUnit;
    use ndarray::arr1;
    use num_dual::Dual64;
    use quantity::si::{METER, MOL, PASCAL};

    #[test]
    fn test_assoc_propanol() {
        let parameters = propanol();
        let contrib = Association {
            parameters: Rc::new(parameters),
        };
        let temperature = 300.0;
        let volume = METER
            .powi(3)
            .to_reduced(EosUnit::reference_volume())
            .unwrap();
        let moles = (1.5 * MOL).to_reduced(EosUnit::reference_moles()).unwrap();
        let state = StateHD::new(
            Dual64::from_re(temperature),
            Dual64::from_re(volume).derive(),
            arr1(&[Dual64::from_re(moles)]),
        );
        let pressure =
            -contrib.helmholtz_energy(&state).eps[0] * temperature * EosUnit::reference_pressure();
        assert_relative_eq!(pressure, -3.6819598891967344 * PASCAL, max_relative = 1e-10);
    }

    #[test]
    fn test_cross_assoc_propanol() {
        let parameters = propanol();
        let contrib = CrossAssociation {
            parameters: Rc::new(parameters),
            max_iter: 50,
            tol: 1e-10,
        };
        let temperature = 300.0;
        let volume = METER
            .powi(3)
            .to_reduced(EosUnit::reference_volume())
            .unwrap();
        let moles = (1.5 * MOL).to_reduced(EosUnit::reference_moles()).unwrap();
        let state = StateHD::new(
            Dual64::from_re(temperature),
            Dual64::from_re(volume).derive(),
            arr1(&[Dual64::from_re(moles)]),
        );
        let pressure =
            -contrib.helmholtz_energy(&state).eps[0] * temperature * EosUnit::reference_pressure();
        assert_relative_eq!(pressure, -3.6819598891967344 * PASCAL, max_relative = 1e-10);
    }

    #[test]
    fn test_cross_assoc_ethanol_propanol() {
        let parameters = ethanol_propanol(false);
        let contrib = CrossAssociation {
            parameters: Rc::new(parameters),
            max_iter: 50,
            tol: 1e-10,
        };
        let temperature = 300.0;
        let volume = METER
            .powi(3)
            .to_reduced(EosUnit::reference_volume())
            .unwrap();
        let moles = (arr1(&[1.5, 2.5]) * MOL)
            .to_reduced(EosUnit::reference_moles())
            .unwrap();
        let state = StateHD::new(
            Dual64::from_re(temperature),
            Dual64::from_re(volume).derive(),
            moles.mapv(Dual64::from_re),
        );
        let pressure =
            -contrib.helmholtz_energy(&state).eps[0] * temperature * EosUnit::reference_pressure();
        assert_relative_eq!(pressure, -26.105606376765632 * PASCAL, max_relative = 1e-10);
    }
}
