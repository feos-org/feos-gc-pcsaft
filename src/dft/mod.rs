use crate::eos::GcPcSaftOptions;
use feos_core::MolarWeight;
use feos_dft::adsorption::FluidParameters;
use feos_dft::fundamental_measure_theory::{FMTContribution, FMTProperties, FMTVersion};
use feos_dft::{FunctionalContribution, HelmholtzEnergyFunctional, MoleculeShape, DFT};
use ndarray::Array1;
use num_dual::DualNum;
use petgraph::graph::UnGraph;
use quantity::si::{SIArray1, SIUnit, GRAM, MOL};
use std::f64::consts::FRAC_PI_6;
use std::rc::Rc;

mod association;
mod dispersion;
mod hard_chain;
mod parameter;
use association::AssociationFunctional;
use dispersion::AttractiveFunctional;
use hard_chain::ChainFunctional;
pub use parameter::GcPcSaftFunctionalParameters;

/// gc-PC-SAFT Helmholtz energy functional.
pub struct GcPcSaftFunctional {
    pub parameters: Rc<GcPcSaftFunctionalParameters>,
    fmt_version: FMTVersion,
    options: GcPcSaftOptions,
    contributions: Vec<Box<dyn FunctionalContribution>>,
}

impl GcPcSaftFunctional {
    pub fn new(parameters: Rc<GcPcSaftFunctionalParameters>) -> DFT<Self> {
        Self::with_options(
            parameters,
            FMTVersion::WhiteBear,
            GcPcSaftOptions::default(),
        )
    }

    pub fn with_options(
        parameters: Rc<GcPcSaftFunctionalParameters>,
        fmt_version: FMTVersion,
        saft_options: GcPcSaftOptions,
    ) -> DFT<Self> {
        let mut contributions: Vec<Box<dyn FunctionalContribution>> = Vec::with_capacity(4);

        // Hard sphere contribution
        let hs = FMTContribution::new(&parameters, fmt_version);
        contributions.push(Box::new(hs));

        // Hard chains
        let chain = ChainFunctional::new(&parameters);
        contributions.push(Box::new(chain));

        // Dispersion
        let att = AttractiveFunctional::new(&parameters);
        contributions.push(Box::new(att));

        // Association
        if !parameters.assoc_segment.is_empty() {
            let assoc = AssociationFunctional::new(
                &parameters,
                saft_options.max_iter_cross_assoc,
                saft_options.tol_cross_assoc,
            );
            contributions.push(Box::new(assoc));
        }

        (Self {
            parameters,
            fmt_version,
            options: saft_options,
            contributions,
        })
        .into()
    }
}

impl HelmholtzEnergyFunctional for GcPcSaftFunctional {
    fn molecule_shape(&self) -> MoleculeShape {
        MoleculeShape::Heterosegmented(&self.parameters.component_index)
    }

    fn subset(&self, component_list: &[usize]) -> DFT<Self> {
        Self::with_options(
            Rc::new(self.parameters.subset(component_list)),
            self.fmt_version,
            self.options,
        )
    }

    fn compute_max_density(&self, moles: &Array1<f64>) -> f64 {
        let p = &self.parameters;
        let moles_segments: Array1<f64> = p.component_index.iter().map(|&i| moles[i]).collect();
        self.options.max_eta * moles.sum()
            / (FRAC_PI_6 * &p.m * p.sigma.mapv(|v| v.powi(3)) * moles_segments).sum()
    }

    fn contributions(&self) -> &[Box<dyn FunctionalContribution>] {
        &self.contributions
    }

    fn bond_lengths(&self, temperature: f64) -> UnGraph<(), f64> {
        // temperature dependent segment diameter
        let d = self.parameters.hs_diameter(temperature);

        self.parameters.bonds.map(
            |_, _| (),
            |e, _| {
                let (i, j) = self.parameters.bonds.edge_endpoints(e).unwrap();
                let di = d[i.index()];
                let dj = d[j.index()];
                0.5 * (di + dj)
            },
        )
    }
}

impl MolarWeight<SIUnit> for GcPcSaftFunctional {
    fn molar_weight(&self) -> SIArray1 {
        self.parameters.molarweight.clone() * GRAM / MOL
    }
}

impl FMTProperties for GcPcSaftFunctionalParameters {
    fn component_index(&self) -> Array1<usize> {
        self.component_index.clone()
    }

    fn chain_length(&self) -> Array1<f64> {
        self.m.clone()
    }

    fn hs_diameter<D: DualNum<f64>>(&self, temperature: D) -> Array1<D> {
        self.hs_diameter(temperature)
    }
}

impl FluidParameters for GcPcSaftFunctional {
    fn epsilon_k_ff(&self) -> Array1<f64> {
        self.parameters.epsilon_k.clone()
    }

    fn sigma_ff(&self) -> &Array1<f64> {
        &self.parameters.sigma
    }
}
