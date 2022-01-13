use serde::{Deserialize, Serialize};

/// gc-PC-SAFT parameter set.
#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct GcPcSaftRecord {
    /// Segment shape factor
    pub m: f64,
    /// Segment diameter in units of Angstrom
    pub sigma: f64,
    /// Energetic parameter in units of Kelvin
    pub epsilon_k: f64,
    /// Dipole moment in units of Debye
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mu: Option<f64>,
    /// Quadrupole moment in units of Debye * Angstrom
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub q: Option<f64>,
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub kappa_ab: Option<f64>,
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub epsilon_k_ab: Option<f64>,
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub na: Option<f64>,
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nb: Option<f64>,
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub psi_dft: Option<f64>,
}

impl GcPcSaftRecord {
    pub fn new(
        m: f64,
        sigma: f64,
        epsilon_k: f64,
        mu: Option<f64>,
        q: Option<f64>,
        kappa_ab: Option<f64>,
        epsilon_k_ab: Option<f64>,
        na: Option<f64>,
        nb: Option<f64>,
        psi_dft: Option<f64>,
    ) -> Self {
        Self {
            m,
            sigma,
            epsilon_k,
            mu,
            q,
            kappa_ab,
            epsilon_k_ab,
            na,
            nb,
            psi_dft,
        }
    }
}

impl std::fmt::Display for GcPcSaftRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GcPcSaftRecord(m={}", self.m)?;
        write!(f, ", sigma={}", self.sigma)?;
        write!(f, ", epsilon_k={}", self.epsilon_k)?;
        if let Some(n) = &self.mu {
            write!(f, ", mu={}", n)?;
        }
        if let Some(n) = &self.q {
            write!(f, ", q={}", n)?;
        }
        if let Some(n) = &self.kappa_ab {
            write!(f, ", kappa_ab={}", n)?;
        }
        if let Some(n) = &self.epsilon_k_ab {
            write!(f, ", epsilon_k_ab={}", n)?;
        }
        if let Some(n) = &self.na {
            write!(f, ", na={}", n)?;
        }
        if let Some(n) = &self.nb {
            write!(f, ", nb={}", n)?;
        }
        write!(f, ")")
    }
}

// pub struct GcPcSaftParameters<B> {
//     pub molarweight: Array1<f64>,
//     pub component_index: Array1<usize>,
//     pub identifiers: Vec<String>,
//     pub m: Array1<f64>,
//     pub sigma: Array1<f64>,
//     pub epsilon_k: Array1<f64>,
//     pub bonds: B,
//     // pub mu: Array1<f64>,
//     // pub q: Array1<f64>,
//     // pub mu2: Array1<f64>,
//     // pub q2: Array1<f64>,
//     pub assoc_segment: Array1<usize>,
//     pub kappa_ab: Array1<f64>,
//     pub epsilon_k_ab: Array1<f64>,
//     pub na: Array1<f64>,
//     pub nb: Array1<f64>,
//     pub k_ij: Array2<f64>,
//     pub sigma_ij: Array2<f64>,
//     pub epsilon_k_ij: Array2<f64>,
//     pub sigma3_kappa_aibj: Array2<f64>,
//     pub epsilon_k_aibj: Array2<f64>,
//     pub pure_records: Vec<GroupContributionRecord>,
//     pub segment_records: Vec<SegmentRecord<GcPcSaftRecord, JobackRecord>>,
//     pub binary_segment_records: Option<Vec<BinaryRecord<String, f64>>>,
//     pub joback_records: Option<Vec<JobackRecord>>,
// }
