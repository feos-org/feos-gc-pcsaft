use crate::record::GcPcSaftRecord;
use feos_core::joback::JobackRecord;
use feos_core::parameter::{
    BinaryRecord, ChemicalRecord, FromSegments, IdentifierOption, ParameterError, SegmentRecord,
};
use indexmap::{IndexMap, IndexSet};
use ndarray::{Array1, Array2};
use quantity::si::{JOULE, KB, KELVIN};
use std::fmt::Write;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Parameter set required for the gc-PC-SAFT equation of state.
pub struct GcPcSaftEosParameters {
    pub molarweight: Array1<f64>,
    pub component_index: Array1<usize>,
    identifiers: Vec<String>,

    pub m: Array1<f64>,
    pub sigma: Array1<f64>,
    pub epsilon_k: Array1<f64>,
    pub bonds: IndexMap<[usize; 2], f64>,

    pub assoc_segment: Array1<usize>,
    pub n: Array1<f64>,
    kappa_ab: Array1<f64>,
    epsilon_k_ab: Array1<f64>,
    pub na: Array1<f64>,
    pub nb: Array1<f64>,

    pub dipole_comp: Array1<usize>,
    mu: Array1<f64>,
    pub mu2: Array1<f64>,
    pub m_mix: Array1<f64>,
    pub s_ij: Array2<f64>,
    pub e_k_ij: Array2<f64>,

    pub k_ij: Array2<f64>,
    pub sigma_ij: Array2<f64>,
    pub epsilon_k_ij: Array2<f64>,
    pub sigma3_kappa_aibj: Array2<f64>,
    pub epsilon_k_aibj: Array2<f64>,

    pub chemical_records: Vec<ChemicalRecord>,
    segment_records: Vec<SegmentRecord<GcPcSaftRecord, JobackRecord>>,
    binary_segment_records: Option<Vec<BinaryRecord<String, f64>>>,
    pub joback_records: Option<Vec<JobackRecord>>,
}

impl GcPcSaftEosParameters {
    pub fn from_segments(
        chemical_records: Vec<ChemicalRecord>,
        segment_records: Vec<SegmentRecord<GcPcSaftRecord, JobackRecord>>,
        binary_segment_records: Option<Vec<BinaryRecord<String, f64>>>,
    ) -> Result<Self, ParameterError> {
        let segment_map: IndexMap<_, _> = segment_records
            .iter()
            .map(|r| (r.identifier.clone(), r.clone()))
            .collect();

        let mut molarweight = Array1::zeros(chemical_records.len());
        let mut component_index = Vec::new();
        let mut identifiers = Vec::new();
        let mut m = Vec::new();
        let mut sigma = Vec::new();
        let mut epsilon_k = Vec::new();
        let mut bonds = IndexMap::with_capacity(segment_records.len());
        let mut assoc_segment = Vec::new();
        let mut n = Vec::new();
        let mut kappa_ab = Vec::new();
        let mut epsilon_k_ab = Vec::new();
        let mut na = Vec::new();
        let mut nb = Vec::new();

        let mut dipole_comp = Vec::new();
        let mut mu = Vec::new();
        let mut mu2 = Vec::new();
        let mut m_mix = Vec::new();
        let mut sigma_mix = Vec::new();
        let mut epsilon_k_mix = Vec::new();

        let mut joback_records = Vec::new();

        for (i, chemical_record) in chemical_records.iter().enumerate() {
            let mut segment_indices = IndexMap::with_capacity(segment_records.len());
            let (segment_counts, bond_counts) = chemical_record.segment_and_bond_count();
            let count: IndexMap<_, _> = segment_counts
                .iter()
                .map(|(id, &count)| {
                    let segment = segment_map
                        .get(id)
                        .ok_or_else(|| ParameterError::ComponentsNotFound(id.clone()))?;
                    Ok((segment, count))
                })
                .collect::<Result<_, ParameterError>>()?;

            let mut m_i = 0.0;
            let mut sigma_i = 0.0;
            let mut epsilon_k_i = 0.0;
            let mut mu2_i = 0.0;

            for (segment, count) in count.iter() {
                segment_indices.insert(segment.identifier.clone(), m.len());

                molarweight[i] += segment.molarweight * count;

                component_index.push(i);
                identifiers.push(segment.identifier.clone());
                m.push(segment.model_record.m * count);
                sigma.push(segment.model_record.sigma);
                epsilon_k.push(segment.model_record.epsilon_k);

                if let (Some(k), Some(e)) = (
                    segment.model_record.kappa_ab,
                    segment.model_record.epsilon_k_ab,
                ) {
                    assoc_segment.push(m.len() - 1);
                    n.push(*count);
                    kappa_ab.push(k);
                    epsilon_k_ab.push(e);
                    na.push(segment.model_record.na.unwrap_or(1.0));
                    nb.push(segment.model_record.nb.unwrap_or(1.0));
                }

                m_i += segment.model_record.m * count;
                sigma_i += segment.model_record.m * segment.model_record.sigma.powi(3) * count;
                epsilon_k_i += segment.model_record.m * segment.model_record.epsilon_k * count;
                if let Some(mu) = segment.model_record.mu {
                    mu2_i += mu.powi(2) * count;
                }
            }

            if mu2_i > 0.0 {
                dipole_comp.push(i);
                mu.push(mu2_i.sqrt());
                mu2.push(mu2_i / m_i * (1e-19 * (JOULE / KELVIN / KB).into_value().unwrap()));
                m_mix.push(m_i);
                sigma_mix.push((sigma_i / m_i).cbrt());
                epsilon_k_mix.push(epsilon_k_i / m_i);
            }

            for (b, &count) in bond_counts.iter() {
                let i1 = segment_indices.get(&b[0]);
                let i2 = segment_indices.get(&b[1]);
                if let (Some(&i1), Some(&i2)) = (i1, i2) {
                    let indices = if i1 > i2 { [i2, i1] } else { [i1, i2] };
                    let bond = bonds.entry(indices).or_insert(0.0);
                    *bond += count;
                }
            }

            let ideal_gas_segments: Option<Vec<_>> = count
                .iter()
                .map(|(s, &n)| s.ideal_gas_record.clone().map(|ig| (ig, n)))
                .collect();

            joback_records.push(
                ideal_gas_segments
                    .as_ref()
                    .map(|s| JobackRecord::from_segments(s)),
            );
        }

        // Binary interaction parameter
        let mut k_ij = Array2::zeros([epsilon_k.len(); 2]);
        if let Some(binary_segment_records) = binary_segment_records.as_ref() {
            let mut binary_segment_records_map = IndexMap::new();
            for binary_record in binary_segment_records {
                binary_segment_records_map.insert(
                    (binary_record.id1.clone(), binary_record.id2.clone()),
                    binary_record.model_record,
                );
                binary_segment_records_map.insert(
                    (binary_record.id2.clone(), binary_record.id1.clone()),
                    binary_record.model_record,
                );
            }
            for (i, id1) in identifiers.iter().enumerate() {
                for (j, id2) in identifiers.iter().cloned().enumerate() {
                    if component_index[i] != component_index[j] {
                        if let Some(k) = binary_segment_records_map.get(&(id1.clone(), id2)) {
                            k_ij[(i, j)] = *k;
                        }
                    }
                }
            }
        }

        // Combining rules dispersion
        let sigma_ij =
            Array2::from_shape_fn([sigma.len(); 2], |(i, j)| 0.5 * (sigma[i] + sigma[j]));
        let epsilon_k_ij = Array2::from_shape_fn([epsilon_k.len(); 2], |(i, j)| {
            (epsilon_k[i] * epsilon_k[j]).sqrt()
        }) * (1.0 - &k_ij);

        // Combining rules polar
        let s_ij = Array2::from_shape_fn([dipole_comp.len(); 2], |(i, j)| {
            0.5 * (sigma_mix[i] + sigma_mix[j])
        });
        let e_k_ij = Array2::from_shape_fn([dipole_comp.len(); 2], |(i, j)| {
            (epsilon_k_mix[i] * epsilon_k_mix[j]).sqrt()
        });

        // Association
        let sigma3_kappa_aibj = Array2::from_shape_fn([kappa_ab.len(); 2], |(i, j)| {
            (sigma[assoc_segment[i]] * sigma[assoc_segment[j]]).powf(1.5)
                * (kappa_ab[i] * kappa_ab[j]).sqrt()
        });
        let epsilon_k_aibj = Array2::from_shape_fn([epsilon_k_ab.len(); 2], |(i, j)| {
            0.5 * (epsilon_k_ab[i] + epsilon_k_ab[j])
        });

        Ok(Self {
            molarweight,
            component_index: Array1::from_vec(component_index),
            identifiers,
            n: Array1::from_vec(n),
            m: Array1::from_vec(m),
            sigma: Array1::from_vec(sigma),
            epsilon_k: Array1::from_vec(epsilon_k),
            bonds,
            assoc_segment: Array1::from_vec(assoc_segment),
            kappa_ab: Array1::from_vec(kappa_ab),
            epsilon_k_ab: Array1::from_vec(epsilon_k_ab),
            na: Array1::from_vec(na),
            nb: Array1::from_vec(nb),
            dipole_comp: Array1::from_vec(dipole_comp),
            mu: Array1::from_vec(mu),
            mu2: Array1::from_vec(mu2),
            m_mix: Array1::from_vec(m_mix),
            s_ij,
            e_k_ij,
            k_ij,
            sigma_ij,
            epsilon_k_ij,
            sigma3_kappa_aibj,
            epsilon_k_aibj,
            chemical_records,
            segment_records,
            binary_segment_records,
            joback_records: joback_records.into_iter().collect(),
        })
    }

    pub fn from_json_segments<P>(
        substances: &[&str],
        file_pure: P,
        file_segments: P,
        file_binary: Option<P>,
        search_option: IdentifierOption,
    ) -> Result<Self, ParameterError>
    where
        P: AsRef<Path>,
    {
        let queried: IndexSet<String> = substances
            .iter()
            .map(|identifier| identifier.to_string())
            .collect();

        let reader = BufReader::new(File::open(file_pure)?);
        let chemical_records: Vec<ChemicalRecord> = serde_json::from_reader(reader)?;
        let mut record_map: IndexMap<_, _> = chemical_records
            .into_iter()
            .filter_map(|record| {
                record
                    .identifier()
                    .as_string(search_option)
                    .map(|i| (i, record))
            })
            .collect();

        // Compare queried components and available components
        let available: IndexSet<String> = record_map
            .keys()
            .map(|identifier| identifier.to_string())
            .collect();
        if !queried.is_subset(&available) {
            let missing: Vec<String> = queried.difference(&available).cloned().collect();
            return Err(ParameterError::ComponentsNotFound(format!("{:?}", missing)));
        };

        // Collect all pure records that were queried
        let chemical_records: Vec<_> = queried
            .iter()
            .filter_map(|identifier| record_map.remove(&identifier.clone()))
            .collect();

        // Read segment records
        let segment_records: Vec<SegmentRecord<GcPcSaftRecord, JobackRecord>> =
            serde_json::from_reader(BufReader::new(File::open(file_segments)?))?;

        // Read binary records
        let binary_records = file_binary
            .map(|file_binary| {
                let reader = BufReader::new(File::open(file_binary)?);
                let binary_records: Result<Vec<BinaryRecord<String, f64>>, ParameterError> =
                    Ok(serde_json::from_reader(reader)?);
                binary_records
            })
            .transpose()?;

        Self::from_segments(chemical_records, segment_records, binary_records)
    }

    pub fn subset(&self, component_list: &[usize]) -> Self {
        let chemical_records: Vec<_> = component_list
            .iter()
            .map(|&i| self.chemical_records[i].clone())
            .collect();
        Self::from_segments(
            chemical_records,
            self.segment_records.clone(),
            self.binary_segment_records.clone(),
        )
        .unwrap()
    }

    pub fn to_markdown(&self) -> String {
        let mut output = String::new();
        let o = &mut output;
        write!(
            o,
            "|component|molarweight|dipole moment|segment|count|$m$|$\\sigma$|$\\varepsilon$|$\\kappa_{{AB}}$|$\\varepsilon_{{AB}}$|$N_A$|$N_B$|$\\mu$|$Q$|\n|-|-|-|-|-|-|-|-|-|-|-|-|-|-|"
        )
        .unwrap();
        for i in 0..self.m.len() {
            let component = if i > 0 && self.component_index[i] == self.component_index[i - 1] {
                "||".to_string()
            } else {
                let pure = self.chemical_records[self.component_index[i]].identifier();
                format!(
                    "{}|{}|{}",
                    pure.name.as_ref().unwrap_or(&pure.cas),
                    self.molarweight[self.component_index[i]],
                    if let Some(d) = self.dipole_comp.iter().position(|&d| d == i) {
                        format!("{}", self.mu[d])
                    } else {
                        "".into()
                    }
                )
            };
            let association = if let Some(a) = self.assoc_segment.iter().position(|&a| a == i) {
                format!(
                    "{}|{}|{}|{}",
                    self.kappa_ab[a], self.epsilon_k_ab[a], self.na[a], self.nb[a]
                )
            } else {
                "|||".to_string()
            };
            write!(
                o,
                "\n|{}|{}|{}|{}|{}|{}|{}|||",
                component,
                self.identifiers[i],
                self.n[i],
                self.m[i],
                self.sigma[i],
                self.epsilon_k[i],
                association
            )
            .unwrap();
        }
        write!(o, "\n\n|component|segment 1|segment 2|bonds|\n|-|-|-|-|").unwrap();

        let mut last_component = None;
        for ([c1, c2], &c) in &self.bonds {
            let pure = self.chemical_records[self.component_index[*c1]].identifier();
            let component = if let Some(last) = last_component {
                if last == pure {
                    ""
                } else {
                    last_component = Some(pure);
                    pure.name.as_ref().unwrap_or(&pure.cas)
                }
            } else {
                last_component = Some(pure);
                pure.name.as_ref().unwrap_or(&pure.cas)
            };
            write!(
                o,
                "\n|{}|{}|{}|{}|",
                component, self.identifiers[*c1], self.identifiers[*c2], c
            )
            .unwrap();
        }

        output
    }
}

impl std::fmt::Display for GcPcSaftEosParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GcPcSaftParameters(")?;
        write!(f, "\n\tmolarweight={}", self.molarweight)?;
        write!(f, "\n\tcomponent_index={}", self.component_index)?;
        write!(f, "\n\tm={}", self.m)?;
        write!(f, "\n\tsigma={}", self.sigma)?;
        write!(f, "\n\tepsilon_k={}", self.epsilon_k)?;
        write!(f, "\n\tbonds={:?}", self.bonds)?;
        if !self.assoc_segment.is_empty() {
            write!(f, "\n\tassoc_segment={}", self.assoc_segment)?;
            write!(f, "\n\tkappa_ab={}", self.kappa_ab)?;
            write!(f, "\n\tepsilon_k_ab={}", self.epsilon_k_ab)?;
            write!(f, "\n\tna={}", self.na)?;
            write!(f, "\n\tnb={}", self.nb)?;
        }
        if !self.dipole_comp.is_empty() {
            write!(f, "\n\tdipole_comp={}", self.dipole_comp)?;
            write!(f, "\n\tmu={}", self.mu)?;
        }
        write!(f, "\n)")
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    use feos_core::parameter::{ChemicalRecord, Identifier};

    fn ch3() -> SegmentRecord<GcPcSaftRecord, JobackRecord> {
        SegmentRecord::new(
            "CH3".into(),
            15.0,
            GcPcSaftRecord::new(0.77247, 3.6937, 181.49, None, None, None, None, None, None),
            None,
        )
    }

    fn ch2() -> SegmentRecord<GcPcSaftRecord, JobackRecord> {
        SegmentRecord::new(
            "CH2".into(),
            14.0,
            GcPcSaftRecord::new(0.7912, 3.0207, 157.23, None, None, None, None, None, None),
            None,
        )
    }

    fn oh() -> SegmentRecord<GcPcSaftRecord, JobackRecord> {
        SegmentRecord::new(
            "OH".into(),
            0.0,
            GcPcSaftRecord::new(
                1.0231,
                2.7702,
                334.29,
                None,
                Some(0.009583),
                Some(2575.9),
                None,
                None,
                None,
            ),
            None,
        )
    }

    pub fn ch3_oh() -> BinaryRecord<String, f64> {
        BinaryRecord::new("CH3".to_string(), "OH".to_string(), -0.0087)
    }

    pub fn propane() -> GcPcSaftEosParameters {
        let pure = ChemicalRecord::new(
            Identifier::new("74-98-6", Some("propane"), None, None, None, None),
            vec!["CH3".into(), "CH2".into(), "CH3".into()],
            None,
        );
        GcPcSaftEosParameters::from_segments(vec![pure], vec![ch3(), ch2()], None).unwrap()
    }

    pub fn propanol() -> GcPcSaftEosParameters {
        let pure = ChemicalRecord::new(
            Identifier::new("71-23-8", Some("1-propanol"), None, None, None, None),
            vec!["CH3".into(), "CH2".into(), "CH2".into(), "OH".into()],
            None,
        );
        GcPcSaftEosParameters::from_segments(vec![pure], vec![ch3(), ch2(), oh()], None).unwrap()
    }

    pub fn ethanol_propanol(binary: bool) -> GcPcSaftEosParameters {
        let ethanol = ChemicalRecord::new(
            Identifier::new("64-17-5", Some("ethanol"), None, None, None, None),
            vec!["CH3".into(), "CH2".into(), "OH".into()],
            None,
        );
        let propanol = ChemicalRecord::new(
            Identifier::new("71-23-8", Some("1-propanol"), None, None, None, None),
            vec!["CH3".into(), "CH2".into(), "CH2".into(), "OH".into()],
            None,
        );
        let binary = if binary { Some(vec![ch3_oh()]) } else { None };
        GcPcSaftEosParameters::from_segments(
            vec![ethanol, propanol],
            vec![ch3(), ch2(), oh()],
            binary,
        )
        .unwrap()
    }

    #[test]
    fn test_kij() {
        let params = ethanol_propanol(true);
        let identifiers: Vec<_> = params.identifiers.iter().enumerate().collect();
        let ch3 = identifiers.iter().find(|&&(_, id)| id == "CH3").unwrap();
        let ch2 = identifiers
            .iter()
            .skip(3)
            .find(|&&(_, id)| id == "CH2")
            .unwrap();
        let oh = identifiers
            .iter()
            .skip(3)
            .find(|&&(_, id)| id == "OH")
            .unwrap();
        println!("{:?}", params.identifiers);
        println!("{}", params.k_ij);
        // CH3 - CH2
        assert_eq!(
            params.epsilon_k_ij[(ch3.0, ch2.0)],
            (181.49f64 * 157.23).sqrt()
        );
        // CH3 - OH
        assert_eq!(
            params.epsilon_k_ij[(ch3.0, oh.0)],
            (181.49f64 * 334.29).sqrt() * 1.0087
        );
    }
}
