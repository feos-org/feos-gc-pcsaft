use crate::parameter::{GcPcSaftParameters, GcPcSaftRecord};
use feos_core::joback::JobackRecord;
use feos_core::parameter::{
    BinaryRecord, GroupContributionRecord, IdentifierOption, ParameterError, SegmentRecord,
};
use indexmap::{IndexMap, IndexSet};
use ndarray::{Array1, Array2};
use petgraph::graph::{Graph, UnGraph};
use std::fmt::Write;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub type GcPcSaftFunctionalParameters = GcPcSaftParameters<UnGraph<(), ()>>;

// pub struct GcPcSaftFunctionalParameters {
//     pub molarweight: Array1<f64>,
//     pub component_index: Array1<usize>,
//     identifiers: Vec<String>,
//     pub m: Array1<f64>,
//     pub sigma: Array1<f64>,
//     pub epsilon_k: Array1<f64>,
//     pub bonds: UnGraph<(), ()>,
//     pub assoc_segment: Array1<usize>,
//     kappa_ab: Array1<f64>,
//     epsilon_k_ab: Array1<f64>,
//     pub na: Array1<f64>,
//     pub nb: Array1<f64>,
//     pub sigma_ij: Array2<f64>,
//     pub epsilon_k_ij: Array2<f64>,
//     pub sigma3_kappa_aibj: Array2<f64>,
//     pub epsilon_k_aibj: Array2<f64>,
//     pub pure_records: Vec<GroupContributionRecord>,
//     segment_records: Vec<SegmentRecord<GcPcSaftRecord, JobackRecord>>,
//     binary_segment_records: Option<Vec<BinaryRecord<String, f64>>>,
//     // pub joback_records: Option<Vec<JobackRecord>>,
// }

impl GcPcSaftFunctionalParameters {
    pub fn from_segments(
        pure_records: Vec<GroupContributionRecord>,
        segment_records: Vec<SegmentRecord<GcPcSaftRecord, JobackRecord>>,
        binary_segment_records: Option<Vec<BinaryRecord<String, f64>>>,
    ) -> Result<Self, ParameterError> {
        let segment_map: IndexMap<_, _> = segment_records
            .iter()
            .map(|r| (r.identifier.clone(), r.clone()))
            .collect();

        let mut molarweight = Array1::zeros(pure_records.len());
        let mut component_index = Vec::new();
        let mut identifiers = Vec::new();
        let mut m = Vec::new();
        let mut sigma = Vec::new();
        let mut epsilon_k = Vec::new();
        let mut bonds = Graph::default();
        let mut assoc_segment = Vec::new();
        let mut kappa_ab = Vec::new();
        let mut epsilon_k_ab = Vec::new();
        let mut na = Vec::new();
        let mut nb = Vec::new();

        // let mut joback_records = Vec::new();

        let mut segment_index = 0;
        for (i, record) in pure_records.iter().enumerate() {
            let chemical_record = record
                .chemical_record
                .as_ref()
                .ok_or(ParameterError::InsufficientInformation)?;

            bonds.extend_with_edges(chemical_record.bonds.iter().map(|x| {
                (
                    (segment_index + x[0]) as u32,
                    (segment_index + x[1]) as u32,
                    (),
                )
            }));

            for id in &chemical_record.segments {
                let segment = segment_map
                    .get(id)
                    .ok_or_else(|| ParameterError::ComponentsNotFound(id.to_string()))?;
                molarweight[i] += segment.molarweight;
                component_index.push(i);
                identifiers.push(id.clone());
                m.push(segment.model_record.m);
                sigma.push(segment.model_record.sigma);
                epsilon_k.push(segment.model_record.epsilon_k);

                if let (Some(k), Some(e)) = (
                    segment.model_record.kappa_ab,
                    segment.model_record.epsilon_k_ab,
                ) {
                    assoc_segment.push(segment_index);
                    kappa_ab.push(k);
                    epsilon_k_ab.push(e);
                    na.push(segment.model_record.na.unwrap_or(1.0));
                    nb.push(segment.model_record.nb.unwrap_or(1.0));
                }

                segment_index += 1;
            }

            // let segment_count = chemical_record.segment_count(&segment_records)?;
            // let ideal_gas_segments: Option<Vec<_>> = segment_count
            //     .iter()
            //     .map(|(s, &n)| s.ideal_gas_record.clone().map(|ig| (ig, n)))
            //     .collect();

            // joback_records.push(
            //     ideal_gas_segments
            //         .as_ref()
            //         .map(|s| JobackRecord::from_segments(s, None))
            //         .transpose()?,
            // );
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
            (epsilon_k[i] * epsilon_k[j]).sqrt() * (1.0 - k_ij[(i, j)])
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
            m: Array1::from_vec(m),
            sigma: Array1::from_vec(sigma),
            epsilon_k: Array1::from_vec(epsilon_k),
            bonds,
            assoc_segment: Array1::from_vec(assoc_segment),
            kappa_ab: Array1::from_vec(kappa_ab),
            epsilon_k_ab: Array1::from_vec(epsilon_k_ab),
            na: Array1::from_vec(na),
            nb: Array1::from_vec(nb),
            sigma_ij,
            epsilon_k_ij,
            sigma3_kappa_aibj,
            epsilon_k_aibj,
            pure_records,
            segment_records,
            binary_segment_records,
            joback_records: None,
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
        let pure_records: Vec<GroupContributionRecord> = serde_json::from_reader(reader)?;
        let mut record_map: IndexMap<_, _> = pure_records
            .into_iter()
            .filter_map(|record| {
                record
                    .identifier
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
        let pure_records: Vec<_> = queried
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

        Self::from_segments(pure_records, segment_records, binary_records)
    }

    pub fn subset(&self, component_list: &[usize]) -> Self {
        let pure_records: Vec<_> = component_list
            .iter()
            .map(|&i| self.pure_records[i].clone())
            .collect();
        Self::from_segments(
            pure_records,
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
            "|component|molarweight|segment|$m$|$\\sigma$|$\\varepsilon$|$\\kappa_{{AB}}$|$\\varepsilon_{{AB}}$|$N_A$|$N_B$|\n|-|-|-|-|-|-|-|-|-|-|"
        )
        .unwrap();
        for i in 0..self.m.len() {
            let component = if i > 0 && self.component_index[i] == self.component_index[i - 1] {
                "|".to_string()
            } else {
                let pure = &self.pure_records[self.component_index[i]].identifier;
                format!(
                    "{}|{}",
                    pure.name.as_ref().unwrap_or(&pure.cas),
                    self.molarweight[self.component_index[i]]
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
                "\n|{}|{}|{}|{}|{}|{}|||",
                component,
                self.identifiers[i],
                self.m[i],
                self.sigma[i],
                self.epsilon_k[i],
                association
            )
            .unwrap();
        }
        // write!(o, "\n\n|component|segment 1|segment 2|bonds|\n|-|-|-|-|").unwrap();

        // let mut last_component = None;
        // for ([c1, c2], &c) in &self.bonds {
        //     let pure = &self.pure_records[self.component_index[*c1]].identifier;
        //     let component = if let Some(last) = last_component {
        //         if last == pure {
        //             ""
        //         } else {
        //             last_component = Some(pure);
        //             pure.name.as_ref().unwrap_or(&pure.cas)
        //         }
        //     } else {
        //         last_component = Some(pure);
        //         pure.name.as_ref().unwrap_or(&pure.cas)
        //     };
        //     write!(
        //         o,
        //         "\n|{}|{}|{}|{}|",
        //         component, self.identifiers[*c1], self.identifiers[*c2], c
        //     )
        //     .unwrap();
        // }

        output
    }
}

impl std::fmt::Display for GcPcSaftFunctionalParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GcPcSaftFunctionalParameters(")?;
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
        write!(f, "\n)")
    }
}

// #[cfg(test)]
// pub mod test {
//     use super::*;
//     use feos_core::parameter::{ChemicalRecord, GroupContributionRecord, Identifier};

//     fn ch3() -> SegmentRecord<GcPcSaftRecord, JobackRecord> {
//         SegmentRecord::new(
//             "CH3".into(),
//             15.0,
//             GcPcSaftRecord::new(0.77247, 3.6937, 181.49, None, None, None, None, None, None),
//             None,
//         )
//     }

//     fn ch2() -> SegmentRecord<GcPcSaftRecord, JobackRecord> {
//         SegmentRecord::new(
//             "CH2".into(),
//             14.0,
//             GcPcSaftRecord::new(0.7912, 3.0207, 157.23, None, None, None, None, None, None),
//             None,
//         )
//     }

//     fn oh() -> SegmentRecord<GcPcSaftRecord, JobackRecord> {
//         SegmentRecord::new(
//             "OH".into(),
//             0.0,
//             GcPcSaftRecord::new(
//                 1.0231,
//                 2.7702,
//                 334.29,
//                 None,
//                 None,
//                 Some(0.009583),
//                 Some(2575.9),
//                 None,
//                 None,
//             ),
//             None,
//         )
//     }

//     pub fn ch3_oh() -> BinaryRecord<String, f64> {
//         BinaryRecord::new("CH3".to_string(), "OH".to_string(), -0.0087)
//     }

//     pub fn propane() -> GcPcSaftFunctionalParameters {
//         let pure = GroupContributionRecord::new(
//             Identifier::new("74-98-6", Some("propane"), None, None, None, None),
//             0.0,
//             Some(ChemicalRecord::new(
//                 vec!["CH3".into(), "CH2".into(), "CH3".into()],
//                 None,
//             )),
//             None,
//             None,
//         );
//         GcPcSaftFunctionalParameters::from_segments(vec![pure], vec![ch3(), ch2()], None).unwrap()
//     }

//     pub fn propanol() -> GcPcSaftFunctionalParameters {
//         let pure = GroupContributionRecord::new(
//             Identifier::new("71-23-8", Some("1-propanol"), None, None, None, None),
//             0.0,
//             Some(ChemicalRecord::new(
//                 vec!["CH3".into(), "CH2".into(), "CH2".into(), "OH".into()],
//                 None,
//             )),
//             None,
//             None,
//         );
//         GcPcSaftFunctionalParameters::from_segments(vec![pure], vec![ch3(), ch2(), oh()], None)
//             .unwrap()
//     }

//     pub fn ethanol_propanol(binary: bool) -> GcPcSaftFunctionalParameters {
//         let ethanol = GroupContributionRecord::new(
//             Identifier::new("64-17-5", Some("ethanol"), None, None, None, None),
//             0.0,
//             Some(ChemicalRecord::new(
//                 vec!["CH3".into(), "CH2".into(), "OH".into()],
//                 None,
//             )),
//             None,
//             None,
//         );
//         let propanol = GroupContributionRecord::new(
//             Identifier::new("71-23-8", Some("1-propanol"), None, None, None, None),
//             0.0,
//             Some(ChemicalRecord::new(
//                 vec!["CH3".into(), "CH2".into(), "CH2".into(), "OH".into()],
//                 None,
//             )),
//             None,
//             None,
//         );
//         let binary = if binary { Some(vec![ch3_oh()]) } else { None };
//         GcPcSaftFunctionalParameters::from_segments(
//             vec![ethanol, propanol],
//             vec![ch3(), ch2(), oh()],
//             binary,
//         )
//         .unwrap()
//     }

//     #[test]
//     fn test_kij() {
//         let params = ethanol_propanol(true);
//         println!("{}", params.epsilon_k_ij);
//         // CH3 - CH2
//         assert_eq!(params.epsilon_k_ij[(0, 4)], (181.49f64 * 157.23).sqrt());
//         // CH3 - OH
//         assert_eq!(
//             params.epsilon_k_ij[(0, 5)],
//             (181.49f64 * 334.29).sqrt() * 1.0087
//         );
//     }
// }
