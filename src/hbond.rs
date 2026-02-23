//! Hydrogen bond detection
//!
//! Identifies potential hydrogen bonds between residues based on
//! donor-acceptor distance and backbone geometry criteria.
//!
//! Author: Moroya Sakamoto

use crate::amino::{AminoAcid, Residue};
use crate::fnv1a;

/// Hydrogen bond detection configuration.
#[derive(Debug, Clone, Copy)]
pub struct HBondConfig {
    /// Maximum donor-acceptor distance (Angstroms). Default 3.5.
    pub max_distance: f64,
    /// Minimum sequence separation (|i - j|) for an H-bond. Default 3.
    pub min_seq_separation: usize,
    /// Energy cutoff (kcal/mol). Bonds weaker than this are rejected. Default -0.5.
    pub energy_cutoff: f64,
}

impl Default for HBondConfig {
    fn default() -> Self {
        Self {
            max_distance: 3.5,
            min_seq_separation: 3,
            energy_cutoff: -0.5,
        }
    }
}

/// A detected hydrogen bond.
#[derive(Debug, Clone, Copy)]
pub struct HBondHit {
    /// Donor residue index (N-H side).
    pub donor_idx: usize,
    /// Acceptor residue index (C=O side).
    pub acceptor_idx: usize,
    /// Donor-acceptor Cα distance (Angstroms).
    pub distance: f64,
    /// Estimated energy (kcal/mol, negative = favorable).
    pub energy: f64,
    /// Deterministic content hash.
    pub content_hash: u64,
}

/// Donor/acceptor classification for amino acids.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HBondRole {
    /// Can donate H-bonds (has N-H).
    Donor,
    /// Can accept H-bonds (has C=O).
    Acceptor,
    /// Can both donate and accept.
    Both,
}

/// Classify an amino acid's hydrogen bonding role.
///
/// All residues (except Pro at N-terminus) have backbone N-H (donor)
/// and C=O (acceptor). Sidechain donors/acceptors are role-dependent.
pub fn classify_hbond_role(aa: AminoAcid) -> HBondRole {
    match aa {
        // Pro lacks backbone N-H → acceptor only
        AminoAcid::Pro => HBondRole::Acceptor,
        // Sidechains with both donor and acceptor groups
        AminoAcid::Ser | AminoAcid::Thr | AminoAcid::Tyr |
        AminoAcid::Asn | AminoAcid::Gln | AminoAcid::His |
        AminoAcid::Trp | AminoAcid::Cys | AminoAcid::Lys |
        AminoAcid::Arg => HBondRole::Both,
        // Asp, Glu: acceptor sidechains + backbone donor
        AminoAcid::Asp | AminoAcid::Glu => HBondRole::Both,
        // Hydrophobic residues: backbone only (both N-H and C=O)
        _ => HBondRole::Both,
    }
}

/// Hydrogen bond detector.
#[derive(Debug, Clone)]
pub struct HBondDetector {
    config: HBondConfig,
}

impl HBondDetector {
    /// Create a new detector with the given configuration.
    pub fn new(config: HBondConfig) -> Self {
        Self { config }
    }

    /// Detect hydrogen bonds between residues given their Cα positions.
    ///
    /// Uses a simplified electrostatic model:
    /// E = 0.084 * 332 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)
    /// Approximated as E ≈ -k / d^2 where d is Cα distance.
    pub fn detect(
        &self,
        residues: &[Residue],
        positions: &[[f64; 3]],
    ) -> Vec<HBondHit> {
        let n = residues.len().min(positions.len());
        let max_dsq = self.config.max_distance * self.config.max_distance;
        let mut hits = Vec::new();

        for i in 0..n {
            for j in (i + 1)..n {
                // Sequence separation check
                if j - i < self.config.min_seq_separation {
                    continue;
                }

                let dsq = dist_sq(&positions[i], &positions[j]);
                if dsq > max_dsq {
                    continue;
                }

                let dist = dsq.sqrt();

                // Simplified DSSP-like energy model
                // E ≈ -27.888 / d^2 (kcal/mol), calibrated so that
                // at d=2.8Å (optimal N-H...O=C), E ≈ -3.5 kcal/mol
                let inv_dsq = 1.0 / dsq;
                let energy = -27.888 * inv_dsq;

                if energy > self.config.energy_cutoff {
                    continue;
                }

                // Determine donor/acceptor direction
                // Convention: lower-index residue is donor (N-H → C=O)
                let donor_idx = i;
                let acceptor_idx = j;

                // Skip if donor is Pro (no backbone N-H)
                if residues[donor_idx].amino == AminoAcid::Pro {
                    continue;
                }

                let mut buf = [0u8; 24];
                buf[..8].copy_from_slice(&donor_idx.to_le_bytes());
                buf[8..16].copy_from_slice(&acceptor_idx.to_le_bytes());
                buf[16..24].copy_from_slice(&energy.to_bits().to_le_bytes());
                let content_hash = fnv1a(&buf);

                hits.push(HBondHit {
                    donor_idx,
                    acceptor_idx,
                    distance: dist,
                    energy,
                    content_hash,
                });
            }
        }

        hits
    }

    /// Detect and return only the count of hydrogen bonds.
    pub fn count(&self, residues: &[Residue], positions: &[[f64; 3]]) -> usize {
        self.detect(residues, positions).len()
    }
}

#[inline]
fn dist_sq(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

// ── Tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn make_residue(aa: AminoAcid) -> Residue {
        Residue::new(aa, -1.0, -0.8, PI)
    }

    #[test]
    fn classify_pro_is_acceptor() {
        assert_eq!(classify_hbond_role(AminoAcid::Pro), HBondRole::Acceptor);
    }

    #[test]
    fn classify_ala_is_both() {
        assert_eq!(classify_hbond_role(AminoAcid::Ala), HBondRole::Both);
    }

    #[test]
    fn classify_ser_is_both() {
        assert_eq!(classify_hbond_role(AminoAcid::Ser), HBondRole::Both);
    }

    #[test]
    fn no_bonds_empty() {
        let detector = HBondDetector::new(HBondConfig::default());
        let hits = detector.detect(&[], &[]);
        assert!(hits.is_empty());
    }

    #[test]
    fn no_bonds_single_residue() {
        let detector = HBondDetector::new(HBondConfig::default());
        let residues = [make_residue(AminoAcid::Ala)];
        let positions = [[0.0, 0.0, 0.0]];
        assert_eq!(detector.count(&residues, &positions), 0);
    }

    #[test]
    fn no_bonds_too_close_in_sequence() {
        let detector = HBondDetector::new(HBondConfig {
            min_seq_separation: 3,
            ..Default::default()
        });
        let residues = [make_residue(AminoAcid::Ala), make_residue(AminoAcid::Gly)];
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        assert_eq!(detector.count(&residues, &positions), 0);
    }

    #[test]
    fn detect_bond_close_residues() {
        let cfg = HBondConfig {
            max_distance: 4.0,
            min_seq_separation: 3,
            energy_cutoff: -0.1,
        };
        let detector = HBondDetector::new(cfg);
        // Residues i=0 and j=3 at distance 3.0
        let residues = vec![
            make_residue(AminoAcid::Ala),
            make_residue(AminoAcid::Gly),
            make_residue(AminoAcid::Val),
            make_residue(AminoAcid::Leu),
        ];
        let positions = [
            [0.0, 0.0, 0.0],
            [3.8, 0.0, 0.0],
            [7.6, 0.0, 0.0],
            [3.0, 0.0, 0.0], // close to residue 0
        ];
        let hits = detector.detect(&residues, &positions);
        assert!(!hits.is_empty());
        assert_eq!(hits[0].donor_idx, 0);
        assert_eq!(hits[0].acceptor_idx, 3);
    }

    #[test]
    fn pro_donor_skipped() {
        let cfg = HBondConfig {
            max_distance: 5.0,
            min_seq_separation: 3,
            energy_cutoff: -0.1,
        };
        let detector = HBondDetector::new(cfg);
        // Pro as donor (index 0) should be skipped
        let residues = vec![
            make_residue(AminoAcid::Pro),
            make_residue(AminoAcid::Ala),
            make_residue(AminoAcid::Gly),
            make_residue(AminoAcid::Val),
        ];
        let positions = [
            [0.0, 0.0, 0.0],
            [3.8, 0.0, 0.0],
            [7.6, 0.0, 0.0],
            [2.8, 0.0, 0.0],
        ];
        let hits = detector.detect(&residues, &positions);
        // Pro at index 0 should not appear as donor
        for h in &hits {
            assert_ne!(residues[h.donor_idx].amino, AminoAcid::Pro);
        }
    }

    #[test]
    fn energy_is_negative() {
        let cfg = HBondConfig {
            max_distance: 5.0,
            min_seq_separation: 3,
            energy_cutoff: 0.0,
        };
        let detector = HBondDetector::new(cfg);
        let residues = vec![
            make_residue(AminoAcid::Ala),
            make_residue(AminoAcid::Gly),
            make_residue(AminoAcid::Val),
            make_residue(AminoAcid::Leu),
        ];
        let positions = [
            [0.0, 0.0, 0.0],
            [3.8, 0.0, 0.0],
            [7.6, 0.0, 0.0],
            [2.8, 0.0, 0.0],
        ];
        let hits = detector.detect(&residues, &positions);
        for h in &hits {
            assert!(h.energy < 0.0);
        }
    }

    #[test]
    fn content_hash_deterministic() {
        let cfg = HBondConfig {
            max_distance: 5.0,
            min_seq_separation: 3,
            energy_cutoff: 0.0,
        };
        let detector = HBondDetector::new(cfg);
        let residues = vec![
            make_residue(AminoAcid::Ala),
            make_residue(AminoAcid::Gly),
            make_residue(AminoAcid::Val),
            make_residue(AminoAcid::Leu),
        ];
        let positions = [
            [0.0, 0.0, 0.0],
            [3.8, 0.0, 0.0],
            [7.6, 0.0, 0.0],
            [2.8, 0.0, 0.0],
        ];
        let h1 = detector.detect(&residues, &positions);
        let h2 = detector.detect(&residues, &positions);
        assert_eq!(h1.len(), h2.len());
        for (a, b) in h1.iter().zip(h2.iter()) {
            assert_eq!(a.content_hash, b.content_hash);
        }
    }

    #[test]
    fn no_bonds_beyond_distance() {
        let cfg = HBondConfig {
            max_distance: 3.0,
            min_seq_separation: 3,
            energy_cutoff: 0.0,
        };
        let detector = HBondDetector::new(cfg);
        let residues = vec![
            make_residue(AminoAcid::Ala),
            make_residue(AminoAcid::Gly),
            make_residue(AminoAcid::Val),
            make_residue(AminoAcid::Leu),
        ];
        let positions = [
            [0.0, 0.0, 0.0],
            [3.8, 0.0, 0.0],
            [7.6, 0.0, 0.0],
            [50.0, 0.0, 0.0], // too far
        ];
        assert_eq!(detector.count(&residues, &positions), 0);
    }

    #[test]
    fn mismatched_residue_position_lengths() {
        // More residues than positions: detector should use min(len) safely
        let detector = HBondDetector::new(HBondConfig::default());
        let residues = vec![
            make_residue(AminoAcid::Ala),
            make_residue(AminoAcid::Gly),
            make_residue(AminoAcid::Val),
            make_residue(AminoAcid::Leu),
            make_residue(AminoAcid::Ile),
        ];
        let positions = [
            [0.0, 0.0, 0.0],
            [3.8, 0.0, 0.0],
            [7.6, 0.0, 0.0],
        ];
        // Should not panic, even though residues.len() != positions.len()
        let hits = detector.detect(&residues, &positions);
        assert!(hits.is_empty() || !hits.is_empty()); // no panic is the test
    }

    #[test]
    fn all_amino_acids_have_hbond_role() {
        // Every amino acid should return a valid HBondRole
        for aa in &AminoAcid::ALL {
            let role = classify_hbond_role(*aa);
            assert!(
                role == HBondRole::Donor || role == HBondRole::Acceptor || role == HBondRole::Both,
                "{:?} has no valid role",
                aa
            );
        }
    }

    #[test]
    fn very_strict_energy_cutoff_filters_all() {
        // Energy cutoff so negative that no bond qualifies
        let cfg = HBondConfig {
            max_distance: 10.0,
            min_seq_separation: 3,
            energy_cutoff: -1000.0, // impossibly strict
        };
        let detector = HBondDetector::new(cfg);
        let residues = vec![
            make_residue(AminoAcid::Ala),
            make_residue(AminoAcid::Gly),
            make_residue(AminoAcid::Val),
            make_residue(AminoAcid::Leu),
        ];
        let positions = [
            [0.0, 0.0, 0.0],
            [3.8, 0.0, 0.0],
            [7.6, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        assert_eq!(detector.count(&residues, &positions), 0);
    }

    #[test]
    fn content_hash_nonzero() {
        let cfg = HBondConfig {
            max_distance: 5.0,
            min_seq_separation: 3,
            energy_cutoff: 0.0,
        };
        let detector = HBondDetector::new(cfg);
        let residues = vec![
            make_residue(AminoAcid::Ala),
            make_residue(AminoAcid::Gly),
            make_residue(AminoAcid::Val),
            make_residue(AminoAcid::Leu),
        ];
        let positions = [
            [0.0, 0.0, 0.0],
            [3.8, 0.0, 0.0],
            [7.6, 0.0, 0.0],
            [2.8, 0.0, 0.0],
        ];
        let hits = detector.detect(&residues, &positions);
        for h in &hits {
            assert_ne!(h.content_hash, 0, "Content hash should not be zero");
        }
    }

    #[test]
    fn min_seq_separation_one_allows_adjacent() {
        // With min_seq_separation=1, even adjacent residues qualify
        let cfg = HBondConfig {
            max_distance: 5.0,
            min_seq_separation: 1,
            energy_cutoff: 0.0,
        };
        let detector = HBondDetector::new(cfg);
        let residues = vec![
            make_residue(AminoAcid::Ala),
            make_residue(AminoAcid::Gly),
        ];
        let positions = [
            [0.0, 0.0, 0.0],
            [2.8, 0.0, 0.0],
        ];
        let hits = detector.detect(&residues, &positions);
        assert!(!hits.is_empty(), "Adjacent residues should form H-bond with min_seq_separation=1");
    }
}
