//! Molecular interaction evaluation — pairwise energies, contact maps.

use crate::amino::Residue;
use crate::cell_list::{CellList, CellListConfig};
use crate::hbond::{HBondConfig, HBondDetector};
use crate::potential::{torsion_potential, TotalEnergy};

/// Lennard-Jones evaluated directly from dist_sq, avoiding sqrt.
/// V = 4ε[(σ/r)^12 - (σ/r)^6] = 4ε[(σ²/r²)^6 - (σ²/r²)^3]
#[inline(always)]
fn lennard_jones_sq(dist_sq: f64, epsilon: f64, sigma: f64) -> f64 {
    let sigma_sq = sigma * sigma;
    let inv_r2 = sigma_sq / dist_sq;
    let inv_r6 = inv_r2 * inv_r2 * inv_r2;
    4.0 * epsilon * (inv_r6 * inv_r6 - inv_r6)
}

/// Evaluate pairwise interaction energy for all residue pairs.
///
/// For small systems (N < 50), uses O(N²) brute-force. For larger systems,
/// uses a cell-linked list for O(N) neighbor search. Hydrogen bonds are
/// evaluated via `HBondDetector`.
pub fn evaluate_pairwise_energy(residues: &[Residue], positions: &[[f64; 3]]) -> TotalEnergy {
    let mut vdw = 0.0;
    let mut torsional = 0.0;

    // Pre-compute per-residue VdW radii to avoid repeated match dispatch in inner loop.
    let radii: Vec<f64> = residues.iter().map(|r| r.amino.van_der_waals_radius()).collect();

    const R_MIN_SQ: f64 = 0.1 * 0.1; // threshold = 0.1 Å, squared
    const CELL_LIST_THRESHOLD: usize = 50;
    const VDW_CUTOFF: f64 = 12.0; // Angstroms

    if positions.len() >= CELL_LIST_THRESHOLD {
        // Cell-linked list path: O(N) for uniform distributions.
        let cfg = CellListConfig { cutoff: VDW_CUTOFF, origin: None };
        let cl = CellList::build(positions, &cfg);
        let pairs = cl.find_pairs(positions);
        for &(i, j, dist_sq) in &pairs {
            if dist_sq > R_MIN_SQ {
                let sigma = (radii[i] + radii[j]) * 0.5;
                vdw += lennard_jones_sq(dist_sq, 0.1, sigma);
            }
        }
    } else {
        // Brute-force path for small systems.
        for i in 0..positions.len() {
            for j in (i + 1)..positions.len() {
                let dx = positions[j][0] - positions[i][0];
                let dy = positions[j][1] - positions[i][1];
                let dz = positions[j][2] - positions[i][2];
                let dist_sq = dx * dx + dy * dy + dz * dz;
                if dist_sq > R_MIN_SQ {
                    let sigma = (radii[i] + radii[j]) * 0.5;
                    vdw += lennard_jones_sq(dist_sq, 0.1, sigma);
                }
            }
        }
    }

    // Torsional energy from backbone angles
    for res in residues {
        torsional += torsion_potential(res.phi, 1.0, 3, 0.0);
        torsional += torsion_potential(res.psi, 0.5, 3, 0.0);
    }

    // Hydrogen bond evaluation via HBondDetector
    let hbond_detector = HBondDetector::new(HBondConfig::default());
    let hbond_hits = hbond_detector.detect(residues, positions);
    let hbond_energy: f64 = hbond_hits.iter().map(|h| h.energy).sum();

    TotalEnergy {
        van_der_waals: vdw,
        electrostatic: 0.0,
        hydrogen_bonds: hbond_energy,
        torsional,
    }
}

/// Returns pairs of residue indices whose Cα atoms are within `threshold` Angstroms.
pub fn contact_map(positions: &[[f64; 3]], threshold: f64) -> Vec<(usize, usize)> {
    let mut contacts = Vec::new();
    let t2 = threshold * threshold;
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            let dx = positions[j][0] - positions[i][0];
            let dy = positions[j][1] - positions[i][1];
            let dz = positions[j][2] - positions[i][2];
            if dx * dx + dy * dy + dz * dz <= t2 {
                contacts.push((i, j));
            }
        }
    }
    contacts
}

/// Radius of gyration: Rg = sqrt(mean(|r_i - r_center|²))
pub fn radius_of_gyration(positions: &[[f64; 3]]) -> f64 {
    if positions.is_empty() {
        return 0.0;
    }
    let n = positions.len() as f64;
    // Pre-compute reciprocal: replace 4 divisions (3 center + 1 mean) with multiplications.
    let inv_n = 1.0 / n;
    let mut center = [0.0; 3];
    for p in positions {
        center[0] += p[0];
        center[1] += p[1];
        center[2] += p[2];
    }
    center[0] *= inv_n;
    center[1] *= inv_n;
    center[2] *= inv_n;

    let mut sum_sq = 0.0;
    for p in positions {
        let dx = p[0] - center[0];
        let dy = p[1] - center[1];
        let dz = p[2] - center[2];
        sum_sq += dx * dx + dy * dy + dz * dz;
    }
    (sum_sq * inv_n).sqrt()
}

/// Euclidean distance between the first and last Cα.
pub fn end_to_end_distance(positions: &[[f64; 3]]) -> f64 {
    if positions.len() < 2 {
        return 0.0;
    }
    let first = &positions[0];
    let last = &positions[positions.len() - 1];
    let dx = last[0] - first[0];
    let dy = last[1] - first[1];
    let dz = last[2] - first[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amino::{AminoAcid, Residue};
    use std::f64::consts::PI;

    #[test]
    fn single_residue_no_pairwise_vdw() {
        let residues = vec![Residue::new(AminoAcid::Ala, -1.0, -0.8, PI)];
        let positions = [[0.0, 0.0, 0.0]];
        let e = evaluate_pairwise_energy(&residues, &positions);
        assert_eq!(e.van_der_waals, 0.0);
    }

    #[test]
    fn two_residues_nonzero_vdw() {
        let residues = vec![
            Residue::new(AminoAcid::Ala, 0.0, 0.0, PI),
            Residue::new(AminoAcid::Gly, 0.0, 0.0, PI),
        ];
        let positions = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let e = evaluate_pairwise_energy(&residues, &positions);
        assert!(e.van_der_waals != 0.0);
    }

    #[test]
    fn contact_map_within_threshold() {
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
        let contacts = contact_map(&positions, 5.0);
        assert!(contacts.contains(&(0, 1)));
        assert!(!contacts.contains(&(0, 2)));
    }

    #[test]
    fn rg_single_point_zero() {
        assert_eq!(radius_of_gyration(&[[0.0, 0.0, 0.0]]), 0.0);
    }

    #[test]
    fn rg_two_points() {
        let positions = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let rg = radius_of_gyration(&positions);
        // Center at (5,0,0), each at distance 5, Rg = 5
        assert!((rg - 5.0).abs() < 1e-10);
    }

    #[test]
    fn end_to_end_two_residues() {
        let positions = [[0.0, 0.0, 0.0], [3.8, 0.0, 0.0]];
        let d = end_to_end_distance(&positions);
        assert!((d - 3.8).abs() < 1e-10);
    }

    #[test]
    fn end_to_end_single_residue() {
        assert_eq!(end_to_end_distance(&[[0.0, 0.0, 0.0]]), 0.0);
    }

    #[test]
    fn hbond_energy_nonzero_with_separated_residues() {
        // 4+ residues with sufficient sequence separation and close distance
        let residues = vec![
            Residue::new(AminoAcid::Ala, -1.0, -0.8, PI),
            Residue::new(AminoAcid::Gly, -1.0, -0.8, PI),
            Residue::new(AminoAcid::Val, -1.0, -0.8, PI),
            Residue::new(AminoAcid::Leu, -1.0, -0.8, PI),
        ];
        // Residue 0 and 3 are within H-bond distance (3.0 Å)
        let positions = [
            [0.0, 0.0, 0.0],
            [3.8, 0.0, 0.0],
            [7.6, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let e = evaluate_pairwise_energy(&residues, &positions);
        assert!(e.hydrogen_bonds < 0.0, "H-bond energy should be negative, got {}", e.hydrogen_bonds);
    }

    #[test]
    fn cell_list_path_large_system() {
        // Create 60 residues to trigger cell-list path (threshold = 50)
        let n = 60;
        let residues: Vec<Residue> = (0..n)
            .map(|_| Residue::new(AminoAcid::Ala, -1.0, -0.8, PI))
            .collect();
        let positions: Vec<[f64; 3]> = (0..n)
            .map(|i| [i as f64 * 3.8, 0.0, 0.0])
            .collect();
        let e = evaluate_pairwise_energy(&residues, &positions);
        assert!(e.van_der_waals.is_finite());
        assert!(e.torsional.is_finite());
    }

    #[test]
    fn contact_map_empty_positions() {
        let contacts = contact_map(&[], 5.0);
        assert!(contacts.is_empty());
    }

    #[test]
    fn contact_map_single_atom() {
        let contacts = contact_map(&[[0.0, 0.0, 0.0]], 5.0);
        assert!(contacts.is_empty());
    }

    #[test]
    fn contact_map_all_within_threshold() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let contacts = contact_map(&positions, 10.0);
        // All 3 pairs should be within threshold
        assert_eq!(contacts.len(), 3);
    }

    #[test]
    fn rg_empty_returns_zero() {
        assert_eq!(radius_of_gyration(&[]), 0.0);
    }

    #[test]
    fn rg_symmetric_cube_vertices() {
        // 8 vertices of a cube centered at origin with side length 2
        let positions = [
            [-1.0, -1.0, -1.0], [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],  [-1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],  [1.0, -1.0, 1.0],
            [1.0, 1.0, -1.0],   [1.0, 1.0, 1.0],
        ];
        let rg = radius_of_gyration(&positions);
        // Each vertex is at distance sqrt(3) from center, Rg = sqrt(3)
        let expected = 3.0_f64.sqrt();
        assert!((rg - expected).abs() < 1e-10, "Rg = {}, expected {}", rg, expected);
    }

    #[test]
    fn end_to_end_empty_returns_zero() {
        assert_eq!(end_to_end_distance(&[]), 0.0);
    }

    #[test]
    fn end_to_end_3d_diagonal() {
        let positions = [[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]];
        let d = end_to_end_distance(&positions);
        assert!((d - 5.0).abs() < 1e-10);
    }

    #[test]
    fn evaluate_pairwise_torsional_nonzero() {
        // Even with widely separated positions, torsional energy from backbone angles is nonzero
        let residues = vec![
            Residue::new(AminoAcid::Ala, -1.0, -0.8, PI),
            Residue::new(AminoAcid::Gly, 0.5, 0.3, PI),
        ];
        let positions = [[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]]; // far apart, no VdW
        let e = evaluate_pairwise_energy(&residues, &positions);
        assert!(e.torsional > 0.0, "Torsional energy should be positive, got {}", e.torsional);
    }

    #[test]
    fn evaluate_pairwise_total_is_sum_of_components() {
        let residues = vec![
            Residue::new(AminoAcid::Ala, -1.0, -0.8, PI),
            Residue::new(AminoAcid::Val, -0.5, -0.3, PI),
        ];
        let positions = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let e = evaluate_pairwise_energy(&residues, &positions);
        let expected = e.van_der_waals + e.electrostatic + e.hydrogen_bonds + e.torsional;
        assert!((e.total() - expected).abs() < 1e-15);
    }
}
