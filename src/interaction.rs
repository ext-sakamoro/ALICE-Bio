//! Molecular interaction evaluation — pairwise energies, contact maps.

use crate::amino::Residue;
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
pub fn evaluate_pairwise_energy(residues: &[Residue], positions: &[[f64; 3]]) -> TotalEnergy {
    let mut vdw = 0.0;
    let mut torsional = 0.0;

    // Pre-compute per-residue VdW radii to avoid repeated match dispatch in inner loop.
    let radii: Vec<f64> = residues.iter().map(|r| r.amino.van_der_waals_radius()).collect();

    // Pairwise van der Waals — use dist_sq to avoid sqrt in O(N²) hot path.
    // (σ/r)^n = (σ²/r²)^(n/2), so LJ can be fully expressed in r² without sqrt.
    const R_MIN_SQ: f64 = 0.1 * 0.1; // threshold = 0.1 Å, squared
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            let dx = positions[j][0] - positions[i][0];
            let dy = positions[j][1] - positions[i][1];
            let dz = positions[j][2] - positions[i][2];
            let dist_sq = dx * dx + dy * dy + dz * dz;
            if dist_sq > R_MIN_SQ {
                // Pre-compute reciprocal of 2 once; sigma = (ri + rj) * 0.5
                let sigma = (radii[i] + radii[j]) * 0.5;
                vdw += lennard_jones_sq(dist_sq, 0.1, sigma);
            }
        }
    }

    // Torsional energy from backbone angles
    for res in residues {
        torsional += torsion_potential(res.phi, 1.0, 3, 0.0);
        torsional += torsion_potential(res.psi, 0.5, 3, 0.0);
    }

    TotalEnergy {
        van_der_waals: vdw,
        electrostatic: 0.0,
        hydrogen_bonds: 0.0,
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
}
