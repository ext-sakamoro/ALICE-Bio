//! Molecular folding as signed distance functions.

use crate::amino::Residue;

/// Cα-Cα distance in Angstroms.
const CA_DISTANCE: f64 = 3.8;

/// Molecular structure represented as an SDF over atom spheres.
#[derive(Debug, Clone)]
pub struct ProteinSdf {
    residues: Vec<Residue>,
    /// Alpha-carbon positions computed from backbone angles.
    positions: Vec<[f64; 3]>,
    /// Pre-computed per-residue VdW radii — avoids repeated match dispatch in eval hot path.
    radii: Vec<f64>,
}

impl ProteinSdf {
    /// Build from residues, computing Cα positions from phi/psi chain geometry.
    pub fn new(residues: Vec<Residue>) -> Self {
        let positions = compute_positions(&residues);
        let radii = residues.iter().map(|r| r.amino.van_der_waals_radius()).collect();
        Self { residues, positions, radii }
    }

    /// Signed distance: minimum over all atom spheres.
    /// Negative = inside an atom VdW radius, positive = outside.
    #[inline]
    pub fn eval(&self, point: &[f64; 3]) -> f64 {
        let mut min_dist = f64::MAX;
        // min_dist is a signed distance; convert to squared threshold for early culling.
        // Any atom whose squared distance exceeds (min_dist + max_vdw)² can never beat
        // the current best. We tighten the bound lazily as min_dist decreases.
        for (i, pos) in self.positions.iter().enumerate() {
            let dx = point[0] - pos[0];
            let dy = point[1] - pos[1];
            let dz = point[2] - pos[2];
            let dist_sq = dx * dx + dy * dy + dz * dz;
            let r = self.radii[i];
            // Early-exit: if the closest possible signed distance from this atom
            // (dist - r) cannot beat min_dist, skip sqrt entirely.
            // dist - r < min_dist  <=>  dist < min_dist + r  <=>  dist_sq < (min_dist + r)²
            // Only valid when min_dist + r > 0 (otherwise everything is inside a sphere).
            let bound = min_dist + r;
            if bound > 0.0 && dist_sq >= bound * bound {
                continue;
            }
            let dist = dist_sq.sqrt();
            let sd = dist - r;
            if sd < min_dist {
                min_dist = sd;
            }
        }
        min_dist
    }

    /// Batch evaluation.
    pub fn eval_batch(&self, points: &[[f64; 3]]) -> Vec<f64> {
        points.iter().map(|p| self.eval(p)).collect()
    }

    /// Axis-aligned bounding box with VdW padding.
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        if self.positions.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let max_r = self.residues.iter()
            .map(|r| r.amino.van_der_waals_radius())
            .fold(0.0_f64, f64::max);

        let mut min = [f64::MAX; 3];
        let mut max = [f64::MIN; 3];
        for pos in &self.positions {
            for i in 0..3 {
                if pos[i] < min[i] { min[i] = pos[i]; }
                if pos[i] > max[i] { max[i] = pos[i]; }
            }
        }
        for i in 0..3 {
            min[i] -= max_r;
            max[i] += max_r;
        }
        (min, max)
    }

    /// Number of residues.
    #[inline]
    pub fn residue_count(&self) -> usize {
        self.residues.len()
    }

    /// Access Cα positions.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// Access residues.
    pub fn residues(&self) -> &[Residue] {
        &self.residues
    }
}

/// Compute Cα positions from backbone angles using simplified chain geometry.
fn compute_positions(residues: &[Residue]) -> Vec<[f64; 3]> {
    let mut positions = Vec::with_capacity(residues.len());
    if residues.is_empty() {
        return positions;
    }

    positions.push([0.0, 0.0, 0.0]);

    let mut dir_theta = 0.0_f64; // accumulated direction in XY plane
    let mut dir_phi = 0.0_f64;   // accumulated direction in Z

    for i in 1..residues.len() {
        // Update direction from previous residue's phi/psi angles
        dir_theta += residues[i - 1].phi * 0.3;
        dir_phi += residues[i - 1].psi * 0.2;

        let prev = &positions[i - 1];
        let dx = CA_DISTANCE * dir_theta.cos() * dir_phi.cos();
        let dy = CA_DISTANCE * dir_theta.sin() * dir_phi.cos();
        let dz = CA_DISTANCE * dir_phi.sin();
        positions.push([prev[0] + dx, prev[1] + dy, prev[2] + dz]);
    }

    positions
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amino::AminoAcid;
    use std::f64::consts::PI;

    fn make_residue(aa: AminoAcid) -> Residue {
        Residue::new(aa, -1.0, -0.8, PI)
    }

    #[test]
    fn single_residue_inside() {
        let sdf = ProteinSdf::new(vec![make_residue(AminoAcid::Ala)]);
        // At origin (center of first atom), distance should be negative (inside)
        let d = sdf.eval(&[0.0, 0.0, 0.0]);
        assert!(d < 0.0);
    }

    #[test]
    fn single_residue_outside() {
        let sdf = ProteinSdf::new(vec![make_residue(AminoAcid::Gly)]);
        // Far away should be positive
        let d = sdf.eval(&[100.0, 100.0, 100.0]);
        assert!(d > 0.0);
    }

    #[test]
    fn multi_residue_spacing() {
        let residues = vec![
            Residue::new(AminoAcid::Ala, 0.0, 0.0, PI),
            Residue::new(AminoAcid::Gly, 0.0, 0.0, PI),
        ];
        let sdf = ProteinSdf::new(residues);
        let pos = sdf.positions();
        let dx = pos[1][0] - pos[0][0];
        let dy = pos[1][1] - pos[0][1];
        let dz = pos[1][2] - pos[0][2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!((dist - CA_DISTANCE).abs() < 0.01);
    }

    #[test]
    fn batch_matches_individual() {
        let sdf = ProteinSdf::new(vec![make_residue(AminoAcid::Ala), make_residue(AminoAcid::Val)]);
        let points = [[0.0, 0.0, 0.0], [5.0, 5.0, 5.0], [10.0, 0.0, 0.0]];
        let batch = sdf.eval_batch(&points);
        for (i, p) in points.iter().enumerate() {
            assert!((batch[i] - sdf.eval(p)).abs() < 1e-12);
        }
    }

    #[test]
    fn bounding_box_contains_positions() {
        let residues = vec![make_residue(AminoAcid::Trp); 5];
        let sdf = ProteinSdf::new(residues);
        let (bb_min, bb_max) = sdf.bounding_box();
        for pos in sdf.positions() {
            for i in 0..3 {
                assert!(pos[i] >= bb_min[i] && pos[i] <= bb_max[i]);
            }
        }
    }

    #[test]
    fn residue_count() {
        let sdf = ProteinSdf::new(vec![make_residue(AminoAcid::Ala); 10]);
        assert_eq!(sdf.residue_count(), 10);
    }

    #[test]
    fn empty_protein() {
        let sdf = ProteinSdf::new(vec![]);
        assert_eq!(sdf.residue_count(), 0);
        assert_eq!(sdf.eval(&[0.0, 0.0, 0.0]), f64::MAX);
    }
}
