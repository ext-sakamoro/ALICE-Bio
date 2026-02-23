//! Simplified DSSP-like secondary structure assignment
//!
//! Assigns secondary structure (alpha-helix, beta-sheet, coil) to each
//! residue based on backbone dihedral angles (phi/psi) and hydrogen bond
//! patterns from `HBondDetector`.
//!
//! Author: Moroya Sakamoto

use crate::amino::Residue;
use crate::hbond::{HBondConfig, HBondDetector, HBondHit};
use crate::fnv1a;

/// Secondary structure type for a single residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SecondaryStructure {
    /// Alpha helix (i → i+4 hydrogen bond pattern).
    Helix,
    /// Beta strand (extended, i → i±≥5 hydrogen bond pattern).
    Sheet,
    /// Coil / loop / turn (neither helix nor sheet).
    Coil,
}

/// Configuration for secondary structure assignment.
#[derive(Debug, Clone, Copy)]
pub struct SecondaryConfig {
    /// Use hydrogen bond patterns for assignment (requires positions).
    pub use_hbonds: bool,
    /// Alpha helix phi center (radians). Default -1.047 (~-60°).
    pub helix_phi_center: f64,
    /// Alpha helix psi center (radians). Default -0.785 (~-45°).
    pub helix_psi_center: f64,
    /// Alpha helix angular tolerance (radians). Default 0.524 (~30°).
    pub helix_tolerance: f64,
    /// Beta sheet phi center (radians). Default -2.094 (~-120°).
    pub sheet_phi_center: f64,
    /// Beta sheet psi center (radians). Default 2.269 (~130°).
    pub sheet_psi_center: f64,
    /// Beta sheet angular tolerance (radians). Default 0.524 (~30°).
    pub sheet_tolerance: f64,
    /// Minimum consecutive helix residues to confirm helix. Default 3.
    pub min_helix_run: usize,
    /// Minimum consecutive sheet residues to confirm sheet. Default 2.
    pub min_sheet_run: usize,
}

impl Default for SecondaryConfig {
    fn default() -> Self {
        Self {
            use_hbonds: false,
            helix_phi_center: -1.047,  // -60°
            helix_psi_center: -0.785,  // -45°
            helix_tolerance: 0.524,    // 30°
            sheet_phi_center: -2.094,  // -120°
            sheet_psi_center: 2.269,   // 130°
            sheet_tolerance: 0.524,    // 30°
            min_helix_run: 3,
            min_sheet_run: 2,
        }
    }
}

/// Result of secondary structure assignment for a chain.
#[derive(Debug, Clone)]
pub struct SecondaryAssignment {
    /// Per-residue structure assignment.
    pub assignments: Vec<SecondaryStructure>,
    /// Number of helix residues.
    pub helix_count: usize,
    /// Number of sheet residues.
    pub sheet_count: usize,
    /// Number of coil residues.
    pub coil_count: usize,
    /// Helix fraction (0.0 to 1.0).
    pub helix_fraction: f64,
    /// Sheet fraction (0.0 to 1.0).
    pub sheet_fraction: f64,
    /// Deterministic content hash.
    pub content_hash: u64,
}

/// Check if a (phi, psi) pair falls within an angular region.
#[inline]
fn in_region(phi: f64, psi: f64, phi_c: f64, psi_c: f64, tol: f64) -> bool {
    let dphi = angle_diff(phi, phi_c);
    let dpsi = angle_diff(psi, psi_c);
    dphi.abs() <= tol && dpsi.abs() <= tol
}

/// Signed angular difference, normalized to [-PI, PI].
#[inline]
fn angle_diff(a: f64, b: f64) -> f64 {
    let mut d = a - b;
    while d > std::f64::consts::PI {
        d -= 2.0 * std::f64::consts::PI;
    }
    while d < -std::f64::consts::PI {
        d += 2.0 * std::f64::consts::PI;
    }
    d
}

/// Classify a single residue by its dihedral angles (phi/psi).
#[inline]
fn classify_by_angles(residue: &Residue, config: &SecondaryConfig) -> SecondaryStructure {
    if in_region(
        residue.phi, residue.psi,
        config.helix_phi_center, config.helix_psi_center,
        config.helix_tolerance,
    ) {
        SecondaryStructure::Helix
    } else if in_region(
        residue.phi, residue.psi,
        config.sheet_phi_center, config.sheet_psi_center,
        config.sheet_tolerance,
    ) {
        SecondaryStructure::Sheet
    } else {
        SecondaryStructure::Coil
    }
}

/// Apply minimum run-length filtering.
///
/// Short isolated runs of helix or sheet below the minimum length
/// are demoted to Coil.
fn apply_run_filter(assignments: &mut [SecondaryStructure], config: &SecondaryConfig) {
    if assignments.is_empty() {
        return;
    }

    let n = assignments.len();
    let mut i = 0;
    while i < n {
        let kind = assignments[i];
        if kind == SecondaryStructure::Coil {
            i += 1;
            continue;
        }

        // Find run end
        let start = i;
        while i < n && assignments[i] == kind {
            i += 1;
        }
        let run_len = i - start;

        let min_run = match kind {
            SecondaryStructure::Helix => config.min_helix_run,
            SecondaryStructure::Sheet => config.min_sheet_run,
            SecondaryStructure::Coil => 0,
        };

        if run_len < min_run {
            for elem in assignments[start..i].iter_mut() {
                *elem = SecondaryStructure::Coil;
            }
        }
    }
}

/// Refine assignments using hydrogen bond patterns.
///
/// Alpha helix: i→i+4 H-bond pattern.
/// Beta sheet: i→j where |i-j| >= 5.
fn refine_with_hbonds(
    assignments: &mut [SecondaryStructure],
    hits: &[HBondHit],
    config: &SecondaryConfig,
) {
    let n = assignments.len();
    if n == 0 || hits.is_empty() {
        return;
    }

    // Score each residue: +1 for helix pattern, +1 for sheet pattern
    let mut helix_score = vec![0i32; n];
    let mut sheet_score = vec![0i32; n];

    for hit in hits {
        let d = hit.donor_idx;
        let a = hit.acceptor_idx;
        if d >= n || a >= n {
            continue;
        }
        let sep = a.abs_diff(d);

        if sep == 4 {
            // Alpha helix pattern: i → i+4
            let lo = d.min(a);
            let hi = d.max(a);
            for elem in helix_score[lo..=hi.min(n - 1)].iter_mut() {
                *elem += 1;
            }
        } else if sep >= 5 {
            // Beta sheet pattern: long-range
            sheet_score[d] += 1;
            sheet_score[a] += 1;
        }
    }

    // Merge H-bond evidence with angle-based assignment
    for i in 0..n {
        if helix_score[i] > 0 && assignments[i] != SecondaryStructure::Sheet {
            assignments[i] = SecondaryStructure::Helix;
        } else if sheet_score[i] > 0 && assignments[i] != SecondaryStructure::Helix {
            assignments[i] = SecondaryStructure::Sheet;
        }
    }

    // Re-apply run filter after H-bond refinement
    apply_run_filter(assignments, config);
}

/// Assign secondary structure to a chain of residues.
///
/// If `positions` is provided and `config.use_hbonds` is true,
/// hydrogen bond patterns are used to refine the assignment.
pub fn assign_secondary_structure(
    residues: &[Residue],
    positions: Option<&[[f64; 3]]>,
    config: &SecondaryConfig,
) -> SecondaryAssignment {
    let n = residues.len();

    if n == 0 {
        return SecondaryAssignment {
            assignments: Vec::new(),
            helix_count: 0,
            sheet_count: 0,
            coil_count: 0,
            helix_fraction: 0.0,
            sheet_fraction: 0.0,
            content_hash: fnv1a(b"empty_secondary"),
        };
    }

    // Phase 1: Angle-based classification
    let mut assignments: Vec<SecondaryStructure> = residues
        .iter()
        .map(|r| classify_by_angles(r, config))
        .collect();

    // Phase 2: Run-length filtering
    apply_run_filter(&mut assignments, config);

    // Phase 3: H-bond refinement (optional)
    if config.use_hbonds {
        if let Some(pos) = positions {
            let hbond_cfg = HBondConfig {
                max_distance: 5.0,
                min_seq_separation: 3,
                energy_cutoff: -0.3,
            };
            let detector = HBondDetector::new(hbond_cfg);
            let hits = detector.detect(residues, pos);
            refine_with_hbonds(&mut assignments, &hits, config);
        }
    }

    // Compute statistics
    let mut helix_count = 0usize;
    let mut sheet_count = 0usize;
    let mut coil_count = 0usize;
    for &a in &assignments {
        match a {
            SecondaryStructure::Helix => helix_count += 1,
            SecondaryStructure::Sheet => sheet_count += 1,
            SecondaryStructure::Coil => coil_count += 1,
        }
    }

    let inv_n = 1.0 / n as f64;
    let helix_fraction = helix_count as f64 * inv_n;
    let sheet_fraction = sheet_count as f64 * inv_n;

    // Content hash: hash the assignment byte sequence
    let mut buf = Vec::with_capacity(n + 8);
    buf.extend_from_slice(&(n as u64).to_le_bytes());
    for &a in &assignments {
        buf.push(match a {
            SecondaryStructure::Helix => 0,
            SecondaryStructure::Sheet => 1,
            SecondaryStructure::Coil => 2,
        });
    }
    let content_hash = fnv1a(&buf);

    SecondaryAssignment {
        assignments,
        helix_count,
        sheet_count,
        coil_count,
        helix_fraction,
        sheet_fraction,
        content_hash,
    }
}

// ── Tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amino::AminoAcid;
    use std::f64::consts::PI;

    fn helix_residue(aa: AminoAcid) -> Residue {
        // phi ≈ -60°, psi ≈ -45° (alpha helix)
        Residue::new(aa, -1.047, -0.785, PI)
    }

    fn sheet_residue(aa: AminoAcid) -> Residue {
        // phi ≈ -120°, psi ≈ 130° (beta sheet)
        Residue::new(aa, -2.094, 2.269, PI)
    }

    fn coil_residue(aa: AminoAcid) -> Residue {
        // phi ≈ 0°, psi ≈ 0° (coil / disallowed)
        Residue::new(aa, 0.0, 0.0, PI)
    }

    #[test]
    fn empty_chain() {
        let result = assign_secondary_structure(&[], None, &SecondaryConfig::default());
        assert!(result.assignments.is_empty());
        assert_eq!(result.helix_count, 0);
        assert_eq!(result.sheet_count, 0);
        assert_eq!(result.coil_count, 0);
        assert_ne!(result.content_hash, 0);
    }

    #[test]
    fn single_residue_coil() {
        let residues = [coil_residue(AminoAcid::Gly)];
        let result = assign_secondary_structure(&residues, None, &SecondaryConfig::default());
        assert_eq!(result.assignments.len(), 1);
        assert_eq!(result.assignments[0], SecondaryStructure::Coil);
        assert_eq!(result.coil_count, 1);
    }

    #[test]
    fn helix_region_classified() {
        // 5 helix residues (above min_helix_run=3)
        let residues: Vec<Residue> = (0..5).map(|_| helix_residue(AminoAcid::Ala)).collect();
        let result = assign_secondary_structure(&residues, None, &SecondaryConfig::default());
        for a in &result.assignments {
            assert_eq!(*a, SecondaryStructure::Helix);
        }
        assert_eq!(result.helix_count, 5);
        assert!((result.helix_fraction - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sheet_region_classified() {
        // 4 sheet residues (above min_sheet_run=2)
        let residues: Vec<Residue> = (0..4).map(|_| sheet_residue(AminoAcid::Val)).collect();
        let result = assign_secondary_structure(&residues, None, &SecondaryConfig::default());
        for a in &result.assignments {
            assert_eq!(*a, SecondaryStructure::Sheet);
        }
        assert_eq!(result.sheet_count, 4);
    }

    #[test]
    fn short_helix_demoted_to_coil() {
        // 2 helix residues (below min_helix_run=3) → coil
        let residues: Vec<Residue> = vec![
            helix_residue(AminoAcid::Ala),
            helix_residue(AminoAcid::Gly),
        ];
        let result = assign_secondary_structure(&residues, None, &SecondaryConfig::default());
        for a in &result.assignments {
            assert_eq!(*a, SecondaryStructure::Coil);
        }
    }

    #[test]
    fn short_sheet_demoted_to_coil() {
        // 1 sheet residue (below min_sheet_run=2) → coil
        let residues = vec![sheet_residue(AminoAcid::Leu)];
        let result = assign_secondary_structure(&residues, None, &SecondaryConfig::default());
        assert_eq!(result.assignments[0], SecondaryStructure::Coil);
    }

    #[test]
    fn mixed_chain() {
        // 4 helix + 3 sheet + 2 coil
        let mut residues = Vec::new();
        for _ in 0..4 { residues.push(helix_residue(AminoAcid::Ala)); }
        for _ in 0..3 { residues.push(sheet_residue(AminoAcid::Val)); }
        for _ in 0..2 { residues.push(coil_residue(AminoAcid::Gly)); }

        let result = assign_secondary_structure(&residues, None, &SecondaryConfig::default());
        assert_eq!(result.helix_count, 4);
        assert_eq!(result.sheet_count, 3);
        assert_eq!(result.coil_count, 2);
        assert_eq!(result.assignments.len(), 9);
    }

    #[test]
    fn fractions_sum_to_one() {
        let residues: Vec<Residue> = vec![
            helix_residue(AminoAcid::Ala),
            helix_residue(AminoAcid::Ala),
            helix_residue(AminoAcid::Ala),
            sheet_residue(AminoAcid::Val),
            sheet_residue(AminoAcid::Val),
            coil_residue(AminoAcid::Gly),
        ];
        let result = assign_secondary_structure(&residues, None, &SecondaryConfig::default());
        let total = result.helix_fraction + result.sheet_fraction
            + (result.coil_count as f64 / residues.len() as f64);
        assert!((total - 1.0).abs() < 1e-10);
    }

    #[test]
    fn content_hash_deterministic() {
        let residues: Vec<Residue> = (0..5).map(|_| helix_residue(AminoAcid::Ala)).collect();
        let cfg = SecondaryConfig::default();
        let r1 = assign_secondary_structure(&residues, None, &cfg);
        let r2 = assign_secondary_structure(&residues, None, &cfg);
        assert_eq!(r1.content_hash, r2.content_hash);
        assert_ne!(r1.content_hash, 0);
    }

    #[test]
    fn content_hash_differs_for_different_input() {
        let helix: Vec<Residue> = (0..5).map(|_| helix_residue(AminoAcid::Ala)).collect();
        let sheet: Vec<Residue> = (0..5).map(|_| sheet_residue(AminoAcid::Val)).collect();
        let cfg = SecondaryConfig::default();
        let r1 = assign_secondary_structure(&helix, None, &cfg);
        let r2 = assign_secondary_structure(&sheet, None, &cfg);
        assert_ne!(r1.content_hash, r2.content_hash);
    }

    #[test]
    fn angle_diff_wraps_correctly() {
        // PI and -PI should have diff ~0
        assert!(angle_diff(PI, -PI).abs() < 1e-10);
        // 0 and 2*PI should have diff ~0
        assert!(angle_diff(0.0, 2.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn classify_boundary_tolerance() {
        // Residue at exact center of helix region
        let r = Residue::new(AminoAcid::Ala, -1.047, -0.785, PI);
        let cfg = SecondaryConfig::default();
        assert_eq!(classify_by_angles(&r, &cfg), SecondaryStructure::Helix);

        // Residue just outside tolerance
        let r2 = Residue::new(AminoAcid::Ala, -1.047 + 0.6, -0.785, PI);
        assert_eq!(classify_by_angles(&r2, &cfg), SecondaryStructure::Coil);
    }

    #[test]
    fn custom_config_relaxed_tolerance() {
        // With relaxed tolerance, more residues become helix
        let residues = vec![
            Residue::new(AminoAcid::Ala, -0.9, -0.6, PI),
            Residue::new(AminoAcid::Ala, -0.9, -0.6, PI),
            Residue::new(AminoAcid::Ala, -0.9, -0.6, PI),
        ];
        let strict = SecondaryConfig { helix_tolerance: 0.1, ..Default::default() };
        let relaxed = SecondaryConfig { helix_tolerance: 0.7, ..Default::default() };

        let r_strict = assign_secondary_structure(&residues, None, &strict);
        let r_relaxed = assign_secondary_structure(&residues, None, &relaxed);

        assert!(r_relaxed.helix_count >= r_strict.helix_count);
    }

    #[test]
    fn run_filter_preserves_long_runs() {
        let mut assignments = vec![
            SecondaryStructure::Helix,
            SecondaryStructure::Helix,
            SecondaryStructure::Helix,
            SecondaryStructure::Helix,
            SecondaryStructure::Sheet,
            SecondaryStructure::Sheet,
        ];
        let cfg = SecondaryConfig::default();
        apply_run_filter(&mut assignments, &cfg);
        // Helix run=4 >= min_helix_run=3 → preserved
        assert_eq!(assignments[0], SecondaryStructure::Helix);
        assert_eq!(assignments[3], SecondaryStructure::Helix);
        // Sheet run=2 >= min_sheet_run=2 → preserved
        assert_eq!(assignments[4], SecondaryStructure::Sheet);
    }

    #[test]
    fn run_filter_demotes_short_runs() {
        let mut assignments = vec![
            SecondaryStructure::Helix,
            SecondaryStructure::Helix,
            SecondaryStructure::Coil,
            SecondaryStructure::Sheet,
        ];
        let cfg = SecondaryConfig::default();
        apply_run_filter(&mut assignments, &cfg);
        // Helix run=2 < min_helix_run=3 → demoted
        assert_eq!(assignments[0], SecondaryStructure::Coil);
        assert_eq!(assignments[1], SecondaryStructure::Coil);
        // Sheet run=1 < min_sheet_run=2 → demoted
        assert_eq!(assignments[3], SecondaryStructure::Coil);
    }

    #[test]
    fn in_region_at_center_is_true() {
        assert!(in_region(-1.047, -0.785, -1.047, -0.785, 0.5));
    }

    #[test]
    fn in_region_outside_tolerance_is_false() {
        // phi = 0 is far from phi_center = -1.047
        assert!(!in_region(0.0, -0.785, -1.047, -0.785, 0.5));
    }

    #[test]
    fn with_hbonds_flag_and_positions() {
        // Test the hbond refinement path: use_hbonds = true, provide positions
        let mut residues = Vec::new();
        for _ in 0..8 {
            residues.push(helix_residue(AminoAcid::Ala));
        }
        // Create positions where residues 0 and 4 are close (alpha helix pattern)
        let positions: Vec<[f64; 3]> = (0..8)
            .map(|i| [i as f64 * 3.8, 0.0, 0.0])
            .collect();
        let cfg = SecondaryConfig {
            use_hbonds: true,
            ..Default::default()
        };
        let result = assign_secondary_structure(&residues, Some(&positions), &cfg);
        assert_eq!(result.assignments.len(), 8);
        // All should be classified (no panic, valid results)
        let total = result.helix_count + result.sheet_count + result.coil_count;
        assert_eq!(total, 8);
    }

    #[test]
    fn all_coil_chain() {
        // All coil residues → all Coil, fractions match
        let residues: Vec<Residue> = (0..6).map(|_| coil_residue(AminoAcid::Gly)).collect();
        let result = assign_secondary_structure(&residues, None, &SecondaryConfig::default());
        assert_eq!(result.coil_count, 6);
        assert_eq!(result.helix_count, 0);
        assert_eq!(result.sheet_count, 0);
        assert!((result.helix_fraction - 0.0).abs() < 1e-15);
        assert!((result.sheet_fraction - 0.0).abs() < 1e-15);
    }

    #[test]
    fn run_filter_on_empty_slice() {
        let mut assignments: Vec<SecondaryStructure> = vec![];
        let cfg = SecondaryConfig::default();
        apply_run_filter(&mut assignments, &cfg);
        assert!(assignments.is_empty());
    }
}
