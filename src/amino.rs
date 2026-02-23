//! Amino acid definitions and residue backbone geometry.

/// The 20 standard amino acids.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AminoAcid {
    Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly,
    His, Ile, Leu, Lys, Met, Phe, Pro, Ser,
    Thr, Trp, Tyr, Val,
}

impl AminoAcid {
    /// Van der Waals radius in Angstroms (representative sidechain).
    #[inline]
    pub fn van_der_waals_radius(&self) -> f64 {
        match self {
            Self::Gly => 1.7,
            Self::Ala => 1.9,
            Self::Val => 2.3,
            Self::Leu => 2.4,
            Self::Ile => 2.4,
            Self::Pro => 2.1,
            Self::Phe => 2.6,
            Self::Trp => 2.8,
            Self::Met => 2.3,
            Self::Ser => 1.9,
            Self::Thr => 2.1,
            Self::Cys => 2.0,
            Self::Tyr => 2.7,
            Self::His => 2.4,
            Self::Asp => 2.1,
            Self::Glu => 2.3,
            Self::Asn => 2.2,
            Self::Gln => 2.4,
            Self::Lys => 2.5,
            Self::Arg => 2.6,
        }
    }

    /// Molecular mass in Daltons.
    #[inline]
    pub fn mass(&self) -> f64 {
        match self {
            Self::Gly => 57.02,  Self::Ala => 71.04,  Self::Val => 99.07,
            Self::Leu => 113.08, Self::Ile => 113.08, Self::Pro => 97.05,
            Self::Phe => 147.07, Self::Trp => 186.08, Self::Met => 131.04,
            Self::Ser => 87.03,  Self::Thr => 101.05, Self::Cys => 103.01,
            Self::Tyr => 163.06, Self::His => 137.06, Self::Asp => 115.03,
            Self::Glu => 129.04, Self::Asn => 114.04, Self::Gln => 128.06,
            Self::Lys => 128.09, Self::Arg => 156.10,
        }
    }

    /// Kyte-Doolittle hydrophobicity scale (-4.5 to 4.5).
    #[inline]
    pub fn hydrophobicity(&self) -> f64 {
        match self {
            Self::Ile =>  4.5, Self::Val =>  4.2, Self::Leu =>  3.8,
            Self::Phe =>  2.8, Self::Cys =>  2.5, Self::Met =>  1.9,
            Self::Ala =>  1.8, Self::Gly => -0.4, Self::Thr => -0.7,
            Self::Ser => -0.8, Self::Trp => -0.9, Self::Tyr => -1.3,
            Self::Pro => -1.6, Self::His => -3.2, Self::Glu => -3.5,
            Self::Gln => -3.5, Self::Asp => -3.5, Self::Asn => -3.5,
            Self::Lys => -3.9, Self::Arg => -4.5,
        }
    }

    /// Single-letter IUPAC code.
    #[inline]
    pub fn one_letter(&self) -> char {
        match self {
            Self::Ala => 'A', Self::Arg => 'R', Self::Asn => 'N', Self::Asp => 'D',
            Self::Cys => 'C', Self::Gln => 'Q', Self::Glu => 'E', Self::Gly => 'G',
            Self::His => 'H', Self::Ile => 'I', Self::Leu => 'L', Self::Lys => 'K',
            Self::Met => 'M', Self::Phe => 'F', Self::Pro => 'P', Self::Ser => 'S',
            Self::Thr => 'T', Self::Trp => 'W', Self::Tyr => 'Y', Self::Val => 'V',
        }
    }

    /// All 20 standard amino acids.
    pub const ALL: [AminoAcid; 20] = [
        Self::Ala, Self::Arg, Self::Asn, Self::Asp, Self::Cys,
        Self::Gln, Self::Glu, Self::Gly, Self::His, Self::Ile,
        Self::Leu, Self::Lys, Self::Met, Self::Phe, Self::Pro,
        Self::Ser, Self::Thr, Self::Trp, Self::Tyr, Self::Val,
    ];
}

/// A single residue in a polypeptide chain.
#[derive(Debug, Clone)]
pub struct Residue {
    pub amino: AminoAcid,
    /// Backbone dihedral angle phi (radians).
    pub phi: f64,
    /// Backbone dihedral angle psi (radians).
    pub psi: f64,
    /// Peptide bond angle omega (radians, typically ~PI).
    pub omega: f64,
}

impl Residue {
    #[inline]
    pub fn new(amino: AminoAcid, phi: f64, psi: f64, omega: f64) -> Self {
        Self { amino, phi, psi, omega }
    }

    /// Basic Ramachandran validation — checks if (phi, psi) falls in
    /// allowed regions (alpha-helix or beta-sheet).
    pub fn is_ramachandran_allowed(&self) -> bool {
        let phi_deg = self.phi.to_degrees();
        let psi_deg = self.psi.to_degrees();

        // Alpha-helix region: phi ∈ [-80, -40], psi ∈ [-60, -30]
        let alpha = (-80.0..=-40.0).contains(&phi_deg)
                 && (-60.0..=-30.0).contains(&psi_deg);

        // Beta-sheet region: phi ∈ [-150, -100], psi ∈ [100, 170]
        let beta = (-150.0..=-100.0).contains(&phi_deg)
                && (100.0..=170.0).contains(&psi_deg);

        // Left-handed helix: phi ∈ [40, 80], psi ∈ [20, 60]
        let left = (40.0..=80.0).contains(&phi_deg)
                && (20.0..=60.0).contains(&psi_deg);

        alpha || beta || left
    }

    /// Classify (phi, psi) into a specific Ramachandran region.
    pub fn ramachandran_region(&self) -> RamachandranRegion {
        let phi_deg = self.phi.to_degrees();
        let psi_deg = self.psi.to_degrees();

        // Alpha-helix: phi ∈ [-80, -40], psi ∈ [-60, -30]
        if (-80.0..=-40.0).contains(&phi_deg) && (-60.0..=-30.0).contains(&psi_deg) {
            return RamachandranRegion::AlphaHelix;
        }
        // Beta-sheet: phi ∈ [-150, -100], psi ∈ [100, 170]
        if (-150.0..=-100.0).contains(&phi_deg) && (100.0..=170.0).contains(&psi_deg) {
            return RamachandranRegion::BetaSheet;
        }
        // Left-handed helix: phi ∈ [40, 80], psi ∈ [20, 60]
        if (40.0..=80.0).contains(&phi_deg) && (20.0..=60.0).contains(&psi_deg) {
            return RamachandranRegion::LeftHandedHelix;
        }
        // Polyproline II helix: phi ∈ [-80, -60], psi ∈ [120, 170]
        if (-80.0..=-60.0).contains(&phi_deg) && (120.0..=170.0).contains(&psi_deg) {
            return RamachandranRegion::PolyProlineII;
        }
        // Beta-turn: phi ∈ [-80, -40], psi ∈ [-10, 40]
        if (-80.0..=-40.0).contains(&phi_deg) && (-10.0..=40.0).contains(&psi_deg) {
            return RamachandranRegion::BetaTurn;
        }
        RamachandranRegion::Disallowed
    }
}

/// Ramachandran region classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RamachandranRegion {
    /// Right-handed alpha-helix.
    AlphaHelix,
    /// Beta-sheet (extended strand).
    BetaSheet,
    /// Left-handed alpha-helix.
    LeftHandedHelix,
    /// Polyproline II helix.
    PolyProlineII,
    /// Beta-turn.
    BetaTurn,
    /// Outside all defined regions.
    Disallowed,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_positive_vdw_radius() {
        for aa in &AminoAcid::ALL {
            assert!(aa.van_der_waals_radius() > 0.0, "{:?}", aa);
        }
    }

    #[test]
    fn all_positive_mass() {
        for aa in &AminoAcid::ALL {
            assert!(aa.mass() > 0.0, "{:?}", aa);
        }
    }

    #[test]
    fn hydrophobicity_range() {
        for aa in &AminoAcid::ALL {
            let h = aa.hydrophobicity();
            assert!(h >= -4.5 && h <= 4.5, "{:?}: {}", aa, h);
        }
    }

    #[test]
    fn one_letter_codes_unique() {
        let codes: Vec<char> = AminoAcid::ALL.iter().map(|a| a.one_letter()).collect();
        let mut dedup = codes.clone();
        dedup.sort();
        dedup.dedup();
        assert_eq!(codes.len(), dedup.len());
    }

    #[test]
    fn ramachandran_alpha_helix_allowed() {
        // phi ≈ -57°, psi ≈ -47° (alpha helix)
        let r = Residue::new(AminoAcid::Ala, -57.0_f64.to_radians(), -47.0_f64.to_radians(), std::f64::consts::PI);
        assert!(r.is_ramachandran_allowed());
    }

    #[test]
    fn ramachandran_disallowed() {
        // phi = 0, psi = 0 — not in any allowed region
        let r = Residue::new(AminoAcid::Gly, 0.0, 0.0, std::f64::consts::PI);
        assert!(!r.is_ramachandran_allowed());
    }

    #[test]
    fn ramachandran_region_alpha_helix() {
        let r = Residue::new(AminoAcid::Ala, -57.0_f64.to_radians(), -47.0_f64.to_radians(), std::f64::consts::PI);
        assert_eq!(r.ramachandran_region(), RamachandranRegion::AlphaHelix);
    }

    #[test]
    fn ramachandran_region_beta_sheet() {
        let r = Residue::new(AminoAcid::Val, -120.0_f64.to_radians(), 130.0_f64.to_radians(), std::f64::consts::PI);
        assert_eq!(r.ramachandran_region(), RamachandranRegion::BetaSheet);
    }

    #[test]
    fn ramachandran_region_ppii() {
        // phi ≈ -75°, psi ≈ 145° → PolyProlineII
        let r = Residue::new(AminoAcid::Pro, -75.0_f64.to_radians(), 145.0_f64.to_radians(), std::f64::consts::PI);
        assert_eq!(r.ramachandran_region(), RamachandranRegion::PolyProlineII);
    }

    #[test]
    fn ramachandran_region_disallowed() {
        let r = Residue::new(AminoAcid::Gly, 0.0, 0.0, std::f64::consts::PI);
        assert_eq!(r.ramachandran_region(), RamachandranRegion::Disallowed);
    }

    #[test]
    fn ramachandran_beta_sheet_allowed() {
        // phi ≈ -120°, psi ≈ 130°
        let r = Residue::new(AminoAcid::Val, -120.0_f64.to_radians(), 130.0_f64.to_radians(), std::f64::consts::PI);
        assert!(r.is_ramachandran_allowed());
    }

    #[test]
    fn ramachandran_left_handed_helix_allowed() {
        // phi ≈ 60°, psi ≈ 40° → left-handed helix region
        let r = Residue::new(AminoAcid::Gly, 60.0_f64.to_radians(), 40.0_f64.to_radians(), std::f64::consts::PI);
        assert!(r.is_ramachandran_allowed());
    }

    #[test]
    fn ramachandran_region_left_handed_helix() {
        // phi ≈ 60°, psi ≈ 40° → LeftHandedHelix
        let r = Residue::new(AminoAcid::Gly, 60.0_f64.to_radians(), 40.0_f64.to_radians(), std::f64::consts::PI);
        assert_eq!(r.ramachandran_region(), RamachandranRegion::LeftHandedHelix);
    }

    #[test]
    fn ramachandran_region_beta_turn() {
        // phi ≈ -60°, psi ≈ 15° → BetaTurn: phi ∈ [-80, -40], psi ∈ [-10, 40]
        let r = Residue::new(AminoAcid::Ala, -60.0_f64.to_radians(), 15.0_f64.to_radians(), std::f64::consts::PI);
        assert_eq!(r.ramachandran_region(), RamachandranRegion::BetaTurn);
    }

    #[test]
    fn ramachandran_alpha_boundary_exact_edge() {
        // phi = -80° (boundary), psi = -60° (boundary) → alpha helix
        let r = Residue::new(AminoAcid::Ala, -80.0_f64.to_radians(), -60.0_f64.to_radians(), std::f64::consts::PI);
        assert!(r.is_ramachandran_allowed());
        assert_eq!(r.ramachandran_region(), RamachandranRegion::AlphaHelix);
    }

    #[test]
    fn ramachandran_just_outside_alpha_region() {
        // phi = -81° (just outside alpha), psi = -47° → should not be alpha
        let r = Residue::new(AminoAcid::Ala, -81.0_f64.to_radians(), -47.0_f64.to_radians(), std::f64::consts::PI);
        assert_ne!(r.ramachandran_region(), RamachandranRegion::AlphaHelix);
    }

    #[test]
    fn all_twenty_amino_acids_in_all() {
        assert_eq!(AminoAcid::ALL.len(), 20);
    }

    #[test]
    fn residue_field_access() {
        let r = Residue::new(AminoAcid::Met, 1.5, -2.3, 3.14);
        assert_eq!(r.amino, AminoAcid::Met);
        assert!((r.phi - 1.5).abs() < 1e-15);
        assert!((r.psi - (-2.3)).abs() < 1e-15);
        assert!((r.omega - 3.14).abs() < 1e-15);
    }

    #[test]
    fn gly_smallest_vdw_radius() {
        // Gly has the smallest sidechain → smallest VdW radius among all 20
        let gly_r = AminoAcid::Gly.van_der_waals_radius();
        for aa in &AminoAcid::ALL {
            assert!(aa.van_der_waals_radius() >= gly_r, "{:?} radius {} < Gly {}", aa, aa.van_der_waals_radius(), gly_r);
        }
    }

    #[test]
    fn trp_heaviest_amino_acid() {
        // Trp (186.08 Da) is the heaviest standard amino acid
        let trp_mass = AminoAcid::Trp.mass();
        for aa in &AminoAcid::ALL {
            assert!(aa.mass() <= trp_mass, "{:?} mass {} > Trp {}", aa, aa.mass(), trp_mass);
        }
    }

    #[test]
    fn ile_most_hydrophobic() {
        // Ile has the highest hydrophobicity (4.5)
        let ile_h = AminoAcid::Ile.hydrophobicity();
        for aa in &AminoAcid::ALL {
            assert!(aa.hydrophobicity() <= ile_h, "{:?}", aa);
        }
    }
}
