//! ALICE-Bio — Molecular Structure as Signed Distance Functions
//!
//! Models molecular structures as SDF primitives instead of raw coordinate
//! point clouds. Amino acid bond angles and energy potentials define compact
//! mathematical representations that can be evaluated at arbitrary resolution.
//!
//! ```
//! use alice_bio::{AminoAcid, Residue, ProteinSdf};
//!
//! let residues = vec![
//!     Residue::new(AminoAcid::Ala, -1.0, -0.8, std::f64::consts::PI),
//!     Residue::new(AminoAcid::Gly, -1.2, -0.5, std::f64::consts::PI),
//! ];
//! let sdf = ProteinSdf::new(residues);
//! assert_eq!(sdf.residue_count(), 2);
//! let d = sdf.eval(&[0.0, 0.0, 0.0]);
//! assert!(d.is_finite());
//! ```

pub mod amino;
pub mod cell_list;
pub mod fold;
pub mod hbond;
pub mod interaction;
pub mod potential;
pub mod secondary;

pub use amino::{AminoAcid, RamachandranRegion, Residue};
pub use cell_list::{CellList, CellListConfig};
pub use fold::ProteinSdf;
pub use hbond::{HBondDetector, HBondHit, HBondConfig};
pub use interaction::{contact_map, end_to_end_distance, evaluate_pairwise_energy, radius_of_gyration};
pub use potential::{coulomb, hydrogen_bond, lennard_jones, torsion_potential, TotalEnergy};
pub use secondary::{assign_secondary_structure, SecondaryStructure};

// ── Shared hash primitive ──────────────────────────────────────────────

#[inline(always)]
pub(crate) fn fnv1a(data: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}
