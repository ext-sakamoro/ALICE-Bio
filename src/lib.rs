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
pub mod potential;
pub mod fold;
pub mod interaction;

pub use amino::{AminoAcid, Residue};
pub use potential::{lennard_jones, coulomb, hydrogen_bond, torsion_potential, TotalEnergy};
pub use fold::ProteinSdf;
pub use interaction::{evaluate_pairwise_energy, contact_map, radius_of_gyration, end_to_end_distance};
