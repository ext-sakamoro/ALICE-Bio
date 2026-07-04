//! Convenience re-export (= `use alice_bio::prelude::*;` で主要 API 一括取得)
//!
//! Molecular biology 13 module のうち core (amino / fold / forcefield / hbond /
//! homology / interaction / md_bio / potential / secondary / sequence /
//! validate / cell_list / fasta v2) から主要型 + 関数を prelude 経由で提供

pub use crate::amino::{AminoAcid, RamachandranRegion, Residue};
pub use crate::cell_list::{CellList, CellListConfig};
pub use crate::fasta::{
    parse as parse_fasta_v2, reverse_complement, FastaError, FastaRecord as FastaRecordV2,
};
pub use crate::fold::ProteinSdf;
pub use crate::forcefield::{AngleParam, AtomType, BondParam, DihedralParam};
pub use crate::hbond::{HBondConfig, HBondDetector, HBondHit};
pub use crate::homology::{build_model, HomologyModel, Template};
pub use crate::interaction::{
    contact_map, end_to_end_distance, evaluate_pairwise_energy, radius_of_gyration,
};
pub use crate::md_bio::{BioParticle, LangevinConfig};
pub use crate::potential::{coulomb, hydrogen_bond, lennard_jones, torsion_potential, TotalEnergy};
pub use crate::secondary::{assign_secondary_structure, SecondaryStructure};
pub use crate::sequence::{parse_fasta, FastaRecord};
pub use crate::validate::{
    detect_clashes, ramachandran_stats, validate_bond_lengths, validation_score,
};
