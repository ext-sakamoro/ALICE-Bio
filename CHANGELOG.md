# Changelog

All notable changes to ALICE-Bio will be documented in this file.

## [0.1.0] - 2026-02-23

### Added
- `amino` — `AminoAcid` (20 standard), `Residue`, `RamachandranRegion` validation
- `potential` — Lennard-Jones, Coulomb, hydrogen bond, torsion potential, `TotalEnergy`
- `fold` — `ProteinSdf` backbone trace and SDF evaluation from residue chain
- `hbond` — `HBondDetector` distance/angle-based hydrogen bond detection
- `interaction` — Pairwise energy, contact map, radius of gyration, end-to-end distance
- `cell_list` — Spatial cell list for O(N) neighbor queries
- `secondary` — Secondary structure assignment (helix, sheet, coil) from phi/psi angles
- FNV-1a shared hash utility
- Zero external dependencies
- 121 unit tests + 1 doc-test
