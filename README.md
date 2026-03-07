# ALICE-Bio

Molecular structure as signed distance functions, amino acid energy potentials, and protein folding SDF.

## Modules

| Module | Description |
|--------|-------------|
| `amino` | Amino acid definitions, Ramachandran regions, residue representation |
| `cell_list` | Spatial cell list for neighbor search |
| `fold` | Protein SDF (signed distance function for folding) |
| `forcefield` | Bond, angle, dihedral force field parameters |
| `hbond` | Hydrogen bond detection |
| `homology` | Homology modeling (template-based structure prediction) |
| `interaction` | Molecular interaction calculations |
| `md_bio` | Langevin dynamics for bio particles |
| `potential` | Lennard-Jones, Coulomb, H-bond, torsion energy potentials |
| `secondary` | Secondary structure assignment (helix, sheet, coil) |
| `sequence` | FASTA sequence parsing |
| `validate` | Clash detection, Ramachandran stats, bond length validation |

## Example

```rust
use alice_bio::{AminoAcid, Residue, lennard_jones, parse_fasta};

// Parse FASTA sequence
let records = parse_fasta(">test\nMKFLV");

// Lennard-Jones potential
let energy = lennard_jones(3.5, 0.01, 3.4);

// Validation
let score = alice_bio::validation_score(&residues);
```

## Quality

| Metric | Value |
|--------|-------|
| clippy (pedantic+nursery) | 0 warnings |
| Tests | 194 |
| fmt | clean |

## License

AGPL-3.0-only
