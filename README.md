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

`AGPL-3.0-only OR LicenseRef-Commercial` — dual-licensed. Pick either.

| Option | Terms | Use it when |
|--------|-------|-------------|
| **AGPL-3.0-only** | [LICENSE-AGPL](LICENSE-AGPL) — free, no reporting obligation | Your project is itself AGPL-compatible open source, or you are only using it internally |
| **Commercial License** | [LICENSE-COMMERCIAL.md](LICENSE-COMMERCIAL.md) — paid, removes the copyleft | Closed-source product, proprietary SaaS, edge / firmware distribution, plugin redistribution, or a platform NDA that forbids source disclosure |

AGPL is a strong copyleft: a product, firmware image, or service that links
`alice-bio` and is distributed or served to users must be released under the AGPL
as well. That is intentional for the open ecosystem, and the Commercial
License exists for the cases where it is not something you are able to do.

Commercial licence enquiries: <contact@extoria.co.jp>
