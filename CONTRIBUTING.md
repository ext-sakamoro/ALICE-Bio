# Contributing to ALICE-Bio

## Build

```bash
cargo build
```

## Test

```bash
cargo test
```

## Lint

```bash
cargo clippy -- -W clippy::all
cargo fmt -- --check
cargo doc --no-deps 2>&1 | grep warning
```

## Design Constraints

- **Molecular SDF**: models protein structures as signed distance functions instead of coordinate point clouds.
- **Residue-level representation**: phi/psi/omega backbone angles per amino acid residue.
- **Physics potentials**: Lennard-Jones, Coulomb, hydrogen bond, and torsion energy evaluations.
- **Cell list spatial index**: O(N) neighbor queries for pairwise interaction computation.
- **Ramachandran validation**: phi/psi angle region classification (alpha-helix, beta-sheet, etc.).
- **Zero external dependencies**: all physics and geometry are self-contained.
