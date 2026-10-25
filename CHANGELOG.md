# Changelog

All notable changes to ALICE-Bio will be documented in this file.

## [Unreleased]

### Changed
- License: `AGPL-3.0-only` → `LicenseRef-Proprietary`; `LICENSE-Proprietary` holds the text
- CI: `ci/**` の branch への送信でも CI を走らせる
- preflight: `.github/preflight-exact` を置いた `scripts/preflight.sh` が生成器の出力 (見出しの `generator rev` の版) と byte で一致しなければ送る前の確認が通らない
- 設定 file の comment を中身の説明だけにし、.gitignore を project の build 出力だけにした
- TRADEMARK_NOTICE が指す licence の file を LICENSE-* の file にした

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
