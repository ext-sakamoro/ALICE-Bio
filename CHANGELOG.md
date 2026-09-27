# Changelog

All notable changes to ALICE-Bio will be documented in this file.

## [Unreleased]

### Changed
- **License: `AGPL-3.0-only` → `AGPL-3.0-only OR LicenseRef-Commercial` (dual-licensed、2026-09-27)** AGPL 側の条件は変更なし (既存 AGPL 利用者への影響ゼロ)、商用という選択肢が追加されただけ SPDX が AGPL 単独だと cargo-deny / FOSSA / SBOM に「商用オプションなし」と見えるため宣言を dual に 変更点: SPDX / `LICENSE` → `LICENSE-AGPL` / `LICENSE-COMMERCIAL.md` (商用トリガー 6 条件 = クローズド製品・商用 SaaS・エッジ / ファームウェア配布・plugin 再配布・プラットフォーム NDA・保証、社内利用は AGPL 側で無償と明記) / README の選択肢表 商用窓口は法人 `contact@extoria.co.jp`

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
