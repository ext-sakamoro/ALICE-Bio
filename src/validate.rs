//! 構造検証 — Ramachandranプロット領域判定、steric clash検出、bond length検証
//!
//! タンパク質構造の妥当性を統計的・幾何学的に評価する。

use crate::amino::Residue;

// ============================================================================
// Ramachandran 統計
// ============================================================================

/// Ramachandranプロット上での残基分類統計。
#[derive(Debug, Clone, Default)]
pub struct RamachandranStats {
    /// 許容領域内の残基数。
    pub favored: usize,
    /// 許容外（Disallowed）の残基数。
    pub outlier: usize,
    /// 全残基数。
    pub total: usize,
}

impl RamachandranStats {
    /// 許容領域内の割合 (0.0–1.0)。
    #[must_use]
    pub fn favored_fraction(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        self.favored as f64 / self.total as f64
    }

    /// Outlier率 (0.0–1.0)。
    #[must_use]
    pub fn outlier_fraction(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        self.outlier as f64 / self.total as f64
    }
}

/// 残基配列のRamachandran統計を算出。
#[must_use]
pub fn ramachandran_stats(residues: &[Residue]) -> RamachandranStats {
    let mut favored = 0usize;
    let mut outlier = 0usize;
    for res in residues {
        if res.is_ramachandran_allowed() {
            favored += 1;
        } else {
            outlier += 1;
        }
    }
    RamachandranStats {
        favored,
        outlier,
        total: residues.len(),
    }
}

// ============================================================================
// Steric Clash 検出
// ============================================================================

/// Steric clash 情報。
#[derive(Debug, Clone, Copy)]
pub struct StericClash {
    /// 残基 i のインデックス。
    pub i: usize,
    /// 残基 j のインデックス。
    pub j: usize,
    /// 実際の距離 (Å)。
    pub distance: f64,
    /// `VdW接触距離` (Å)。
    pub contact_distance: f64,
}

/// Steric clash 検出設定。
#[derive(Debug, Clone, Copy)]
pub struct ClashConfig {
    /// `VdW` 半径の合計に対する許容率。0.7 = VdW合計の70%未満で clash。
    pub overlap_factor: f64,
    /// 配列上の最小距離。近接残基は除外。
    pub min_seq_separation: usize,
}

impl Default for ClashConfig {
    fn default() -> Self {
        Self {
            overlap_factor: 0.7,
            min_seq_separation: 2,
        }
    }
}

/// Cα座標に基づくsteric clash検出。
///
/// `二残基のCα距離がVdW半径合計` × `overlap_factor` 未満の場合 clash と判定。
#[must_use]
pub fn detect_clashes(
    residues: &[Residue],
    positions: &[[f64; 3]],
    config: &ClashConfig,
) -> Vec<StericClash> {
    let n = residues.len().min(positions.len());
    let mut clashes = Vec::new();

    for i in 0..n {
        for j in (i + 1)..n {
            if j - i < config.min_seq_separation {
                continue;
            }
            let dx = positions[j][0] - positions[i][0];
            let dy = positions[j][1] - positions[i][1];
            let dz = positions[j][2] - positions[i][2];
            let dist = dz.mul_add(dz, dx.mul_add(dx, dy * dy)).sqrt();

            let contact = (residues[i].amino.van_der_waals_radius()
                + residues[j].amino.van_der_waals_radius())
                * config.overlap_factor;

            if dist < contact {
                clashes.push(StericClash {
                    i,
                    j,
                    distance: dist,
                    contact_distance: contact,
                });
            }
        }
    }

    clashes
}

// ============================================================================
// Bond Length 検証
// ============================================================================

/// Cα-Cα 結合長の理想値 (3.8 Å)。
const IDEAL_CA_DISTANCE: f64 = 3.8;

/// Bond length 逸脱情報。
#[derive(Debug, Clone, Copy)]
pub struct BondDeviation {
    /// 残基 i のインデックス。
    pub i: usize,
    /// 残基 i+1 のインデックス。
    pub j: usize,
    /// 実測距離 (Å)。
    pub distance: f64,
    /// 理想距離からの偏差 (Å)。
    pub deviation: f64,
}

/// Cα-Cα結合長を検証し、逸脱の大きい結合を返す。
///
/// `tolerance` は許容偏差 (Å)。デフォルト 0.5 Å。
#[must_use]
pub fn validate_bond_lengths(positions: &[[f64; 3]], tolerance: f64) -> Vec<BondDeviation> {
    let mut deviations = Vec::new();
    if positions.len() < 2 {
        return deviations;
    }

    for i in 0..(positions.len() - 1) {
        let j = i + 1;
        let dx = positions[j][0] - positions[i][0];
        let dy = positions[j][1] - positions[i][1];
        let dz = positions[j][2] - positions[i][2];
        let dist = dz.mul_add(dz, dx.mul_add(dx, dy * dy)).sqrt();
        let deviation = (dist - IDEAL_CA_DISTANCE).abs();
        if deviation > tolerance {
            deviations.push(BondDeviation {
                i,
                j,
                distance: dist,
                deviation,
            });
        }
    }

    deviations
}

/// 全体の構造検証スコア (0.0–1.0)。
///
/// Ramachandran favored率 × (1 - clash率) × bond length正常率 の統合指標。
#[must_use]
pub fn validation_score(residues: &[Residue], positions: &[[f64; 3]]) -> f64 {
    if residues.is_empty() || positions.is_empty() {
        return 0.0;
    }
    let n = residues.len().min(positions.len());

    // Ramachandran
    let rama = ramachandran_stats(residues);
    let rama_score = rama.favored_fraction();

    // Steric clashes
    let max_pairs = n * (n - 1) / 2;
    let clashes = detect_clashes(residues, positions, &ClashConfig::default());
    let clash_score = if max_pairs > 0 {
        1.0 - (clashes.len() as f64 / max_pairs as f64).min(1.0)
    } else {
        1.0
    };

    // Bond lengths
    let bond_devs = validate_bond_lengths(positions, 0.5);
    let n_bonds = if positions.len() > 1 {
        positions.len() - 1
    } else {
        1
    };
    let bond_score = 1.0 - (bond_devs.len() as f64 / n_bonds as f64).min(1.0);

    rama_score * clash_score * bond_score
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amino::AminoAcid;
    use std::f64::consts::PI;

    fn alpha_helix_residue() -> Residue {
        Residue::new(
            AminoAcid::Ala,
            -57.0_f64.to_radians(),
            -47.0_f64.to_radians(),
            PI,
        )
    }

    fn disallowed_residue() -> Residue {
        Residue::new(AminoAcid::Ala, 0.0, 0.0, PI)
    }

    #[test]
    fn rama_stats_all_favored() {
        let residues: Vec<Residue> = (0..5).map(|_| alpha_helix_residue()).collect();
        let stats = ramachandran_stats(&residues);
        assert_eq!(stats.favored, 5);
        assert_eq!(stats.outlier, 0);
        assert!((stats.favored_fraction() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rama_stats_all_outlier() {
        let residues: Vec<Residue> = (0..3).map(|_| disallowed_residue()).collect();
        let stats = ramachandran_stats(&residues);
        assert_eq!(stats.outlier, 3);
        assert!((stats.outlier_fraction() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rama_stats_empty() {
        let stats = ramachandran_stats(&[]);
        assert_eq!(stats.total, 0);
        assert!((stats.favored_fraction() - 0.0).abs() < 1e-10);
        assert!((stats.outlier_fraction() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn rama_stats_mixed() {
        let residues = vec![alpha_helix_residue(), disallowed_residue()];
        let stats = ramachandran_stats(&residues);
        assert_eq!(stats.favored, 1);
        assert_eq!(stats.outlier, 1);
        assert!((stats.favored_fraction() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn no_clashes_normal_spacing() {
        let residues = vec![
            Residue::new(AminoAcid::Ala, 0.0, 0.0, PI),
            Residue::new(AminoAcid::Gly, 0.0, 0.0, PI),
            Residue::new(AminoAcid::Val, 0.0, 0.0, PI),
        ];
        let positions = [[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [7.6, 0.0, 0.0]];
        let clashes = detect_clashes(&residues, &positions, &ClashConfig::default());
        assert!(clashes.is_empty());
    }

    #[test]
    fn detect_clash_overlapping() {
        let residues = vec![
            Residue::new(AminoAcid::Ala, 0.0, 0.0, PI),
            Residue::new(AminoAcid::Val, 0.0, 0.0, PI),
            Residue::new(AminoAcid::Leu, 0.0, 0.0, PI),
        ];
        // 残基0と残基2が非常に近い
        let positions = [[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let clashes = detect_clashes(&residues, &positions, &ClashConfig::default());
        assert!(!clashes.is_empty());
        assert_eq!(clashes[0].i, 0);
        assert_eq!(clashes[0].j, 2);
    }

    #[test]
    fn clash_respects_seq_separation() {
        let residues = vec![
            Residue::new(AminoAcid::Ala, 0.0, 0.0, PI),
            Residue::new(AminoAcid::Ala, 0.0, 0.0, PI),
        ];
        // 非常に近いが seq_separation=2 で隣接は除外
        let positions = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]];
        let clashes = detect_clashes(&residues, &positions, &ClashConfig::default());
        assert!(clashes.is_empty());
    }

    #[test]
    fn bond_length_normal() {
        let positions = [[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [7.6, 0.0, 0.0]];
        let devs = validate_bond_lengths(&positions, 0.5);
        assert!(devs.is_empty());
    }

    #[test]
    fn bond_length_deviation_detected() {
        let positions = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let devs = validate_bond_lengths(&positions, 0.5);
        assert_eq!(devs.len(), 1);
        assert!((devs[0].distance - 10.0).abs() < 1e-10);
    }

    #[test]
    fn bond_length_empty() {
        let devs = validate_bond_lengths(&[], 0.5);
        assert!(devs.is_empty());
    }

    #[test]
    fn bond_length_single_position() {
        let devs = validate_bond_lengths(&[[0.0, 0.0, 0.0]], 0.5);
        assert!(devs.is_empty());
    }

    #[test]
    fn validation_score_perfect_structure() {
        let residues: Vec<Residue> = (0..5).map(|_| alpha_helix_residue()).collect();
        let positions: Vec<[f64; 3]> = (0..5).map(|i| [i as f64 * 3.8, 0.0, 0.0]).collect();
        let score = validation_score(&residues, &positions);
        assert!(
            score > 0.5,
            "Perfect structure should score > 0.5, got {score}"
        );
    }

    #[test]
    fn validation_score_empty() {
        assert!((validation_score(&[], &[]) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn validation_score_in_range() {
        let residues = vec![
            alpha_helix_residue(),
            disallowed_residue(),
            alpha_helix_residue(),
        ];
        let positions = [[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [7.6, 0.0, 0.0]];
        let score = validation_score(&residues, &positions);
        assert!((0.0..=1.0).contains(&score));
    }

    #[test]
    fn clash_config_default() {
        let cfg = ClashConfig::default();
        assert!((cfg.overlap_factor - 0.7).abs() < 1e-10);
        assert_eq!(cfg.min_seq_separation, 2);
    }
}
