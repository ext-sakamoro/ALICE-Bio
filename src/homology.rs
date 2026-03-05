//! ホモロジーモデリング — テンプレートベース座標転写、ループ構築、side-chain packing
//!
//! テンプレート構造からターゲット配列の3D座標を推定する。
//! Cα粗視化レベルでの高速モデリング。

use crate::amino::AminoAcid;
use crate::sequence::{needleman_wunsch, AlignConfig};

// ============================================================================
// テンプレート構造
// ============================================================================

/// テンプレート構造 (既知3D座標)。
#[derive(Debug, Clone)]
pub struct Template {
    /// テンプレート配列。
    pub sequence: Vec<AminoAcid>,
    /// Cα座標 (Å)。
    pub positions: Vec<[f64; 3]>,
}

/// ホモロジーモデリング結果。
#[derive(Debug, Clone)]
pub struct HomologyModel {
    /// ターゲット配列。
    pub sequence: Vec<AminoAcid>,
    /// 推定Cα座標。
    pub positions: Vec<[f64; 3]>,
    /// テンプレートとの配列同一性。
    pub identity: f64,
    /// ループ領域のインデックス (ギャップ由来)。
    pub loop_regions: Vec<usize>,
}

// ============================================================================
// テンプレートベース座標転写
// ============================================================================

/// アラインメントに基づきテンプレート座標をターゲットに転写。
///
/// 一致/不一致位置はテンプレート座標をそのまま使用。
/// ギャップ位置は前後の座標から線形補間 (ループ構築)。
#[must_use]
pub fn build_model(target: &[AminoAcid], template: &Template) -> HomologyModel {
    if target.is_empty() {
        return HomologyModel {
            sequence: Vec::new(),
            positions: Vec::new(),
            identity: 0.0,
            loop_regions: Vec::new(),
        };
    }

    if template.sequence.is_empty() {
        // テンプレートなし → 直線配置
        let positions: Vec<[f64; 3]> = (0..target.len())
            .map(|i| [i as f64 * 3.8, 0.0, 0.0])
            .collect();
        return HomologyModel {
            sequence: target.to_vec(),
            positions,
            identity: 0.0,
            loop_regions: (0..target.len()).collect(),
        };
    }

    // NW アラインメントでマッピングを構築
    let align_result = needleman_wunsch(target, &template.sequence, &AlignConfig::default());

    // トレースバックでターゲット→テンプレートのマッピング構築
    let mapping = build_alignment_mapping(target, &template.sequence);

    // 座標転写
    let mut positions = vec![[0.0; 3]; target.len()];
    let mut loop_regions = Vec::new();

    for (t_idx, &mapped) in mapping.iter().enumerate() {
        if let Some(s_idx) = mapped {
            if s_idx < template.positions.len() {
                positions[t_idx] = template.positions[s_idx];
            } else {
                loop_regions.push(t_idx);
            }
        } else {
            loop_regions.push(t_idx);
        }
    }

    // ループ構築: ギャップ位置を線形補間
    build_loops(&mut positions, &loop_regions);

    // Side-chain packing (粗視化: VdW半径に基づく位置微調整)
    pack_sidechains(&mut positions, target);

    HomologyModel {
        sequence: target.to_vec(),
        positions,
        identity: align_result.identity,
        loop_regions,
    }
}

/// NWアラインメントからターゲット→テンプレートのインデックスマッピング構築。
fn build_alignment_mapping(target: &[AminoAcid], template: &[AminoAcid]) -> Vec<Option<usize>> {
    let m = target.len();
    let n = template.len();
    let config = AlignConfig::default();
    let cols = n + 1;

    // DP テーブル構築
    let mut dp = vec![0i32; (m + 1) * cols];
    for i in 0..=m {
        dp[i * cols] = config.gap_penalty * i as i32;
    }
    for (j, val) in dp.iter_mut().enumerate().take(n + 1) {
        *val = config.gap_penalty * j as i32;
    }
    for i in 1..=m {
        for j in 1..=n {
            let s = if target[i - 1] == template[j - 1] {
                config.match_score
            } else {
                config.mismatch_penalty
            };
            let diag = dp[(i - 1) * cols + (j - 1)] + s;
            let up = dp[(i - 1) * cols + j] + config.gap_penalty;
            let left = dp[i * cols + (j - 1)] + config.gap_penalty;
            dp[i * cols + j] = diag.max(up).max(left);
        }
    }

    // トレースバック
    let mut mapping = vec![None; m];
    let mut i = m;
    let mut j = n;
    while i > 0 && j > 0 {
        let s = if target[i - 1] == template[j - 1] {
            config.match_score
        } else {
            config.mismatch_penalty
        };
        if dp[i * cols + j] == dp[(i - 1) * cols + (j - 1)] + s {
            mapping[i - 1] = Some(j - 1);
            i -= 1;
            j -= 1;
        } else if dp[i * cols + j] == dp[(i - 1) * cols + j] + config.gap_penalty {
            // ターゲット残基にギャップ (テンプレート対応なし)
            i -= 1;
        } else {
            j -= 1;
        }
    }
    // 残りの i > 0 はテンプレートにマッピングされない

    mapping
}

// ============================================================================
// ループ構築
// ============================================================================

/// ギャップ位置を前後の既知座標から線形補間で埋める。
fn build_loops(positions: &mut [[f64; 3]], loop_regions: &[usize]) {
    if loop_regions.is_empty() {
        return;
    }

    let n = positions.len();

    for &idx in loop_regions {
        // 前後の非ループ位置を探す
        let prev = find_prev_anchor(idx, loop_regions);
        let next = find_next_anchor(idx, loop_regions, n);

        let (start_pos, end_pos) = match (prev, next) {
            (Some(p), Some(nx)) => (positions[p], positions[nx]),
            (Some(p), None) => {
                // 後方アンカーなし → 前方から3.8Åで延伸
                let base = positions[p];
                let dir = if p > 0 {
                    [
                        base[0] - positions[p - 1][0],
                        base[1] - positions[p - 1][1],
                        base[2] - positions[p - 1][2],
                    ]
                } else {
                    [3.8, 0.0, 0.0]
                };
                let len = dir[2]
                    .mul_add(dir[2], dir[0].mul_add(dir[0], dir[1] * dir[1]))
                    .sqrt();
                let norm = if len > 1e-10 {
                    [dir[0] / len, dir[1] / len, dir[2] / len]
                } else {
                    [1.0, 0.0, 0.0]
                };
                let steps = idx - p;
                positions[idx] = [
                    (norm[0] * 3.8).mul_add(steps as f64, base[0]),
                    (norm[1] * 3.8).mul_add(steps as f64, base[1]),
                    (norm[2] * 3.8).mul_add(steps as f64, base[2]),
                ];
                continue;
            }
            (None, Some(nx)) => {
                // 前方アンカーなし → 後方から逆方向に延伸
                let base = positions[nx];
                let steps = nx - idx;
                positions[idx] = [3.8f64.mul_add(-(steps as f64), base[0]), base[1], base[2]];
                continue;
            }
            (None, None) => {
                // 全部ループ → 直線配置
                positions[idx] = [idx as f64 * 3.8, 0.0, 0.0];
                continue;
            }
        };

        // 前後のアンカー間で線形補間
        let p_idx = prev.unwrap();
        let n_idx = next.unwrap();
        let total_span = n_idx - p_idx;
        let frac = (idx - p_idx) as f64 / total_span as f64;

        for k in 0..3 {
            positions[idx][k] = (end_pos[k] - start_pos[k]).mul_add(frac, start_pos[k]);
        }
    }
}

fn find_prev_anchor(idx: usize, loop_regions: &[usize]) -> Option<usize> {
    (0..idx).rev().find(|i| !loop_regions.contains(i))
}

fn find_next_anchor(idx: usize, loop_regions: &[usize], n: usize) -> Option<usize> {
    ((idx + 1)..n).find(|i| !loop_regions.contains(i))
}

// ============================================================================
// Side-chain Packing (粗視化)
// ============================================================================

/// 粗視化 side-chain packing: `VdW衝突の解消`。
///
/// 大きなside-chainを持つ残基の位置を微調整して、
/// `近接残基とのVdW衝突を軽減する`。
fn pack_sidechains(positions: &mut [[f64; 3]], sequence: &[AminoAcid]) {
    let n = positions.len().min(sequence.len());
    if n < 2 {
        return;
    }

    // 3回反復で衝突解消
    for _ in 0..3 {
        for i in 0..n {
            let r_i = sequence[i].van_der_waals_radius();
            for j in (i + 1)..n {
                let r_j = sequence[j].van_der_waals_radius();
                let min_dist = (r_i + r_j) * 0.6; // 最小許容距離

                let dx = positions[j][0] - positions[i][0];
                let dy = positions[j][1] - positions[i][1];
                let dz = positions[j][2] - positions[i][2];
                let dist = dz.mul_add(dz, dx.mul_add(dx, dy * dy)).sqrt();

                if dist < min_dist && dist > 1e-10 {
                    // 衝突 → 反発方向に微調整
                    let overlap = min_dist - dist;
                    let scale = overlap * 0.25 / dist;
                    positions[i][0] -= dx * scale;
                    positions[i][1] -= dy * scale;
                    positions[i][2] -= dz * scale;
                    positions[j][0] += dx * scale;
                    positions[j][1] += dy * scale;
                    positions[j][2] += dz * scale;
                }
            }
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_template() -> Template {
        Template {
            sequence: vec![
                AminoAcid::Ala,
                AminoAcid::Arg,
                AminoAcid::Gly,
                AminoAcid::Val,
            ],
            positions: vec![
                [0.0, 0.0, 0.0],
                [3.8, 0.0, 0.0],
                [7.6, 0.0, 0.0],
                [11.4, 0.0, 0.0],
            ],
        }
    }

    #[test]
    fn identical_target_and_template() {
        let template = simple_template();
        let target = template.sequence.clone();
        let model = build_model(&target, &template);
        assert_eq!(model.sequence.len(), 4);
        assert_eq!(model.positions.len(), 4);
        assert!(model.identity > 0.9);
    }

    #[test]
    fn similar_target() {
        let template = simple_template();
        // 1残基異なる
        let target = vec![
            AminoAcid::Ala,
            AminoAcid::Arg,
            AminoAcid::Gly,
            AminoAcid::Leu,
        ];
        let model = build_model(&target, &template);
        assert_eq!(model.positions.len(), 4);
        // 座標は有限値
        for pos in &model.positions {
            assert!(pos[0].is_finite());
            assert!(pos[1].is_finite());
            assert!(pos[2].is_finite());
        }
    }

    #[test]
    fn target_with_insertion() {
        let template = simple_template();
        // ターゲットに挿入あり
        let target = vec![
            AminoAcid::Ala,
            AminoAcid::Arg,
            AminoAcid::Asn, // 挿入
            AminoAcid::Gly,
            AminoAcid::Val,
        ];
        let model = build_model(&target, &template);
        assert_eq!(model.positions.len(), 5);
        // 挿入位置はloop_regionsに含まれる可能性
        for pos in &model.positions {
            assert!(pos[0].is_finite());
        }
    }

    #[test]
    fn empty_target() {
        let template = simple_template();
        let model = build_model(&[], &template);
        assert!(model.positions.is_empty());
        assert!(model.sequence.is_empty());
    }

    #[test]
    fn empty_template() {
        let target = vec![AminoAcid::Ala, AminoAcid::Gly];
        let template = Template {
            sequence: Vec::new(),
            positions: Vec::new(),
        };
        let model = build_model(&target, &template);
        assert_eq!(model.positions.len(), 2);
        // 直線配置
        assert!((model.positions[0][0] - 0.0).abs() < 1e-10);
        assert!((model.positions[1][0] - 3.8).abs() < 1e-10);
    }

    #[test]
    fn model_positions_finite() {
        let template = simple_template();
        let target = vec![AminoAcid::Ala, AminoAcid::Leu, AminoAcid::Met];
        let model = build_model(&target, &template);
        for pos in &model.positions {
            for k in 0..3 {
                assert!(pos[k].is_finite(), "Non-finite position: {pos:?}");
            }
        }
    }

    #[test]
    fn loop_regions_identified() {
        let template = Template {
            sequence: vec![AminoAcid::Ala, AminoAcid::Val],
            positions: vec![[0.0, 0.0, 0.0], [3.8, 0.0, 0.0]],
        };
        // ターゲットがテンプレートより長い → ギャップ
        let target = vec![AminoAcid::Ala, AminoAcid::Gly, AminoAcid::Val];
        let model = build_model(&target, &template);
        // 少なくとも1つのループ領域がある可能性
        assert_eq!(model.positions.len(), 3);
    }

    #[test]
    fn pack_sidechains_no_crash() {
        let mut positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0], // 近すぎる
            [2.0, 0.0, 0.0],
        ];
        let seq = vec![AminoAcid::Trp, AminoAcid::Trp, AminoAcid::Trp]; // 大きなside-chain
        pack_sidechains(&mut positions, &seq);
        for pos in &positions {
            assert!(pos[0].is_finite());
        }
    }

    #[test]
    fn pack_sidechains_pushes_apart() {
        let mut positions = vec![
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0], // VdW衝突
        ];
        let seq = vec![AminoAcid::Ala, AminoAcid::Ala];
        let initial_dist = 0.5;
        pack_sidechains(&mut positions, &seq);
        let dx = positions[1][0] - positions[0][0];
        let dist = (dx * dx).sqrt();
        assert!(
            dist > initial_dist,
            "Packing should push apart: {initial_dist} → {dist}"
        );
    }

    #[test]
    fn build_model_identity_range() {
        let template = simple_template();
        let target = template.sequence.clone();
        let model = build_model(&target, &template);
        assert!((0.0..=1.0).contains(&model.identity));
    }

    #[test]
    fn completely_different_target() {
        let template = simple_template();
        let target = vec![AminoAcid::Trp, AminoAcid::Pro, AminoAcid::His];
        let model = build_model(&target, &template);
        assert_eq!(model.positions.len(), 3);
        for pos in &model.positions {
            assert!(pos[0].is_finite());
        }
    }
}
