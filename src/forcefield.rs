//! AMBER力場パラメータ — bond stretching, angle bending, dihedral, LJ 14
//!
//! AMBER ff14SB 力場の主要パラメータを提供する。
//! 結合伸縮・角度変形・二面角・1-4 LJ相互作用のエネルギー関数。

// ============================================================================
// 原子タイプ
// ============================================================================

/// AMBER原子タイプ (主要骨格原子)。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomType {
    /// Cα (sp3 carbon)。
    CA,
    /// カルボニル炭素 C。
    C,
    /// アミド窒素 N。
    N,
    /// カルボニル酸素 O。
    O,
    /// 水素 H (N-H)。
    H,
    /// Cβ (aliphatic carbon)。
    CB,
}

/// `原子タイプのVdW` パラメータ (σ: Å, ε: kcal/mol)。
#[must_use]
pub const fn lj_params(atom: AtomType) -> (f64, f64) {
    match atom {
        #[allow(clippy::match_same_arms)] // CA and C share physical VdW values
        AtomType::CA => (3.40, 0.086),
        AtomType::C => (3.40, 0.086),
        AtomType::N => (3.25, 0.170),
        AtomType::O => (2.96, 0.210),
        AtomType::H => (1.07, 0.016),
        AtomType::CB => (3.40, 0.109),
    }
}

// ============================================================================
// Bond Stretching: E = K_b * (r - r_eq)^2
// ============================================================================

/// 結合伸縮パラメータ。
#[derive(Debug, Clone, Copy)]
pub struct BondParam {
    /// 力定数 (kcal/mol/Å²)。
    pub k: f64,
    /// 平衡距離 (Å)。
    pub r_eq: f64,
}

/// 主要結合タイプのパラメータ (AMBER近似値)。
#[must_use]
pub const fn bond_param(a: AtomType, b: AtomType) -> BondParam {
    // 対称化: (a,b) と (b,a) を同じに扱う
    let pair = normalize_pair(a, b);
    match pair {
        (AtomType::CA, AtomType::C) | (AtomType::C, AtomType::CA) => BondParam {
            k: 317.0,
            r_eq: 1.522,
        },
        (AtomType::CA, AtomType::N) | (AtomType::N, AtomType::CA) => BondParam {
            k: 337.0,
            r_eq: 1.449,
        },
        (AtomType::C, AtomType::N) | (AtomType::N, AtomType::C) => BondParam {
            k: 490.0,
            r_eq: 1.335,
        },
        (AtomType::C, AtomType::O) | (AtomType::O, AtomType::C) => BondParam {
            k: 570.0,
            r_eq: 1.229,
        },
        (AtomType::N, AtomType::H) | (AtomType::H, AtomType::N) => BondParam {
            k: 434.0,
            r_eq: 1.010,
        },
        (AtomType::CA, AtomType::CB) | (AtomType::CB, AtomType::CA) => BondParam {
            k: 310.0,
            r_eq: 1.526,
        },
        _ => BondParam {
            k: 300.0,
            r_eq: 1.50,
        }, // デフォルト近似
    }
}

/// 結合伸縮エネルギー: E = `K_b` * (r - `r_eq)²`
#[inline]
#[must_use]
pub fn bond_energy(r: f64, param: &BondParam) -> f64 {
    let dr = r - param.r_eq;
    param.k * dr * dr
}

// ============================================================================
// Angle Bending: E = K_a * (θ - θ_eq)^2
// ============================================================================

/// 角度変形パラメータ。
#[derive(Debug, Clone, Copy)]
pub struct AngleParam {
    /// 力定数 (kcal/mol/rad²)。
    pub k: f64,
    /// 平衡角度 (rad)。
    pub theta_eq: f64,
}

/// 主要角度タイプのパラメータ (AMBER近似値)。
#[must_use]
pub const fn angle_param(a: AtomType, b: AtomType, c: AtomType) -> AngleParam {
    match (a, b, c) {
        // N-CA-C (backbone)
        (AtomType::N, AtomType::CA, AtomType::C) | (AtomType::C, AtomType::CA, AtomType::N) => {
            AngleParam {
                k: 63.0,
                theta_eq: 1.9146,
            } // 109.7°
        }
        // CA-C-N (peptide bond)
        (AtomType::CA, AtomType::C, AtomType::N) | (AtomType::N, AtomType::C, AtomType::CA) => {
            AngleParam {
                k: 70.0,
                theta_eq: 2.0350,
            } // 116.6°
        }
        // C-N-CA
        (AtomType::C, AtomType::N, AtomType::CA) | (AtomType::CA, AtomType::N, AtomType::C) => {
            AngleParam {
                k: 50.0,
                theta_eq: 2.1240,
            } // 121.7°
        }
        // CA-C-O
        (AtomType::CA, AtomType::C, AtomType::O) | (AtomType::O, AtomType::C, AtomType::CA) => {
            AngleParam {
                k: 80.0,
                theta_eq: 2.1014,
            } // 120.4°
        }
        _ => AngleParam {
            k: 50.0,
            theta_eq: 1.9111,
        }, // デフォルト ~109.5°
    }
}

/// 角度変形エネルギー: E = `K_a` * (θ - `θ_eq)²`
#[inline]
#[must_use]
pub fn angle_energy(theta: f64, param: &AngleParam) -> f64 {
    let dt = theta - param.theta_eq;
    param.k * dt * dt
}

// ============================================================================
// Dihedral: E = (V_n / 2) * [1 + cos(n*φ - γ)]
// ============================================================================

/// 二面角パラメータ。
#[derive(Debug, Clone, Copy)]
pub struct DihedralParam {
    /// バリア高さ `V_n` (kcal/mol)。
    pub barrier: f64,
    /// 多重度 n。
    pub periodicity: u32,
    /// 位相角 γ (rad)。
    pub phase: f64,
}

/// 主要二面角パラメータ (AMBER backbone)。
#[must_use]
pub const fn backbone_dihedral_params() -> [DihedralParam; 3] {
    [
        // phi (C-N-CA-C)
        DihedralParam {
            barrier: 0.80,
            periodicity: 1,
            phase: 0.0,
        },
        // psi (N-CA-C-N)
        DihedralParam {
            barrier: 0.65,
            periodicity: 2,
            phase: std::f64::consts::PI,
        },
        // omega (CA-C-N-CA) — ペプチド結合平面性
        DihedralParam {
            barrier: 10.5,
            periodicity: 2,
            phase: std::f64::consts::PI,
        },
    ]
}

/// 二面角エネルギー: E = (`V_n` / 2) * [1 + cos(n*φ - γ)]
#[inline]
#[must_use]
pub fn dihedral_energy(phi: f64, param: &DihedralParam) -> f64 {
    (param.barrier * 0.5) * (1.0 + (param.periodicity as f64).mul_add(phi, -param.phase).cos())
}

// ============================================================================
// 1-4 LJ 相互作用
// ============================================================================

/// 1-4 LJ 相互作用エネルギー (AMBER: scee=1.2, scnb=2.0)。
///
/// E = (1/scnb) * 4ε[(σ/r)^12 - (σ/r)^6]
#[inline]
#[must_use]
pub fn lj_14_energy(r: f64, atom_i: AtomType, atom_j: AtomType) -> f64 {
    const SCNB: f64 = 2.0; // AMBER 1-4 VdW scaling

    let (sigma_i, eps_i) = lj_params(atom_i);
    let (sigma_j, eps_j) = lj_params(atom_j);

    // Lorentz-Berthelot混合則
    let sigma = (sigma_i + sigma_j) * 0.5;
    let epsilon = (eps_i * eps_j).sqrt();

    let inv_r = 1.0 / r;
    let s_over_r = sigma * inv_r;
    let s6 = s_over_r * s_over_r * s_over_r * s_over_r * s_over_r * s_over_r;
    (1.0 / SCNB) * 4.0 * epsilon * s6.mul_add(s6, -s6)
}

// ============================================================================
// 内部ユーティリティ
// ============================================================================

/// 原子ペアの正規化 (辞書順)。
const fn normalize_pair(a: AtomType, b: AtomType) -> (AtomType, AtomType) {
    if (a as u8) <= (b as u8) {
        (a, b)
    } else {
        (b, a)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn bond_energy_at_equilibrium_is_zero() {
        let p = bond_param(AtomType::CA, AtomType::C);
        let e = bond_energy(p.r_eq, &p);
        assert!(e.abs() < 1e-15);
    }

    #[test]
    fn bond_energy_positive_away_from_eq() {
        let p = bond_param(AtomType::CA, AtomType::C);
        let e = bond_energy(p.r_eq + 0.1, &p);
        assert!(e > 0.0);
    }

    #[test]
    fn bond_energy_symmetric() {
        let p = bond_param(AtomType::CA, AtomType::C);
        let e_plus = bond_energy(p.r_eq + 0.05, &p);
        let e_minus = bond_energy(p.r_eq - 0.05, &p);
        assert!((e_plus - e_minus).abs() < 1e-10);
    }

    #[test]
    fn angle_energy_at_equilibrium_is_zero() {
        let p = angle_param(AtomType::N, AtomType::CA, AtomType::C);
        let e = angle_energy(p.theta_eq, &p);
        assert!(e.abs() < 1e-15);
    }

    #[test]
    fn angle_energy_positive_away_from_eq() {
        let p = angle_param(AtomType::N, AtomType::CA, AtomType::C);
        let e = angle_energy(p.theta_eq + 0.1, &p);
        assert!(e > 0.0);
    }

    #[test]
    fn dihedral_energy_range() {
        let params = backbone_dihedral_params();
        for p in &params {
            for angle in [0.0, PI / 2.0, PI, -PI / 2.0] {
                let e = dihedral_energy(angle, p);
                assert!(e >= 0.0, "Dihedral energy should be >= 0, got {e}");
                assert!(
                    e <= p.barrier,
                    "Dihedral energy should be <= barrier, got {e}"
                );
            }
        }
    }

    #[test]
    fn dihedral_energy_at_minimum() {
        // cos(n*phi - phase) = -1 → E = 0
        let p = DihedralParam {
            barrier: 2.0,
            periodicity: 1,
            phase: 0.0,
        };
        let e = dihedral_energy(PI, &p);
        assert!(e.abs() < 1e-10);
    }

    #[test]
    fn lj_14_repulsive_at_short_range() {
        let e = lj_14_energy(1.0, AtomType::CA, AtomType::CA);
        assert!(e > 0.0);
    }

    #[test]
    fn lj_14_small_at_long_range() {
        let e = lj_14_energy(50.0, AtomType::CA, AtomType::CA);
        assert!(e.abs() < 1e-4);
    }

    #[test]
    fn lj_params_all_positive() {
        let types = [
            AtomType::CA,
            AtomType::C,
            AtomType::N,
            AtomType::O,
            AtomType::H,
            AtomType::CB,
        ];
        for t in &types {
            let (sigma, eps) = lj_params(*t);
            assert!(sigma > 0.0, "{t:?} sigma = {sigma}");
            assert!(eps > 0.0, "{t:?} eps = {eps}");
        }
    }

    #[test]
    fn bond_param_symmetric() {
        let p1 = bond_param(AtomType::CA, AtomType::C);
        let p2 = bond_param(AtomType::C, AtomType::CA);
        assert!((p1.k - p2.k).abs() < 1e-10);
        assert!((p1.r_eq - p2.r_eq).abs() < 1e-10);
    }

    #[test]
    fn backbone_dihedral_count() {
        let params = backbone_dihedral_params();
        assert_eq!(params.len(), 3);
    }

    #[test]
    fn lj_14_scaling_factor() {
        // 1-4 energy は full LJの半分 (scnb=2.0)
        let r = 4.0;
        let e14 = lj_14_energy(r, AtomType::CA, AtomType::N);
        // Full LJ
        let (si, ei) = lj_params(AtomType::CA);
        let (sj, ej) = lj_params(AtomType::N);
        let sigma = (si + sj) * 0.5;
        let epsilon = (ei * ej).sqrt();
        let inv_r = 1.0 / r;
        let s_over_r = sigma * inv_r;
        let s6 = s_over_r * s_over_r * s_over_r * s_over_r * s_over_r * s_over_r;
        let e_full = 4.0 * epsilon * s6.mul_add(s6, -s6);
        assert!((e14 - e_full * 0.5).abs() < 1e-10);
    }

    #[test]
    fn normalize_pair_symmetric() {
        let (a, b) = normalize_pair(AtomType::N, AtomType::CA);
        let (c, d) = normalize_pair(AtomType::CA, AtomType::N);
        assert_eq!(a, c);
        assert_eq!(b, d);
    }
}
