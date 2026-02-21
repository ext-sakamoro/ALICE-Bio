//! Energy potential functions for molecular interactions.

use std::f64::consts::PI;

/// Lennard-Jones potential: V = 4ε[(σ/r)^12 - (σ/r)^6]
///
/// Uses a pre-computed reciprocal `inv_r = 1/r` to replace division with multiplication,
/// reducing latency for callers that already have r available.
#[inline(always)]
pub fn lennard_jones(r: f64, epsilon: f64, sigma: f64) -> f64 {
    let inv_r = 1.0 / r; // single division; subsequent ops are multiplications
    let s_over_r = sigma * inv_r;
    let s6 = s_over_r * s_over_r * s_over_r * s_over_r * s_over_r * s_over_r;
    4.0 * epsilon * (s6 * s6 - s6)
}

/// Coulomb potential: V = (q1*q2) / (4π*ε*r)
///
/// Pre-computes `inv_r` to make the denominator a multiply chain.
#[inline(always)]
pub fn coulomb(q1: f64, q2: f64, r: f64, dielectric: f64) -> f64 {
    // Pre-computed constant reciprocal: 1/(4π) ≈ 0.07957747…
    const INV_4PI: f64 = 1.0 / (4.0 * PI);
    let inv_r = 1.0 / r;
    q1 * q2 * INV_4PI * (1.0 / dielectric) * inv_r
}

/// Simplified hydrogen bond potential: depth * [(r_eq/r)^12 - 2*(r_eq/r)^6]
///
/// Uses `inv_r` to replace division with multiplication.
#[inline(always)]
pub fn hydrogen_bond(r: f64, r_eq: f64, depth: f64) -> f64 {
    let inv_r = 1.0 / r; // single division
    let ratio = r_eq * inv_r;
    let r6 = ratio * ratio * ratio * ratio * ratio * ratio;
    depth * (r6 * r6 - 2.0 * r6)
}

/// Torsion potential: V = (barrier/2) * [1 + cos(n*angle - phase)]
///
/// Uses multiplication by 0.5 instead of division by 2.
#[inline(always)]
pub fn torsion_potential(angle: f64, barrier: f64, n: u32, phase: f64) -> f64 {
    (barrier * 0.5) * (1.0 + (n as f64 * angle - phase).cos())
}

/// Aggregated energy from different interaction types.
#[derive(Debug, Clone, Default)]
pub struct TotalEnergy {
    pub van_der_waals: f64,
    pub electrostatic: f64,
    pub hydrogen_bonds: f64,
    pub torsional: f64,
}

impl TotalEnergy {
    #[inline]
    pub fn total(&self) -> f64 {
        self.van_der_waals + self.electrostatic + self.hydrogen_bonds + self.torsional
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lj_minimum_at_equilibrium() {
        // At r = sigma * 2^(1/6), LJ reaches minimum = -epsilon
        let sigma = 3.4;
        let epsilon = 1.0;
        let r_min = sigma * 2.0_f64.powf(1.0 / 6.0);
        let v = lennard_jones(r_min, epsilon, sigma);
        assert!((v - (-epsilon)).abs() < 1e-10);
    }

    #[test]
    fn lj_large_at_small_r() {
        let v = lennard_jones(0.1, 1.0, 3.4);
        assert!(v > 1e10);
    }

    #[test]
    fn coulomb_same_sign_positive() {
        let v = coulomb(1.0, 1.0, 5.0, 1.0);
        assert!(v > 0.0);
    }

    #[test]
    fn coulomb_opposite_sign_negative() {
        let v = coulomb(1.0, -1.0, 5.0, 1.0);
        assert!(v < 0.0);
    }

    #[test]
    fn hbond_minimum_at_equilibrium() {
        // At r = r_eq, hydrogen_bond = depth * (1 - 2) = -depth
        let v = hydrogen_bond(2.8, 2.8, 5.0);
        assert!((v - (-5.0)).abs() < 1e-10);
    }

    #[test]
    fn torsion_periodicity() {
        let v1 = torsion_potential(0.0, 2.0, 3, 0.0);
        let v2 = torsion_potential(2.0 * PI / 3.0, 2.0, 3, 0.0);
        assert!((v1 - v2).abs() < 1e-10);
    }

    #[test]
    fn total_energy_sum() {
        let e = TotalEnergy {
            van_der_waals: 1.0,
            electrostatic: 2.0,
            hydrogen_bonds: -3.0,
            torsional: 0.5,
        };
        assert!((e.total() - 0.5).abs() < 1e-10);
    }
}
