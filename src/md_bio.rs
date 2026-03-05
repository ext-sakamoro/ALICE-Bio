//! 生体分子MD — Velocity Verlet, Langevin thermostat, SHAKE拘束
//!
//! タンパク質系の分子動力学シミュレーション基盤。
//! Cα粗視化モデル上での NVT (Langevin) シミュレーションを実装。

use crate::amino::AminoAcid;

/// Boltzmann定数 (kcal/mol/K)。
const KB: f64 = 0.001_987_204;

// ============================================================================
// 粒子
// ============================================================================

/// Cα 粗視化粒子。
#[derive(Debug, Clone)]
pub struct BioParticle {
    /// アミノ酸タイプ。
    pub amino: AminoAcid,
    /// 位置 (Å)。
    pub pos: [f64; 3],
    /// 速度 (Å/fs)。
    pub vel: [f64; 3],
    /// 力 (kcal/mol/Å)。
    pub force: [f64; 3],
    /// 質量 (Da)。
    pub mass: f64,
}

impl BioParticle {
    /// 新しい粒子を生成。
    #[must_use]
    pub const fn new(amino: AminoAcid, pos: [f64; 3]) -> Self {
        Self {
            amino,
            pos,
            vel: [0.0; 3],
            force: [0.0; 3],
            mass: amino.mass(),
        }
    }
}

// ============================================================================
// LJ力場計算 (Cα粗視化)
// ============================================================================

/// Cα 間のLJ力を計算して particles.force に蓄積。
///
/// σ はVdW半径の算術平均、ε は固定値。
fn compute_forces(particles: &mut [BioParticle], epsilon: f64) {
    // 力リセット
    for p in particles.iter_mut() {
        p.force = [0.0; 3];
    }

    let n = particles.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = particles[j].pos[0] - particles[i].pos[0];
            let dy = particles[j].pos[1] - particles[i].pos[1];
            let dz = particles[j].pos[2] - particles[i].pos[2];
            let dist_sq = dz.mul_add(dz, dx.mul_add(dx, dy * dy));

            if dist_sq < 0.01 {
                continue; // 距離が近すぎる場合はスキップ
            }

            let sigma = (particles[i].amino.van_der_waals_radius()
                + particles[j].amino.van_der_waals_radius())
                * 0.5;

            let inv_r2 = 1.0 / dist_sq;
            let sigma2 = sigma * sigma;
            let s2_over_r2 = sigma2 * inv_r2;
            let s6 = s2_over_r2 * s2_over_r2 * s2_over_r2;
            let s12 = s6 * s6;

            // F = -dV/dr * (r_vec/r) = 24ε/r² * [2(σ/r)^12 - (σ/r)^6] * r_vec
            let f_mag = 24.0 * epsilon * inv_r2 * 2.0f64.mul_add(s12, -s6);

            let fx = f_mag * dx;
            let fy = f_mag * dy;
            let fz = f_mag * dz;

            particles[i].force[0] += fx;
            particles[i].force[1] += fy;
            particles[i].force[2] += fz;
            particles[j].force[0] -= fx;
            particles[j].force[1] -= fy;
            particles[j].force[2] -= fz;
        }
    }
}

// ============================================================================
// Velocity Verlet 積分
// ============================================================================

/// Velocity Verlet 1ステップ (NVE)。
///
/// `dt` はタイムステップ (fs)、`epsilon` はLJ ε (kcal/mol)。
pub fn velocity_verlet_step(particles: &mut [BioParticle], dt: f64, epsilon: f64) {
    let half_dt = dt * 0.5;

    // 位置更新: x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt²
    for p in particles.iter_mut() {
        let inv_mass = 1.0 / p.mass;
        for k in 0..3 {
            p.vel[k] += half_dt * p.force[k] * inv_mass;
            p.pos[k] += dt * p.vel[k];
        }
    }

    // 力の再計算: F(t+dt)
    compute_forces(particles, epsilon);

    // 速度更新の後半: v(t+dt) = v(t+dt/2) + 0.5*a(t+dt)*dt
    for p in particles.iter_mut() {
        let inv_mass = 1.0 / p.mass;
        for k in 0..3 {
            p.vel[k] += half_dt * p.force[k] * inv_mass;
        }
    }
}

// ============================================================================
// Langevin Thermostat
// ============================================================================

/// Langevin thermostat 設定。
#[derive(Debug, Clone, Copy)]
pub struct LangevinConfig {
    /// 摩擦係数 γ (1/fs)。
    pub gamma: f64,
    /// 目標温度 (K)。
    pub temperature: f64,
    /// タイムステップ (fs)。
    pub dt: f64,
    /// LJ ε (kcal/mol)。
    pub epsilon: f64,
}

impl Default for LangevinConfig {
    fn default() -> Self {
        Self {
            gamma: 0.01,
            temperature: 300.0,
            dt: 1.0,
            epsilon: 0.1,
        }
    }
}

/// Langevin dynamics 1ステップ。
///
/// 乱数シードから簡易擬似乱数を生成し、確率的力を適用。
pub fn langevin_step(particles: &mut [BioParticle], config: &LangevinConfig, seed: u64) {
    let dt = config.dt;
    let half_dt = dt * 0.5;
    let gamma = config.gamma;

    let sigma_force = (2.0 * gamma * KB * config.temperature / dt).sqrt();

    // 力の計算
    compute_forces(particles, config.epsilon);

    let mut rng_state = seed;

    for p in particles.iter_mut() {
        let inv_mass = 1.0 / p.mass;
        let mass_sqrt = p.mass.sqrt();

        for k in 0..3 {
            // 擬似ガウスノイズ (Box-Muller近似: 2つの一様乱数から)
            rng_state = rng_state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            let u1 = (rng_state >> 32) as f64 / u32::MAX as f64;
            rng_state = rng_state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            let u2 = (rng_state >> 32) as f64 / u32::MAX as f64;

            // Box-Muller
            let u1_clamped = u1.max(1e-15);
            let gauss = (-2.0 * u1_clamped.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();

            let random_force = sigma_force * mass_sqrt * gauss;
            let friction_force = -gamma * p.mass * p.vel[k];
            let total_force = p.force[k] + friction_force + random_force;

            // Velocity Verlet with Langevin
            p.vel[k] += half_dt * total_force * inv_mass;
            p.pos[k] += dt * p.vel[k];
        }
    }

    // 力再計算
    compute_forces(particles, config.epsilon);

    rng_state = rng_state.wrapping_add(42);

    for p in particles.iter_mut() {
        let inv_mass = 1.0 / p.mass;
        let mass_sqrt = p.mass.sqrt();
        for k in 0..3 {
            rng_state = rng_state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            let u1 = (rng_state >> 32) as f64 / u32::MAX as f64;
            rng_state = rng_state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            let u2 = (rng_state >> 32) as f64 / u32::MAX as f64;
            let u1_clamped = u1.max(1e-15);
            let gauss = (-2.0 * u1_clamped.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();

            let random_force = sigma_force * mass_sqrt * gauss;
            let friction_force = -gamma * p.mass * p.vel[k];
            let total_force = p.force[k] + friction_force + random_force;
            p.vel[k] += half_dt * total_force * inv_mass;
        }
    }
}

// ============================================================================
// SHAKE 拘束
// ============================================================================

/// SHAKE拘束: Cα-Cα結合長を固定距離に拘束。
///
/// `target_dist` は拘束距離 (Å)、`max_iter` は反復回数上限。
/// `tolerance` は収束判定閾値 (Å²)。
pub fn shake_constraints(
    particles: &mut [BioParticle],
    target_dist: f64,
    max_iter: usize,
    tolerance: f64,
) {
    if particles.len() < 2 {
        return;
    }

    let target_sq = target_dist * target_dist;

    for _ in 0..max_iter {
        let mut converged = true;

        for i in 0..(particles.len() - 1) {
            let j = i + 1;
            let dx = particles[j].pos[0] - particles[i].pos[0];
            let dy = particles[j].pos[1] - particles[i].pos[1];
            let dz = particles[j].pos[2] - particles[i].pos[2];
            let dist_sq = dz.mul_add(dz, dx.mul_add(dx, dy * dy));

            let diff = dist_sq - target_sq;
            if diff.abs() < tolerance {
                continue;
            }
            converged = false;

            let inv_mass_i = 1.0 / particles[i].mass;
            let inv_mass_j = 1.0 / particles[j].mass;
            let lambda = diff / (2.0 * (inv_mass_i + inv_mass_j) * dist_sq);

            let corr_x = lambda * dx;
            let corr_y = lambda * dy;
            let corr_z = lambda * dz;

            particles[i].pos[0] += inv_mass_i * corr_x;
            particles[i].pos[1] += inv_mass_i * corr_y;
            particles[i].pos[2] += inv_mass_i * corr_z;
            particles[j].pos[0] -= inv_mass_j * corr_x;
            particles[j].pos[1] -= inv_mass_j * corr_y;
            particles[j].pos[2] -= inv_mass_j * corr_z;
        }

        if converged {
            break;
        }
    }
}

// ============================================================================
// ユーティリティ
// ============================================================================

/// 系の運動エネルギーを計算 (kcal/mol)。
#[must_use]
pub fn kinetic_energy(particles: &[BioParticle]) -> f64 {
    let mut ke = 0.0;
    for p in particles {
        let v2 = p.vel[2].mul_add(p.vel[2], p.vel[0].mul_add(p.vel[0], p.vel[1] * p.vel[1]));
        ke += 0.5 * p.mass * v2;
    }
    ke
}

/// 瞬間温度 (K)。
///
/// T = 2*KE / (3*N*kB)
#[must_use]
pub fn instantaneous_temperature(particles: &[BioParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let ke = kinetic_energy(particles);
    let dof = 3 * particles.len();
    2.0 * ke / (dof as f64 * KB)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn two_particle_system() -> Vec<BioParticle> {
        vec![
            BioParticle::new(AminoAcid::Ala, [0.0, 0.0, 0.0]),
            BioParticle::new(AminoAcid::Gly, [5.0, 0.0, 0.0]),
        ]
    }

    #[test]
    fn particle_new_default_vel_zero() {
        let p = BioParticle::new(AminoAcid::Ala, [1.0, 2.0, 3.0]);
        assert!(p.vel[0].abs() < 1e-15);
        assert!(p.vel[1].abs() < 1e-15);
        assert!(p.vel[2].abs() < 1e-15);
        assert!((p.mass - AminoAcid::Ala.mass()).abs() < 1e-10);
    }

    #[test]
    fn compute_forces_two_particles() {
        let mut particles = two_particle_system();
        compute_forces(&mut particles, 0.1);
        // 粒子0は粒子1方向に引っ張られるか反発
        // 5.0ÅでLJ minimum ~2.15σ → σ=(1.9+1.7)/2=1.8 → 2^(1/6)*1.8≈2.02
        // 5.0 > 2.02 → 引力領域
        assert!(particles[0].force[0] != 0.0);
        // ニュートン第三法則
        assert!((particles[0].force[0] + particles[1].force[0]).abs() < 1e-10);
    }

    #[test]
    fn velocity_verlet_conserves_momentum() {
        let mut particles = two_particle_system();
        // 初期運動量 = 0 (全速度ゼロ)
        compute_forces(&mut particles, 0.1);
        for _ in 0..10 {
            velocity_verlet_step(&mut particles, 1.0, 0.1);
        }
        let total_px: f64 = particles.iter().map(|p| p.mass * p.vel[0]).sum();
        let total_py: f64 = particles.iter().map(|p| p.mass * p.vel[1]).sum();
        let total_pz: f64 = particles.iter().map(|p| p.mass * p.vel[2]).sum();
        assert!(total_px.abs() < 1e-6, "px = {total_px}");
        assert!(total_py.abs() < 1e-6, "py = {total_py}");
        assert!(total_pz.abs() < 1e-6, "pz = {total_pz}");
    }

    #[test]
    fn velocity_verlet_particles_move() {
        let mut particles = two_particle_system();
        let initial_pos = particles[0].pos;
        compute_forces(&mut particles, 0.1);
        for _ in 0..10 {
            velocity_verlet_step(&mut particles, 1.0, 0.1);
        }
        assert!(
            (particles[0].pos[0] - initial_pos[0]).abs() > 1e-10,
            "Particle should have moved"
        );
    }

    #[test]
    fn langevin_step_runs_without_panic() {
        let mut particles = two_particle_system();
        let config = LangevinConfig::default();
        for step in 0..10 {
            langevin_step(&mut particles, &config, step);
        }
        // 位置は有限値
        for p in &particles {
            assert!(p.pos[0].is_finite());
            assert!(p.pos[1].is_finite());
            assert!(p.pos[2].is_finite());
        }
    }

    #[test]
    fn langevin_different_seeds_different_trajectories() {
        let mut p1 = two_particle_system();
        let mut p2 = two_particle_system();
        let config = LangevinConfig::default();
        langevin_step(&mut p1, &config, 42);
        langevin_step(&mut p2, &config, 99);
        // 異なるシードで異なる位置
        assert!(
            (p1[0].pos[0] - p2[0].pos[0]).abs() > 1e-15
                || (p1[0].pos[1] - p2[0].pos[1]).abs() > 1e-15
        );
    }

    #[test]
    fn shake_constrains_distance() {
        let mut particles = vec![
            BioParticle::new(AminoAcid::Ala, [0.0, 0.0, 0.0]),
            BioParticle::new(AminoAcid::Gly, [4.5, 0.0, 0.0]), // 少しずれた初期位置
        ];
        let target = 3.8;
        shake_constraints(&mut particles, target, 100, 1e-10);

        let dx = particles[1].pos[0] - particles[0].pos[0];
        let dy = particles[1].pos[1] - particles[0].pos[1];
        let dz = particles[1].pos[2] - particles[0].pos[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!(
            (dist - target).abs() < 0.01,
            "SHAKE should constrain to {target}, got {dist}"
        );
    }

    #[test]
    fn shake_single_particle_no_panic() {
        let mut particles = vec![BioParticle::new(AminoAcid::Ala, [0.0, 0.0, 0.0])];
        shake_constraints(&mut particles, 3.8, 100, 1e-10);
    }

    #[test]
    fn shake_chain_all_constrained() {
        let mut particles: Vec<BioParticle> = (0..5)
            .map(|i| BioParticle::new(AminoAcid::Ala, [i as f64 * 4.0, 0.0, 0.0]))
            .collect();
        let target = 3.8;
        shake_constraints(&mut particles, target, 200, 1e-10);

        for i in 0..4 {
            let dx = particles[i + 1].pos[0] - particles[i].pos[0];
            let dy = particles[i + 1].pos[1] - particles[i].pos[1];
            let dz = particles[i + 1].pos[2] - particles[i].pos[2];
            let dist = dz.mul_add(dz, dx.mul_add(dx, dy * dy)).sqrt();
            assert!(
                (dist - target).abs() < 0.05,
                "Bond {i}-{}: dist = {dist}, target = {target}",
                i + 1
            );
        }
    }

    #[test]
    fn kinetic_energy_zero_at_rest() {
        let particles = two_particle_system();
        assert!((kinetic_energy(&particles) - 0.0).abs() < 1e-15);
    }

    #[test]
    fn kinetic_energy_positive_with_velocity() {
        let mut particles = two_particle_system();
        particles[0].vel = [1.0, 0.0, 0.0];
        assert!(kinetic_energy(&particles) > 0.0);
    }

    #[test]
    fn temperature_zero_at_rest() {
        let particles = two_particle_system();
        assert!((instantaneous_temperature(&particles) - 0.0).abs() < 1e-15);
    }

    #[test]
    fn temperature_empty_system() {
        assert!((instantaneous_temperature(&[]) - 0.0).abs() < 1e-15);
    }

    #[test]
    fn temperature_positive_with_velocity() {
        let mut particles = two_particle_system();
        particles[0].vel = [0.1, 0.0, 0.0];
        particles[1].vel = [-0.1, 0.0, 0.0];
        assert!(instantaneous_temperature(&particles) > 0.0);
    }

    #[test]
    fn langevin_config_default() {
        let cfg = LangevinConfig::default();
        assert!((cfg.gamma - 0.01).abs() < 1e-10);
        assert!((cfg.temperature - 300.0).abs() < 1e-10);
    }
}
