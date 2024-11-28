use super::{eos::EoS};
use self::eq::euler::conserved;
use self::{eq::RHS, bc::BCHD};

mod eq;
mod bc;

const NF: usize = 4;

pub trait Marching {
  fn calc_cfl(&self, cfl_coeff: &f64) -> f64;
  fn march_ssprk3(&mut self, dt: &f64, bc_type: &str, reconstruction: &str, flux_scheme: &str);
}


impl<const NI: usize, const NJ: usize, const NB: usize> Marching
for super::IdealGas::<NI,NJ,NB> {
  /// calc dt thet meets CFL condition
  #[inline]
  fn calc_cfl(&self, cfl_coeff: &f64) -> f64 {
    let mut nu: f64 = 1.0e+10;
    for i in 0..(NI-2*NB) {
      for j in 0..(NJ-2*NB) {
        let tmp: f64 = self.coord.dx[i][j] /
        (self.eos.calc_cs(&self.basic.rho[i+NB][j+NB], &self.basic.u[i+NB][j+NB], &self.basic.v[i+NB][j+NB], &self.basic.e[i+NB][j+NB]) + (self.basic.u[i+NB][j+NB] * self.basic.u[i+NB][j+NB] + self.basic.v[i+NB][j+NB] * self.basic.v[i+NB][j+NB]).sqrt());
        nu = nu.min(tmp);
      }
    }
    cfl_coeff * nu
  }


  /// marching dt with 3rd order SSP Rungr-Kutta method
  fn march_ssprk3(&mut self, dt: &f64, bc_type: &str, reconstruction: &str, flux_scheme: &str) {
    // construct conservative var from basic var
    for i in 0..(NI-2*NB) {
      for j in 0..(NJ-2*NB) {
        let s_a = 0.25 * (self.coord.s[NB+i-1][NB+j-1] + self.coord.s[NB+i][NB+j-1] + self.coord.s[NB+i-1][NB+j] + self.coord.s[NB+i][NB+j]);
        conserved::calc_conservative(&self.basic.rho[NB+i][NB+j], &self.basic.u[NB+i][NB+j], &self.basic.v[NB+i][NB+j], &self.basic.e[NB+i][NB+j], &s_a, &mut self.arr_q0[i][j]);
      }
    }
    // 1st stage
    self.eq.calc_rhs(reconstruction, flux_scheme, &self.basic.rho, &self.basic.u, &self.basic.v, &self.basic.e, &self.coord, &self.eos, &mut self.arr_q1);
    for i in 0..(NI-2*NB) {
      for j in 0..(NJ-2*NB) {
        for k in 0..NF {
          self.arr_q1[i][j][k] = self.arr_q0[i][j][k] + dt * self.arr_q1[i][j][k];
        }
        let s_a = 0.25 * (self.coord.s[NB+i-1][NB+j-1] + self.coord.s[NB+i][NB+j-1] + self.coord.s[NB+i-1][NB+j] + self.coord.s[NB+i][NB+j]);
        conserved::calc_basic(&self.arr_q1[i][j], &s_a, &mut self.basic.rho[NB+i][NB+j], &mut self.basic.u[NB+i][NB+j], &mut self.basic.v[NB+i][NB+j], &mut self.basic.e[NB+i][NB+j]);
      }
    }
    self.bound.reflect_bc(bc_type, &self.eos, &mut self.basic);

    // 2nd stage
    self.eq.calc_rhs(reconstruction, flux_scheme, &self.basic.rho, &self.basic.u, &self.basic.v, &self.basic.e, &self.coord, &self.eos, &mut self.arr_q2);
    for i in 0..(NI-2*NB) {
      for j in 0..(NJ-2*NB) {
        for k in 0..NF {
          self.arr_q2[i][j][k] = 0.75 * self.arr_q0[i][j][k] + 0.25 * (self.arr_q1[i][j][k] + dt * self.arr_q2[i][j][k]);
        }
        let s_a = 0.25 * (self.coord.s[NB+i-1][NB+j-1] + self.coord.s[NB+i][NB+j-1] + self.coord.s[NB+i-1][NB+j] + self.coord.s[NB+i][NB+j]);
        conserved::calc_basic(&self.arr_q2[i][j], &s_a, &mut self.basic.rho[NB+i][NB+j], &mut self.basic.u[NB+i][NB+j], &mut self.basic.v[NB+i][NB+j], &mut self.basic.e[NB+i][NB+j]);
      }
    }
    self.bound.reflect_bc(bc_type, &self.eos, &mut self.basic);

    // 3rd stage
    self.eq.calc_rhs(reconstruction, flux_scheme, &self.basic.rho, &self.basic.u, &self.basic.v, &self.basic.e, &self.coord, &self.eos, &mut self.arr_q1);
    for i in 0..(NI-2*NB) {
      for j in 0..(NJ-2*NB) {
        for k in 0..NF {
          self.arr_q0[i][j][k] = self.arr_q0[i][j][k] / 3.0 + 2.0 / 3.0 * (self.arr_q2[i][j][k] + dt * self.arr_q1[i][j][k]);
        }
        let s_a = 0.25 * (self.coord.s[NB+i-1][NB+j-1] + self.coord.s[NB+i][NB+j-1] + self.coord.s[NB+i-1][NB+j] + self.coord.s[NB+i][NB+j]);
        conserved::calc_basic(&self.arr_q0[i][j], &s_a, &mut self.basic.rho[NB+i][NB+j], &mut self.basic.u[NB+i][NB+j], &mut self.basic.v[NB+i][NB+j], &mut self.basic.e[NB+i][NB+j]);
      }
    }
    self.bound.reflect_bc(bc_type, &self.eos, &mut self.basic);
  }
}