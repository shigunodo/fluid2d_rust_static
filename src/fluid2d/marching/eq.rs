pub mod euler;
mod fnd;
use super::super::{eos::EoS, GenStructCoord};

const NF: usize = 4;

pub trait RHS<const NI: usize, const NJ: usize, const NB: usize> {


  fn calc_rhs<T: EoS>(&mut self, reconstruction: &str, flux_scheme: &str, 
  rho: &[[f64; NJ]; NI], u: &[[f64; NJ]; NI], v: &[[f64; NJ]; NI],
  e: &[[f64; NJ]; NI], coord: &GenStructCoord::<NI,NJ,NB>, eos: &T, arr_q: &mut [[[f64; NF]; NJ]; NI]);

}




impl<const NI: usize, const NJ: usize, const NB: usize> RHS::<NI,NJ,NB>
for super::super::EulerEq::<NI,NJ,NB> {


  /// calc RHS of equation with selected reconstruction/flux scheme
  fn calc_rhs<T: EoS>(&mut self, reconstruction: &str, flux_scheme: &str, 
  rho: &[[f64; NJ]; NI], u: &[[f64; NJ]; NI], v: &[[f64; NJ]; NI],
  e: &[[f64; NJ]; NI], coord: &GenStructCoord::<NI,NJ,NB>, eos: &T, arr_q: &mut [[[f64; NF]; NJ]; NI]) {

    //--------------------i-direction---------------------
    // evaluating numerical flux at (i+0.5,j)
    for i in 0..(NI-2*NB+1) {
      for j in 0..(NJ-2*NB) {
        // inverse of Jacobian
        // evaluated at(i + 0.5, j)
        // by averaging the values at(i + 0.5, j + 0.5) and (i + 0.5, j - 0.5)
        let s_a = 0.5 * (coord.s[NB+i-1][NB+j-1] + coord.s[NB+i-1][NB+j]);
        // metrices devided by Jacobian
        // evaluated at (i+0.5,j)
        // by averaging the values at (i,j) and (i+1,j)
        let ixs_a = 0.5 * (coord.ixs[i][j+1] + coord.ixs[i+1][j+1]);
        let iys_a = 0.5 * (coord.iys[i][j+1] + coord.iys[i+1][j+1]);
        //let jxs_a = 0.5 * (coord.jxs[i][j+1] + coord.jxs[i+1][j+1]);
        //let jys_a = 0.5 * (coord.jys[i][j+1] + coord.jys[i+1][j+1]);
        // reconstruction
        let mut rho_l = 0.0;
        let mut u_l = 0.0;
        let mut v_l = 0.0;
        let mut e_l = 0.0;
        let mut rho_r = 0.0;
        let mut u_r = 0.0;
        let mut v_r = 0.0;
        let mut e_r = 0.0;
        match reconstruction {
          "MUSCL_minmod_basic" => euler::flux_scheme::reconst_by_basic_muscl(
            &rho[NB+i-2][NB+j], &rho[NB+i-1][NB+j], 
            &rho[NB+i][NB+j], &rho[NB+i+1][NB+j],
            &u[NB+i-2][NB+j], &u[NB+i-1][NB+j], 
            &u[NB+i][NB+j], &u[NB+i+1][NB+j],
            &v[NB+i-2][NB+j], &v[NB+i-1][NB+j], 
            &v[NB+i][NB+j], &v[NB+i+1][NB+j],
            &e[NB+i-2][NB+j], &e[NB+i-1][NB+j], 
            &e[NB+i][NB+j], &e[NB+i+1][NB+j],
            &mut rho_l, &mut u_l, &mut v_l, &mut e_l, 
            &mut rho_r, &mut u_r, &mut v_r, &mut e_r),
          "MP5_basic" => euler::flux_scheme::reconst_by_basic_mp5(
            &rho[NB+i-3][NB+j], &rho[NB+i-2][NB+j],
					&rho[NB+i-1][NB+j], &rho[NB+i][NB+j],
					&rho[NB+i+1][NB+j], &rho[NB+i+2][NB+j],
					&u[NB+i-3][NB+j], &u[NB+i-2][NB+j],
					&u[NB+i-1][NB+j], &u[NB+i][NB+j],
					&u[NB+i+1][NB+j], &u[NB+i+2][NB+j],
					&v[NB+i-3][NB+j], &v[NB+i-2][NB+j],
					&v[NB+i-1][NB+j], &v[NB+i][NB+j],
					&v[NB+i+1][NB+j], &v[NB+i+2][NB+j],
					&e[NB+i-3][NB+j], &e[NB+i-2][NB+j],
					&e[NB+i-1][NB+j], &e[NB+i][NB+j],
					&e[NB+i+1][NB+j], &e[NB+i+2][NB+j],
					&mut rho_l, &mut u_l, &mut v_l, &mut e_l, 
          &mut rho_r, &mut u_r, &mut v_r, &mut e_r),
          _ => panic!("Reconstruction method not specified."),
        }
        // evaluating flux using flux scheme
        euler::flux_scheme::calc_num_flux(flux_scheme, &rho_l, &u_l, &v_l, &e_l, &rho_r, &u_r, &v_r, &e_r, &ixs_a, &iys_a, &s_a, eos, &mut self.arr_fi[i][j])
      }
    }


    //--------------------j-direction---------------------
    // evaluating numerical flux at (i,j+0.5)
    for i in 0..(NI-2*NB) {
      for j in 0..(NJ-2*NB+1) {
        // inverse of Jacobian
        // evaluated at(i, j+0.5)
        // by averaging the values at(i+0.5, j+0.5) and (i-0.5, j+0.5)
        let s_a = 0.5 * (coord.s[NB+i-1][NB+j-1] + coord.s[NB+i][NB+j-1]);
        // metrices devided by Jacobian
        // evaluated at (i,j+0.5)
        // by averaging the values at (i,j) and (i,j+1)
        //let ixs_a = 0.5 * (coord.ixs[i+1][j] + coord.ixs[i+1][j+1]);
        //let iys_a = 0.5 * (coord.iys[i+1][j] + coord.iys[i+1][j+1]);
        let jxs_a = 0.5 * (coord.jxs[i+1][j] + coord.jxs[i+1][j+1]);
        let jys_a = 0.5 * (coord.jys[i+1][j] + coord.jys[i+1][j+1]);
        // reconstruction
        let mut rho_l = 0.0;
        let mut u_l = 0.0;
        let mut v_l = 0.0;
        let mut e_l = 0.0;
        let mut rho_r = 0.0;
        let mut u_r = 0.0;
        let mut v_r = 0.0;
        let mut e_r = 0.0;
        match reconstruction {
          "MUSCL_minmod_basic" => euler::flux_scheme::reconst_by_basic_muscl(
            &rho[NB+i][NB+j-2], &rho[NB+i][NB+j-1], 
            &rho[NB+i][NB+j], &rho[NB+i][NB+j+1],
            &u[NB+i][NB+j-2], &u[NB+i][NB+j-1], 
            &u[NB+i][NB+j], &u[NB+i][NB+j+1],
            &v[NB+i][NB+j-2], &v[NB+i][NB+j-1], 
            &v[NB+i][NB+j], &v[NB+i][NB+j+1],
            &e[NB+i][NB+j-2], &e[NB+i][NB+j-1], 
            &e[NB+i][NB+j], &e[NB+i][NB+j+1],
            &mut rho_l, &mut u_l, &mut v_l, &mut e_l, 
            &mut rho_r, &mut u_r, &mut v_r, &mut e_r),
          "MP5_basic" => euler::flux_scheme::reconst_by_basic_mp5(
            &rho[NB+i][NB+j-3], &rho[NB+i][NB+j-2],
					&rho[NB+i][NB+j-1], &rho[NB+i][NB+j],
					&rho[NB+i][NB+j+1], &rho[NB+i][NB+j+2],
					&u[NB+i][NB+j-3], &u[NB+i][NB+j-2],
					&u[NB+i][NB+j-1], &u[NB+i][NB+j],
					&u[NB+i][NB+j+1], &u[NB+i][NB+j+2],
					&v[NB+i][NB+j-3], &v[NB+i][NB+j-2],
					&v[NB+i][NB+j-1], &v[NB+i][NB+j],
					&v[NB+i][NB+j+1], &v[NB+i][NB+j+2],
					&e[NB+i][NB+j-3], &e[NB+i][NB+j-2],
					&e[NB+i][NB+j-1], &e[NB+i][NB+j],
					&e[NB+i][NB+j+1], &e[NB+i][NB+j+2],
					&mut rho_l, &mut u_l, &mut v_l, &mut e_l, 
          &mut rho_r, &mut u_r, &mut v_r, &mut e_r),
          _ => panic!("Reconstruction method not specified."),
        }
        // evaluating flux using flux scheme
        euler::flux_scheme::calc_num_flux(flux_scheme, &rho_l, &u_l, &v_l, &e_l, &rho_r, &u_r, &v_r, &e_r, &jxs_a, &jys_a, &s_a, eos, &mut self.arr_fj[i][j]);
      }
    }


    //---------------------calculation og RHS---------------------
    for i in 0..(NI-2*NB) {
      for j in 0..(NJ-2*NB) {
        for k in 0..NF {
          arr_q[i][j][k] =
          self.arr_fi[i][j][k] - self.arr_fi[i+1][j][k] +
					self.arr_fj[i][j][k] - self.arr_fj[i][j+1][k];
        }
      }
    }
  }


}

