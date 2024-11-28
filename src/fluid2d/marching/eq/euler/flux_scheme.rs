const NF: usize = 4;
use super::super::super::super::eos::EoS;

use super::super::fnd;



/// calc numerical flux with scheme selected as the first arg
#[inline]
pub fn calc_num_flux<T: EoS>(flux_scheme: &str, rho_l: &f64, u_l: &f64, v_l: &f64, e_l: &f64,
  rho_r: &f64, u_r: &f64, v_r: &f64, e_r: &f64, 
  ixs: &f64, iys: &f64, s: &f64, eos: &T, vec_fc: &mut [f64; NF]) {
    match flux_scheme {
      "Roe_FDS" => roe_fds(rho_l, u_l, v_l, e_l, rho_r, u_r, v_r, e_r, ixs, iys, s, eos, vec_fc),
      _ => panic!("Flux scheme not specified."),
    }
  }



/// reconstruct the cell-boundary values
/// with MUSCL-minmod adapted for basic variables
#[inline]
pub fn reconst_by_basic_muscl(
rho0: &f64, rho1: &f64, rho2: &f64, rho3: &f64,
u0: &f64, u1: &f64, u2: &f64, u3: &f64,
v0: &f64, v1: &f64, v2: &f64, v3: &f64,
e0: &f64, e1: &f64, e2: &f64, e3: &f64,
rho_l: &mut f64, u_l: &mut f64, v_l: &mut f64, e_l: &mut f64,
rho_r: &mut f64, u_r: &mut f64, v_r: &mut f64, e_r: &mut f64) {
  fnd::muscl_minmod(rho0, rho1, rho2, rho3, rho_l, rho_r);
  fnd::muscl_minmod(u0, u1, u2, u3, u_l, u_r);
  fnd::muscl_minmod(v0, v1, v2, v3, v_l, v_r);
  fnd::muscl_minmod(e0, e1, e2, e3, e_l, e_r);
}



/// reconstruct the cell-boundary values
/// with MP5 adapted for basic variables
#[inline]
pub fn reconst_by_basic_mp5(
rho0: &f64, rho1: &f64, rho2: &f64, rho3: &f64, rho4: &f64, rho5: &f64,
u0: &f64, u1: &f64, u2: &f64, u3: &f64, u4: &f64, u5: &f64,
v0: &f64, v1: &f64, v2: &f64, v3: &f64, v4: &f64, v5: &f64,
e0: &f64, e1: &f64, e2: &f64, e3: &f64, e4: &f64, e5: &f64,
rho_l: &mut f64, u_l: &mut f64, v_l: &mut f64, e_l: &mut f64,
rho_r: &mut f64, u_r: &mut f64, v_r: &mut f64, e_r: &mut f64) {
  fnd::mp5(rho0, rho1, rho2, rho3, rho4, rho5, rho_l, rho_r);
  fnd::mp5(u0, u1, u2, u3, u4, u5, u_l, u_r);
  fnd::mp5(v0, v1, v2, v3, v4, v5, v_l, v_r);
  fnd::mp5(e0, e1, e2, e3, e4, e5, e_l, e_r);
}

/// calc Roe average
#[inline]
fn roe_average<T: EoS>(rho_l: &f64, u_l: &f64, v_l: &f64, e_l: &f64,
  rho_r: &f64, u_r: &f64, v_r: &f64, e_r: &f64, eos: &T,
  rho_a: &mut f64, u_a: &mut f64, v_a: &mut f64, e_a: &mut f64) {
  let p_l = eos.calc_p(rho_l, u_l, v_l, e_l);
  let p_r = eos.calc_p(rho_r, u_r, v_r, e_r);
  let h_l = (e_l + p_l) / rho_l;
  let h_r = (e_r + p_r) / rho_r;
  let sqr_l = rho_l.sqrt();
  let sqr_r = rho_r.sqrt();

  *rho_a = sqr_l * sqr_r;
  *u_a = (sqr_l * u_l + sqr_r * u_r) / (sqr_l + sqr_r);
  *v_a = (sqr_l * v_l + sqr_r * v_r) / (sqr_l + sqr_r);
  let h_a = (sqr_l * h_l + sqr_r * h_r) / (sqr_l + sqr_r);
  *e_a = eos.calc_e(rho_a, u_a, v_a, &h_a);
}
/// calc convective numerical flux using classical Roe-type FDS scheme
fn roe_fds<T: EoS>(rho_l: &f64, u_l: &f64, v_l: &f64, e_l: &f64,
  rho_r: &f64, u_r: &f64, v_r: &f64, e_r: &f64,
  ixs: &f64, iys: &f64, s: &f64, eos: &T, vec_fc: &mut [f64; NF]) {
  let ix = ixs / s;
  let iy = iys / s;
  let mut rho_a = 0.0;
  let mut u_a = 0.0;
  let mut v_a = 0.0;
  let mut e_a = 0.0;
  roe_average(rho_l, u_l, v_l, e_l, rho_r, u_r, v_r, e_r, eos,
    &mut rho_a, &mut u_a, &mut v_a, &mut e_a);
  let mut vec_q_l = [0.0; NF];
  super::conserved::calc_conservative(rho_l, u_l, v_l, e_l, s, &mut vec_q_l);
  let mut vec_q_r = [0.0; NF];
 super::conserved::calc_conservative(rho_r, u_r, v_r, e_r, s, &mut vec_q_r);
  let mut vec_f_l = [0.0; NF];
  super::conserved::calc_flux_conv::<T>(rho_l, u_l, v_l, e_l, ixs, iys, eos, &mut vec_f_l);
  let mut vec_f_r = [0.0; NF];
  super::conserved::calc_flux_conv::<T>(rho_r, u_r, v_r, e_r, ixs, iys, eos, &mut vec_f_r);
  let mut dia_lam = [0.0; NF];
  let mut mat_r = [[0.0; NF]; NF];
  let mut mat_rinv = [[0.0; NF]; NF];
  eos.calc_eigen(&rho_a, &u_a, &v_a, &e_a, &ix, &iy, &mut dia_lam, &mut mat_r, &mut mat_rinv);

  for i in 0..NF {
    vec_fc[i] = vec_f_r[i] + vec_f_l[i];
    for j in 0..NF {
      for k in 0.. NF {
        vec_fc[i] -= mat_r[i][j] * (dia_lam[j].abs()) * mat_rinv[j][k] * (vec_q_r[k] - vec_q_l[k]);
      }
    }
    vec_fc[i] *= 0.5;
  }
}