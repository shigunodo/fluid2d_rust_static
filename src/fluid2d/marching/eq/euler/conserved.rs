const NF: usize = 4;
use super::super::super::super::eos::EoS;



/// construct conservative var from basic var
#[inline]
pub fn calc_conservative(rho: &f64, u: &f64, v: &f64, e: &f64, 
  s: &f64, vec_q: &mut [f64; NF]) {
    vec_q[0] = rho * s;
    vec_q[1] = rho * u * s;
    vec_q[2] = rho * v * s;
    vec_q[3] = e * s;
  }



/// construct convection flux from basic var
#[inline]
pub fn calc_flux_conv<T: EoS>(rho: &f64, u: &f64, v: &f64, e:&f64,
  ixs: &f64, iys: &f64, eos: &T, vec_fc: &mut[f64; NF]) {
    let bigu = ixs * u + iys * v;
    let p = eos.calc_p(rho, u, v, e);
    vec_fc[0] = rho * bigu;
    vec_fc[1] = rho * u * bigu + ixs * p;
    vec_fc[2] = rho * v * bigu + iys * p;
    vec_fc[3] = (e + p) * bigu;
  }



  
/// construct basic var from conservative var
#[inline]
pub fn calc_basic(vec_q: &[f64; NF], s: &f64, 
  rho: &mut f64, u: &mut f64, v: &mut f64, e: &mut f64) {
    *rho = vec_q[0] / s;
    *u = vec_q[1] / *rho / s;
    *v = vec_q[2] / *rho / s;
    *e = vec_q[3] / s;
  }