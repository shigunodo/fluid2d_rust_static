/// minmod for MUSCL and MP5 procedures
#[inline]
fn minmod(x: &f64, y: &f64) -> f64 {
  x.signum() * (x.abs().min(x.signum() * y).max(0.0))
}

/// minmod of 4 numbers, for MP5 procedure
#[inline]
fn minmod4(x1: &f64, x2: &f64, x3: &f64, x4: &f64) -> f64 {
  0.5 * (x1.signum() + x2.signum())
  * (0.5 * (x1.signum() + x3.signum())
  * 0.5 * (x1.signum() + x4.signum())).abs()
  * x1.abs().min(x2.abs()).min(x3.abs()).min(x4.abs())
}

/// median for MP5 procedure
#[inline]
fn median(x: &f64, y:&f64, z: &f64) -> f64 {
  x + minmod(&(y - x), &(z - x))
}

/// MUSCL-minmod reconstruction
pub fn muscl_minmod(qm: &f64, q: &f64, qp: &f64, q2p: &f64,
  q_l: &mut f64, q_r: &mut f64) {
  let k = 1.0 / 3.0;
  let b = (3.0 - k) / (1.0 - k);
  // calculation of L
  let dp = qp - q;
  let dm = q - qm;
  let ddp = minmod(&dp, &(b * dm));
  let ddm = minmod(&dm, &(b * dp));
  *q_l = q + 0.25 * ((1.0 - k) * ddm + (1.0 + k) * ddp);
  // calculation of R
  let dp = q2p - qp;
  let dm = qp - q;
  let ddp = minmod(&dp, &(b * dm));
  let ddm = minmod(&dm, &(b * dp));
  *q_r = qp - 0.25 * ((1.0 - k) * ddp + (1.0 + k) * ddm);
}

/// MP5 reconstruction
pub fn mp5_sub(q2m: &f64, qm: &f64, q: &f64, qp: &f64, q2p: &f64,
  q_l: &mut f64) {
  let alp = 2.0;
  *q_l = (2.0 * q2m - 13.0 * qm + 47.0 * q + 27.0 * qp - 3.0 * q2p) / 60.0;
  let q_mp = q + minmod(&(qp - q), &(alp * (q - qm)));
  if (*q_l - q) * (*q_l - q_mp) > 1.0e-10 {
    let dm = q2m + q - 2.0 * qm;
    let d = qm + qp - 2.0 * q;
    let dp = q + q2p - 2.0 * qp;
    let d_mm = minmod4(&(4.0 * dm - d), &(4.0 * d - dm), &dm, &d);
    let d_mp = minmod4(&(4.0 * d - dp), &(4.0 * dp - d), &d, &dp);
    let q_ul = q + alp * (q - qm);
    let q_md = 0.5 * (q + qp) - 0.5 * d_mp;
    let q_lc = q + 0.5 * (q - qm) + 4.0 / 3.0 * d_mm;
    let qmin = q.min(*qp).min(q_md).max(q.min(q_ul).min(q_lc));
    let qmax = q.max(*qp).max(q_md).min(q.max(q_ul).max(q_lc));
    *q_l = median(&q_l, &qmin, &qmax);
  }
}

pub fn mp5(q2m: &f64, qm: &f64, q: &f64, qp: &f64, q2p: &f64, q3p: &f64,
  q_l: &mut f64, q_r: &mut f64) {
    mp5_sub(q2m, qm, q, qp, q2p, q_l);
    mp5_sub(q3p, q2p, qp, q, qm, q_r);
  }

/// 3rd order central difference
/// returning first dericative at i=1.5 from values at i=0,1,2,3
pub fn central_diff3(q0: &f64, q1: &f64, q2: &f64, q3: &f64) -> f64 {
  q0 / 24.0 - 9.0 / 8.0 * q1 + 9.0 / 8.0 * q2 - q3 / 24.0
}

/// 4th order central difference
/// returning first derivative at i=2 from values at i=0,1,3,4
pub fn central_diff4(q0: &f64, q1: &f64, q3: &f64, q4: &f64) -> f64 {
  q0 / 12.0 - 2.0 / 3.0 * q1 + 2.0 / 3.0 * q3 - q4 / 12.0
} 

#[test]
fn test() {
  println!("{}", minmod4(&1.0,&2.0,&3.0,&2.0));
}