const NF: usize = 4;

pub trait EoS {
  /// calc pressure from density, velocity, total energy per volume
  fn calc_p(&self, rho: &f64, u: &f64, v: &f64, e: &f64) -> f64;
  /// calc temperature from density, velocity, total energy per volume
  fn calc_temp(&self, rho: &f64, u: &f64, v: &f64, e: &f64) -> f64;
  /// calc sound-speed from density, velocity, total energy per volume
  fn calc_cs(&self, rho: &f64, u: &f64, v: &f64, e: &f64) -> f64;
  /// calc total energy per volume
  /// from density, velocity, total specific enthalpy
  fn calc_e(&self, rho: &f64, u: &f64, v: &f64, h: &f64) -> f64;
  /// calc density and total energy per volume
  /// from pressure, temperature, velocity
  fn calc_rho_e(&self, p: &f64, temp: &f64, u: &f64, v: &f64,
    rho: &mut f64, e: &mut f64);
  /// calc total energy per volume from density, velocity, pressure
  fn calc_e_wp(&self, rho: &f64, u: &f64, v: &f64, p: &f64) -> f64;
  /// eigen values/vectors of flux Jacobian
  fn calc_eigen(&self, rho: &f64, u: &f64, v: &f64, e: &f64,
    ix: &f64, iy: &f64, dia_lam: &mut [f64; NF], 
    mat_r: &mut [[f64; NF]; NF], mat_rinv: &mut [[f64; NF]; NF]);
}

impl EoS for super::IdealEoS {
  /// calc pressure from density, velocity, total energy per volume
  #[inline]
  fn calc_p(&self, rho: &f64, u: &f64, v: &f64, e: &f64) -> f64 {
    (self.gamma - 1.0) * (e - 0.5 * rho * (u * u + v * v))
  }
  /// calc temperature from density, velocity, total energy per volume
  #[inline]
  fn calc_temp(&self, rho: &f64, u: &f64, v: &f64, e: &f64) -> f64 {
    e / rho - 0.5 * (u * u + v * v)
  }
  /// calc sound-speed from density, velocity, total energy per volume
  #[inline]
  fn calc_cs(&self, rho: &f64, u: &f64, v: &f64, e: &f64) -> f64 {
    (self.gamma * (self.gamma - 1.0) * (e / rho - 0.5 * (u * u + v * v))).sqrt()
  }
  /// calc total energy per volume
  /// from density, velocity, total specific enthalpy
  #[inline]
  fn calc_e(&self, rho: &f64, u: &f64, v: &f64, h: &f64) -> f64 {
    rho * (h + (self.gamma - 1.0) * 0.5 * (u * u + v * v)) / self.gamma
  }
  /// calc density and total energy per volume
  /// from pressure, temperature, velocity
  #[inline]
  fn calc_rho_e(&self, p: &f64, temp: &f64, u: &f64, v: &f64,
    rho: &mut f64, e: &mut f64) {
    *rho = p / temp / (self.gamma - 1.0);
    *e = *rho * (temp + 0.5 * (u * u + v * v));
  }
  /// calc total energy per volume from density, velocity, pressure
  #[inline]
  fn calc_e_wp(&self, rho: &f64, u: &f64, v: &f64, p: &f64) -> f64 {
    p / (self.gamma - 1.0) + 0.5 * rho * (u * u + v * v)
  }
  /// eigen values/vectors of flux Jacobian
  #[inline]
  fn calc_eigen(&self, rho: &f64, u: &f64, v: &f64, e: &f64,
    ix: &f64, iy: &f64, dia_lam: &mut [f64; NF], 
    mat_r: &mut [[f64; NF]; NF], mat_rinv: &mut [[f64; NF]; NF]) {
      // preperation
      let cs = self.calc_cs(rho, u, v, e);
      let p = self.calc_p(rho, u, v, e);
      let h = (e + p) / rho;
      let sqr = (ix * ix + iy * iy).sqrt();
      let ixb = ix / sqr;
      let iyb = iy / sqr;
      let bigu = ix * u + iy * v;
      let bigub = bigu / sqr;
      let b1 = 0.5 * (u * u + v * v) * (self.gamma - 1.0) / cs / cs;
      let b2 = (self.gamma - 1.0) / cs / cs;
      // diagonal matrix of eigen values
      dia_lam[0] = bigu - cs * sqr;
      dia_lam[1] = bigu;
      dia_lam[2] = bigu + cs * sqr;
      dia_lam[3] = bigu;
      // right eigen matrix
      mat_r[0][0] = 1.0;
      mat_r[0][1] = 1.0;
      mat_r[0][2] = 1.0;
      mat_r[0][3] = 0.0;
      mat_r[1][0] = u - ixb * cs;
      mat_r[1][1] = *u;
      mat_r[1][2] = u + ixb * cs;
      mat_r[1][3] = -iyb;
      mat_r[2][0] = v - iyb * cs;
      mat_r[2][1] = *v;
      mat_r[2][2] = v + iyb * cs;
      mat_r[2][3] = ixb;
      mat_r[3][0] = h - cs * bigub;
      mat_r[3][1] = 0.5 * (u * u + v * v);
      mat_r[3][2] = h + cs * bigub;
      mat_r[3][3] = -(iyb * u - ixb * v);
      // inverse of righr eigen matrix
      mat_rinv[0][0] = 0.5 * (b1 + bigub / cs);
      mat_rinv[0][1] = -0.5 * (ixb / cs + b2 * u);
      mat_rinv[0][2] = -0.5 * (iyb / cs + b2 * v);
      mat_rinv[0][3] = 0.5 * b2;
      mat_rinv[1][0] = 1.0 - b1;
      mat_rinv[1][1] = b2 * u;
      mat_rinv[1][2] = b2 * v;
      mat_rinv[1][3] = -b2;
      mat_rinv[2][0] = 0.5 * (b1 - bigub / cs);
      mat_rinv[2][1] = 0.5 * (ixb / cs - b2 * u);
      mat_rinv[2][2] = 0.5 * (iyb / cs - b2 * v);
      mat_rinv[2][3] = 0.5 * b2;
      mat_rinv[3][0] = iyb * u - ixb * v;
      mat_rinv[3][1] = -iyb;
      mat_rinv[3][2] = ixb;
      mat_rinv[3][3] = 0.0;
    }
}