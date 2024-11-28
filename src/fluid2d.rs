use self::{basic_var::BasicVar, coordinate::Coord};

mod coordinate;
pub mod basic_var;
pub mod marching;
mod eos;
pub mod settings;

const NF: usize = 4;

/// struct representing non-viscous ideal gas
/// should have trait Marching
pub struct IdealGas<const NI: usize, const NJ: usize, const NB: usize> {
  /// directory for data-output
  //dir_o: String,

  /// file name for coordinate
  //f_coordinate: String,

  /// file name for setting file
  //f_settings: String,

  /// settings
  pub settings: SetStructEuler::<NI,NJ,NB>,

  /// general coordinate
  coord: GenStructCoord::<NI,NJ,NB>,

  /// basic variables
  pub basic: BasicVarHD::<NI,NJ,NB>,

  /// Euler equation
  eq: EulerEq::<NI,NJ,NB>,

  /// equation of state
  pub eos: IdealEoS,

  /// boundary condition
  bound: Boundary::<NI,NJ,NB>,

  // arrays needed for procedures
  /// for march_ssprk3
  /// only the partial array [0:NI-2*NB][0:NJ-2*NB][NF] is used actually
  arr_q0: [[[f64; NF]; NJ]; NI],
  /// for march_ssprk3
  /// only the partial array [0:NI-2*NB][0:NJ-2*NB][NF] is used actually
  arr_q1: [[[f64; NF]; NJ]; NI],
  /// for march_ssprk3
  /// only the partial array [0:NI-2*NB][0:NJ-2*NB][NF] is used actually
  arr_q2: [[[f64; NF]; NJ]; NI],
}

impl<const NI: usize, const NJ: usize, const NB: usize> IdealGas::<NI,NJ,NB> {
  /// Constructor of IdealGas
  pub const fn new(gamma: &f64) -> Self {
    IdealGas::<NI,NJ,NB> {
      settings: SetStructEuler::<NI,NJ,NB>::new(),
      coord: GenStructCoord::<NI,NJ,NB>::new(),
      basic: BasicVarHD::<NI,NJ,NB>::new(),
      eq: EulerEq::<NI,NJ,NB>::new(),
      eos: IdealEoS::new(gamma),
      bound: Boundary::<NI,NJ,NB>::new(),
      arr_q0: [[[0.0; NF]; NJ]; NI],
      arr_q1: [[[0.0; NF]; NJ]; NI],
      arr_q2: [[[0.0; NF]; NJ]; NI],
    }
  }

  pub fn initialize(&mut self, dir_o: &str, f_coordinate: &str) {
    self.coord.input(&f_coordinate);
    self.coord.calc_metrices_dx();
    let f_initial = dir_o.to_string() + "b0000000.dat";
    self.basic.input(&f_initial);
  }
}





/// struct representing settings of calculation
/// should have trait Settings
pub struct SetStructEuler<const NI: usize, const NJ: usize, const NB: usize> {}

impl<const NI: usize, const NJ: usize, const NB: usize> SetStructEuler::<NI,NJ,NB> {
  /// Constructor of SetStructEuler
  const fn new() -> Self {
    SetStructEuler::<NI,NJ,NB> {}
  }
}






/// struct representing basic variables
/// should have trait BasicVar
pub struct BasicVarHD<const NI: usize, const NJ: usize, const NB: usize> {
  // basic variables
  /// density
  rho: [[f64; NJ]; NI],
  /// x-component of velocity
  u: [[f64; NJ]; NI],
  /// y-component of velocity
  v: [[f64; NJ]; NI],
  /// total energy per volume
  e: [[f64; NJ]; NI],
}

impl<const NI: usize, const NJ: usize, const NB: usize> BasicVarHD::<NI,NJ,NB> {
  /// constructor of BasicVarHD
  const fn new() -> Self {
    BasicVarHD::<NI,NJ,NB> {
      rho: [[0.0; NJ]; NI],
      u: [[0.0; NJ]; NI],
      v: [[0.0; NJ]; NI],
      e: [[0.0; NJ]; NI],
    }
  }
}





/// struct representing general coordinate
/// should have trait Coord
pub struct GenStructCoord<const NI: usize, const NJ: usize, const NB: usize> {
  // coordinate
  /// x-list of grids
  x: [[f64; NJ]; NI],
  /// y-list of grids
  y: [[f64; NJ]; NI],

  // quantities derived from coordinate at first
  /// metrix ix devided by Jacobian
  /// only the partial array [0:NI-2*NB+2][0:NJ-2*NB+2] is used actually
  ixs: [[f64; NJ]; NI],
  /// metrix iy devided by Jacobian
  /// only the partial array [0:NI-2*NB+2][0:NJ-2*NB+2] is used actually
  iys: [[f64; NJ]; NI],
  /// metrix jx devided by Jacobian
  /// only the partial array [0:NI-2*NB+2][0:NJ-2*NB+2] is used actually
  jxs: [[f64; NJ]; NI],
  /// metrix jy devided by Jacobian
  /// only the partial array [0:NI-2*NB+2][0:NJ-2*NB+2] is used actually
  jys: [[f64; NJ]; NI],
  /// inverse of Jacobian
  /// only the partial array [0:NI-1][0:NJ-1] is used actually
  s: [[f64; NJ]; NI],
  /// dx for CFL condition
  /// only the partial array [0:NI-2*NB][0:NJ-2*NB] is used actually
  dx: [[f64; NJ]; NI],
}

impl<const NI: usize, const NJ: usize, const NB: usize>  GenStructCoord::<NI,NJ,NB> {
  /// constructor of GenStructCoord
  const fn new() -> Self {
    GenStructCoord::<NI,NJ,NB> {
      x: [[0.0; NJ]; NI],
      y: [[0.0; NJ]; NI],
      ixs: [[0.0; NJ]; NI],
      iys: [[0.0; NJ]; NI],
      jxs: [[0.0; NJ]; NI],
      jys: [[0.0; NJ]; NI],
      s: [[0.0; NJ]; NI],
      dx: [[0.0; NJ]; NI],
    }
  }
}





/// struct representing rhs of Euler equation
/// should have trait RHS
struct EulerEq<const NI: usize, const NJ: usize, const NB: usize> {
  /// for rhs
  /// only the partial array [0:NI-2*NB+1][0:NJ-2*NB][NF] is used actually
  arr_fi: [[[f64; NF]; NJ]; NI],
  /// for rhs
  /// only the partial array [0:NI-2*NB][0:NJ-2*NB+1][NF] is used actually
  arr_fj: [[[f64; NF]; NJ]; NI],
}

impl<const NI: usize, const NJ: usize, const NB: usize> EulerEq::<NI,NJ,NB> {
  /// constructor of EulerEq
  const fn new() -> Self {
    EulerEq::<NI,NJ,NB> {
      arr_fi: [[[0.0; NF]; NJ]; NI],
      arr_fj: [[[0.0; NF]; NJ]; NI],
    }
  }
}





/// struct representing ideal equation of state
/// should have trait EoS
pub struct IdealEoS {
  /// specific heat
  gamma: f64,
}

impl IdealEoS {
  /// constructor of IdealEoS
  const fn new(gamma: &f64) -> Self {
    IdealEoS {
      gamma: *gamma,
    }
  }
}





/// struct representing boundary condition
/// should have trait BC
struct Boundary<const NI: usize, const NJ: usize, const NB: usize> {}

impl<const NI: usize, const NJ: usize, const NB: usize> Boundary::<NI,NJ,NB> {
  /// constructor of Boundary
  const fn new() -> Self {
    Boundary::<NI,NJ,NB> {}
  }
}