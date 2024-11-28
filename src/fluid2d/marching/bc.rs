use super::super::{eos::EoS, BasicVarHD};

pub trait BCHD<const NI: usize, const NJ: usize, const NB: usize> {
  fn reflect_bc<T: EoS>(&self, bc_type: &str, eos: &T, basic: &mut BasicVarHD::<NI,NJ,NB>);
  fn bc_periodical_in_i(&self, basic: &mut BasicVarHD::<NI,NJ,NB>);
}

impl<const NI: usize, const NJ:usize, const NB: usize> BCHD::<NI,NJ,NB>
for super::super::Boundary::<NI,NJ,NB> {
  fn reflect_bc<T: EoS>(&self, bc_type: &str, eos: &T, basic: &mut BasicVarHD::<NI,NJ,NB>) {
    match bc_type {
      "periodical_in_i" => self.bc_periodical_in_i(basic),
      _ => panic!("BC not specified."),
    }
  }

  fn bc_periodical_in_i(&self, basic: &mut BasicVarHD::<NI,NJ,NB>)
  {
    // preiodical in i-direction
    // for i=0
    for i in 0..NB {
      for j in 0..NJ {
        basic.rho[i][j] = basic.rho[NI-2*NB+i][j];
        basic.u[i][j] = basic.u[NI-2*NB+i][j];
        basic.v[i][j] = basic.v[NI-2*NB+i][j];
        basic.e[i][j] = basic.e[NI-2*NB+i][j];
      }
    }
    // for i=Ni
    for i in 0..NB {
      for j in 0..NJ {
        basic.rho[NI-NB+i][j] = basic.rho[NB+i][j];
        basic.u[NI-NB+i][j] = basic.u[NB+i][j];
        basic.v[NI-NB+i][j] = basic.v[NB+i][j];
        basic.e[NI-NB+i][j] = basic.e[NB+i][j];
      }
    }
    // No updates for j-boundaries since Dirichlet
  } 
}