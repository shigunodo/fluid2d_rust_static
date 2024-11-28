use std::fs::File;
use std::io::Write;
use super::IdealEoS;


pub trait Settings {
  fn output(&self, f_settings: &str, t_max: &f64, n_out: &u64, eos: &IdealEoS);
}


impl<const NI: usize, const NJ: usize, const NB: usize> Settings 
for super::SetStructEuler::<NI,NJ,NB> {
  fn output(&self, f_settings: &str, t_max: &f64, n_out: &u64, eos: &IdealEoS) {
    let mut file = File::create(f_settings).unwrap();
    let line = format!("Solving Euler eq. with ideal gas eos.
    gamma = {}
    NI = {}
    NJ = {}
    NB = {}
    Tmax = {}
    Nout = {}
    
    ", eos.gamma, NI, NJ, NB, t_max, n_out);
    write!(file, "{}", line).unwrap();
  }
}