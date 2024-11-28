use std::fs::{File, OpenOptions};
use std::io::{Write, BufReader, BufWriter, BufRead, stdout};

pub trait BasicVar {
  fn input(&mut self, f_name: &str);
  fn output(&self, dir_o: &str, f_settings: &str, t: &f64, tstep: &usize, iter: &usize, cpu_time: &u64, rest_time: &u64);
}


impl<const NI: usize, const NJ: usize, const NB: usize> BasicVar
for super::BasicVarHD::<NI,NJ,NB> {


  fn input(&mut self, f_name: &str) {
    let f = File::open(f_name).unwrap();
    let buf = BufReader::new(f);
    let mut n_line: usize = 0;
    for line in buf.lines() {
        let l: &str = &line.unwrap();
        let c = n_line / NI / NJ;
        let n_line_sub = n_line - NI * NJ * c;
        let i = n_line_sub / NI;
        let j = n_line_sub - NI * i;
        match c {
            0 => self.rho[i][j] = l.trim().parse::<f64>().unwrap(),
            1 => self.u[i][j] = l.trim().parse::<f64>().unwrap(),
            2 => self.v[i][j] = l.trim().parse::<f64>().unwrap(),
            3 => self.e[i][j] = l.trim().parse::<f64>().unwrap(),
            _ => panic!("Data size of basic var may be wrong."),
        }
        n_line += 1;
    }
  }



  fn output(&self, dir_o: &str, f_settings: &str, t: &f64, tstep: &usize, iter: &usize, cpu_time: &u64, rest_time: &u64) {
    // output into file
    let mut fo_name = format!("b{:07}.dat", tstep);
    fo_name = dir_o.to_string() + &fo_name;
    let mut file = BufWriter::new(File::create(fo_name).unwrap());
    for i in 0..NI {
        for j in 0..NJ {
            write!(file, "{:.18e}\n", self.rho[i][j]).unwrap();
        }
    }
    for i in 0..NI {
        for j in 0..NJ {
            write!(file, "{:.18e}\n", self.u[i][j]).unwrap();
        }
    }
    for i in 0..NI {
      for j in 0..NJ {
          write!(file, "{:.18e}\n", self.v[i][j]).unwrap();
      }
    }
    for i in 0..NI {
      for j in 0..NJ {
          write!(file, "{:.18e}\n", self.e[i][j]).unwrap();
      }
    }

    // convert cpu_time and rest_time into h/m/s
    const SECS: u64 = 60;
    let h_c = cpu_time / SECS / SECS;
    let mut s_c = cpu_time - h_c * SECS * SECS;
    let m_c = s_c / SECS;
    s_c -= m_c * SECS;

    let h_r = rest_time / SECS / SECS;
    let mut s_r = rest_time - h_r * SECS * SECS;
    let m_r = s_r / SECS;
    s_r -= m_r * SECS;

    // output status
    let line = format!("elapsed: {h_c:3} h {m_c:02} m {s_c:02} s | tstep = {tstep:5} | t = {t:10.4} | iter = {iter:7} | rest: {h_r:3} h {m_r:02} m {s_r:02} s");

    let mut file = BufWriter::new(OpenOptions::new().append(true).open(f_settings).unwrap());
    writeln!(file, "{}", line).unwrap();

    let out = stdout();
    let mut out = BufWriter::new(out.lock());
    writeln!(out, "{}", line).unwrap();
  }



}

