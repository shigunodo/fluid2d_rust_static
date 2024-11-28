use std::fs::File;
use std::io::{BufReader, BufRead};


pub trait Coord {
  fn calc_metrices_dx(&mut self);
  fn input(&mut self, f_name: &str);
}


impl<const NI: usize, const NJ: usize, const NB: usize> Coord
for super::GenStructCoord::<NI,NJ,NB> {
  fn calc_metrices_dx(&mut self) {
    // calculate metrices devided by Jacobian
    // evaluated at cell-points (i,j)
    // returning (Ni-2*Nb+2,Nj-2*Nb+2) arrays
    // using 2nd order central difference
    for i in 0..(NI-2*NB+2) {
      for j in 0..(NJ-2*NB+2) {
        self.ixs[i][j] = 0.5 * (self.y[NB-1+i][NB+j] - self.y[NB-1+i][NB-2+j]);
        self.iys[i][j] = -0.5 * (self.x[NB-1+i][NB+j] - self.x[NB-1+i][NB-2+j]);
        self.jxs[i][j] = -0.5 * (self.y[NB+i][NB-1+j] - self.y[NB-2+i][NB-1+j]);
        self.jys[i][j] = 0.5 * (self.x[NB+i][NB-1+j] - self.x[NB-2+i][NB-1+j]);
      }
    }
    // calculate inverse of Jacobian
    // evaluated at cell-verteces (i+0.5,j+0.5)
    // returning (Ni-1,Nj-1) array
    for i in 0..(NI-1) {
      for j in 0..(NJ-1) {
        self.s[i][j] = ((self.x[1+i][1+j] - self.x[i][j]) * (self.y[i][1+j] - self.y[1+i][j]) - (self.y[1+i][1+j] - self.y[i][j]) * (self.x[i][1+j] - self.x[1+i][j]))  * 0.5;
      }
    }
    // calculate dx for CFL condition
    // evaluated at cell-points (i,j)
    // returning (Ni-2*Nb,Nj-2*Nb) array
    // first, calculate the adjacent 4 values
    // then, minimize the 4 values
    for i in 0..(NI-2*NB) {
      for j in 0..(NJ-2*NB) {
        let dx1 = self.x[NB+1+i][NB+j] - self.x[NB+i][NB+j];
        let dy1 = self.y[NB+1+i][NB+j] - self.y[NB+i][NB+j];
        self.dx[i][j] = dx1 * dx1 + dy1 * dy1;
        let dx1 = self.x[NB+i][NB+1+j] - self.x[NB+i][NB+j];
        let dy1 = self.y[NB+i][NB+1+j] - self.y[NB+i][NB+j];
        let dd = dx1 * dx1 + dy1 * dy1;
        self.dx[i][j] = self.dx[i][j].min(dd);
        let dx1 = self.x[NB-1+i][NB+j] - self.x[NB+i][NB+j];
        let dy1 = self.y[NB-1+i][NB+j] - self.y[NB+i][NB+j];
        let dd = dx1 * dx1 + dy1 * dy1;
        self.dx[i][j] = self.dx[i][j].min(dd);
        let dx1 = self.x[NB+i][NB-1+j] - self.x[NB+i][NB+j];
        let dy1 = self.y[NB+i][NB-1+j] - self.y[NB+i][NB+j];
        let dd = dx1 * dx1 + dy1 * dy1;
        self.dx[i][j] = self.dx[i][j].min(dd);
        self.dx[i][j] = self.dx[i][j].sqrt();
      }
    }
  }


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
            0 => self.x[i][j] = l.trim().parse::<f64>().unwrap(),
            1 => self.y[i][j] = l.trim().parse::<f64>().unwrap(),
            _ => panic!("Data size of coordinate may be wrong."),
        }
        n_line += 1;
    }
  }


}