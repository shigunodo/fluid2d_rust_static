use crate::fluid2d::marching::Marching;

mod fluid2d;

use fluid2d::{basic_var::BasicVar, settings::Settings};

fn main() {
    const NI: usize = 408;
    const NJ: usize = 408;
    const NB: usize = 4;
    // !!!!!!!!!!ワーキングディレクトリを設定してください
    const DIR: &str = "data/";
    let dir_o: &str = &(DIR.to_string() + "");
    let f_coordinate: &str = &(DIR.to_string() + "coordinate.dat");
    let f_settings: &str = &(DIR.to_string() + "settings.dat");
    const GAMMA: f64 = 1.4;

    let start = std::time::Instant::now();

    unsafe {
        static mut FLUID: fluid2d::IdealGas<NI, NJ, NB> = fluid2d::IdealGas::<NI, NJ, NB>::new(&GAMMA);
        FLUID.initialize(dir_o, f_coordinate);

        let t_max = 3.0;
        let n_out = 100;
        let dt_out = t_max / n_out as f64;
        let cfl_coeff = 0.7;
        let mut t = 0.0;
        FLUID
            .settings
            .output(f_settings, &t_max, &n_out, &FLUID.eos);

        for tstep in 1..=n_out {
            let mut iter = 0;

            while t < dt_out * tstep as f64 {
                let dt = FLUID.calc_cfl(&cfl_coeff);

                FLUID.march_ssprk3(&dt, "periodical_in_i", "MP5_basic", "Roe_FDS");

                t += dt;
                iter += 1;
            }

            let cpu_time = start.elapsed().as_secs();
            let rest_time = cpu_time * (n_out - tstep) / tstep;
            FLUID.basic.output(
                dir_o,
                f_settings,
                &t,
                &(tstep as usize),
                &iter,
                &cpu_time,
                &rest_time,
            );
        }
    } // unsafe

    println!("Program ended.");
}
