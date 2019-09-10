#![warn(clippy::all)]

mod app_timers;
mod checksum;
mod fortran_array_2d;
mod kernels;
mod model_parameters;
mod timer;
mod simulation_variables;

use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use app_timers::AppTimers;
use fortran_array_2d::FortranArray2D;
use kernels::{boundary_conditions_kernel, continuity_kernel, momentum_kernel, next_kernel};
use model_parameters::ModelParameters;
use simulation_variables::SimulationVariables;

type WorkingPrecision = f64;

pub struct GridConstants {
    pub e1t: FortranArray2D<WorkingPrecision>,
    pub e2t: FortranArray2D<WorkingPrecision>,

    pub e1u: FortranArray2D<WorkingPrecision>,
    pub e2u: FortranArray2D<WorkingPrecision>,

    pub e1v: FortranArray2D<WorkingPrecision>,
    pub e2v: FortranArray2D<WorkingPrecision>,

    pub e1f: FortranArray2D<WorkingPrecision>,
    pub e2f: FortranArray2D<WorkingPrecision>,

    pub e12t: FortranArray2D<WorkingPrecision>,
    pub e12u: FortranArray2D<WorkingPrecision>,
    pub e12v: FortranArray2D<WorkingPrecision>,

    pub gphiu: FortranArray2D<WorkingPrecision>,
    pub gphiv: FortranArray2D<WorkingPrecision>,
    pub gphif: FortranArray2D<WorkingPrecision>,

    pub xt: FortranArray2D<WorkingPrecision>,
    pub yt: FortranArray2D<WorkingPrecision>,

    pub ht: FortranArray2D<WorkingPrecision>,
    pub hu: FortranArray2D<WorkingPrecision>,
    pub hv: FortranArray2D<WorkingPrecision>,

    // -1 = Water cell outside computational domain
    //  0 = Land cell
    //  1 = Water cell inside computational domain
    pub pt: FortranArray2D<i8>,
}

impl GridConstants {
    pub fn new(model_params: &ModelParameters) -> Self {
        let mut grid_constants = GridConstants {
            e1t: FortranArray2D::new(1, 1, model_params.jpi, model_params.jpj),
            e2t: FortranArray2D::new(1, 1, model_params.jpi, model_params.jpj),

            e1u: FortranArray2D::new(0, 1, model_params.jpi, model_params.jpj),
            e2u: FortranArray2D::new(0, 1, model_params.jpi, model_params.jpj),

            e1v: FortranArray2D::new(1, 0, model_params.jpi, model_params.jpj),
            e2v: FortranArray2D::new(1, 0, model_params.jpi, model_params.jpj),

            e1f: FortranArray2D::new(0, 0, model_params.jpi, model_params.jpj),
            e2f: FortranArray2D::new(0, 0, model_params.jpi, model_params.jpj),

            e12t: FortranArray2D::new(1, 1, model_params.jpi, model_params.jpj),
            e12u: FortranArray2D::new(0, 1, model_params.jpi, model_params.jpj),
            e12v: FortranArray2D::new(1, 0, model_params.jpi, model_params.jpj),

            gphiu: FortranArray2D::new(0, 1, model_params.jpi, model_params.jpj),
            gphiv: FortranArray2D::new(1, 0, model_params.jpi, model_params.jpj),
            gphif: FortranArray2D::new(0, 0, model_params.jpi, model_params.jpj),

            xt: FortranArray2D::new(1, 1, model_params.jpi, model_params.jpj),
            yt: FortranArray2D::new(1, 1, model_params.jpi, model_params.jpj),

            ht: FortranArray2D::new(1, 1, model_params.jpi, model_params.jpj),
            hu: FortranArray2D::new(0, 1, model_params.jpi, model_params.jpj),
            hv: FortranArray2D::new(1, 0, model_params.jpi, model_params.jpj),

            pt: FortranArray2D::new(0, 0, model_params.jpi + 1, model_params.jpj + 1),
        };

        grid_constants.initialise(model_params);
        grid_constants
    }

    fn initialise(&mut self, model_params: &ModelParameters) {
        let jpi = model_params.jpi;
        let jpj = model_params.jpj;

        // Water inside computational domain in inner cells
        for jj in 0..=jpj + 1 {
            for ji in 0..=jpi + 1 {
                self.pt.set(ji, jj, 1);
            }
        }

        for jj in 0..=jpj + 1 {
            self.pt.set(0, jj, 0); // Solid boundary on West
            self.pt.set(jpi + 1, jj, 0); // Solid boundary on East
        }

        // Solid boundary on North
        for ji in 0..=jpi + 1 {
            self.pt.set(ji, jpj + 1, 0);
        }

        // Open boundary (water outside computational domain) on South
        for ji in 1..=jpi {
            self.pt.set(ji, 0, -1);
        }

        self.e1t.set_all(model_params.dx);
        self.e2t.set_all(model_params.dy);

        self.e1u.set_all(model_params.dx);
        self.e2u.set_all(model_params.dy);

        self.e1v.set_all(model_params.dx);
        self.e2v.set_all(model_params.dy);

        self.e1f.set_all(model_params.dx);
        self.e2f.set_all(model_params.dy);

        self.gphiu.set_all(50.0);
        self.gphiv.set_all(50.0);
        self.gphif.set_all(50.0);

        self.ht.set_all(model_params.dep_const);
        self.hu.set_all(model_params.dep_const);
        self.hv.set_all(model_params.dep_const);

        for jj in 1..=jpj {
            for ji in 1..=jpi {
                self.e12t
                    .set(ji, jj, self.e1t.get(ji, jj) * self.e2t.get(ji, jj));

                // NOTE: The NEMOLite2D Fortran code was designed to handle a dx that
                // varies, indicating a non-linear physical grid size (different cells
                // have different sizes). Here we assume that the dx and dy are fixed and
                // not variant on the grid cell. This makes the calculation much easier
                // and makes parallelising the below xt, yt initilisation possible.
                self.xt.set(
                    ji,
                    jj,
                    self.e1t.get(ji, jj) * ((ji as WorkingPrecision) - 0.5),
                );
                self.yt.set(
                    ji,
                    jj,
                    self.e2t.get(ji, jj) * ((jj as WorkingPrecision) - 0.5),
                );
            }
        }

        for jj in 1..=jpj {
            for ji in 0..=jpi {
                self.e12u
                    .set(ji, jj, self.e1u.get(ji, jj) * self.e2u.get(ji, jj));
            }
        }

        for jj in 0..=jpj {
            for ji in 1..=jpi {
                self.e12v
                    .set(ji, jj, self.e1v.get(ji, jj) * self.e2v.get(ji, jj));
            }
        }
    }
}

pub fn run_simulation() {
    let model_params = ModelParameters::new("Config");
    println!("Model params: {:?}", model_params);

    let grid_constants = GridConstants::new(&model_params);
    let mut simulation_vars = SimulationVariables::new(&model_params, &grid_constants);

    println!("Initialised grid constants and simulation variables.");

    output_values(&model_params, &grid_constants, &simulation_vars, 0);

    let initial_step_index = model_params.initial_step_index;
    let final_step_index = model_params.final_step_index;

    let num_iters = final_step_index - initial_step_index + 1;
    let mut app_timers = AppTimers::new(num_iters as usize);

    for step_idx in initial_step_index..=final_step_index {
        step_simulation(
            &model_params,
            &grid_constants,
            &mut simulation_vars,
            step_idx,
            &mut app_timers,
        );
    }

    let ua_checksum = checksum::field_checksum(&simulation_vars.ua);
    let va_checksum = checksum::field_checksum(&simulation_vars.va);

    println!("ua checksum = {:.8E}", ua_checksum);
    println!("va checksum = {:.8E}", va_checksum);

    println!("Kernel timing report:\n{}", app_timers.generate_report());

    write_timings_csv(&app_timers);
}

fn step_simulation(
    model_params: &ModelParameters,
    grid_constants: &GridConstants,
    simulation_vars: &mut SimulationVariables,
    step_idx: u32,
    app_timers: &mut AppTimers,
) {
    let current_time = WorkingPrecision::from(step_idx) * model_params.rdt;

    app_timers.step.start();

    continuity_kernel(
        model_params,
        grid_constants,
        simulation_vars,
        &mut app_timers.continuity,
    );
    momentum_kernel(
        model_params,
        grid_constants,
        simulation_vars,
        &mut app_timers.momentum,
    );
    boundary_conditions_kernel(
        model_params,
        grid_constants,
        simulation_vars,
        current_time,
        &mut app_timers.boundary_conditions,
    );
    next_kernel(
        model_params,
        grid_constants,
        simulation_vars,
        &mut app_timers.next,
    );

    app_timers.step.stop();

    if step_idx % model_params.output_interval == 0 {
        output_values(model_params, grid_constants, simulation_vars, step_idx);
    }
}

fn output_values(
    model_params: &ModelParameters,
    grid_constants: &GridConstants,
    simulation_vars: &SimulationVariables,
    step_index: u32,
) {
    let jpi = model_params.jpi;
    let jpj = model_params.jpj;

    let field_separater = ", ";
    let line_separator = "\n";

    let headers = [
        "xt".to_owned(),
        "yt".to_owned(),
        "ht".to_owned(),
        "sshn".to_owned(),
        "un-variant".to_owned(),
        "vn-variant".to_owned(),
    ];
    let mut rows = vec![headers];
    for jj in 1..=jpj {
        for ji in 1..=jpi {
            let un_variant =
                0.5 * simulation_vars.un.get(ji - 1, jj) + simulation_vars.un.get(ji, jj);
            let vn_variant =
                0.5 * simulation_vars.vn.get(ji, jj - 1) + simulation_vars.vn.get(ji, jj);

            rows.push([
                format!("{:6.8E}", grid_constants.xt.get(ji, jj)),
                format!("{:6.8E}", grid_constants.yt.get(ji, jj)),
                format!("{:6.8E}", grid_constants.ht.get(ji, jj)),
                format!("{:6.8E}", simulation_vars.sshn.get(ji, jj)),
                format!("{:6.8E}", un_variant),
                format!("{:6.8E}", vn_variant),
            ]);
        }
    }

    let out_csv = rows
        .iter()
        .map(|row| row.join(field_separater))
        .collect::<Vec<String>>()
        .join(line_separator);

    let output_fname = format!("nemolite2d_output_{:0>5}.csv", step_index);
    let output_path = Path::new(&output_fname);
    let path_display = output_path.display();

    let mut file = match File::create(&output_path) {
        Err(why) => panic!(
            "Couldn't open output file {}: {}",
            path_display,
            why.description()
        ),
        Ok(file) => file,
    };

    match file.write_all(out_csv.as_bytes()) {
        Err(why) => panic!("Couldn't output CSV data to {}: {}", path_display, why),
        Ok(_) => println!(
            "Wrote data for step {} to file {}.",
            step_index, path_display
        ),
    }
}

fn write_timings_csv(app_timers: &AppTimers) {
    let output_fname = "nemolite2d_timings.csv";
    let output_path = Path::new(&output_fname);
    let path_display = output_path.display();

    let mut file = match File::create(&output_path) {
        Err(why) => panic!(
            "Couldn't open timings output file {}: {}",
            path_display,
            why.description(),
        ),
        Ok(file) => file,
    };

    match file.write_all(app_timers.generate_timings_csv().as_bytes()) {
        Err(why) => panic!("Couldn't output timings CSV to {}: {}", path_display, why),
        Ok(_) => println!("Wrote kernel timings data to file {}.", path_display),
    };
}
