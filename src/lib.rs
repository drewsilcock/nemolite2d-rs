#![warn(clippy::all)]

mod app_timers;
mod checksum;
mod fortran_array_2d;
mod grid_constants;
mod kernels;
mod model_parameters;
mod simulation_variables;
mod timer;

use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use app_timers::AppTimers;
use fortran_array_2d::FortranArray2D;
use grid_constants::GridConstants;
use kernels::{boundary_conditions_kernel, continuity_kernel, momentum_kernel, next_kernel};
use model_parameters::ModelParameters;
use simulation_variables::SimulationVariables;

type WorkingPrecision = f64;

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
