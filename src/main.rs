#![warn(clippy::all)]

mod checksum;
mod fortran_array_2d;
mod model_parameters;

use std::cmp::{max, min};
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use fortran_array_2d::FortranArray2D;
use model_parameters::ModelParameters;

type WorkingPrecision = f64;

const AMP_TIDE: WorkingPrecision = 0.2;
const OMEGA_TIDE: WorkingPrecision = 2.0;
const GRAVITY_FORCE: WorkingPrecision = 9.80665;

struct GridConstants {
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
        for jj in 0..jpj + 1 {
            for ji in 0..jpi + 1 {
                self.pt.set(ji, jj, 1);
            }
        }

        for jj in 0..jpj + 1 {
            self.pt.set(0, jj, 0); // Solid boundary on West
            self.pt.set(jpi + 1, jj, 0); // Solid boundary on East
        }

        // Solid boundary on North
        for ji in 0..jpi + 1 {
            self.pt.set(ji, jpj + 1, 0);
        }

        // Open boundary (water outside computational domain) on South
        for ji in 1..jpi {
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

        for jj in 1..jpj {
            for ji in 1..jpi {
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

        for jj in 1..jpj {
            for ji in 0..jpi {
                self.e12u
                    .set(ji, jj, self.e1u.get(ji, jj) * self.e2u.get(ji, jj));
            }
        }

        for jj in 0..jpj {
            for ji in 1..jpi {
                self.e12v
                    .set(ji, jj, self.e1v.get(ji, jj) * self.e2v.get(ji, jj));
            }
        }
    }
}

struct SimulationVariables {
    // Sea surface height - current values
    sshn: FortranArray2D<WorkingPrecision>,
    sshn_u: FortranArray2D<WorkingPrecision>,
    sshn_v: FortranArray2D<WorkingPrecision>,

    // Sea surface height - next step's values
    ssha: FortranArray2D<WorkingPrecision>,
    ssha_u: FortranArray2D<WorkingPrecision>,
    ssha_v: FortranArray2D<WorkingPrecision>,

    // Velocities - current values
    un: FortranArray2D<WorkingPrecision>,
    vn: FortranArray2D<WorkingPrecision>,

    // Velocities - next step's values
    ua: FortranArray2D<WorkingPrecision>,
    va: FortranArray2D<WorkingPrecision>,

    // We need to double buffer the ua and va due to possible race conditions in
    // the Flather boundary conditions.
    ua_buffer: FortranArray2D<WorkingPrecision>,
    va_buffer: FortranArray2D<WorkingPrecision>,
}

impl SimulationVariables {
    pub fn new(model_params: &ModelParameters, grid_constants: &GridConstants) -> Self {
        let mut simulation_vars = SimulationVariables {
            sshn: FortranArray2D::new(1, 1, model_params.jpi, model_params.jpj),
            sshn_u: FortranArray2D::new(0, 1, model_params.jpi, model_params.jpj),
            sshn_v: FortranArray2D::new(1, 0, model_params.jpi, model_params.jpj),

            ssha: FortranArray2D::new(1, 1, model_params.jpi, model_params.jpj),
            ssha_u: FortranArray2D::new(0, 1, model_params.jpi, model_params.jpj),
            ssha_v: FortranArray2D::new(1, 0, model_params.jpi, model_params.jpj),

            un: FortranArray2D::new(0, 1, model_params.jpi, model_params.jpj),
            vn: FortranArray2D::new(1, 0, model_params.jpi, model_params.jpj),

            ua: FortranArray2D::new(0, 1, model_params.jpi, model_params.jpj),
            va: FortranArray2D::new(1, 0, model_params.jpi, model_params.jpj),

            ua_buffer: FortranArray2D::new(0, 1, model_params.jpi, model_params.jpj),
            va_buffer: FortranArray2D::new(1, 0, model_params.jpi, model_params.jpj),
        };

        simulation_vars.initialise(model_params, grid_constants);
        simulation_vars
    }

    fn initialise(&mut self, model_params: &ModelParameters, grid_constants: &GridConstants) {
        let jpi = model_params.jpi;
        let jpj = model_params.jpj;

        for jj in 1..jpj {
            for ji in 0..jpi {
                let itmp1 = min(ji + 1, jpi);
                let itmp2 = max(ji, 1);
                let rtmp1 = grid_constants.e12t.get(itmp1, jj) * self.sshn.get(itmp1, jj)
                    + grid_constants.e12t.get(itmp2, jj) * self.sshn.get(itmp2, jj);
                self.sshn_u
                    .set(ji, jj, 0.5 * rtmp1 / grid_constants.e12u.get(ji, jj));
            }
        }

        for jj in 0..jpj {
            for ji in 1..jpi {
                let itmp1 = min(jj + 1, jpj);
                let itmp2 = max(jj, 1);
                let rtmp1 = grid_constants.e12t.get(ji, itmp1) * self.sshn.get(ji, itmp1)
                    + grid_constants.e12t.get(ji, itmp2) * self.sshn.get(ji, itmp2);
                self.sshn_v
                    .set(ji, jj, 0.5 * rtmp1 / grid_constants.e12v.get(ji, jj));
            }
        }
    }
}

fn main() {
    let model_params = ModelParameters::new("Config").unwrap();
    println!("Model params: {:?}", model_params);

    let grid_constants = GridConstants::new(&model_params);
    let mut simulation_vars = SimulationVariables::new(&model_params, &grid_constants);

    println!("Initialised grid constants and simulation variables.");

    for step_idx in model_params.initial_step_index..model_params.final_step_index {
        step_simulation(
            &model_params,
            &grid_constants,
            &mut simulation_vars,
            step_idx,
        );
    }

    let ua_checksum = checksum::field_checksum(&simulation_vars.ua);
    let va_checksum = checksum::field_checksum(&simulation_vars.va);

    println!("ua checksum = {}", ua_checksum);
    println!("va checksum = {}", va_checksum);
}

fn step_simulation(
    model_params: &ModelParameters,
    grid_constants: &GridConstants,
    simulation_vars: &mut SimulationVariables,
    step_idx: u32,
) {
    let current_time = WorkingPrecision::from(step_idx) * model_params.rdt;

    continuity_kernel(model_params, grid_constants, simulation_vars);
    momentum_kernel(model_params, grid_constants, simulation_vars);
    boundary_conditions_kernel(model_params, grid_constants, simulation_vars, current_time);
    next_kernel(model_params, grid_constants, simulation_vars);

    if step_idx % model_params.output_interval == 0 {
        output_values(model_params, grid_constants, simulation_vars, step_idx);
    }
}

fn continuity_kernel(
    model_params: &ModelParameters,
    grid_constants: &GridConstants,
    simulation_vars: &mut SimulationVariables,
) {
    let jpi = model_params.jpi;
    let jpj = model_params.jpj;

    for jj in 1..jpj {
        for ji in 1..jpi {
            let rtmp1 = (simulation_vars.sshn_u.get(ji, jj) + grid_constants.hu.get(ji, jj))
                * simulation_vars.un.get(ji, jj);

            let rtmp2 = (simulation_vars.sshn_u.get(ji - 1, jj)
                + grid_constants.hu.get(ji - 1, jj))
                * simulation_vars.un.get(ji - 1, jj);

            let rtmp3 = (simulation_vars.sshn_v.get(ji, jj) + grid_constants.hv.get(ji, jj))
                * simulation_vars.vn.get(ji, jj);

            let rtmp4 = (simulation_vars.sshn_v.get(ji, jj - 1)
                + grid_constants.hv.get(ji, jj - 1))
                * simulation_vars.vn.get(ji, jj - 1);

            simulation_vars.ssha.set(
                ji,
                jj,
                simulation_vars.sshn.get(ji, jj)
                    + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * model_params.rdt
                        / grid_constants.e12t.get(ji, jj),
            );
        }
    }
}

fn momentum_kernel(
    model_params: &ModelParameters,
    grid_constants: &GridConstants,
    simulation_vars: &mut SimulationVariables,
) {
    let jpi = model_params.jpi;
    let jpj = model_params.jpj;
}

fn boundary_conditions_kernel(
    model_params: &ModelParameters,
    grid_constants: &GridConstants,
    simulation_vars: &mut SimulationVariables,
    current_time: WorkingPrecision,
) {
    let jpi = model_params.jpi;
    let jpj = model_params.jpj;

    let pt = &grid_constants.pt;

    // Apply open boundary conditions of clamped sea surface height
    for jj in 1..jpj {
        for ji in 1..jpi {
            if pt.get(ji, jj) <= 0 {
                continue;
            }

            if pt.get(ji, jj - 1) < 0
                || pt.get(ji, jj + 1) < 0
                || pt.get(ji - 1, jj) < 0
                || pt.get(ji + 1, jj) < 0
            {
                simulation_vars
                    .ssha
                    .set(ji, jj, AMP_TIDE * (OMEGA_TIDE * current_time).sin());
            }
        }
    }

    // Apply solid boundary conditions for u-velocity
    for jj in 1..jpj {
        for ji in 0..jpi {
            if pt.get(ji, jj) * pt.get(ji + 1, jj) == 0 {
                simulation_vars.ua.set(ji, jj, 0.0);
            }
        }
    }

    // Apply solid boundary conditions for v-velocity
    for jj in 0..jpj {
        for ji in 1..jpi {
            if pt.get(ji, jj) * pt.get(ji, jj + 1) == 0 {
                simulation_vars.va.set(ji, jj, 0.0);
            }
        }
    }

    // Flather open boundary condition for u
    for jj in 1..jpj {
        for ji in 0..jpi {
            if pt.get(ji, jj) + pt.get(ji + 1, jj) <= -1 {
                continue;
            }

            if pt.get(ji, jj) < 0 {
                let jiu = ji + 1;
                simulation_vars.ua_buffer.set(
                    ji,
                    jj,
                    simulation_vars.ua.get(jiu, jj)
                        + (GRAVITY_FORCE / grid_constants.hu.get(ji, jj)).sqrt()
                            * (simulation_vars.sshn_u.get(ji, jj)
                                - simulation_vars.sshn_u.get(jiu, jj)),
                );
            } else if pt.get(ji + 1, jj) < 0 {
                let jiu = ji - 1;
                simulation_vars.ua_buffer.set(
                    ji,
                    jj,
                    simulation_vars.ua.get(jiu, jj)
                        + (GRAVITY_FORCE / grid_constants.hu.get(ji, jj)).sqrt()
                            * (simulation_vars.sshn_u.get(ji, jj)
                                - simulation_vars.sshn_u.get(jiu, jj)),
                );
            }
        }
    }

    // Flather open boundary condition for v
    for jj in 0..jpj {
        for ji in 1..jpi {
            if pt.get(ji, jj) + pt.get(ji, jj + 1) <= -1 {
                continue;
            }

            if pt.get(ji, jj) < 0 {
                let jiv = jj + 1;
                simulation_vars.va_buffer.set(
                    ji,
                    jj,
                    simulation_vars.va.get(ji, jiv)
                        + (GRAVITY_FORCE / grid_constants.hv.get(ji, jj).sqrt()
                            * (simulation_vars.sshn_v.get(ji, jj)
                                - simulation_vars.sshn_v.get(ji, jiv))),
                );
            } else if pt.get(ji, jj + 1) < 0 {
                let jiv = jj - 1;
                simulation_vars.va_buffer.set(
                    ji,
                    jj,
                    simulation_vars.va.get(ji, jiv)
                        + (GRAVITY_FORCE / grid_constants.hv.get(ji, jj).sqrt()
                            * (simulation_vars.sshn_v.get(ji, jj)
                                - simulation_vars.sshn_v.get(ji, jiv))),
                );
            }
        }
    }
}

fn next_kernel(
    model_params: &ModelParameters,
    grid_constants: &GridConstants,
    simulation_vars: &mut SimulationVariables,
) {
    let jpi = model_params.jpi;
    let jpj = model_params.jpj;

    for jj in 1..jpj {
        for ji in 0..jpi {
            simulation_vars
                .un
                .set(ji, jj, simulation_vars.ua.get(ji, jj));
        }
    }

    for jj in 0..jpj {
        for ji in 1..jpi {
            simulation_vars
                .vn
                .set(ji, jj, simulation_vars.va.get(ji, jj));
        }
    }

    for jj in 1..jpj {
        for ji in 1..jpi {
            simulation_vars
                .sshn
                .set(ji, jj, simulation_vars.ssha.get(ji, jj));
        }
    }

    for jj in 1..jpj {
        for ji in 0..jpi {
            let this_cell_type = grid_constants.pt.get(ji, jj);
            let next_cell_type = grid_constants.pt.get(ji + 1, jj);

            if this_cell_type + next_cell_type <= 0 {
                continue;
            }

            if this_cell_type == 1 && next_cell_type == 1 {
                let rtmp1 = grid_constants.e12t.get(ji, jj) * simulation_vars.sshn.get(ji, jj)
                    + grid_constants.e12t.get(ji + 1, jj) * simulation_vars.sshn.get(ji + 1, jj);

                simulation_vars
                    .sshn_u
                    .set(ji, jj, 0.5 * rtmp1 / grid_constants.e12u.get(ji, jj));
            } else if this_cell_type <= 0 {
                simulation_vars
                    .ssha_u
                    .set(ji, jj, simulation_vars.sshn.get(ji + 1, jj));
            } else if next_cell_type <= 0 {
                simulation_vars
                    .sshn_u
                    .set(ji, jj, simulation_vars.sshn.get(ji, jj));
            }
        }
    }

    for jj in 0..jpj {
        for ji in 1..jpi {
            let this_cell_type = grid_constants.pt.get(ji, jj);
            let next_cell_type = grid_constants.pt.get(ji, jj + 1);

            if this_cell_type + next_cell_type <= 0 {
                continue;
            }

            if this_cell_type == 1 && next_cell_type == 1 {
                let rtmp1 = grid_constants.e12t.get(ji, jj) * simulation_vars.sshn.get(ji, jj)
                    + grid_constants.e12t.get(ji, jj + 1) * simulation_vars.sshn.get(ji, jj + 1);

                simulation_vars
                    .sshn_u
                    .set(ji, jj, 0.5 * rtmp1 / grid_constants.e12v.get(ji, jj));
            } else if this_cell_type <= 0 {
                simulation_vars
                    .ssha_u
                    .set(ji, jj, simulation_vars.sshn.get(ji, jj + 1));
            } else if next_cell_type <= 0 {
                simulation_vars
                    .sshn_u
                    .set(ji, jj, simulation_vars.sshn.get(ji, jj));
            }
        }
    }
}

fn output_values(
    model_params: &ModelParameters,
    grid_constants: &GridConstants,
    simulation_vars: &SimulationVariables,
    step_index: u32,
) {
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
    for jj in 1..model_params.jpj {
        for ji in 1..model_params.jpi {
            let un_variant = 0.5 * simulation_vars.un.get(ji, jj) + simulation_vars.un.get(ji, jj);
            let vn_variant =
                0.5 * simulation_vars.vn.get(ji, jj - 1) + simulation_vars.vn.get(ji, jj);

            rows.push([
                format!("{}", grid_constants.xt.get(ji, jj)),
                format!("{}", grid_constants.yt.get(ji, jj)),
                format!("{}", grid_constants.ht.get(ji, jj)),
                format!("{}", simulation_vars.sshn.get(ji, jj)),
                format!("{}", un_variant),
                format!("{}", vn_variant),
            ]);
        }
    }

    let out_csv = rows
        .iter()
        .map(|row| row.join(field_separater))
        .collect::<Vec<String>>()
        .join(line_separator);

    let output_fname = format!("nemolite2d_output_{}.csv", step_index);
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
