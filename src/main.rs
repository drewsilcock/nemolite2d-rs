#![warn(clippy::all)]

mod checksum;
mod model_parameters;
mod fortran_array_2d;

use std::cmp::{min, max};

use nalgebra::DMatrix;

use fortran_array_2d::FortranArray2D;
use model_parameters::ModelParameters;

type WorkingPrecision = f64;

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
        for ji in 0.. jpi + 1 {
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

        for ji in 1..jpi {
            for jj in 1..jpj {
                self.e12t.set(ji, jj, self.e1t.get(ji, jj) * self.e2t.get(ji, jj));

                // NOTE: The NEMOLite2D Fortran code was designed to handle a dx that
                // varies, indicating a non-linear physical grid size (different cells
                // have different sizes). Here we assume that the dx and dy are fixed and
                // not variant on the grid cell. This makes the calculation much easier
                // and makes parallelising the below xt, yt initilisation possible.
                self.xt.set(ji, jj, self.e1t.get(ji, jj) * ((ji as WorkingPrecision) - 0.5));
                self.yt.set(ji, jj, self.e2t.get(ji, jj) * ((jj as WorkingPrecision) - 0.5));
            }
        }

        for ji in 0..jpi {
            for jj in 1..jpj {
                self.e12u.set(ji, jj, self.e1u.get(ji, jj) * self.e2u.get(ji, jj));
            }
        }

        for ji in 1..jpi {
            for jj in 0..jpj {
                self.e12v.set(ji, jj, self.e1v.get(ji, jj) * self.e2v.get(ji, jj));
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

        for ji in 0..jpi {
            for jj in 1..jpj {
                let itmp1 = min(ji + 1, jpi);
                let itmp2 = max(ji, 1);
                let rtmp1 = grid_constants.e12t.get(itmp1, jj) * self.sshn.get(itmp1, jj) + grid_constants.e12t.get(itmp2, jj) * self.sshn.get(itmp2, jj);
                self.sshn_u.set(ji, jj, 0.5 * rtmp1 / grid_constants.e12u.get(ji, jj));
            }
        }

        for ji in 1..jpi {
            for jj in 0..jpj {
                let itmp1 = min(jj + 1, jpj);
                let itmp2 = max(jj, 1);
                let rtmp1 = grid_constants.e12t.get(ji, itmp1) * self.sshn.get(ji, itmp1) + grid_constants.e12t.get(ji, itmp2) * self.sshn.get(ji, itmp2);
                self.sshn_v.set(ji, jj, 0.5 * rtmp1 / grid_constants.e12v.get(ji, jj));
            }
        }
    }
}

fn main() {
    let model_params = ModelParameters::new("Config").unwrap();
    println!("Model params: {:?}", model_params);

    let grid_constants = GridConstants::new(&model_params);
    let simulation_vars = SimulationVariables::new(&model_params, &grid_constants);

    println!("Initialised grid constants and simulation variables.");

    for step_idx in model_params.initial_step_index..model_params.final_step_index {
        step_simulation(&model_params, &grid_constants, &mut simulation_vars, step_idx);
    }

    let ua_checksum = checksum::field_checksum(&simulation_vars.ua);
    let va_checksum = checksum::field_checksum(&simulation_vars.va);

    println!("ua checksum = {}", ua_checksum);
    println!("va checksum = {}", va_checksum);
}

fn step_simulation(model_params: &ModelParameters, grid_constants: &GridConstants, simulation_vars: &mut SimulationVariables, step_idx: u32) {
    let current_time = WorkingPrecision::from(step_idx) * model_params.rdt;

    continuity_kernel(model_params, grid_constants, simulation_vars);
    momentum_kernel(model_params, grid_constants, simulation_vars);
    boundary_conditions_kernel(model_params, grid_constants, simulation_vars);
    next_kernel(model_params, grid_constants, simulation_vars);

    if (step_idx % model_params.output_interval == 0) {
        output_values(model_params, grid_constants, simulation_vars);
    }
}

fn continuity_kernel(model_params: &ModelParameters, grid_constants: &GridConstants, simulation_vars: &mut SimulationVariables) {
    //
}

fn momentum_kernel(model_params: &ModelParameters, grid_constants: &GridConstants, simulation_vars: &mut SimulationVariables) {
    //
}

fn boundary_conditions_kernel(model_params: &ModelParameters, grid_constants: &GridConstants, simulation_vars: &mut SimulationVariables) {
    //
}

fn next_kernel(model_params: &ModelParameters, grid_constants: &GridConstants, simulation_vars: &mut SimulationVariables) {
    //
}

fn output_values(model_params: &ModelParameters, grid_constants: &GridConstants, simulation_vars: &SimulationVariables) {
    //
}
