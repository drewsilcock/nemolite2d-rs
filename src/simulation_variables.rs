use std::cmp::{max, min};

use crate::{ModelParameters, GridConstants, FortranArray2D, WorkingPrecision};

pub struct SimulationVariables {
    // Sea surface height - current values
    pub sshn: FortranArray2D<WorkingPrecision>,
    pub sshn_u: FortranArray2D<WorkingPrecision>,
    pub sshn_v: FortranArray2D<WorkingPrecision>,

    // Sea surface height - next step's values
    pub ssha: FortranArray2D<WorkingPrecision>,
    pub ssha_u: FortranArray2D<WorkingPrecision>,
    pub ssha_v: FortranArray2D<WorkingPrecision>,

    // Velocities - current values
    pub un: FortranArray2D<WorkingPrecision>,
    pub vn: FortranArray2D<WorkingPrecision>,

    // Velocities - next step's values
    pub ua: FortranArray2D<WorkingPrecision>,
    pub va: FortranArray2D<WorkingPrecision>,
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
        };

        simulation_vars.initialise(model_params, grid_constants);
        simulation_vars
    }

    fn initialise(&mut self, model_params: &ModelParameters, grid_constants: &GridConstants) {
        let jpi = model_params.jpi;
        let jpj = model_params.jpj;

        for jj in 1..=jpj {
            for ji in 0..=jpi {
                let itmp1 = min(ji + 1, jpi);
                let itmp2 = max(ji, 1);
                let rtmp1 = grid_constants.e12t.get(itmp1, jj) * self.sshn.get(itmp1, jj)
                    + grid_constants.e12t.get(itmp2, jj) * self.sshn.get(itmp2, jj);
                self.sshn_u
                    .set(ji, jj, 0.5 * rtmp1 / grid_constants.e12u.get(ji, jj));
            }
        }

        for jj in 0..=jpj {
            for ji in 1..=jpi {
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
