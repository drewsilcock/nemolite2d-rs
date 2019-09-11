use crate::{FortranArray2D, ModelParameters, WorkingPrecision};

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
