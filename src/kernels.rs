use crate::timer::Timer;
use crate::{GridConstants, ModelParameters, SimulationVariables, WorkingPrecision};

const PI: WorkingPrecision = 3.141_592_653_589_793_2;
const AMP_TIDE: WorkingPrecision = 0.2;
const OMEGA_TIDE: WorkingPrecision = 2.0 * 3.14159 / (12.42 * 3600.0);
const GRAVITY_FORCE: WorkingPrecision = 9.80665;
const EARTH_ROTATIONAL_SPEED: WorkingPrecision = 7.292_116E-05; // Earth rotation speed (s^(-1)) a.k.a. omega
const DEGREES_TO_RADIANS: WorkingPrecision = PI / 180.0;

pub fn continuity_kernel(
    model_params: &ModelParameters,
    grid_constants: &GridConstants,
    simulation_vars: &mut SimulationVariables,
    kernel_timer: &mut Timer,
) {
    let jpi = model_params.jpi;
    let jpj = model_params.jpj;
    let rdt = model_params.rdt;

    let e12t = &grid_constants.e12t;
    let hu = &grid_constants.hu;
    let hv = &grid_constants.hv;

    let sshn = &simulation_vars.sshn;
    let sshn_u = &simulation_vars.sshn_u;
    let sshn_v = &simulation_vars.sshn_v;
    let un = &simulation_vars.un;
    let vn = &simulation_vars.vn;

    kernel_timer.start();

    for jj in 1..=jpj {
        for ji in 1..=jpi {
            let rtmp1 = (sshn_u.get(ji, jj) + hu.get(ji, jj)) * un.get(ji, jj);
            let rtmp2 = (sshn_u.get(ji - 1, jj) + hu.get(ji - 1, jj)) * un.get(ji - 1, jj);
            let rtmp3 = (sshn_v.get(ji, jj) + hv.get(ji, jj)) * vn.get(ji, jj);
            let rtmp4 = (sshn_v.get(ji, jj - 1) + hv.get(ji, jj - 1)) * vn.get(ji, jj - 1);

            simulation_vars.ssha.set(
                ji,
                jj,
                sshn.get(ji, jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * rdt / e12t.get(ji, jj),
            );
        }
    }

    kernel_timer.stop();
}

pub fn momentum_kernel(
    model_params: &ModelParameters,
    grid_constants: &GridConstants,
    simulation_vars: &mut SimulationVariables,
    kernel_timer: &mut Timer,
) {
    let jpi = model_params.jpi;
    let jpj = model_params.jpj;
    let rdt = model_params.rdt;
    let cbfr = model_params.cbfr;
    let visc = model_params.visc;

    let e1u = &grid_constants.e1u;
    let e2u = &grid_constants.e2u;
    let e1v = &grid_constants.e1v;
    let e2v = &grid_constants.e2v;
    let e1t = &grid_constants.e1t;
    let e2t = &grid_constants.e2t;
    let e12u = &grid_constants.e12u;
    let e12v = &grid_constants.e12v;
    let ht = &grid_constants.ht;
    let hu = &grid_constants.hu;
    let hv = &grid_constants.hv;
    let gphiu = &grid_constants.gphiu;
    let gphiv = &grid_constants.gphiv;
    let pt = &grid_constants.pt;

    let sshn = &simulation_vars.sshn;
    let sshn_u = &simulation_vars.sshn_u;
    let sshn_v = &simulation_vars.sshn_v;
    let ssha_u = &simulation_vars.ssha_u;
    let ssha_v = &simulation_vars.ssha_v;
    let un = &simulation_vars.un;
    let vn = &simulation_vars.vn;

    kernel_timer.start();

    for jj in 1..=jpj {
        for ji in 1..=jpi - 1 {
            // Jump over non-computational domain
            if pt.get(ji, jj) + pt.get(ji + 1, jj) <= 0 {
                continue;
            }

            // Jump over u boundary
            if pt.get(ji, jj) <= 0 || pt.get(ji + 1, jj) <= 0 {
                continue;
            }

            let u_e = 0.5 * (un.get(ji, jj) + un.get(ji + 1, jj)) * e2t.get(ji + 1, jj);
            let depe = ht.get(ji + 1, jj) + sshn.get(ji + 1, jj);

            let u_w = 0.5 * (un.get(ji, jj) + un.get(ji - 1, jj)) * e2t.get(ji, jj);
            let depw = ht.get(ji, jj) + sshn.get(ji, jj);

            let v_sc = 0.5 * (vn.get(ji, jj - 1) + vn.get(ji + 1, jj - 1));
            let v_s = 0.5 * v_sc * (e1v.get(ji, jj - 1) + e1v.get(ji + 1, jj - 1));
            let deps = 0.5
                * (hv.get(ji, jj - 1)
                    + sshn_v.get(ji, jj - 1)
                    + hv.get(ji + 1, jj - 1)
                    + sshn_v.get(ji + 1, jj - 1));

            let v_nc = 0.5 * (vn.get(ji, jj) + vn.get(ji + 1, jj));
            let v_n = 0.5 * v_nc * (e1v.get(ji, jj) + e1v.get(ji + 1, jj));
            let depn = 0.5
                * (hv.get(ji, jj)
                    + sshn_v.get(ji, jj)
                    + hv.get(ji + 1, jj)
                    + sshn_v.get(ji + 1, jj));

            // Advection (currently first order upwind)
            let uu_w = (0.5 - (0.5 as WorkingPrecision).copysign(u_w)) * un.get(ji, jj)
                + (0.5 + (0.5 as WorkingPrecision).copysign(u_w)) * un.get(ji - 1, jj);
            let uu_e = (0.5 + (0.5 as WorkingPrecision).copysign(u_e)) * un.get(ji, jj)
                + (0.5 - (0.5 as WorkingPrecision).copysign(u_e)) * un.get(ji + 1, jj);

            let uu_s = if pt.get(ji, jj - 1) <= 0 || pt.get(ji + 1, jj - 1) <= 0 {
                (0.5 - (0.5 as WorkingPrecision).copysign(v_s)) * un.get(ji, jj)
            } else {
                (0.5 - (0.5 as WorkingPrecision).copysign(v_s)) * un.get(ji, jj)
                    + (0.5 + (0.5 as WorkingPrecision).copysign(v_s)) * un.get(ji, jj - 1)
            };

            let uu_n = if pt.get(ji, jj + 1) <= 0 || pt.get(ji + 1, jj + 1) <= 0 {
                (0.5 + (0.5 as WorkingPrecision).copysign(v_n)) * un.get(ji, jj)
            } else {
                (0.5 + (0.5 as WorkingPrecision).copysign(v_n)) * un.get(ji, jj)
                    + (0.5 - (0.5 as WorkingPrecision).copysign(v_n)) * un.get(ji, jj + 1)
            };

            let adv = uu_w * u_w * depw - uu_e * u_e * depe + uu_s * v_s * depe + uu_s * v_s * deps
                - uu_n * v_n * depn;

            // Viscosity
            let dudx_e = (un.get(ji + 1, jj) - un.get(ji, jj)) / e1t.get(ji + 1, jj)
                * (ht.get(ji + 1, jj) + sshn.get(ji + 1, jj));
            let dudx_w = (un.get(ji, jj) - un.get(ji - 1, jj)) / e1t.get(ji, jj)
                * (ht.get(ji, jj) + sshn.get(ji, jj));

            let dudy_s = if pt.get(ji, jj - 1) <= 0 || pt.get(ji + 1, jj - 1) <= 0 {
                0.0 // Slip boundary
            } else {
                (un.get(ji, jj) - un.get(ji, jj - 1)) / (e2u.get(ji, jj) + e2u.get(ji, jj - 1))
                    * (hu.get(ji, jj)
                        + sshn_u.get(ji, jj)
                        + hu.get(ji, jj - 1)
                        + sshn_u.get(ji, jj - 1))
            };

            let dudy_n = if pt.get(ji, jj + 1) <= 0 || pt.get(ji + 1, jj + 1) <= 0 {
                0.0 // Slip boundary
            } else {
                (un.get(ji, jj + 1) - un.get(ji, jj)) / (e2u.get(ji, jj) + e2u.get(ji, jj + 1))
                    * (hu.get(ji, jj)
                        + sshn_u.get(ji, jj)
                        + hu.get(ji, jj + 1)
                        + sshn_u.get(ji, jj + 1))
            };

            let vis = visc * (dudx_e - dudx_w) * e2u.get(ji, jj)
                + (dudy_n - dudy_s) * e1u.get(ji, jj) * 0.5;

            // Coriolis' force (can be implemented implicitly)
            let cor = 0.5
                * (2.0
                    * EARTH_ROTATIONAL_SPEED
                    * (gphiu.get(ji, jj) * DEGREES_TO_RADIANS).sin()
                    * (v_sc + v_nc))
                * e12u.get(ji, jj)
                * (hu.get(ji, jj) + sshn_u.get(ji, jj));

            // Pressure gradient
            let hpg = -GRAVITY_FORCE
                * (hu.get(ji, jj) + sshn_u.get(ji, jj))
                * e2u.get(ji, jj)
                * (sshn.get(ji + 1, jj) - sshn.get(ji, jj));

            // Linear bottom friction (implemented implicitly.

            // Final ua calculation based on combining all the other factors
            let ua_value = (un.get(ji, jj) * (hu.get(ji, jj) + sshn_u.get(ji, jj))
                + rdt * (adv + vis + cor + hpg) / e12u.get(ji, jj))
                / (hu.get(ji, jj) + ssha_u.get(ji, jj))
                / (1.0 + cbfr * rdt);
            simulation_vars.ua.set(ji, jj, ua_value);
        }
    }

    for jj in 1..=jpj - 1 {
        for ji in 1..=jpi {
            if pt.get(ji, jj) + pt.get(ji + 1, jj) <= 0 {
                continue; // Jump over non-computatinal domain
            }

            if pt.get(ji, jj) <= 0 || pt.get(ji, jj + 1) <= 0 {
                continue; // Jump over v boundary cells
            }

            // Advection
            let v_n = 0.5 * (vn.get(ji, jj) + vn.get(ji, jj + 1)) * e1t.get(ji, jj + 1); // Add length scale.
            let depn = ht.get(ji, jj + 1) + sshn.get(ji, jj + 1);

            let v_s = 0.5 * (vn.get(ji, jj) + vn.get(ji, jj - 1)) * e1t.get(ji, jj); // Add length scale
            let deps = ht.get(ji, jj) + sshn.get(ji, jj);

            let u_wc = 0.5 * (un.get(ji - 1, jj) + un.get(ji - 1, jj + 1));
            let u_w = 0.5 * u_wc * (e2u.get(ji - 1, jj) + e2u.get(ji - 1, jj + 1));
            let depw = 0.5
                * (hu.get(ji - 1, jj)
                    + sshn_u.get(ji - 1, jj)
                    + hu.get(ji - 1, jj + 1)
                    + sshn_u.get(ji - 1, jj + 1));

            let u_ec = 0.5 * (un.get(ji, jj) + un.get(ji, jj + 1));
            let u_e = 0.5 * u_ec * (e2u.get(ji, jj) + e2u.get(ji, jj + 1));
            let depe = 0.5
                * (hu.get(ji, jj)
                    + sshn_u.get(ji, jj)
                    + hu.get(ji, jj + 1)
                    + sshn_u.get(ji, jj + 1));

            // Advection (currently first order upwind)
            let vv_s = (0.5 - (0.5 as WorkingPrecision).copysign(v_s)) * vn.get(ji, jj)
                + (0.5 + (0.5 as WorkingPrecision).copysign(v_s)) * vn.get(ji, jj - 1);
            let vv_n = (0.5 + (0.5 as WorkingPrecision).copysign(v_n)) * vn.get(ji, jj)
                + (0.5 - (0.5 as WorkingPrecision).copysign(v_n)) * vn.get(ji, jj + 1);

            let vv_w = if pt.get(ji - 1, jj) <= 0 || pt.get(ji - 1, jj + 1) <= 0 {
                (0.5 - (0.5 as WorkingPrecision).copysign(u_w)) * vn.get(ji, jj)
            } else {
                (0.5 - (0.5 as WorkingPrecision).copysign(u_w)) * vn.get(ji, jj)
                    + (0.5 + (0.5 as WorkingPrecision).copysign(u_w)) * vn.get(ji - 1, jj)
            };

            let vv_e = if pt.get(ji + 1, jj) <= 0 || pt.get(ji + 1, jj + 1) <= 0 {
                (0.5 + (0.5 as WorkingPrecision).copysign(u_e)) * vn.get(ji, jj)
            } else {
                (0.5 + (0.5 as WorkingPrecision).copysign(u_e)) * vn.get(ji, jj)
                    + (0.5 - (0.5 as WorkingPrecision).copysign(u_e)) * vn.get(ji + 1, jj)
            };

            let adv = vv_w * u_w * depw - vv_e * u_e * depe + vv_s * v_s * deps - vv_n * v_n * depn;

            // Viscosity
            let dvdy_n = (vn.get(ji, jj + 1) - vn.get(ji, jj)) / e2t.get(ji, jj + 1)
                * (ht.get(ji, jj + 1) + sshn.get(ji, jj + 1));
            let dvdy_s = (vn.get(ji, jj) - vn.get(ji, jj - 1)) / e2t.get(ji, jj)
                * (ht.get(ji, jj) + sshn.get(ji, jj));

            let dvdx_w = if pt.get(ji - 1, jj) <= 0 || pt.get(ji - 1, jj + 1) <= 0 {
                0.0 // Slip boundary
            } else {
                (vn.get(ji, jj) - vn.get(ji - 1, jj)) / (e1v.get(ji, jj) + e1v.get(ji - 1, jj))
                    * (hv.get(ji, jj)
                        + sshn_v.get(ji, jj)
                        + hv.get(ji - 1, jj)
                        + sshn_v.get(ji - 1, jj))
            };

            let dvdx_e = if pt.get(ji + 1, jj) <= 0 || pt.get(ji + 1, jj + 1) <= 0 {
                0.0 // Slip boundary
            } else {
                (vn.get(ji + 1, jj) - vn.get(ji, jj)) / (e1v.get(ji, jj) + e1v.get(ji + 1, jj))
                    * (hv.get(ji, jj)
                        + sshn_v.get(ji, jj)
                        + hv.get(ji + 1, jj)
                        + sshn_v.get(ji + 1, jj))
            };

            let vis = visc * (dvdy_n - dvdy_s) * e1v.get(ji, jj)
                + (dvdx_e - dvdx_w) * e2v.get(ji, jj) * 0.5;

            // Coriolis' force (can be implemented implicitly)
            let cor = -0.5
                * (2.0
                    * EARTH_ROTATIONAL_SPEED
                    * (gphiv.get(ji, jj) * DEGREES_TO_RADIANS).sin()
                    * (u_ec + u_wc))
                * e12v.get(ji, jj)
                * (hv.get(ji, jj) + sshn_v.get(ji, jj));

            // Pressure gradient
            let hpg = -GRAVITY_FORCE
                * (hv.get(ji, jj) + sshn_v.get(ji, jj))
                * e1v.get(ji, jj)
                * (sshn.get(ji, jj + 1) - sshn.get(ji, jj));

            // Linear bottom friction (implemented implicitly.

            // Final va calculation based on combining all the other factors
            let va_value = (vn.get(ji, jj) * (hv.get(ji, jj) + sshn_v.get(ji, jj))
                + rdt * (adv + vis + cor + hpg) / e12v.get(ji, jj))
                / (hv.get(ji, jj) + ssha_v.get(ji, jj))
                / (1.0 + cbfr * rdt);
            simulation_vars.va.set(ji, jj, va_value);
        }
    }

    kernel_timer.stop();
}

pub fn boundary_conditions_kernel(
    model_params: &ModelParameters,
    grid_constants: &GridConstants,
    simulation_vars: &mut SimulationVariables,
    current_time: WorkingPrecision,
    kernel_timer: &mut Timer,
) {
    let jpi = model_params.jpi;
    let jpj = model_params.jpj;

    let hu = &grid_constants.hu;
    let hv = &grid_constants.hv;
    let pt = &grid_constants.pt;

    let sshn_u = &simulation_vars.sshn_u;
    let sshn_v = &simulation_vars.sshn_v;

    kernel_timer.start();

    // Apply open boundary conditions of clamped sea surface height
    for jj in 1..=jpj {
        for ji in 1..=jpi {
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
    for jj in 1..=jpj {
        for ji in 0..=jpi {
            if pt.get(ji, jj) * pt.get(ji + 1, jj) == 0 {
                simulation_vars.ua.set(ji, jj, 0.0);
            }
        }
    }

    // Apply solid boundary conditions for v-velocity
    for jj in 0..=jpj {
        for ji in 1..=jpi {
            if pt.get(ji, jj) * pt.get(ji, jj + 1) == 0 {
                simulation_vars.va.set(ji, jj, 0.0);
            }
        }
    }

    // Flather open boundary condition for u
    for jj in 1..=jpj {
        for ji in 0..=jpi {
            if pt.get(ji, jj) + pt.get(ji + 1, jj) <= -1 {
                continue;
            }

            if pt.get(ji, jj) < 0 {
                let jiu = ji + 1;
                let flather_ua = simulation_vars.ua.get(jiu, jj)
                    + (GRAVITY_FORCE / hu.get(ji, jj)).sqrt()
                        * (sshn_u.get(ji, jj) - sshn_u.get(jiu, jj));
                simulation_vars.ua.set(ji, jj, flather_ua);
            } else if pt.get(ji + 1, jj) < 0 {
                let jiu = ji - 1;
                let flather_ua = simulation_vars.ua.get(jiu, jj)
                    + (GRAVITY_FORCE / hu.get(ji, jj)).sqrt()
                        * (sshn_u.get(ji, jj) - sshn_u.get(jiu, jj));
                simulation_vars.ua.set(ji, jj, flather_ua);
            }
        }
    }

    // Flather open boundary condition for v
    for jj in 0..=jpj {
        for ji in 1..=jpi {
            if pt.get(ji, jj) + pt.get(ji, jj + 1) <= -1 {
                continue;
            }

            if pt.get(ji, jj) < 0 {
                let jiv = jj + 1;
                let flather_va = simulation_vars.va.get(ji, jiv)
                    + (GRAVITY_FORCE / hv.get(ji, jj)).sqrt()
                        * (sshn_v.get(ji, jj) - sshn_v.get(ji, jiv));
                simulation_vars.va.set(ji, jj, flather_va);
            } else if pt.get(ji, jj + 1) < 0 {
                let jiv = jj - 1;
                let flather_va = simulation_vars.va.get(ji, jiv)
                    + (GRAVITY_FORCE / hv.get(ji, jj)).sqrt()
                        * (sshn_v.get(ji, jj) - sshn_v.get(ji, jiv));
                simulation_vars.va.set(ji, jj, flather_va);
            }
        }
    }

    kernel_timer.stop();
}

pub fn next_kernel(
    model_params: &ModelParameters,
    grid_constants: &GridConstants,
    simulation_vars: &mut SimulationVariables,
    kernel_timer: &mut Timer,
) {
    let jpi = model_params.jpi;
    let jpj = model_params.jpj;

    let e12t = &grid_constants.e12t;
    let e12u = &grid_constants.e12u;
    let e12v = &grid_constants.e12v;
    let pt = &grid_constants.pt;

    kernel_timer.start();

    for jj in 1..=jpj {
        for ji in 0..=jpi {
            simulation_vars
                .un
                .set(ji, jj, simulation_vars.ua.get(ji, jj));
        }
    }

    for jj in 0..=jpj {
        for ji in 1..=jpi {
            simulation_vars
                .vn
                .set(ji, jj, simulation_vars.va.get(ji, jj));
        }
    }

    for jj in 1..=jpj {
        for ji in 1..=jpi {
            simulation_vars
                .sshn
                .set(ji, jj, simulation_vars.ssha.get(ji, jj));
        }
    }

    for jj in 1..=jpj {
        for ji in 0..=jpi {
            let this_cell_type = pt.get(ji, jj);
            let next_cell_type = pt.get(ji + 1, jj);

            if this_cell_type + next_cell_type <= 0 {
                continue;
            }

            if this_cell_type == 1 && next_cell_type == 1 {
                let rtmp1 = e12t.get(ji, jj) * simulation_vars.sshn.get(ji, jj)
                    + e12t.get(ji + 1, jj) * simulation_vars.sshn.get(ji + 1, jj);

                simulation_vars
                    .sshn_u
                    .set(ji, jj, 0.5 * rtmp1 / e12u.get(ji, jj));
            } else if this_cell_type <= 0 {
                simulation_vars
                    .sshn_u
                    .set(ji, jj, simulation_vars.sshn.get(ji + 1, jj));
            } else if next_cell_type <= 0 {
                simulation_vars
                    .sshn_u
                    .set(ji, jj, simulation_vars.sshn.get(ji, jj));
            }
        }
    }

    for jj in 0..=jpj {
        for ji in 1..=jpi {
            let this_cell_type = pt.get(ji, jj);
            let next_cell_type = pt.get(ji, jj + 1);

            if this_cell_type + next_cell_type <= 0 {
                continue;
            }

            if this_cell_type == 1 && next_cell_type == 1 {
                let rtmp1 = e12t.get(ji, jj) * simulation_vars.sshn.get(ji, jj)
                    + e12t.get(ji, jj + 1) * simulation_vars.sshn.get(ji, jj + 1);

                simulation_vars
                    .sshn_v
                    .set(ji, jj, 0.5 * rtmp1 / e12v.get(ji, jj));
            } else if this_cell_type <= 0 {
                simulation_vars
                    .sshn_v
                    .set(ji, jj, simulation_vars.sshn.get(ji, jj + 1));
            } else if next_cell_type <= 0 {
                simulation_vars
                    .sshn_v
                    .set(ji, jj, simulation_vars.sshn.get(ji, jj));
            }
        }
    }

    kernel_timer.stop();
}
