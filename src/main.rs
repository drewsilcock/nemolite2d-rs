mod checksum;

use nalgebra::DMatrix;

#[macro_use]
extern crate serde_derive;

type WorkingPrecision = f64;

#[derive(Debug, Deserialize)]
struct ModelParameters {
    /// Number of columns in model grid
    pub jpi: u32,

    // Number of rows in model grid
    pub jpj: u32,

    // Grid size in x direction (in meters)
    pub dx: WorkingPrecision,

    // Grid size in y direciton (in meters)
    pub dy: WorkingPrecision,

    /// Constant depth of simulation (in meters)
    pub dep_const: WorkingPrecision,

    /// Initial time step
    pub initial_step_index: u32,

    /// Final time step
    pub final_step_index: u32,

    /// Interval to record output
    pub output_interval: u32,

    /// Size of time step (in seconds)
    pub rdt: WorkingPrecision,

    /// Bottom friction coefficient
    pub cbfr: WorkingPrecision,

    /// Horizontal kinematic viscosity coefficient
    pub visc: WorkingPrecision,
}

impl Default for ModelParameters {
    fn default() -> Self {
        ModelParameters {
            jpi: 256,
            jpj: 256,
            dx: 1000.0,
            dy: 1000.0,
            dep_const: 100.0,
            initial_step_index: 1,
            final_step_index: 1000,
            output_interval: 1,
            rdt: 10.0,
            cbfr: 0.001,
            visc: 100.0,
        }
    }
}

impl ModelParameters {
    pub fn new(config_fname: &str) -> Result<Self, config::ConfigError> {
        let mut settings = config::Config::default();

        settings
            .merge(config::File::with_name(config_fname)).unwrap()
            .merge(config::Environment::with_prefix("NEMOLITE2D")).unwrap();

        settings.try_into()
    }
}

fn main() {
    let model_params = ModelParameters::new("Config").unwrap();
    println!("Model params: {:?}", model_params);

    let x = -3.0f64;
    println!("Abs x is: {}", x.abs());

    let my_matrix = DMatrix::from_vec(2, 2, vec![1.0f64, 2.0f64, 3.0f64, 4.0f64]);

    let chksm = checksum::field_checksum(&my_matrix);
    println!("Checksum is: {}", chksm);
}
