use nalgebra::DMatrix;

use crate::WorkingPrecision;

pub fn field_checksum(field: &DMatrix<WorkingPrecision>) -> WorkingPrecision {
    field.iter().map(|value| value.abs()).sum()
}
