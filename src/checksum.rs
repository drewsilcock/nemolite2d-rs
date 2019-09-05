use crate::fortran_array_2d::FortranArray2D;
use crate::WorkingPrecision;

pub fn field_checksum(field: &FortranArray2D<WorkingPrecision>) -> WorkingPrecision {
    field.iter().map(|value| value.abs()).sum()
}
