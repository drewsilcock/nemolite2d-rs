use std::any::Any;
use std::fmt::Debug;

pub trait Scalar: Copy + PartialEq + Debug + Any {}
impl<T: Copy + PartialEq + Debug + Any> Scalar for T {}

pub struct FortranArray2D<T: Scalar> {
    values: Vec<T>,

    start_row_idx: usize,
    start_column_idx: usize,

    num_rows: usize,
}

impl<T: Scalar> FortranArray2D<T>
where
    T: Default,
{
    pub fn new(
        start_row_idx: usize,
        start_column_idx: usize,
        end_row_idx: usize,
        end_column_idx: usize,
    ) -> Self {
        let num_rows = end_row_idx - start_row_idx + 1;
        let num_columns = end_column_idx - start_column_idx + 1;

        FortranArray2D {
            start_row_idx,
            start_column_idx,

            num_rows,
            values: vec![T::default(); num_rows * num_columns],
        }
    }

    #[inline(always)]
    pub fn get(&self, row_idx: usize, column_idx: usize) -> T {
        let idx = self.index_from_row_and_column(row_idx, column_idx);
        self.values[idx]
    }

    #[inline(always)]
    pub fn set(&mut self, row_idx: usize, column_idx: usize, value: T) {
        let idx = self.index_from_row_and_column(row_idx, column_idx);
        self.values[idx] = value
    }

    pub fn set_all(&mut self, value: T) {
        self.values.iter_mut().map(|x| *x = value).count();
    }

    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.values.iter()
    }

    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, T> {
        self.values.iter_mut()
    }

    #[inline(always)]
    fn index_from_row_and_column(&self, row_idx: usize, column_idx: usize) -> usize {
        (row_idx - self.start_row_idx) + (column_idx - self.start_column_idx) * (self.num_rows)
    }
}
