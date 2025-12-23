use std::{
    fmt,
    ops::{Add, Mul},
};

use crate::linalg::{number::Real, vector::Vector};

#[derive(Debug, thiserror::Error)]
pub enum MatrixError {
    #[error("Dimension mismatch: {curr:?} {other:?}")]
    DimensionMismatch {
        curr: (usize, usize),
        other: (usize, usize),
    },

    #[error("Row-Column Mismatch: {0} : {1}")]
    RowColMismatch(usize, usize),

    #[error("Inconsistent system")]
    Inconsistent,
}

#[derive(Debug, Clone)]
pub struct Matrix {
    entries: Vec<Vec<f64>>,
    augmented: bool,
}

pub enum ElementaryOp {
    I { i: usize, j: usize },
    II { i: usize, c: f64 },
    III { i: usize, j: usize, c: f64 },
}

pub struct Pivot {
    pub pivot: f64,
    pub pos: usize,
}

impl Matrix {
    pub fn from_cols(columns: Vec<Vector>) -> Self {
        let c = columns.len();
        let r = columns[0].dimension();

        let mut matrix = Matrix {
            entries: vec![vec![0f64; c]; r],
            augmented: false,
        };

        (0..r).for_each(|i| {
            (0..c).for_each(|j| {
                if let Some(comp) = columns[j].component(i + 1) {
                    matrix.entries[i][j] = comp;
                }
            })
        });

        matrix
    }

    pub fn from_rows(rows: Vec<Vector>) -> Self {
        Self {
            entries: rows
                .into_iter()
                .map(|row| row.components().clone())
                .collect(),
            augmented: false,
        }
    }

    pub fn from_size(m: usize, n: usize) -> Self {
        Self {
            entries: vec![vec![0f64; m]; n],
            augmented: false,
        }
    }

    pub fn identity<const N: usize>() -> Self {
        Matrix::from_cols((0..N).map(|i| Vector::e::<N>(i + 1)).collect())
    }

    pub fn row_size(&self) -> usize {
        self.entries.len()
    }

    pub fn column_size(&self) -> usize {
        self.entries[0].len() - (self.augmented as usize)
    }

    /**
     * NOTE: `E = true` will exclude the `b` vector components from the row
     */
    pub fn row<const E: bool>(&self, i: usize) -> Option<Vector> {
        if i <= 0 || i > self.entries.len() {
            return None;
        }

        if self.augmented {}

        Some(Vector::from_borrowed(
            &self.entries[i - 1][0..(self.column_size() - ((self.augmented && E) as usize))]
                .to_vec(),
        ))
    }

    pub fn rows<const E: bool>(&self) -> Vec<Vector> {
        (1..(self.row_size() + 1))
            .filter_map(|i| self.row::<E>(i))
            .collect()
    }

    pub fn column(&self, j: usize) -> Option<Vector> {
        let n = self.entries[0].len();
        if j <= 0 || j > n {
            return None;
        }

        let mut vec = Vec::new();
        vec.reserve(n);

        self.entries.iter().for_each(|row| vec.push(row[j - 1]));

        Some(Vector::from(vec))
    }

    pub fn columns(&self) -> Vec<Vector> {
        (1..(self.column_size() + 1))
            .filter_map(|j| self.column(j))
            .collect()
    }

    pub fn entry(&self, i: usize, j: usize) -> Option<f64> {
        if i <= 0 || i > self.row_size() || j <= 0 || j > self.column_size() {
            return None;
        }

        Some(self.entries[i - 1][j - 1])
    }

    pub fn is_square(&self) -> bool {
        self.row_size() == self.column_size()
    }

    pub fn transpose(&self) -> Matrix {
        Matrix::from_rows(self.columns())
    }

    pub fn elementary(&mut self, op: ElementaryOp) {
        let _ = 'op: {
            match op {
                ElementaryOp::I { i, j } => self.entries.swap(i - 1, j - 1),
                ElementaryOp::II { i, c } => self.entries[i - 1].iter_mut().for_each(|x| *x *= c),
                ElementaryOp::III { i, j, c } => {
                    let row = self.row::<false>(j);
                    if row.is_none() {
                        break 'op;
                    }

                    self.entries[i - 1]
                        .iter_mut()
                        .zip((&row.unwrap() * c).components().iter())
                        .for_each(|(x, y)| *x += y);
                }
            }
        };
    }

    pub fn augment(&self, b: Vector) -> Matrix {
        let mut cols = self.columns();
        cols.push(b);

        let mut matrix = Matrix::from_cols(cols);
        matrix.augmented = true;

        matrix
    }

    /**
     * NOTE: Modifies the matrix in place
     */
    pub fn gaussian(&mut self) -> Result<(Vec<Pivot>, usize), MatrixError> {
        let (mut pivots, mut swaps) = (Vec::<Pivot>::new(), 0usize);

        let mut p_idx = 0;
        let (m, n) = (self.row_size(), self.column_size());

        (0..m).for_each(|i| {
            let mut max_row = i;

            // Partial Pivoting (also aids in shifting zero rows to the bottom)
            ((i + 1)..m).for_each(|j| {
                if self.entries[j][p_idx].abs() > self.entries[max_row][p_idx].abs() {
                    max_row = j;
                }
            });

            if i != max_row {
                self.elementary(ElementaryOp::I {
                    i: i + 1,
                    j: max_row + 1,
                });
                swaps += 1
            }

            if p_idx >= n {
                return;
            }

            // Proper pivot searching as pivots are not directly
            // next to each other in terms of columns and rows
            let mut pivot = self.entries[i][p_idx];
            while pivot == 0f64 && p_idx < n {
                p_idx += 1;
                pivot = if p_idx < n {
                    self.entries[i][p_idx]
                } else {
                    0_f64
                };
            }

            if p_idx >= n {
                return;
            }

            // Eliminate non-zero entries below the current pivot
            ((i + 1)..m).for_each(|j| {
                let below = self.entries[j][p_idx];
                if below == 0f64 {
                    return;
                }

                self.elementary(ElementaryOp::III {
                    i: j + 1,
                    j: i + 1,
                    c: (1f64 / pivot) * -below,
                });
            });

            pivots.push(Pivot { pivot, pos: p_idx });
            p_idx += 1
        });

        if self.augmented {
            let n = self.column_size();

            for i in 0..self.row_size() {
                let row = self.row::<true>(i + 1).unwrap();
                if row.is_zero() && self.entries[i][n - 1] != 0f64 {
                    return Err(MatrixError::Inconsistent);
                }
            }
        }

        Ok((pivots, swaps))
    }

    pub fn gauss_jordan(&mut self) -> Result<Vec<Pivot>, MatrixError> {
        let (mut pivots, _) = self.gaussian()?;

        (0..pivots.len()).rev().for_each(|i| {
            let Pivot { pivot: piv, pos } = &mut pivots[i];

            self.elementary(ElementaryOp::II {
                i: i + 1,
                c: 1f64 / *piv,
            });

            *piv = self.entries[i][*pos];

            (0..(i)).rev().for_each(|j| {
                let above = self.entries[j][*pos];
                if above == 0f64 {
                    return;
                }

                self.elementary(ElementaryOp::III {
                    i: j + 1,
                    j: i + 1,
                    c: -above,
                });
            });
        });

        Ok(pivots)
    }
}

impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let n = self.entries.len();
        if n == 0 {
            write!(f, "[] {}x{}", n, n)
        } else {
            write!(f, "[ {}x{}\n", self.entries.len(), self.entries[0].len())?;

            for row in self.entries.iter() {
                for entry in row {
                    write!(f, "\t{},", entry.fraction())?;
                }

                write!(f, "\n")?;
            }

            write!(f, "\n]")
        }
    }
}

impl<'a, 'b> Add<&'b Matrix> for &'a Matrix {
    type Output = Result<Matrix, MatrixError>;

    fn add(self, rhs: &'b Matrix) -> Self::Output {
        let (sr, rr, sc, rc) = (
            self.row_size(),
            rhs.row_size(),
            self.column_size(),
            rhs.column_size(),
        );

        if sr != rr || sc != rc {
            return Err(MatrixError::DimensionMismatch {
                curr: (sr, sc),
                other: (rr, rc),
            });
        }

        let mut matrix = self.clone();
        (1..sr).for_each(|i| {
            (1..sc).for_each(|j| {
                if let Some(entry) = rhs.entry(i, j) {
                    matrix.entries[i][j] += entry;
                }
            })
        });

        Ok(matrix)
    }
}

impl<'a, T> Mul<T> for &'a Matrix
where
    T: Into<f64> + Copy,
{
    type Output = Matrix;

    fn mul(self, rhs: T) -> Self::Output {
        let mut matrix = self.clone();
        (1..self.row_size()).for_each(|i| {
            (1..self.column_size()).for_each(|j| {
                matrix.entries[i][j] *= rhs.into();
            })
        });

        matrix
    }
}

impl<'a, 'b> Mul<&'b Matrix> for &'a Matrix {
    type Output = Result<Matrix, MatrixError>;

    fn mul(self, rhs: &'b Matrix) -> Self::Output {
        let (sc, rr) = (self.column_size(), rhs.row_size());
        if sc != rr {
            return Err(MatrixError::RowColMismatch(sc, rr));
        }

        // (m x p) * (p x n) = (m x n)
        let (m, n) = (self.row_size(), rhs.column_size());
        let mut matrix = Matrix::from_size(m, n);

        (1..sc).for_each(|i| {
            (1..sc).for_each(|j| {
                if let Ok(prod) = &(self.row::<false>(i).unwrap()) * &(rhs.column(j).unwrap()) {
                    matrix.entries[i - 1][j - 1] = prod;
                }
            })
        });

        Ok(matrix)
    }
}

/// Constructs a matrix with the provided vectors as rows
#[macro_export]
macro_rules! matrix {
    [$($vec:expr),* $(,)?] => {
        Matrix::from_rows(vec![$($vec),*])
    };
}

/// Constructs a matrix with the provided vectors as columns
#[macro_export]
macro_rules! matrix_c {
    [$($vec:expr),* $(,)?] => {
        Matrix::from_cols(vec![$($vec),*])
    };
}
