use std::{
    collections::HashMap,
    fmt,
    ops::{Add, Mul},
};

use crate::linalg::{number::Real, vector::Vector};

#[derive(Debug, thiserror::Error)]
pub enum MatrixError {
    #[error("dimension mismatch: {curr:?} {other:?}")]
    DimensionMismatch {
        curr: (usize, usize),
        other: (usize, usize),
    },

    #[error("row-column mismatch: {0} : {1}")]
    RowColMismatch(usize, usize),

    #[error("inconsistent system")]
    Inconsistent,

    #[error("non-square matrix")]
    NonSquareMatrix,

    #[error("matrix is singular")]
    Singular,

    #[error("invalid vector dimension: {0:?}")]
    InvalidVectorDimension(usize),
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

    #[allow(non_snake_case)]
    pub fn I<const N: usize>() -> Self {
        Matrix::from_cols((0..N).map(|i| Vector::e::<N>(i + 1)).collect())
    }

    pub fn identity(n: usize) -> Self {
        Matrix::from_cols((0..n).map(|i| Vector::standard(n, i + 1)).collect())
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
        if i == 0 || i > self.entries.len() {
            None
        } else {
            Some(Vector::from_borrowed(
                &self.entries[i - 1][0..((self.column_size() + (self.augmented as usize))
                    - ((self.augmented && E) as usize))]
                    .to_vec(),
            ))
        }
    }

    pub fn rows<const E: bool>(&self) -> Vec<Vector> {
        (1..(self.row_size() + 1))
            .filter_map(|i| self.row::<E>(i))
            .collect()
    }

    pub fn column(&self, j: usize) -> Option<Vector> {
        let n = self.entries[0].len();
        if j == 0 || j > n {
            None
        } else {
            let mut vec = Vec::new();
            vec.reserve(n);

            self.entries.iter().for_each(|row| vec.push(row[j - 1]));

            Some(Vector::from(vec))
        }
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

    pub fn is_singular(&self) -> Result<bool, MatrixError> {
        Ok(self.gaussian_diagonal_det()? == 0_f64)
    }

    pub fn transpose(&self) -> Matrix {
        Matrix::from_rows(self.columns())
    }

    pub fn elementary(&mut self, op: ElementaryOp) {
        match op {
            ElementaryOp::I { i, j } => self.entries.swap(i - 1, j - 1),
            ElementaryOp::II { i, c } => self.entries[i - 1].iter_mut().for_each(|x| *x *= c),
            ElementaryOp::III { i, j, c } => {
                if let Some(row) = self.row::<false>(j) {
                    self.entries[i - 1]
                        .iter_mut()
                        .zip((&row * c).components().iter())
                        .for_each(|(x, y)| *x += y);
                }
            }
        }
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

        println!("{self}");

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
            while pivot == 0_f64 && p_idx < n {
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
                    c: -(below / pivot),
                });
            });

            pivots.push(Pivot { pivot, pos: p_idx });
            p_idx += 1;
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

        println!("{self}");

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

    pub fn invert(&self) -> Result<Matrix, MatrixError> {
        if self.is_singular()? {
            return Err(MatrixError::Singular);
        }

        let n = self.row_size();

        let mut cols = self.columns(); // cols(A | I^n)
        cols.extend(Matrix::identity(n).columns());

        // [A | I^n]
        let mut inv = Matrix::from_cols(cols);
        inv.gauss_jordan()?;

        // Take A^-1 from [I^n | A^-1] after Gauss-Jordan
        Ok(Matrix::from_cols(inv.columns().drain(n..(n * 2)).collect()))
    }

    pub fn cofactor_expansion(&self) -> Result<f64, MatrixError> {
        if !self.is_square() {
            return Err(MatrixError::NonSquareMatrix);
        }

        let n = self.row_size();

        if n == 0 {
            Ok(1_f64)
        } else if n == 1 {
            Ok(self.entries[0][0])
        } else if n == 2 {
            Ok((self.entries[0][0] * self.entries[1][1])
                - (self.entries[0][1] * self.entries[1][0]))
        } else {
            let mut det = 0_f64;
            let top = &self.entries[0];

            for i in 0..n {
                let cofactor = (if i % 2 == 0 { 1 } else { -1 }) as f64 * top[i];

                let mut columns: Vec<Vector> = Vec::new();
                (0..n).for_each(|j| {
                    if i == j {
                        return;
                    }

                    if let Some(col) = self.column(j + 1) {
                        columns.push(Vector::from(col.components().clone().drain(1..n).collect()));
                    }
                });

                det += cofactor * Matrix::from_cols(columns).cofactor_expansion()?;
            }

            Ok(det)
        }
    }

    pub fn gaussian_diagonal_det(&self) -> Result<f64, MatrixError> {
        if !self.is_square() {
            Err(MatrixError::NonSquareMatrix)
        } else {
            let mut clone = self.clone();
            let (_, swaps) = clone.gaussian()?;

            println!("{clone} {swaps}");

            Ok((1..=self.row_size())
                .filter_map(|i| clone.entry(i, i))
                .reduce(|acc, x| acc * x)
                .unwrap_or(1_f64)
                * (-1_i64).pow(swaps as u32) as f64)
        }
    }

    pub fn solve_linear(&self, b: Vector) -> Result<Vec<Vector>, MatrixError> {
        let b_dim = b.dimension();
        if b_dim != self.row_size() {
            return Err(MatrixError::InvalidVectorDimension(b_dim));
        }

        let n = self.column_size();
        let is_homogeneous = b.is_zero();

        let mut augmented = self.augment(b);

        let pivs = augmented.gauss_jordan()?;
        let pivots =
            HashMap::<usize, &Pivot>::from_iter(pivs.iter().map(|pivot| (pivot.pos, pivot)));

        println!("{augmented}");

        let mut solutions = Vec::<Vector>::new();
        solutions.reserve((n - pivots.len()) + (!is_homogeneous) as usize);

        // Particular Solution
        if !is_homogeneous {
            let mut soln = vec![0_f64; n];
            if let Some(xp) = augmented.column(augmented.column_size() + 1) {
                pivs.iter()
                    .zip(xp.components().iter())
                    .for_each(|(piv, xc)| {
                        soln[piv.pos] = *xc;
                    });
            }

            solutions.push(Vector::from(soln));
        }

        (0..n).for_each(|i| {
            if pivots.get(&i).is_some() {
                return;
            }

            let mut soln = vec![0_f64; n];
            if let Some(col) = augmented.column(i + 1) {
                pivs.iter()
                    .zip(col.components().iter())
                    .for_each(|(piv, comp)| soln[piv.pos] = -*comp);
            }

            soln[i] = 1_f64;
            solutions.push(Vector::from(soln));
        });

        Ok(solutions)
    }
}

impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let n = self.entries.len();
        if n == 0 {
            write!(f, "[] {}x{}", n, n)
        } else {
            write!(f, "[ {}x{}\n", self.entries.len(), self.entries[0].len())?;

            for row in &self.entries {
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

impl PartialEq for Matrix {
    fn eq(&self, other: &Self) -> bool {
        self.entries
            .iter()
            .zip(other.entries.iter())
            .all(|(a, b)| a.iter().zip(b.iter()).all(|(x, y)| x == y))
    }

    fn ne(&self, other: &Self) -> bool {
        !(self == other)
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
