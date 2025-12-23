use std::{
    fmt,
    ops::{Add, Div, Mul},
};

use crate::linalg::number::Real;

#[derive(Debug, thiserror::Error)]
pub enum VectorError {
    #[error("Dimension mismatch")]
    DimensionMismatch(usize, usize),

    #[error("Zero vector")]
    ZeroVector,
}

#[derive(Debug, Clone)]
pub struct Vector {
    components: Vec<f64>,
}

impl Vector {
    pub fn new(size: usize) -> Self {
        let mut inst = Self {
            components: Vec::new(),
        };

        inst.components.resize(size, 0 as f64);
        inst
    }

    pub fn from_borrowed(vec: &Vec<f64>) -> Self {
        Self {
            components: vec.clone(),
        }
    }

    pub fn from<T>(vec: Vec<T>) -> Self
    where
        T: Into<f64> + Copy,
    {
        Self {
            components: vec.into_iter().map(|x| x.into()).collect(),
        }
    }

    pub fn zero<const N: usize>() -> Self {
        Vector::new(N)
    }

    pub fn e<const N: usize>(i: usize) -> Self {
        let mut inst = Vector::new(N);
        inst.components[i - 1] = 1 as f64;

        inst
    }

    pub fn is_orthogonal_set(set: &Vec<Vector>) -> bool {
        let n = set.len();

        for i in 0..n {
            let current = &set[i];

            for j in (i + 1)..n {
                let next = &set[j];
                if !(current.is_orthogonal(next)) {
                    return false;
                }
            }
        }

        true
    }

    pub fn is_orthonormal_set(set: &Vec<Vector>) -> bool {
        let n = set.len();

        for i in 0..n {
            let current = &set[i];

            for j in (i + 1)..n {
                let next = &set[j];
                if !(current.is_orthonormal(next)) {
                    return false;
                }
            }
        }

        true
    }

    pub fn gram_schmidt<const NM: bool>(set: &Vec<Vector>) -> Vec<Vector> {
        let n = set.len();

        let mut vec = Vec::<Vector>::new();
        vec.reserve(n);

        todo!();
    }

    pub fn component(&self, i: usize) -> Option<f64> {
        if i <= self.components.len() {
            return Some(self.components[i - 1]);
        }

        None
    }

    pub fn components(&self) -> &Vec<f64> {
        &self.components
    }

    pub fn is_zero(&self) -> bool {
        self.components.iter().all(|x| *x == 0f64)
    }

    pub fn dimension(&self) -> usize {
        self.components.len()
    }

    pub fn magnitude(&self) -> f64 {
        self.components.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    pub fn norm(&self) -> f64 {
        self.magnitude()
    }

    pub fn is_unit(&self) -> bool {
        self.norm() == 1f64
    }

    pub fn unit(&self) -> Result<Vector, VectorError> {
        if self.is_zero() {
            return Err(VectorError::ZeroVector);
        }

        Ok(self / self.norm())
    }

    pub fn normalize(&mut self) -> Result<(), VectorError> {
        if self.is_zero() {
            return Err(VectorError::ZeroVector);
        }

        *self = self as &Self / self.norm();
        Ok(())
    }

    pub fn angle_between(&self, other: &Self) -> Result<f64, VectorError> {
        let prod = &self.unit()? * &other.unit()?;
        Ok(prod?.clamp(-1f64, 1f64).acos())
    }

    pub fn is_orthogonal(&self, other: &Self) -> bool {
        match self * other {
            Ok(prod) => prod.abs() < f64::EPSILON,
            Err(_) => false,
        }
    }

    pub fn is_orthonormal(&self, other: &Self) -> bool {
        self.is_orthogonal(other) && self.is_unit()
    }
}

impl fmt::Display for Vector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[\n")?;

        for component in self.components.iter() {
            write!(f, "\t{},\n", component.fraction())?;
            // write!(f, "\t{},\n", component)?;
        }

        write!(f, "\n]")
    }
}

impl<'a, 'b> Add<&'b Vector> for &'a Vector {
    type Output = Result<Vector, VectorError>;

    fn add(self, rhs: &'b Vector) -> Self::Output {
        let mut vec = self.clone();

        let n = vec.dimension();
        let m = rhs.dimension();

        if n != m {
            return Err(VectorError::DimensionMismatch(n, m));
        }

        vec.components
            .iter_mut()
            .zip(rhs.components.iter())
            .for_each(|(a, b)| *a += *b);

        Ok(vec)
    }
}

// Vector * Vector (Dot Product)
impl<'a, 'b> Mul<&'b Vector> for &'a Vector {
    type Output = Result<f64, VectorError>;

    fn mul(self, rhs: &'b Vector) -> Self::Output {
        let n = self.dimension();
        let m = rhs.dimension();

        if n != m {
            return Err(VectorError::DimensionMismatch(n, m));
        }

        Ok(self
            .components
            .iter()
            .zip(rhs.components.iter())
            .map(|(a, b)| a * b)
            .sum())
    }
}

// Scalar * Vector
impl<'a, T> Mul<T> for &'a Vector
where
    T: Into<f64> + Copy,
{
    type Output = Vector;

    fn mul(self, rhs: T) -> Self::Output {
        let mut vec = self.clone();
        vec.components.iter_mut().for_each(|a| *a *= rhs.into());

        vec
    }
}

// Vector / Scalar
impl<'a, T> Div<T> for &'a Vector
where
    T: Into<f64> + Copy,
{
    type Output = Vector;

    fn div(self, rhs: T) -> Self::Output {
        self * (1.0 / rhs.into())
    }
}

#[macro_export]
macro_rules! vector {
    [$($component:expr),* $(,)?] => {
        Vector::from(vec![$($component),*])
    };
}
