mod linalg;
use std::error::Error;

#[allow(unused)]
use crate::linalg::{
    matrix::{Matrix, MatrixError},
    number::Real,
    vector::{Vector, VectorError},
};

#[allow(unused)]
fn matrix_test() -> Result<(), MatrixError> {
    let sample = matrix_c![
        //
        vector![1, 2, 3],
        vector![4, 5, 6],
    ];

    println!("{sample}");

    Ok(())
}

#[allow(unused)]
fn vector_test() -> Result<(), VectorError> {
    let a = vector![1, 2, 3];
    let b = vector![4, 5, 6];

    if let Ok(sum) = &a + &b {
        println!("{sum}");
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    vector_test()?;

    Ok(())
}
