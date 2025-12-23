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
    let mut sample = matrix_c![
        //
        vector![1, 2, 3],
        vector![4, 5, 6],
    ];

    println!("{sample}");

    sample.gauss_jordan()?;

    println!("{sample}");

    Ok(())
}

#[allow(unused)]
fn vector_test() -> Result<(), VectorError> {
    let a = vector![1, 2, 3];
    let b = vector![0, 1, 1];

    if let Ok(unit) = a.unit() {
        println!("Unit (A): {unit}");
    }

    if let Ok(unit) = b.unit() {
        println!("Unit (B): {unit}");
    }

    if let Ok(angle) = a.angle_between(&b) {
        println!("Angle: {angle}", angle = angle.to_degrees());
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    matrix_test()?;

    Ok(())
}
