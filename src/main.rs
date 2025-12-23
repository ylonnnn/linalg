mod linalg;

use linalg::{matrix::Matrix, vector::Vector};

use crate::linalg::matrix::MatrixError;

fn main() -> Result<(), MatrixError> {
    #[allow(unused)]
    use linalg::matrix::ElementaryOp;

    let mut sample = Matrix::from_rows(vec![
        Vector::from(vec![1, 2, -1, 1]),  // 1
        Vector::from(vec![3, 3, 0, 4]),   // 2
        Vector::from(vec![-1, 0, 2, -3]), // 3
    ]);

    println!("{sample}");

    sample.gauss_jordan()?;

    println!("{sample}");

    Ok(())
}
