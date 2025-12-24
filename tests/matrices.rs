#[cfg(test)]
mod tests {
    use linalg::*;

    #[test]
    fn matrix_from_cols() {
        // [1, -1, 4]
        // [2, -2, 5]
        // [3, -3, 6]
        let a = Matrix::from_cols(vec![
            Vector::from(vec![1, 2, 3]),
            Vector::from(vec![-1, -2, -3]),
            Vector::from(vec![4, 5, 6]),
        ]);

        println!("{a}");
    }

    #[test]
    fn matrix_from_rows() {
        // [ 1,  2,  3]
        // [-1, -2, -3]
        // [ 4,  5,  6]
        let b = Matrix::from_rows(vec![
            Vector::from(vec![1, 2, 3]),
            Vector::from(vec![-1, -2, -3]),
            Vector::from(vec![4, 5, 6]),
        ]);

        println!("{b}");
    }

    #[test]
    fn matrix_from_size() {
        // Matrix::from_size() creates a zero matrix
        let mat = Matrix::from_size(3, 3);

        // Every column and row are zero vectors
        assert!(mat.columns().iter().all(|col| col.is_zero()));
        assert!(mat.rows::<false>().iter().all(|row| row.is_zero()));
    }

    #[test]
    fn matrix_identity() {
        // Creates a 3x3 Identity Matrix I^n (I^3)
        #[allow(non_snake_case)]
        let I = Matrix::I::<3>();

        println!("{I}");

        // To clarify, each row vector must be orthonormal of/from each other.
        // The same goes for column vectors
        assert!(Vector::is_orthogonal_set(&I.rows::<false>()));
        assert!(Vector::is_orthogonal_set(&I.columns()));
    }

    #[test]
    fn matrix_macro() {
        let x = matrix![vector![1, 0, 1], vector![0, 1, 0], vector![1, 0, 1],];

        println!("{x}");
    }

    #[test]
    fn matrix_square_check() {
        let nonsquare = matrix![
            //
            vector![1, 2, 3],
            vector![4, 5, 6]
        ];

        assert_eq!(nonsquare.is_square(), false);

        let square = matrix![
            //
            vector![1, 2],
            vector![2, 1]
        ];

        assert!(square.is_square());
    }

    #[test]
    fn matrix_transpose() {
        let a = matrix![
            //
            vector![6, 5, 4],
            vector![9, 8, 7],
        ];

        let a_t = a.transpose();

        assert_ne!(a, a_t);
        assert_eq!(a, a_t.transpose());
    }

    #[test]
    fn matrix_elementary_ops() {
        let mut x = matrix![
            //
            vector![1, 2, 3],
            vector![4, 5, 6]
        ];

        // NOTE: .unwrap() is unsafe unless logically possible
        let row_2 = x.row::<false>(2).unwrap();

        // Type I - Row Swapping
        x.elementary(matrix::ElementaryOp::I { i: 1, j: 2 });

        assert_eq!(x.row::<false>(1).unwrap(), row_2);

        // Type II - Row Scaling
        x.elementary(matrix::ElementaryOp::II { i: 1, c: 5_f64 });

        assert_eq!(x.row::<false>(1).unwrap(), &row_2 * 5_f64);

        let row_2 = x.row::<false>(2).unwrap(); // Shadow row_2

        // Type III - Row Product Addition
        x.elementary(matrix::ElementaryOp::III {
            i: 2,
            j: 1,
            c: 3_f64,
        });

        assert_ne!(x.row::<false>(2).unwrap(), row_2);
    }

    #[test]
    fn matrix_augment() {
        #[allow(non_snake_case)]
        let A = matrix![
            //
            vector![3, 2, 4],
            vector![9, 4, 2]
        ];
        let b = Vector::zero::<2>();

        let augmented = A.augment(b.clone());
        let b_vec = augmented.column(augmented.column_size() + 1).unwrap();

        assert_eq!(b, b_vec);
    }

    #[test]
    fn matrix_det_cofactor_expansion_nonsquare() {
        let non_square = matrix![vector![1, 2]];

        assert!(matches!(
            non_square.cofactor_expansion(),
            Err(MatrixError::NonSquareMatrix)
        ));

        let valid = matrix![
            // ad - bc
            vector![1, 1],
            vector![2, 0]
        ];

        assert!(matches!(valid.cofactor_expansion(), Ok(-2_f64)));
    }

    #[test]
    fn matrix_gaussian_ref_diagonal_det() {
        let sample = matrix![
            //
            vector![1, 1, 2],
            vector![0, -1, 3],
            vector![3, -2, 1]
        ];

        if let Ok(det) = sample.gaussian_diagonal_det() {
            dbg!(&det, sample.cofactor_expansion().unwrap());
            assert_eq!(det, sample.cofactor_expansion().unwrap());
        }
    }

    #[test]
    fn matrix_solve_linear_systems() {
        let sample = matrix_c![
            vector![1, 2, 3],
            vector![1, 2, 3],
            vector![0, 1, 0],
            vector![0, 1, 1],
            vector![1, 0, 2],
        ];

        // Expected RREF
        // [1, 1, 0, 0,  1] =     [-1] [-1]
        // [0, 0, 1, 0, -1] =  =  [ 0] [ 1]
        // [0, 0, 0, 1,  1] =     [ 0] [-1]

        // x1 = [-x2 - x5]
        // x2 = [ x2 +  0]
        // x3 = [  0 + x5]
        // x4 = [  0 - x5]
        // x5 = [  0 + x5]

        let sample_b = vector![8, 5, 2];
        if let Ok(solutions) = sample.solve_linear(sample_b) {
            solutions.iter().for_each(|sol| {
                dbg!(&sol);
            });

            assert!(false);
        }
    }
}
