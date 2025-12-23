#[cfg(test)]
mod tests {
    use linalg::*;

    #[test]
    fn matrix_det_cofactor_expansion_nonsquare() {
        let non_square = matrix![vector![1, 2]];

        assert!(matches!(
            non_square.cofactor_expansion(),
            Err(MatrixError::NonSquareMatrix)
        ));
    }

    #[test]
    fn matrix_det_cofactor_expansion_valid() {
        let valid = matrix![
            // ad - bc
            vector![1, 1],
            vector![2, 0]
        ];

        assert!(matches!(valid.cofactor_expansion(), Ok(-2_f64)));
    }
}
