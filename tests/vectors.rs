#[cfg(test)]
mod tests {
    use linalg::*;

    #[test]
    fn vector_new() {
        let vec = Vector::new(3); // 3 is the dimension

        // Expect that a new vector will always be a
        // zero vector unless components are specified
        assert_eq!(vec, Vector::new(3));
    }

    #[test]
    fn vector_macro() {
        let a = vector![1, 2, 3]; // For convenience
        let hardcoded = Vector::from(vec![1, 2, 3]);

        assert_eq!(a, hardcoded);
    }

    #[test]
    fn vector_from_borrowed() {
        let a = vector![1, 2, 3];
        let b = Vector::from_borrowed(a.components());

        assert_eq!(a, b); // Equal and a is not moved
    }

    #[test]
    fn vector_new_is_always_zero() {
        let vec = Vector::new(3);

        // Dimension of the zero vector is provided as a constant
        // generic parameter
        let zero = Vector::zero::<3>();

        assert_eq!(vec, zero);
    }

    #[test]
    fn vector_standard_basis_vectors() {
        // Dimension - constant generic parameter
        // Index - Index of the nonzero (1) scalar or the index of the vector
        let e1 = Vector::e::<3>(1);

        // Both parameters if dimension is not constant
        let e2 = Vector::standard(3, 2);

        assert_eq!(e1, Vector::standard(3, 1));
        assert_eq!(e2, Vector::e::<3>(2));
    }

    #[test]
    fn vector_mag_norm() {
        let a = vector![1, 2, 3];
        let (magnitude, norm) = (a.magnitude(), a.norm());

        // Magnitude and norm will always be equal
        assert_eq!(magnitude, norm);
    }

    #[test]
    fn vector_unit_vectors() {
        let sample = vector![4, 1, 3];

        // Vector::unit() returns a Result as not all vectors can be unit vectors
        // Specifically, zero vectors
        let sample_unit = sample.unit();

        if let Ok(unit) = sample_unit {
            dbg!(&unit);

            assert_eq!(unit.is_unit(), true);

            // Constant (EPSILON) proximity check as floating-point values are represented in
            // binary in computers which may lead to floating-point errors such as being under or
            // over the expected value
            let norm = unit.norm();
            assert!((1_f64 - norm).abs() < f64::EPSILON);
        }

        // In-place creation/convertion of vectors to unit vectors
        // NOTE: This will require the vector to be mutable
        let mut sample_copy = sample.clone();
        // Normalization Result simply holds an error and a
        // unit (`()`) value
        let norm_res = sample_copy.normalize();

        if norm_res.is_ok() {
            // The sample_copy that is normalized can now be
            // expected to be a unit vector
            assert!(sample_copy.is_unit());
        }
    }

    #[test]
    fn vector_angle_between_two_vectors() {
        let a = vector![1, 0];
        let b = vector![0, 1];

        // May result to an error as angle calculation uses the inverse of the cosine ofthe angle
        // between two vectors to calculate the angle. And with that, vectors may be vectors that
        // cannot be normalized (zero vectors, etc)
        match a.angle_between(&b) {
            Ok(angle) => {
                dbg!("Angle (Radians): {}", angle);
                dbg!("Angle (Degree): {}", angle.to_degrees());
            }

            Err(err) => {
                dbg!(&err);
            }
        }
    }

    #[test]
    fn vector_orthogonality_two_vectors() {
        let a = vector![3, 0, 0];
        let b = vector![0, 4, 0];

        assert!(a.is_orthogonal(&b));
    }

    #[test]
    fn vector_orthonormality_two_vectors() {
        let a = vector![1, 0, 0];
        let b = vector![0, 1, 0];

        assert!(a.is_orthogonal(&b));
    }

    #[test]
    fn vector_orthogonality_vectors_in_set() {
        let set = vec![
            vector![2, -1, 4],
            vector![-3, 5, 0],
            vector![1, 2, -2],
            vector![4, 1, 3],
        ];

        assert_eq!(Vector::is_orthogonal_set(&set), false);

        let orthogonal_set: Vec<Vector> = (1..=3).map(|i| &Vector::e::<3>(i) * 3).collect();
        assert!(Vector::is_orthogonal_set(&orthogonal_set));
    }

    #[test]
    fn vector_orthonormality_vectors_in_set() {
        let set = vec![
            vector![2, -1, 4],
            vector![-3, 5, 0],
            vector![1, 2, -2],
            vector![4, 1, 3],
        ];

        assert_eq!(Vector::is_orthonormal_set(&set), false);

        let orthogonal_set: Vec<Vector> = (1..=3).map(|i| &Vector::e::<3>(i) * 3).collect();
        // Orthogonal, but not orthonormal
        assert_eq!(Vector::is_orthonormal_set(&orthogonal_set), false);

        let orthonormal_set: Vec<Vector> = orthogonal_set
            .iter()
            .filter_map(|vec| match vec.unit() {
                Ok(vec) => Some(vec),
                Err(_) => None,
            })
            .collect();

        assert!(Vector::is_orthonormal_set(orthonormal_set.as_slice()));
    }
}
