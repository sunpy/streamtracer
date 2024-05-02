#[cfg(test)]
mod field_tests {
    use float_eq::assert_float_eq;
    use numpy::ndarray::{array, s, Array, Array4};

    use super::super::field::VectorField;

    #[test]
    fn test_grid_idx() {
        let xgrid = array![0., 0.2, 0.3];
        let ygrid = array![0., 1.1, 1.2, 1.3];
        let zgrid = array![0.0, 1.0, 50.0, 56.0, 100.0];
        let mut field: Array4<f64> = Array::zeros((3, 4, 5, 3));
        field.slice_mut(s![2, 2.., 2.., 0]).fill(1.);
        let cyclic = array![true, false, true];
        let f = VectorField::new(
            xgrid.view(),
            ygrid.view(),
            zgrid.view(),
            field.view(),
            cyclic.view(),
        );

        let x = array![0.15, 1.05, 0.05];
        assert_eq!(f.grid_idx(x.view()), array![0, 0, 0]);

        let x = array![0.25, 1.05, 0.05];
        assert_eq!(f.grid_idx(x.view()), array![1, 0, 0]);

        let x = array![0.25, 1.15, 53.0];
        assert_eq!(f.grid_idx(x.view()), array![1, 1, 2]);
    }

    #[test]
    fn test_vector_at_position() {
        let xgrid = array![0., 0.2, 0.3];
        let ygrid = array![0., 1.1, 1.2, 1.3];
        let zgrid = array![0.0, 1.0, 50.0, 56.0, 100.0];
        let mut field: Array4<f64> = Array::zeros((3, 4, 5, 3));
        field.slice_mut(s![2, 2.., 2.., 0]).fill(1.);
        let cyclic = array![true, false, true];
        let f = VectorField::new(
            xgrid.view(),
            ygrid.view(),
            zgrid.view(),
            field.view(),
            cyclic.view(),
        );

        let mut x = array![0.15, 1.05, 0.05];

        assert_eq!(f.vector_at_position(x.view()), array![0., 0., 0.]);

        // Exactly on a grid boundary in {y, z}, halfway along in {x}
        x = array![0.25, 1.2, 50.0];
        assert_eq!(f.vector_at_position(x.view()), array![0.5, 0., 0.]);

        // Exactly on a grid boundary in {y, z}, 3/4 along in {x}
        x = array![0.275, 1.2, 50.0];
        let vec = f.vector_at_position(x.view());
        let expected = array![0.75, 0., 0.];
        for i in 0..2 {
            assert_float_eq!(vec[i], expected[i], abs <= 0.000_000_1);
        }
    }
}
