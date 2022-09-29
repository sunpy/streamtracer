#[cfg(test)]
mod tests{
    use ndarray::{array, Array, s};
    use float_eq::assert_float_eq;

    use streamtracer::field::VectorField;


    /// Get a simple vector field
    /// This implicitly tests `new()`.
    fn get_field() -> VectorField {
        let xcoords = array![0.1, 0.2, 0.3];
        let ycoords = array![1.0, 1.1, 1.2, 1.3];
        let zcoords = array![-0.1, 0.0, 1.0, 50.0, 56.0];
        let mut field = Array::zeros((3, 4, 5, 3));

        field.slice_mut(s![2, 2.., 2.., 0]).fill(1.);
        let cyclic = array![true, false, true];

        return VectorField::new(xcoords, ycoords, zcoords, field, cyclic);
    }

    #[test]
    fn test_grid_idx() {
        let f = get_field();

        let x = array![0.15, 1.05, -0.05];
        assert_eq!(f.grid_idx(x.view()), array![0, 0, 0]);

        let x = array![0.25, 1.05, -0.05];
        assert_eq!(f.grid_idx(x.view()), array![1, 0, 0]);

        let x = array![0.25, 1.15, 53.0];
        assert_eq!(f.grid_idx(x.view()), array![1, 1, 3]);
    }

    #[test]
    fn test_vector_at_position() {
        let f = get_field();

        let mut x = array![0.15, 1.05, -0.05];

        assert_eq!(f.vector_at_position(x.view()), array![0., 0., 0.]);

        // Exactly on a grid boundary in {y, z}, halfway along in {x}
        x = array![0.25, 1.2, 50.0];
        assert_eq!(f.vector_at_position(x.view()), array![0.5, 0., 0.]);

        // Exactly on a grid boundary in {y, z}, 3/4 along in {x}
        x = array![0.275, 1.2, 50.0];
        let vec = f.vector_at_position(x.view());
        let expected = array![0.75, 0., 0.];
        for i in 0..2{
            assert_float_eq!(vec[i], expected[i], abs <= 0.000_000_1);
        }
    }
}
