#[cfg(test)]
mod tests{
    use streamtracer::field::VectorField;
    use ndarray::{array, Array};

    #[test]
    fn test_field() {
        let xcoords = array![0.1, 0.2, 0.3];
        let ycoords = array![1.0, 1.1, 1.2, 1.3];
        let zcoords = array![-0.1, 0.0, 1.0, 50.0, 56.0];
        let field = Array::zeros((3, 4, 5, 3));
        let cyclic = array![true, false, true];

        let f = VectorField::new(&xcoords, &ycoords, &zcoords, &field, &cyclic);

        // Test grid_idx
        let x = array![0.15, 1.05, -0.05];
        assert_eq!(f.grid_idx(&x), array![0, 0, 0]);

        let x = array![0.25, 1.05, -0.05];
        assert_eq!(f.grid_idx(&x), array![1, 0, 0]);

        let x = array![0.25, 1.15, 53.0];
        assert_eq!(f.grid_idx(&x), array![1, 1, 3]);
    }
}
