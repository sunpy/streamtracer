#[cfg(test)]
mod tests{
    use streamtracer::field::VectorField;
    use ndarray::{array, Array};

    #[test]
    fn test_field() {
        let xcoords = array![0.1, 0.2, 0.3];
        let ycoords = array![1.0, 1.1, 1.2, 1.3];
        let zcoords = array![-0.1, -0.2];
        let field = Array::zeros((3, 4, 2, 3));
        let cyclic = array![true, false, true];

        let _f = VectorField::new(&xcoords, &ycoords, &zcoords, &field, &cyclic);
    }
}
