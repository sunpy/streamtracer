#[cfg(test)]
mod tests{
    use numpy::ndarray::{array, Array, Array2, s};

    use super::super::field::VectorField;
    use super::super::trace::{TracerStatus, trace_streamline};

    #[test]
    fn test_uniform_field() {
        let xgrid = Array::range(0., 10.1, 0.5);
        let ygrid = Array::range(0., 10.1, 0.5);
        let zgrid = Array::range(0., 10.1, 0.5);
        // Create a vector field pointing in the x direction
        let mut field = Array::zeros((xgrid.len(), ygrid.len(), zgrid.len(), 3));
        field.slice_mut(s![.., .., .., 0]).fill(1.);

        let cyclic = array![false, false, false];
        let f = VectorField::new(xgrid.view(), ygrid.view(), zgrid.view(), field.view(), cyclic.view());

        let seed = array![5., 5., 5.];
        let direction = 1;
        let step_size = 0.1;

        // Should take 50 steps before going out of bounds
        let mut xs: Array2<f64> = Array::zeros((100, 3));
        let status = trace_streamline(seed.view(), &f, &direction, &step_size, xs.view_mut());
        assert_eq![status.n_points, 51];
        assert_eq![status.rot, TracerStatus::OutOfBounds];

        let mut xs: Array2<f64> = Array::zeros((10, 3));
        let status = trace_streamline(seed.view(), &f, &direction, &step_size, xs.view_mut());
        assert_eq![status.n_points, 10];
        assert_eq![status.rot, TracerStatus::RanOutOfSteps];
    }
}
