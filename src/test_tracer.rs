#[cfg(test)]
mod tests {
    use numpy::ndarray::{array, s, Array};

    use super::super::field::VectorField;
    use super::super::trace::{trace_streamline, TracerStatus};

    #[test]
    fn test_uniform_field() {
        let xgrid = Array::range(0., 10.1, 0.5);
        let ygrid = Array::range(0., 10.1, 0.5);
        let zgrid = Array::range(0., 10.1, 0.5);
        // Create a vector field pointing in the x direction
        let mut field = Array::zeros((xgrid.len(), ygrid.len(), zgrid.len(), 3));
        field.slice_mut(s![.., .., .., 0]).fill(1.);

        let cyclic = array![false, false, false];
        let f = VectorField::new(
            xgrid.view(),
            ygrid.view(),
            zgrid.view(),
            field.view(),
            cyclic.view(),
        );

        let seed = array![5., 5., 5.];
        let direction = 1;
        let step_size = 0.1;

        // Should take 50 steps before going out of bounds
        let max_steps = 100;
        let result = trace_streamline(seed.view(), &f, &direction, &step_size, max_steps);
        assert_eq![result.status.n_points, 51];
        assert_eq![result.status.rot, TracerStatus::OutOfBounds];

        let max_steps = 10;
        let result = trace_streamline(seed.view(), &f, &direction, &step_size, max_steps);
        assert_eq![result.status.n_points, max_steps];
        assert_eq![result.status.rot, TracerStatus::RanOutOfSteps];
    }
}
