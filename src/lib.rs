//! This crate contains code for tracing streamlines through a 3D vector field defined
//! on rectilinear grids.
#![warn(missing_docs)]
pub mod field;
pub mod interp;
pub mod trace;

#[cfg(test)]
mod test_field;
mod test_interp;
mod test_tracer;

use numpy::{
    ndarray::Array, IntoPyArray, PyArray1, PyArray3, PyReadonlyArray1, PyReadonlyArray2,
    PyReadonlyArray4,
};
use pyo3::prelude::{pymodule, Bound, PyModule, PyResult, Python};

#[pymodule]
#[pyo3(name = "_streamtracer_rust")]
fn streamtracer(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::type_complexity)]
    fn trace_streamlines<'py>(
        py: Python<'py>,
        seeds: PyReadonlyArray2<f64>,
        xgrid: PyReadonlyArray1<f64>,
        ygrid: PyReadonlyArray1<f64>,
        zgrid: PyReadonlyArray1<f64>,
        values: PyReadonlyArray4<f64>,
        cyclic: PyReadonlyArray1<bool>,
        direction: i32,
        step_size: f64,
        max_steps: usize,
    ) -> (
        Bound<'py, PyArray3<f64>>,
        Bound<'py, PyArray1<i64>>,
        Bound<'py, PyArray1<i64>>,
    ) {
        let (statuses, xs) = trace::trace_streamlines(
            seeds.as_array(),
            xgrid.as_array(),
            ygrid.as_array(),
            zgrid.as_array(),
            values.as_array(),
            cyclic.as_array(),
            direction,
            step_size,
            max_steps,
        );

        let mut termination_reasons = Array::zeros(statuses.len());
        let mut n_points = Array::zeros(statuses.len());
        for (i, status) in statuses.iter().enumerate() {
            termination_reasons[[i]] = status.rot as i64;
            n_points[[i]] = status.n_points as i64;
        }

        return (
            xs.into_pyarray_bound(py),
            n_points.into_pyarray_bound(py),
            termination_reasons.into_pyarray_bound(py),
        );
    }

    return Ok(());
}
