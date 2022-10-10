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

use numpy::{PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray4, PyArray3, IntoPyArray};
use pyo3::prelude::{pymodule, PyModule, PyResult, Python};


#[pymodule]
#[pyo3(name = "_streamtracer_rust")]
fn streamtracer(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    #[pyfn(m)]
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
    ) -> &'py PyArray3<f64> {
        let (_statuses, xs) = trace::trace_streamlines(
            seeds.as_array(),
            xgrid.as_array(),
            ygrid.as_array(),
            zgrid.as_array(),
            values.as_array(),
            cyclic.as_array(),
            direction,
            step_size,
            max_steps
        );

        xs.into_pyarray(py)
    }

    Ok(())
}
