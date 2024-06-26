//! Helper functions for interpolation.

use numpy::ndarray::{Array1, ArrayBase, Data, Ix3};

/// Trilinear-interpolation of a scalar defined on
/// the eight corners of a cuboid.
///
/// # Arguments
///
/// * `values` - Values on the eight cube corners. Must be shape `(2, 2, 2)`.
/// * `x` - Coordinate to interpolate at. Components must be `>= 0` and `<=1`. Must be shape `(3,)`.
pub fn interp_trilinear<S>(values: &ArrayBase<S, Ix3>, x: &Array1<f64>) -> f64
where
    S: Data<Elem = f64>,
{
    if values.dim() != (2, 2, 2) {
        panic!("Interp values are not the right shape {:?}", values.shape());
    }
    let m_x = 1. - x;

    let mut c: [f64; 4] = [0.0; 4];
    // Interpolate over x
    for iy in 0..2 {
        for iz in 0..2 {
            let iix = (2 * iy) + iz;
            c[iix] = values[[0, iy, iz]] * m_x[[0]] + values[[1, iy, iz]] * x[[0]];
        }
    }

    // Interpolate over y
    let mut c1: [f64; 2] = [0.0; 2];
    for iz in 0..2 {
        c1[iz] = c[iz] * m_x[[1]] + c[iz + 2] * x[[1]];
    }

    // Interpolate over z
    return c1[0] * m_x[[2]] + c1[1] * x[[2]];
}
