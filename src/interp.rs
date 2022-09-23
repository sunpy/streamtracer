use ndarray::{Array3, Array1, Array};

/// Trilinear-interpolation of a scalar defined on
/// the eight corners of a cuboid.
///
/// # Arguments
///
/// * `values` - Values on the eight cube corners
/// * `x` - Coordinate to interpolate at. Components must be >= 0 and <=1
pub fn interp_trilinear(values: &Array3<f64>, x: &Array1<f64>) -> f64 {
    let m_x = 1. - x;

    let mut c = Array::zeros((2, 2));
    // Interpolate over x
    for iy in 0..2 {
        for iz in 0..2 {
            c[[iy, iz]] = values[[0, iy, iz]] * m_x[[0]] + values[[1, iy, iz]] * x[[0]];
        }
    }

    // Interpolate over y
    let mut c1 = Array::zeros(2);
    for iz in 0..2 {
        c1[[iz]] = c[[0, iz]] * m_x[[1]] + c[[1, iz]] * x[[1]];
    }

    // Interpolate over z
    return c1[[0]] * m_x[[2]] + c1[[1]] * x[[2]];
}
