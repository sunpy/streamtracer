//! Structure for representing a 3D vector field defined on the corners
//! of a rectilinear grid.
use numpy::ndarray::{array, s, Array1, ArrayView1, ArrayView4};

use crate::interp::interp_trilinear;

/// Enum denoting whether a point is in or out of the bounds
/// of a VectorField grid.
pub enum Bounds {
    /// A position vector is in bounds.
    In,
    /// A position vector is out of bounds.
    Out,
}

/// A 3D vector field defined at grid corners.
pub struct VectorField<'a> {
    /// Grid points along x dimension. Must start at 0.
    pub xgrid: ArrayView1<'a, f64>,
    /// Grid points along y dimension. Must start at 0.
    pub ygrid: ArrayView1<'a, f64>,
    /// Grid points along z dimension. Must start at 0.
    pub zgrid: ArrayView1<'a, f64>,
    /// Vector values at each grid point. Must be shape
    /// (nx, ny, ny, 3), where (nx, ny, nz) are the number
    /// of coordinates along dimension.
    pub values: ArrayView4<'a, f64>,
    /// Whether each dimension should be treated as cyclic
    /// or not. Must be shape (3,).
    cyclic: ArrayView1<'a, bool>,
    /// Number of x coordinates.
    nx: usize,
    /// Number of y coordinates.
    ny: usize,
    /// Number of z coordinates.
    nz: usize,

    /// Upper boundaries
    upper_bounds: Array1<f64>,
}

impl VectorField<'_> {
    /// Create a new VectorField, checking for appropriate array shapes.
    pub fn new<'a>(
        xgrid: ArrayView1<'a, f64>,
        ygrid: ArrayView1<'a, f64>,
        zgrid: ArrayView1<'a, f64>,
        values: ArrayView4<'a, f64>,
        cyclic: ArrayView1<'a, bool>,
    ) -> VectorField<'a> {
        // Do some shape checking
        let nx = xgrid.len();
        let ny = ygrid.len();
        let nz = zgrid.len();
        let field_shape = values.shape();

        assert_eq!(field_shape[0], nx);
        assert_eq!(field_shape[1], ny);
        assert_eq!(field_shape[2], nz);
        assert_eq!(field_shape[3], 3);

        assert_eq!(cyclic.shape()[0], 3);

        // Check first coordinates are zero.
        assert_eq!(xgrid[0], 0.);
        assert_eq!(ygrid[0], 0.);
        assert_eq!(zgrid[0], 0.);

        let upper_bounds = array![xgrid[nx - 1], ygrid[ny - 1], zgrid[nz - 1]];

        return VectorField {
            xgrid,
            ygrid,
            zgrid,
            values,
            cyclic,
            nx,
            ny,
            nz,
            upper_bounds,
        };
    }

    /// Return grid index of the cell containing `x`.
    pub fn grid_idx(&self, x: ArrayView1<f64>) -> Array1<usize> {
        // Output array
        let mut grid_idx = array![self.nx - 2, self.ny - 2, self.nz - 2];

        // x
        for i in 0..self.nx - 1 {
            if x[0] >= self.xgrid[[i]] && x[0] < self.xgrid[[i + 1]] {
                grid_idx[0] = i;
                break;
            }
        }
        // y
        for i in 0..self.ny - 1 {
            if x[1] >= self.ygrid[[i]] && x[1] < self.ygrid[[i + 1]] {
                grid_idx[1] = i;
                break;
            }
        }
        // z
        for i in 0..self.nz - 1 {
            if x[2] >= self.zgrid[[i]] && x[2] < self.zgrid[[i + 1]] {
                grid_idx[2] = i;
                break;
            }
        }

        return grid_idx;
    }

    /// Get vector at position `x` using tri-linear interpolation.
    pub fn vector_at_position(&self, x0: ArrayView1<f64>) -> Array1<f64> {
        let cell_idx = self.grid_idx(x0);
        let cell_origin = array![
            self.xgrid[cell_idx[0]],
            self.ygrid[cell_idx[1]],
            self.zgrid[cell_idx[2]]
        ];
        let cell_size = array![
            self.xgrid[cell_idx[0] + 1] - cell_origin[0],
            self.ygrid[cell_idx[1] + 1] - cell_origin[1],
            self.zgrid[cell_idx[2] + 1] - cell_origin[2]
        ];

        // Distance along each cell edge in normalised units
        let cell_dist: Array1<f64> = (x0.to_owned() - cell_origin) / cell_size;

        // Eight corners of the cube that the position vector is
        // currently in
        let vec_cube = self.values.slice(s![
            cell_idx[0]..(cell_idx[0] + 2),
            cell_idx[1]..(cell_idx[1] + 2),
            cell_idx[2]..(cell_idx[2] + 2),
            ..
        ]);

        let mut vector_at_pos = array![0., 0., 0.];
        // Loop over vector components
        for i in 0..3 {
            vector_at_pos[[i]] = interp_trilinear(&vec_cube.slice(s![.., .., .., i]), &cell_dist);
        }
        return vector_at_pos;
    }

    /// If any of the dimensions of the grid are cyclic, wrap a coordinate.
    pub fn wrap_cyclic(&self, mut x: Array1<f64>) -> Array1<f64> {
        if self.cyclic[0] {
            x[0] = (x[0] + self.xgrid[self.nx - 1]) % self.upper_bounds[0];
        }
        if self.cyclic[1] {
            x[1] = (x[1] + self.ygrid[self.ny - 1]) % self.upper_bounds[1];
        }
        if self.cyclic[2] {
            x[2] = (x[2] + self.zgrid[self.nz - 1]) % self.upper_bounds[2];
        }

        return x;
    }

    /// Check whether a coordinate is in bounds of the grid.
    pub fn check_bounds(&self, x: ArrayView1<f64>) -> Bounds {
        if (x[0] < self.xgrid[0])
            || (x[0] > self.xgrid[self.nx - 1])
            || (x[1] < self.ygrid[0])
            || (x[1] > self.ygrid[self.ny - 1])
            || (x[2] < self.zgrid[0])
            || (x[2] > self.zgrid[self.nz - 1])
        {
            return Bounds::Out;
        }
        return Bounds::In;
    }
}
