//! Structure for representing a 3D vector field defined on the corners
//! of a rectilinear grid.

use ndarray::{array, Array1, Array4, s};

use crate::interp::interp_trilinear;

/// A 3D vector field defined at grid corners.
pub struct VectorField<'a> {
    /// Grid points along x dimension.
    pub xgrid: &'a Array1<f64>,
    /// Grid points along y dimension.
    pub ygrid: &'a Array1<f64>,
    /// Grid points along z dimension.
    pub zgrid: &'a Array1<f64>,
    /// Vector values at each grid point. Must be shape
    /// (nx, ny, ny, 3), where (nx, ny, nz) are the number
    /// of coordinates along dimension.
    pub values: &'a Array4<f64>,
    /// Whether each dimension should be treated as cyclic
    /// or not. Must be shape (3,).
    cyclic: &'a Array1<bool>,
    /// Number of x coordinates.
    nx: usize,
    /// Number of y coordinates.
    ny: usize,
    /// Number of z coordinates.
    nz: usize
}

impl<'a> VectorField<'a> {
    /// Create a new VectorField, checking for appropriate array shapes.
    pub fn new(
        xgrid: &'a Array1<f64>,
        ygrid: &'a Array1<f64>,
        zgrid: &'a Array1<f64>,
        values: &'a Array4<f64>,
        cyclic: &'a Array1<bool>
    ) -> Self{
        // Do some shape checking
        let nx = xgrid.shape()[0];
        let ny = ygrid.shape()[0];
        let nz = zgrid.shape()[0];
        let field_shape = values.shape();

        assert_eq!(field_shape[0], nx);
        assert_eq!(field_shape[1], ny);
        assert_eq!(field_shape[2], nz);
        assert_eq!(field_shape[3], 3);

        assert_eq!(cyclic.shape()[0], 3);

        Self{xgrid, ygrid, zgrid, values, cyclic, nx, ny, nz}
    }

    /// Return grid index of the cell containing `x`.
    pub fn grid_idx(
        &self,
        x: &Array1<f64>
    ) -> Array1<usize>{

        let mut grid_idx = array![self.nx, self.ny, self.nz];

        // x
        for i in 0..self.nx{
            if x[0] >= self.xgrid[[i]] && x[0] < self.xgrid[[i+1]]{
                grid_idx[0] = i;
                break;
            }
        }
        // y
        for i in 0..self.ny{
            if x[1] >= self.ygrid[[i]] && x[1] < self.ygrid[[i+1]]{
                grid_idx[1] = i;
                break;
            }
        }
        // z
        for i in 0..self.nz{
            if x[2] >= self.zgrid[[i]] && x[2] < self.zgrid[[i+1]]{
                grid_idx[2] = i;
                break;
            }
        }

        return grid_idx
    }

    /// Get vector at position `x` using tri-linear interpolation.
    pub fn vector_at_position(
        &self,
        x0: &Array1<f64>
    ) -> Array1<f64>{
        let cell_idx = self.grid_idx(x0);
        let cell_origin = array![
            self.xgrid[cell_idx[0]],
            self.ygrid[cell_idx[1]],
            self.zgrid[cell_idx[2]]
        ];
        let cell_size = array![
            self.xgrid[cell_idx[0]+1] - cell_origin[0],
            self.ygrid[cell_idx[1]+1] - cell_origin[1],
            self.zgrid[cell_idx[2]+1] - cell_origin[2]
        ];

        // Distance along each cell edge in normalised units
        let cell_dist: Array1<f64> = (x0 - cell_origin) / cell_size;

        // Eight corners of the cube that the position vector is
        // currently in
        let vec_cube = self.values.slice(s![
            cell_idx[0]..(cell_idx[0] + 1),
            cell_idx[1]..(cell_idx[1] + 1),
            cell_idx[2]..(cell_idx[3] + 1),
            ..
        ]);


        let mut vector_at_pos = array![0., 0., 0.];
        // Loop over vector components
        for i in 0..3{
            vector_at_pos[[i]] = interp_trilinear(
                &vec_cube.slice(s![.., .., .., i]),
                &cell_dist);
        }
        return vector_at_pos;
    }
}
