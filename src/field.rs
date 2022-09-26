/// Structure for holding a 3D vector field
use ndarray::{array, Array1, Array4};


pub struct VectorField<'a> {
    xgrid: &'a Array1<f64>,
    ygrid: &'a Array1<f64>,
    zgrid: &'a Array1<f64>,
    field: &'a Array4<f64>,
    cyclic: &'a Array1<bool>,
    nx: usize,
    ny: usize,
    nz: usize
}

impl<'a> VectorField<'a> {
    // Create a new VectorField, checking for appropriate array shapes
    pub fn new(
        xgrid: &'a Array1<f64>,
        ygrid: &'a Array1<f64>,
        zgrid: &'a Array1<f64>,
        field: &'a Array4<f64>,
        cyclic: &'a Array1<bool>
    ) -> Self{
        // Do some shape checking
        let nx = xgrid.shape()[0];
        let ny = ygrid.shape()[0];
        let nz = zgrid.shape()[0];
        let field_shape = field.shape();

        assert_eq!(field_shape[0], nx);
        assert_eq!(field_shape[1], ny);
        assert_eq!(field_shape[2], nz);
        assert_eq!(field_shape[3], 3);

        assert_eq!(cyclic.shape()[0], 3);

        Self{xgrid, ygrid, zgrid, field, cyclic, nx, ny, nz}
    }

    // Return grid index of the cell containing `x`.
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
}
