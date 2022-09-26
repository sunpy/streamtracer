/// Structure for holding a 3D vector field
use ndarray::{Array1, Array4};


pub struct VectorField<'a> {
    xgrid: &'a Array1<f64>,
    ygrid: &'a Array1<f64>,
    zgrid: &'a Array1<f64>,
    field: &'a Array4<f64>,
    cyclic: &'a Array1<bool>,
}

impl<'a> VectorField<'a> {
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

        Self{xgrid, ygrid, zgrid, field, cyclic}
    }
}
