import numpy as np
from streamtracer.fortran.streamtracer import streamtracer
from scipy.interpolate import RegularGridInterpolator as interpolate
import matplotlib.pyplot as plt


__all__ = ['StreamTracer']


class StreamTracer:
    """
    Parameters
    ----------
    n_steps : int
        Number of steps.
    step_size : float
        Step size in fraction of cell size.
    direction : int, optional
        Integration direction. ``0`` for both directions, ``1`` for forward, or
        ``-1`` for backwards.
    inner_boundary : bool, optional
        Whether to terminate calculation at the inner boundary.
    r_IB : float, optional
        Radius of the inner boundary.

    Attributes
    ----------
    xs : array of (*, 3) arrays
        An array of the streamlines, which in general have different numbers
        of points.
    """
    def __init__(self, n_steps, step_size, direction=0,
                 inner_boundary=True, r_IB=1.):
        self.ns = n_steps
        self.ns0 = n_steps  # Save original number
        self.ds = step_size
        self.dir = direction

        streamtracer.inner_boundary = inner_boundary
        streamtracer.r_IB = 1.

        self._ROT_reasons = ['Uncalculated',
                             'Out of steps',
                             'Out of domain',
                             'Isnan']
        self._dir_str = {-1: 'Reverse',
                         0: 'Both',
                         1: 'Forward'}

        self.var = {}
        self.var_names = []
        self.cell_data = {}

    def reset(self, ns=None, ds=None):
        del self.xs
        del self.ROT

        if(ns is None):
            ns = self.ns0
        if(ds is None):
            ds = self.ds

        self.__init__(ns, ds)

    # Calculate the streamline from a vector array
    def calc(self, x0, v, d, xc, v_name='v'):
        """
        Parameters
        ----------
        x0 : (n, 3) array
            Seed points.
        v : (nx, ny, nz, 3) array
            Box of magnetic field vectors.
        d : (3,) array
            Box gridpoint spacing in (x, y, z) directions.
        xc : (3,) array
            Coordinate of the box center.
        """
        self.x0 = x0.copy()
        self.n_lines = x0.shape[0]
        streamtracer.ds = self.ds
        streamtracer.xc = xc.copy()

        # Validate shapes
        if len(x0.shape) != 2:
            raise ValueError('x0 must be a 2D array')
        if x0.shape[1] != 3:
            raise ValueError(f'x0 must have shape (n, 3), got {x0.shape}')

        if len(v.shape) != 4:
            raise ValueError('field must be a 4D array')
        if v.shape[-1] != 3:
            raise ValueError(f'field must have shape (nx, ny, nz, 3), got {v.shape}')

        self.x0 = np.array([xi + xc for xi in self.x0])

        if self.dir == 1 or self.dir == -1:
            # Calculate streamlines
            self.xs, vs, ROT, self.ns = streamtracer.streamline_array(
                self.x0, v, d, self.dir, self.ns)

            # Reduce the size of the array
            self.xs = np.array([xi[:ni, :] for xi, ni in zip(self.xs, self.ns)])
            vs = np.array([vi[:ni, :] for vi, ni in zip(vs, self.ns)])

            # Save the Reason of Termination
            self.ROT = ROT

        elif self.dir == 0:
            # Calculate forward streamline
            xs_f, vs_f, ROT_f, ns_f = streamtracer.streamline_array(
                self.x0, v, d, 1, self.ns)
            # Calculate backward streamline
            xs_r, vs_r, ROT_r, ns_r = streamtracer.streamline_array(
                self.x0, v, d, -1, self.ns)

            # Reduce the size of the arrays, and flip the reverse streamlines
            xs_f = np.array([xi[:ni, :] for xi, ni in zip(xs_f, ns_f)])
            vs_f = np.array([vi[:ni, :] for vi, ni in zip(vs_f, ns_f)])

            xs_r = np.array([xi[ni - 1:0:-1, :] for xi, ni in zip(xs_r, ns_r)])
            vs_r = np.array([vi[ni - 1:0:-1, :] for vi, ni in zip(vs_r, ns_r)])

            # Stack the forward and reverse arrays
            self.xs = np.array([np.vstack([xri, xfi]) for xri, xfi in zip(xs_r, xs_f)])
            vs = np.array([np.vstack([vri, vfi]) for vri, vfi in zip(vs_r, vs_f)])
            self.ns = np.fromiter([len(xsi) for xsi in self.xs], int)

            self.ROT = np.vstack([ROT_f, ROT_r]).T

        # Remove streamlines with zero size
        el = self.ns > 1
        self.ROT = self.ROT[el]  # , :
        self.ns = self.ns[el]

        self.xs = np.array([xi - xc for xi in self.xs])

        self.var[v_name] = vs.copy()
        self.var_names = np.array([s for s in self.var])
        del vs

        for s in self.cell_data:
            self.cell_data[s] = self.cell_data[s][el]

    def interp(self, x, y, z, v, var_name):
        """
        Interpolate for other quantities.
        """

        if len(v.shape) == 3:
            vi = self._interp_scalar(x, y, z, v)
        elif len(v.shape) == 4:
            vi = self._interp_vector(x, y, z, v)

        self.var[var_name] = vi

        self.var_names = np.array([s for s in self.var])

    def _interp_scalar(self, x, y, z, f):
        I = interpolate((x, y, z), f, bounds_error=False)

        xI = np.vstack(self.xs)
        fI = I(xI)

        fI = np.array(np.split(fI, np.cumsum(self.ns)))[:-1]

        return fI

    def _interp_vector(self, x, y, z, v):

        Ix = interpolate((x, y, z), v[:,:,:,0], bounds_error=False, fill_value=None)
        Iy = interpolate((x, y, z), v[:,:,:,1], bounds_error=False, fill_value=None)
        Iz = interpolate((x, y, z), v[:,:,:,2], bounds_error=False, fill_value=None)

        xI = np.vstack(self.xs)
        vI = np.array([Ix(xI), Iy(xI), Iz(xI)]).T

        vI = np.array(np.vsplit(vI, np.cumsum(self.ns)))[:-1]

        return vI

        #  Write to vtp

    def _write_vtp(self, fname, pts_step=10):
        import vtk
        from vtk.util import numpy_support as vtk_np

        if(pts_step is None):
            pts_step = int(max(1, 1./self.ds))
        print(pts_step)

        # Points
        pts = np.vstack([xi[::pts_step] for xi in self.xs]).ravel()
        doubleArray = vtk_np.numpy_to_vtk(pts)
        doubleArray.SetNumberOfComponents(3)

        points = vtk.vtkPoints()
        points.SetData(doubleArray)

        # Cells

        n_pts_in_cell = np.array([len(xi[::pts_step]) for xi in self.xs])

        i = np.arange(np.sum(n_pts_in_cell), dtype=np.int64)
        i = np.array(np.split(i, n_pts_in_cell.cumsum())[:-1])

        id_array = np.array([np.hstack([ni, ii]) for ni, ii in zip(n_pts_in_cell, i)])

        id_array = np.hstack([ii for ii in id_array])

        cellArray = vtk.vtkCellArray()
        idArray = vtk_np.numpy_to_vtkIdTypeArray(id_array)
        cellArray.SetCells(len(self.ns), idArray)

        # Pointdata

        point_arrays = self._vtk_pointData_arrays(pts_step)

        # Cell Data

        cell_arrays = self._vtk_cellData_arrays(pts_step)

        # Polydata

        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.SetLines(cellArray)

        for Arr in point_arrays:
            polyData.GetPointData().AddArray(Arr)

        for Arr in cell_arrays:
            polyData.GetCellData().AddArray(Arr)

        # Write to file

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fname)
        writer.SetInputData(polyData)
        writer.Write()

    def _vtk_pointData_arrays(self, pts_step):

        def create_pointDataArrays(arr, name):

            if(len(arr[0].shape)==1):
                data = np.hstack([fi[::pts_step] for fi in arr])
            else:
                data = np.vstack([fi[::pts_step] for fi in arr]).ravel()

            data = data.astype(np.float64)

            doubleArray = vtk_np.numpy_to_vtk(data, deep=1)
            doubleArray.SetName(name)

            if(len(arr[0].shape)>1):
                doubleArray.SetNumberOfComponents(arr[0].shape[1])

            return doubleArray

        return [create_pointDataArrays(self.var[name], name) for name in self.var]

    def _vtk_cellData_arrays(self, pts_step):

        def create_cellDataArrays(arr, name):
            from vtk.util import numpy_support as vtk_np

            data = arr.ravel()

            data = data.astype(np.float64)

            doubleArray = vtk_np.numpy_to_vtk(data, deep=1)
            doubleArray.SetName(name)

            if(len(arr.shape)>1):
                doubleArray.SetNumberOfComponents(arr.shape[1])

            return doubleArray

        cell_arrays = [self.ns, self.ROT]
        cell_names = ['ns', 'ROT']

        for s in self.cell_data:
            cell_arrays.append(self.cell_data[s])
            cell_names.append(s)

        return [create_cellDataArrays(arr, name) for arr, name in zip(cell_arrays, cell_names)]
