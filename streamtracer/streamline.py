import numpy as np
from streamtracer.fortran.streamtracer import streamtracer
from scipy.interpolate import RegularGridInterpolator as interpolate


__all__ = ['StreamTracer']


class StreamTracer:
    """
    A streamline tracing class.

    Parameters
    ----------
    n_steps : int
        Number of steps available for each line.
    step_size : float
        Step size in fraction of cell size.
    inner_boundary : bool, optional
        Whether to terminate calculation at a spherical inner boundary.
    r_IB : float, optional
        Radius of the inner boundary.
    outer_boundary : bool, optional
        Whether to terminate calculate at a spherical outer boundary.
    r_OB : float, optional
        Radius of the outer boundary.

    Attributes
    ----------
    xs : array of (n, 3) arrays
        An array of the streamlines, which in general have different numbers
        of points.
    """
    def __init__(self, n_steps, step_size,
                 inner_boundary=False, r_IB=None,
                 outer_boundary=False, r_OB=None):
        self.ns = n_steps
        self.ns0 = n_steps  # Save original number
        self.ds = step_size

        streamtracer.inner_boundary = inner_boundary
        streamtracer.r_IB = r_IB

        self._ROT_reasons = ['Uncalculated',
                             'Out of steps',
                             'Out of domain',
                             'Isnan']
        self._dir_str = {-1: 'Reverse',
                         0: 'Both',
                         1: 'Forward'}

    # Calculate the streamline from a vector array
    def trace(self, seeds, field, grid_spacing, box_center, direction=0):
        """
        Trace streamlines.

        This traces streamlines from a series of seeds, through a vector
        field.

        Parameters
        ----------
        seeds : (n, 3) array
            Seed points.
        field : (nx, ny, nz, 3) array
            Box of field vectors.
        grid_spacing : (3,) array
            Box gridpoint spacing in (x, y, z) directions.
        box_center : (3,) array
            Coordinate of the box center.
        direction : int, optional
            Integration direction. ``0`` for both directions, ``1`` for forward, or
            ``-1`` for backwards.
        """
        self.x0 = seeds.copy()
        self.n_lines = seeds.shape[0]
        streamtracer.ds = self.ds
        streamtracer.box_center = box_center.copy()

        grid_spacing = np.array(grid_spacing)
        seeds = np.atleast_2d(seeds)

        # Validate shapes
        if len(seeds.shape) != 2:
            raise ValueError('seeds must be a 2D array')
        if seeds.shape[1] != 3:
            raise ValueError(f'seeds must have shape (n, 3), got {seeds.shape}')

        if len(field.shape) != 4:
            raise ValueError('field must be a 4D array')
        if field.shape[-1] != 3:
            raise ValueError(f'field must have shape (nx, ny, nz, 3), got {field.shape}')

        if grid_spacing.shape != (3,):
            raise ValueError(f'grid spacing must have shape (3,), got {grid_spacing.shape}')

        self.x0 = np.array([xi + box_center for xi in self.x0])

        if direction == 1 or direction == -1:
            # Calculate streamlines
            self.xs, vs, ROT, self.ns = streamtracer.streamline_array(
                self.x0, field, grid_spacing, direction, self.ns)

            # Reduce the size of the array
            self.xs = np.array([xi[:ni, :] for xi, ni in zip(self.xs, self.ns)])
            vs = np.array([vi[:ni, :] for vi, ni in zip(vs, self.ns)])

            # Save the Reason of Termination
            self.ROT = ROT

        elif direction == 0:
            # Calculate forward streamline
            xs_f, vs_f, ROT_f, ns_f = streamtracer.streamline_array(
                self.x0, field, grid_spacing, 1, self.ns)
            # Calculate backward streamline
            xs_r, vs_r, ROT_r, ns_r = streamtracer.streamline_array(
                self.x0, field, grid_spacing, -1, self.ns)

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
        self.ROT = self.ROT[el]
        self.ns = self.ns[el]

        self.xs = np.array([xi - box_center for xi in self.xs])
        # Filter out nans
        nanrows = np.any(np.isnan(self.xs), axis=1)
        self.xs = self.xs[~nanrows, :]

        del vs
