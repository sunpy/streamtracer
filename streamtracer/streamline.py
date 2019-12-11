import numpy as np
from streamtracer.fortran.streamtracer import streamtracer
from scipy.interpolate import RegularGridInterpolator as interpolate


__all__ = ['StreamTracer', 'VectorGrid']


class VectorGrid:
    """
    A grid of vectors.

    Parameters
    ----------
    vectors : array
        A (nx, ny, nz, 3) shaped array. The three values at (i, j, k, :)
        specify the (x, y, z) components of the vector at index (i, j, k).
    grid_spacing : array
        A (3,) shaped array, that contains the grid spacings in the (x, y, z)
        directions.
    cyclic : [bool, bool, bool], optional
        Whether to have cyclic boundary conditions in each of the (x, y, z)
        directions.

    Notes
    -----
    If any of *cyclic* are ``True``, then the grid values on each side of the
    cyclic dimension **must** match, e.g. if ``cyclic=[False, True, False]``,
    ``vectors[:, 0, :, :]`` must equal ``vectors[:, -1, :, :]``.
    """
    def __init__(self, vectors, grid_spacing, cyclic=[False, False, False]):
        grid_spacing = np.array(grid_spacing)
        self._validate_vectors(vectors)
        self._validate_spacing(grid_spacing)
        self._validate_cyclic(vectors, cyclic)

        self.vectors = vectors
        self.grid_spacing = grid_spacing
        self.cyclic = cyclic

    @staticmethod
    def _validate_vectors(vectors):
        if len(vectors.shape) != 4:
            raise ValueError('vectors must be a 4D array')
        if vectors.shape[-1] != 3:
            raise ValueError('vectors must have shape (nx, ny, nz, 3), '
                             f'got {vectors.shape}')

    @staticmethod
    def _validate_spacing(grid_spacing):
        if grid_spacing.shape != (3,):
            raise ValueError(f'grid spacing must have shape (3,), got '
                             f'{grid_spacing.shape}')

    @staticmethod
    def _validate_cyclic(vectors, cyclic):
        s = [slice(None)] * 4
        for i, c in enumerate(cyclic):
            if c:
                slc = s.copy()
                slc[i] = slice(0)
                side1 = vectors[tuple(slc)]
                slc[i] = slice(-1)
                side2 = vectors[tuple(slc)]

                np.testing.assert_equal(
                    side1, side2,
                    err_msg=f'grid values in dimension {i} '
                    'do not match on each side of the cube')

    @property
    def cyclic(self):
        return self._cyclic

    @cyclic.setter
    def cyclic(self, val):
        self._cyclic = np.array(val, dtype=int)


class StreamTracer:
    """
    A streamline tracing class.

    Parameters
    ----------
    max_steps : int
        Number of steps available for each line.
    step_size : float
        Step size as a the fraction of cell size.
    cyclic : [bool, bool, bool], optional
        Whether to have cyclic boundary conditions in each dimension.

    Attributes
    ----------
    xs : array of (n, 3) arrays
        An array of the streamlines, which in general can have varying
        numbers of points.
    """
    def __init__(self, max_steps, step_size):
        self.max_steps = max_steps
        self.ds = step_size

    # Calculate the streamline from a vector array
    def trace(self, seeds, grid, direction=0):
        """
        Trace streamlines.

        This traces streamlines from a series of seeds, through a vector
        field.

        Parameters
        ----------
        seeds : (n, 3) array
            Seed points.
        grid : VectorGrid
            Grid of field vectors.
        direction : int, optional
            Integration direction. ``0`` for both directions, ``1`` for
            forward, or ``-1`` for backwards.
        """
        if not isinstance(grid, VectorGrid):
            raise ValueError('grid must be an instance of StreamTracer')
        self.grid = grid
        self.x0 = seeds.copy()
        self.n_lines = seeds.shape[0]

        # Set the step size (this is a module level variable for streamtracer)
        streamtracer.ds = self.ds

        field = grid.vectors
        grid_spacing = grid.grid_spacing
        cyclic = grid.cyclic
        seeds = np.atleast_2d(seeds)

        # Validate shapes
        if len(seeds.shape) != 2:
            raise ValueError('seeds must be a 2D array')
        if seeds.shape[1] != 3:
            raise ValueError(f'seeds must have shape (n, 3), got {seeds.shape}')

        if direction == 1 or direction == -1:
            # Calculate streamlines
            self.xs, vs, ROT, self.ns = streamtracer.streamline_array(
                seeds, field, grid_spacing, direction, self.max_steps, cyclic)

            # Reduce the size of the array
            self.xs = np.array([xi[:ni, :] for xi, ni in zip(self.xs, self.ns)])
            vs = np.array([vi[:ni, :] for vi, ni in zip(vs, self.ns)])

            # Save the Reason of Termination
            self.ROT = ROT

        elif direction == 0:
            # Calculate forward streamline
            xs_f, vs_f, ROT_f, ns_f = streamtracer.streamline_array(
                seeds, field, grid_spacing, 1, self.max_steps, cyclic)
            # Calculate backward streamline
            xs_r, vs_r, ROT_r, ns_r = streamtracer.streamline_array(
                seeds, field, grid_spacing, -1, self.max_steps, cyclic)

            # Reduce the size of the arrays, and flip the reverse streamlines
            xs_f = np.array([xi[:ni, :] for xi, ni in zip(xs_f, ns_f)])
            vs_f = np.array([vi[:ni, :] for vi, ni in zip(vs_f, ns_f)])

            xs_r = np.array([xi[ni - 1:0:-1, :] for xi, ni in zip(xs_r, ns_r)])
            vs_r = np.array([vi[ni - 1:0:-1, :] for vi, ni in zip(vs_r, ns_r)])

            # Stack the forward and reverse arrays
            self.xs = np.array([np.vstack([xri, xfi]) for xri, xfi in zip(xs_r, xs_f)])
            vs = np.array([np.vstack([vri, vfi]) for vri, vfi in zip(vs_r, vs_f)])
            self.n_lines = np.fromiter([len(xsi) for xsi in self.xs], int)

            self.ROT = np.vstack([ROT_f, ROT_r]).T
        else:
            raise ValueError(f'Direction must be -1, 1 or 0 (got {direction})')

        # Filter out nans
        xi = self.xs[0]
        self.xs = [xi[~np.any(np.isnan(xi), axis=1), :] for xi in self.xs]

        del vs
