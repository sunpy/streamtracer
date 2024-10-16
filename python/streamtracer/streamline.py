import numpy as np

from streamtracer._streamtracer_rust import trace_streamlines

__all__ = ["StreamTracer", "VectorGrid"]


class VectorGrid:
    """
    A grid of vectors.

    .. note::

        If any of *cyclic* are ``True``, then the grid values on each side of the
        cyclic dimension **must** match, e.g. if ``cyclic=[False, True, False]``,
        ``vectors[:, 0, :, :]`` must equal ``vectors[:, -1, :, :]``.

    Parameters
    ----------
    vectors : array-like
        A (nx, ny, nz, 3) shaped array. The three values at (i, j, k, :)
        specify the (x, y, z) components of the vector at index (i, j, k).
    grid_spacing : array-like, optional
        A (3,) shaped array, that contains the grid spacings in the (x, y, z)
        directions. If not specified ``grid_coords`` must be specified.
    origin_coord : [`float`, `float`, `float`], optional
        The coordinate of the ``vectors[0, 0, 0, :]`` vector at the corner of
        the box. Defaults to ``[0, 0, 0]``. This is not used if ``grid_coords``
        is specified.
    cyclic : [`bool`, `bool`, `bool`], optional
        Whether to have cyclic boundary conditions in each of the (x, y, z)
        directions. Defaults to ``[False, False, False]``.
    grid_coords : list[array], optional
        A list of length 3 storing the (x, y, z) coordinates of the grid. If not
        specified ``grid_spacing`` must be specified.
    """

    def __init__(
        self,
        vectors,
        grid_spacing=None,
        origin_coord=None,
        cyclic=None,
        *,
        grid_coords=None,
    ):
        if grid_spacing is not None and grid_coords is not None:
            raise ValueError(
                'Only one of "grid_spacing" and "grid_coords" can be specified.'
            )
        if grid_spacing is None and grid_coords is None:
            raise ValueError(
                'One of "grid_spacing" and "grid_coords" must be specified.'
            )
        if grid_coords is not None and origin_coord is not None:
            raise ValueError(
                'Specifying both "grid_coords" and "origin_coord" is ambiguous.'
            )
        self.grid_spacing = grid_spacing
        self.vectors = vectors
        self.cyclic = cyclic
        self.coords = grid_coords
        self.origin_coord = origin_coord

    @property
    def grid_spacing(self):
        "Physical spacing between grid points along each axis."
        return self._grid_spacing

    @grid_spacing.setter
    def grid_spacing(self, val):
        if val is not None:
            val = np.array(val)
            if val.shape != (3,):
                raise ValueError(
                    f"grid spacing must have shape (3,), got " f"{val.shape}"
                )
        self._grid_spacing = val

    @property
    def vectors(self):
        "Three-dimensional vector field through which the streamlines will be traced."
        return self._vectors

    @vectors.setter
    def vectors(self, val):
        if len(val.shape) != 4:
            raise ValueError("vectors must be a 4D array")
        if val.shape[-1] != 3:
            raise ValueError(
                "vectors must have shape (nx, ny, nz, 3), " f"got {val.shape}"
            )
        self._vectors = val

    @property
    def coords(self):
        "The physical coordinates along each axis of the grid."
        return self._coords

    @coords.setter
    def coords(self, val):
        if val is not None:
            if len(val) != 3:
                raise ValueError("coords must be len(3)")
            for i, dim in zip(range(3), ["x", "y", "z"]):
                shape = np.array(val[i]).shape
                if shape != (self.vectors.shape[i],):
                    raise ValueError(
                        f"Expected {self.vectors.shape[i]} {dim} "
                        f"coordinates but got {shape}"
                    )
        self._coords = val

    @property
    def cyclic(self):
        """
        Boolean describing whether to have cyclic boundary conditions in each of the (x, y, z)
        directions.
        """
        return self._cyclic

    @cyclic.setter
    def cyclic(self, val):
        if val is None:
            val = [False, False, False]
        dims = {0: "x", 1: "y", 2: "z"}
        s = [slice(None)] * 4
        for i, c in enumerate(val):
            if c:
                slc = s.copy()
                slc[i] = slice(0, 1)
                side1 = self.vectors[tuple(slc)]
                slc[i] = slice(-1, None)
                side2 = self.vectors[tuple(slc)]

                np.testing.assert_equal(
                    side1,
                    side2,
                    err_msg=f"grid values in dimension {dims[i]} (size {self.vectors.shape[i]}) "
                    "do not match on each side of the cube",
                )
        self._cyclic = np.array(val, dtype=bool)

    @property
    def origin_coord(self):
        """
        The physical coordinate corresponding to the index at ``(0,0,0)``.
        """
        return self._origin_coord

    @origin_coord.setter
    def origin_coord(self, val):
        if val is None:
            if self.grid_spacing is not None:
                self._origin_coord = np.array([0, 0, 0])
            else:
                self._origin_coord = np.array(
                    [self.xcoords[0], self.ycoords[0], self.zcoords[0]]
                )
        else:
            self._origin_coord = np.array(val)

    def _get_coords(self, i):
        if self.grid_spacing is not None:
            return (
                self.grid_spacing[i] * np.arange(self.vectors.shape[i])
                + self.origin_coord[i]
            )
        else:
            return self.coords[i]

    @property
    def xcoords(self):
        """
        Coordinates of the x grid points.
        """
        return self._get_coords(0)

    @property
    def ycoords(self):
        """
        Coordinates of the y grid points.
        """
        return self._get_coords(1)

    @property
    def zcoords(self):
        """
        Coordinates of the z grid points.
        """
        return self._get_coords(2)


class StreamTracer:
    """
    A streamline tracing class.

    Parameters
    ----------
    max_steps : `int`
        Number of steps available for each line. The maximum number of points
        on a single stream line is ``max_steps``.
    step_size : `float`
        Step size as a the fraction of cell size.
    """

    def __init__(self, max_steps, step_size):
        self.max_steps = max_steps
        self.ds = step_size
        self.xs = None

    @property
    def xs(self):
        """
        An array of the streamlines each with shape ``(n,3)``, which in general can have varying
        numbers of points.
        """
        return self._xs

    @xs.setter
    def xs(self, val):
        self._xs = val

    @property
    def ROT(self):
        """
        Reason(s) of termination.

        Integer array with shape ``len(xs)`` if traced in one direction,
        or ``(len(xs), 2)`` if traced in both directions.
        Can take the following values:

        - -1: Encountered a NaN
        - 1: Reached maximum available steps
        - 2: Out of bounds
        """
        return self._ROT

    @ROT.setter
    def ROT(self, val):
        self._ROT = val

    @property
    def max_steps(self):
        """
        Number of steps available for each line.

        The maximum number of points on a single stream line is ``max_steps``.
        """
        return self._max_steps

    @max_steps.setter
    def max_steps(self, val):
        if not isinstance(val, int):
            raise ValueError(f"max_steps must be an integer (got {type(val)})")

        if not val > 0:
            raise ValueError("max_steps must be greater than zero " f"(got {val})")

        self._max_steps = val

    def trace(self, seeds, grid, direction=0):
        """
        Trace streamlines.

        This traces streamlines from a series of seeds, through a vector
        field.

        Parameters
        ----------
        seeds : array-like with shape ``(n, 3)``
            Seed points.
        grid : `VectorGrid`
            Grid of field vectors.
        direction : `int`, optional
            Integration direction. ``0`` for both directions, ``1`` for
            forward, or ``-1`` for backwards.
        """
        if not isinstance(grid, VectorGrid):
            raise ValueError("grid must be an instance of StreamTracer")
        self.grid = grid
        self.x0 = seeds.copy()
        self.n_lines = seeds.shape[0]

        field = grid.vectors

        cyclic = grid.cyclic
        seeds = np.atleast_2d(seeds)

        # Validate shapes
        if len(seeds.shape) != 2:
            raise ValueError("seeds must be a 2D array")
        if seeds.shape[1] != 3:
            raise ValueError(f"seeds must have shape (n, 3), got {seeds.shape}")

        seeds = (seeds - grid.origin_coord).astype(np.float64)
        xcoords = (grid.xcoords - grid.origin_coord[0]).astype(np.float64)
        ycoords = (grid.ycoords - grid.origin_coord[1]).astype(np.float64)
        zcoords = (grid.zcoords - grid.origin_coord[2]).astype(np.float64)

        if direction == 1 or direction == -1:
            # Calculate streamlines
            self.xs, self.ns, self.ROT = trace_streamlines(
                seeds,
                xcoords,
                ycoords,
                zcoords,
                field,
                cyclic,
                direction,
                self.ds,
                self.max_steps,
            )

            self.xs += grid.origin_coord
            # Reduce the size of the arrays
            self.xs = [xi[:ni, :] for xi, ni in zip(self.xs, self.ns)]

        elif direction == 0:
            # Calculate forward streamline
            xs_f, ns_f, ROT_f = trace_streamlines(
                seeds,
                xcoords,
                ycoords,
                zcoords,
                field,
                cyclic,
                1,
                self.ds,
                self.max_steps,
            )
            # Calculate backward streamline
            xs_r, ns_r, ROT_r = trace_streamlines(
                seeds,
                xcoords,
                ycoords,
                zcoords,
                field,
                cyclic,
                -1,
                self.ds,
                self.max_steps,
            )

            xs_f += grid.origin_coord
            xs_r += grid.origin_coord

            # Stack the forward and reverse arrays
            self.xs = [
                np.vstack([xri[int(nr) - 1 : 0 : -1, :], xfi[: int(nf)]])
                for xri, xfi, nr, nf in zip(xs_r, xs_f, ns_r, ns_f)
            ]
            self.n_lines = np.fromiter([len(xsi) for xsi in self.xs], int)

            self.ROT = np.vstack([ROT_f, ROT_r]).T
        else:
            raise ValueError(f"Direction must be -1, 1 or 0 (got {direction})")

        # Filter out nans
        self.xs = [xi[~np.any(np.isnan(xi), axis=1), :] for xi in self.xs]
