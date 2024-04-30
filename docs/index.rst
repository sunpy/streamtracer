************
streamtracer
************

streamtracer is a python package for rapid streamline tracing on regularly spaced grids.
The actual streamline tracing is done at a low level in Rust, with a nice Python API provided on top.

To use, create a :class:`streamtracer.StreamTracer` object

.. jupyter-execute::

  import numpy as np
  from streamtracer import StreamTracer, VectorGrid

  nsteps = 10000
  step_size = 0.1
  tracer = StreamTracer(nsteps, step_size)

and a :class:`streamtracer.VectorGrid`

.. jupyter-execute::

  field = np.ones((10, 10, 10, 3))
  grid_spacing = [1, 2, 1]
  grid = VectorGrid(field, grid_spacing)

This can then be used to trace lines through a 3D cartesian vector field

.. jupyter-execute::

  seeds = np.array([[0, 0, 0], [0, 0, 1]])
  tracer.trace(seeds, grid)

and the traced field lines can be accessed via. the ``.xs`` attribute

.. jupyter-execute::

  print(f'Number of traced lines: {len(tracer.xs)}')
  line_lengths = [len(x) for x in tracer.xs]
  print(f'Line lengths: {line_lengths}')

For more information see the :mod:`streamtracer` API docs.

Boundary handling
=================
When the stream tracer steps outside the boundary of the grid, the first point
outside the grid is saved in the traced stream line.

Installing
==========

It should be possible to build and install streamtracer in one go with::

  pip install streamtracer

This requires a fortran compiler; currently known to work is gfortran. On
macOS this can be installed with `brew <https://brew.sh/>`_  using
``brew install gfortran``


If you have problems installing, please open an issue at
https://github.com/dstansby/streamtracer/issues

Code reference
==============

.. toctree::
   :maxdepth: 1

   streamtracer

Changelog
=========

2.0.1
-----
streamtracer now includes wheels for Python 3.11, and these have now been uploaded to PyPI as a release (version 2.0.0 was only ever uploaded as an alpha).
streamtracer still does not work in parallel, if you require parallel stream tracing then for now either:

- Downgrade to version 1.2.0
- Manually run several instances of the streamtracer in parallel

2.0.0
-----
The low level streamline tracing code has been ported from FORTRAN to Rust.
This has several benefits on the development side, and from a user perspective brings the first built wheels for Windows 🪟🎉

The new code is **not** *yet* parallelised, so when tracing multiple streamlines on multiple cores runs slower than version 1.
Version 2.1 (hopefully coming soon!) should implement parallel line tracing.

The minimum supported Python version is now 3.8.

1.2.0
-----

New features
~~~~~~~~~~~~
Added the ability to trace field lines through non-uniformly spaced grids.
To do this pass a list of grid coordinates to the new ``grid_coords``
argument in `~streamtracer.Grid`.

Bug fixes
~~~~~~~~~
Fixed coordinate values returned by `streamtracer.Grid.xcoords` etc. Previously
the origin coordinate was interpreted negative what it should be. This doesn't
affect any traced streamlines.

1.1.2
-----
Fixed the example code listed above.

1.1.1
-----
Fixed wheel building to use the oldest supported version of numpy for each
major version of python. This fixes errors like
"module compiled against API version 0xe but this version of numpy is 0xd"
that were in the 1.1.0 wheels.

1.1.0
-----
- Fixed handling of steps going out of bounds. Previously, in the forward
  direction a single step would be saved out of bounds, but in the backwards
  direction the streamline ended before going out of bounds. Now both the forward
  and negative directions both save a single out of bounds step in each stream
  line.
- Linux and macOS binary distributions (wheels) are now automatically built and
  uploaded, meaning you should no longer need a FORTRAN compiler installed
  locally to use streamtracer.
- Minor performance and memory improvements have been made.

1.0.1
-----
- Fix compilation issues on MacOS Big Sur.

1.0.0
-----
- Nothing major, just re-versioning to a stable release number.

0.1.2
-----
- Make sure to install numpy before trying to install streamtracer.

0.1.1
-----
- Added validation for the ``max_steps`` argument to
  :class:`~streamtracer.StreamTracer`.

0.1.0
-----
First streamtracer release.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
