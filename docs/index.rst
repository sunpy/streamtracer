************
streamtracer
************

streamtracer is a Python package for rapid streamline tracing on regularly spaced grids.
The actual streamline tracing is done at a low level in Rust, with a nice Python API provided on top.

.. toctree::
    :maxdepth: 1

    streamtracer
    whatsnew/index

Installing
==========

It is possible to install streamtracer in one go with:

.. code-block:: bash

  pip install streamtracer

or using conda with:

.. code-block:: bash

  conda install -c conda-forge streamtracer

There are wheels available for Linux, macOS and Windows.
If you need to compile from source, you will need to have a Rust compiler installed.

If you have problems installing, please open an issue at https://github.com/sunpy/streamtracer/issues

Usage
=====

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

When the stream tracer steps outside the boundary of the grid, the first point outside the grid is saved in the traced stream line.
