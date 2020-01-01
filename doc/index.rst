streamtracer
============

streamtracer is a python package for rapid streamline tracing. It is a wrapper
to compiled fortran code that does the heavy lifting, and is therefore
relatively fast.

To use, create a :class:`streamtracer.StreamTracer` object::

  nsteps = 10000
  step_size = 0.1
  tracer = StreamTracer(nsteps, step_size)

and a :class:`streamtracer.VectorGrid`::

  field = np.ones((10, 10, 10, 3))
  grid_spacing = [1, 2, 1]
  grid = VectorGrid(field, grid_spacing)

This can then be used to trace lines through a 3D cartesian vector field::

  seeds = np.array([[0, 0, 0], [0, 0, 1]])
  streamlines = StreamTracer.trace(seeds, grid)

For more information see the :mod:`streamtracer` API docs.

Installing
==========

In theory, it should be possible to build and install streamtracer in one
go with::

  pip install streamtracer

Note that this requires a fortran compiler; currently known to work is gfortran.
If you have problems installing, please open an issue at
https://github.com/dstansby/streamtracer/issues

Code reference
==============

.. toctree::
   :maxdepth: 1

   streamtracer

Changelog
=========

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
