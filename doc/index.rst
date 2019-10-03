streamtracer
============

streamtracer is a python package for rapid streamline tracing. It is a wrapper
to compiled fortran code that does the heavy lifting.

To use, create a :class:`streamtracer.StreamTracer` object::

  nsteps = 10000
  step_size = 0.1
  tracer = StreamTracer(nsteps, step_size)

This can then be used to trace lines through a 3D cartesian vector field::

  seeds = np.array([[0, 0, 0], [0, 0, 1]])
  field = np.ones((10, 10, 10, 3))
  grid_spacing = 1
  box_center = np.array([0, 0, 0])
  streamlines = StreamTracer.trace(seeds, field, grid_spacing, box_center)

For more imformation see the :mod:`streamtracer` API docs.

Caveats
=======

Code reference
==============

.. toctree::
   :maxdepth: 1

   streamtracer


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
