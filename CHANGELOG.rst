2.5.0 (2025-09-01)
==================

New Features
------------

- Bump minimum supported numpy to 1.26, and add support for Python 3.14. (`#211 <https://github.com/sunpy/streamtracer/pull/211>`__)


Internal Changes
----------------

- Bump various rust build dependencies for Python 3.14 support. (`#211 <https://github.com/sunpy/streamtracer/pull/211>`__)


2.4.1 (2025-04-09)
==================

Breaking Changes
----------------

- Increased version of pyo3 and numpy to account for "Risk of buffer overflow" exploit (`#201 <https://github.com/sunpy/streamtracer/pull/201>`__)


2.4.0 (2025-03-13)
==================

Breaking Changes
----------------

- The minimum supported version for macos binaries is now Sierra (10.12). (`#192 <https://github.com/sunpy/streamtracer/pull/192>`__)


2.3.0
=====

* Add support for Python 3.13 `#176 <https://github.com/sunpy/streamtracer/pull/176>`__
* Minor docstring formatting fixes, refactor input argument validation, and fix bug when tracing in +/-1 direction `#172 <https://github.com/sunpy/streamtracer/pull/172>`__
* Build binaries for aarch64 on native runners `#167 <https://github.com/sunpy/streamtracer/pull/167>`__

2.2.0
=====

* Improved performance of the Rust streamtracing algorithm, so it is now close to or faster than the old FORTRAN version.

2.1.0
=====

* Wheels for Python 3.12 are now built and published on PyPI.
* Rust code is now parallelized and while still not as fast as the old FORTRAN version, is near parity.

2.0.1
=====

* streamtracer now includes wheels for Python 3.11, and these have now been uploaded to PyPI as a release (version 2.0.0 was only ever uploaded as an alpha).
* streamtracer is now available on conda-forge.
* streamtracer still does not work in parallel, if you require parallel stream tracing then for now either:

    * Downgrade to version 1.2.0
    * Manually run several instances of the streamtracer in parallel

2.0.0
=====

* The low level streamline tracing code has been ported from FORTRAN to Rust.
  This has several benefits on the development side, and from a user perspective brings the first built wheels for Windows 🪟🎉
  The new code is **not** *yet* parallelized, so when tracing multiple streamlines on multiple cores runs slower than version 1.
  Version 2.1 (hopefully coming soon!) should implement parallel line tracing.
* The minimum supported Python version is now 3.8.

1.2.0
=====

New features
------------

* Added the ability to trace field lines through non-uniformly spaced grids.
  To do this pass a list of grid coordinates to the new ``grid_coords`` argument in `-streamtracer.Grid`.

Bug fixes
---------

* Fixed coordinate values returned by `streamtracer.Grid.xcoords` etc.
  Previously the origin coordinate was interpreted negative what it should be.
  This doesn't affect any traced streamlines.

1.1.2
=====

* Fixed the example code listed above.

1.1.1
=====

* Fixed wheel building to use the oldest supported version of numpy for each major version of python.
  This fixes errors like "module compiled against API version 0xe but this version of numpy is 0xd" that were in the 1.1.0 wheels.

1.1.0
=====

* Fixed handling of steps going out of bounds.
  Previously, in the forward direction a single step would be saved out of bounds, but in the backwards direction the streamline ended before going out of bounds.
  Now both the forward and negative directions both save a single out of bounds step in each stream line.
* Linux and macOS binary distributions (wheels) are now automatically built and uploaded, meaning you should no longer need a FORTRAN compiler installed locally to use streamtracer.
* Minor performance and memory improvements have been made.

1.0.1
=====

* Fix compilation issues on MacOS Big Sur.

1.0.0
=====

* Nothing major, just re-versioning to a stable release number.

0.1.2
=====

* Make sure to install numpy before trying to install streamtracer.

0.1.1
=====

* Added validation for the ``max_steps`` argument to :class:`-streamtracer.StreamTracer`.

0.1.0
=====

* First streamtracer release.
