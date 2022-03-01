import numpy as np
import pytest

from streamtracer.fortran.streamtracer import streamtracer
from streamtracer import StreamTracer, VectorGrid


@pytest.fixture
def tracer():
    return StreamTracer(2000, 0.1)


@pytest.fixture
def uniform_x_field():
    # A uniform field pointing in the x direction
    v = np.zeros((101, 101, 101, 3))
    # Make all vectors point in the x-direction
    v[:, :, :, 0] = 1
    spacing = [1, 1, 1]
    return VectorGrid(v, spacing)


@pytest.mark.parametrize('direction, ROTs', [[1, np.array([2, 2, 2])],
                                             [-1, np.array([2, 2, 2])],
                                             [0, np.array([[2, 2],
                                                           [2, 2],
                                                           [2, 2]])]
                                             ])
def test_rot(tracer, uniform_x_field, direction, ROTs):
    seeds = np.array([[50, 50, 50], [50, 50, 50], [50, 50, 50]])
    tracer.trace(seeds, uniform_x_field, direction=direction)
    np.testing.assert_equal(tracer.ROT, ROTs)


@pytest.mark.parametrize('x0', [0, 50])
def test_different_seeds(tracer, uniform_x_field, x0):
    # Check that different seed points give sensible results
    seed = np.array([0, x0, x0])
    tracer.trace(seed, uniform_x_field)

    sline = tracer.xs[0]
    # Check that y, z coordinates are all zero
    np.testing.assert_equal(sline[:, 1], x0)
    np.testing.assert_equal(sline[:, 2], x0)


@pytest.mark.parametrize('vec_len', [1, 2])
def test_uniform_field(tracer, uniform_x_field, vec_len):
    # Check that tracing thought a uniform field gives sensible results
    seed = np.array([0, 0, 0])
    uniform_x_field.vectors *= vec_len
    tracer.trace(seed, uniform_x_field)
    assert isinstance(tracer.xs[0], np.ndarray)
    assert len(tracer.xs) == len(tracer.ROT)

    sline = tracer.xs[0]
    # Check that y, z coordinates are all zero
    np.testing.assert_almost_equal(sline[:, 0], np.linspace(-0.1, 100.1, 1003))
    np.testing.assert_equal(sline[:, 1], 0)
    np.testing.assert_equal(sline[:, 2], 0)

    # Check that streamline always goes in a positive direction
    assert np.all(np.diff(sline[:, 0]) > 0)
    # Check that there are 100 * 0.1 = 1000 steps in the streamline
    assert sline.shape[0] == 1003


def test_trace_direction(tracer, uniform_x_field):
    # Check that the direction keyword argument works
    seed = np.array([0, 0, 0])
    tracer.trace(seed, uniform_x_field, direction=1)
    sline = tracer.xs[0]
    assert np.all(sline[:, 0] >= 0)

    tracer.trace(seed, uniform_x_field, direction=-1)
    sline = tracer.xs[0]

    assert np.all(sline[:, 0] <= 0)


def test_cyclic(uniform_x_field):
    # Check the cyclic option
    maxsteps = 4
    tracer = StreamTracer(maxsteps, 0.1)
    seed = np.array([99.95, 50, 50])

    uniform_x_field.cyclic = [True, False, False]
    tracer.trace(seed, uniform_x_field, direction=1)
    fline = tracer.xs[0]
    # Check that the tracer uses all the steps available
    assert len(fline) == maxsteps
    assert tracer.max_steps == maxsteps
    # Check that the cyclic boundary works properly
    np.testing.assert_equal(fline[0, :], seed)
    np.testing.assert_almost_equal(fline[1, :], np.array([0.05, 50, 50]))
    np.testing.assert_almost_equal(fline[2, :], np.array([0.15, 50, 50]))

    # Check that going the other way across the boundary (through zero) works
    seed = np.array([0.1, 50, 50])
    tracer.trace(seed, uniform_x_field, direction=-1)
    fline = tracer.xs[0]
    np.testing.assert_equal(fline[0, :], seed)
    np.testing.assert_almost_equal(fline[1, :], np.array([0, 50, 50]))
    np.testing.assert_almost_equal(fline[2, :], np.array([99.9, 50, 50]))

    # Check that turning cyclic off interrupts the tracing
    uniform_x_field.cyclic = [False, False, False]
    # Going forwards
    seed = np.array([99.9, 50, 50])
    tracer.trace(seed, uniform_x_field, direction=1)
    assert len(tracer.xs[0]) == 3

    # Going backwards
    seed = np.array([0.1, 50, 50])
    tracer.trace(seed, uniform_x_field, direction=-1)
    assert len(tracer.xs[0]) == 3


@pytest.mark.parametrize('origin_coord', [[0, 0, 0], [1, 1, 1]])
def test_origin(uniform_x_field, tracer, origin_coord):
    uniform_x_field.origin_coord = origin_coord
    seed = np.array([50, 50, 50])
    tracer.trace(seed, uniform_x_field, direction=1)

    np.testing.assert_equal(tracer.xs[0][0, :], seed)
    end_coord = np.array([100 + origin_coord[0] + 0.1, 50, 50])
    np.testing.assert_almost_equal(tracer.xs[0][-1, :], end_coord)


def test_cyclic_field():
    # A uniform field pointing in the x direction
    v = np.zeros((100, 100, 100, 3))
    # Make all vectors point in the x-direction
    v[:, :, :, 0] = 1
    # Mismatch vectors on the x-faces
    v[0, :, :, 0] = -1

    spacing = [1, 1, 1]
    cyclic = [True, False, False]
    with pytest.raises(AssertionError, match='Arrays are not equal'):
        VectorGrid(v, spacing, cyclic=cyclic)

    # Check that with cyclic off no error is thrown
    cyclic = [False, False, False]
    VectorGrid(v, spacing, cyclic=cyclic)


def test_grid_points():
    # A uniform field pointing in the x direction
    v = np.zeros((100, 100, 100, 3))
    # Make all vectors point in the x-direction
    spacing = [2, 3, 4]
    origin_coord = [4, 9, 16]
    grid = VectorGrid(v, spacing, origin_coord=origin_coord)
    assert grid.xcoords[0] == 4
    assert grid.xcoords[1] == 6
    assert grid.xcoords[-1] == 4 + 99 * 2

    assert grid.ycoords[0] == 9
    assert grid.ycoords[1] == 12
    assert grid.ycoords[-1] == 9 + 99 * 3

    assert grid.zcoords[0] == 16
    assert grid.zcoords[1] == 20
    assert grid.zcoords[-1] == 16 + 99 * 4


def test_bad_input(tracer, uniform_x_field):
    # Check input validation
    seed = np.array([0, 0, 0])
    with pytest.raises(ValueError, match='seeds must be a 2D array'):
        tracer.trace(np.array([[[1], [1]], [[1], [1]]]), uniform_x_field)

    with pytest.raises(ValueError, match='seeds must have shape'):
        tracer.trace(np.array([1, 1]), uniform_x_field)

    with pytest.raises(ValueError, match='grid must be an instance of StreamTracer'):
        tracer.trace(seed, 1)

    with pytest.raises(ValueError, match='Direction must be -1, 1 or 0'):
        tracer.trace(seed, uniform_x_field, direction=2)

    with pytest.raises(ValueError, match='vectors must be a 4D array'):
        VectorGrid(np.array([1]), [1, 1, 1])

    with pytest.raises(ValueError, match='vectors must have shape'):
        VectorGrid(np.zeros((1, 1, 1, 2)), [1, 1, 1])

    with pytest.raises(ValueError, match='grid spacing must have shape'):
        VectorGrid(np.zeros((1, 1, 1, 3)), [1, 1])


@pytest.mark.parametrize('val, errstr',
                         [(0.1, 'max_steps must be an integer'),
                          (0, 'max_steps must be greater than zero'),
                          (-1, 'max_steps must be greater than zero')])
def test_invalid_max_steps(val, errstr):
    with pytest.raises(ValueError, match=errstr):
        StreamTracer(val, 0.1)


# Paramatrize to make sure behaviour is same in x,y,z directions
@pytest.mark.parametrize('dir', [0, 1, 2])
def test_bounds(dir):
    # A uniform field pointing in the x direction
    v = np.zeros((2, 2, 2, 3))
    # Make all vectors point in the x-direction
    v[:, :, :, dir] = 1
    spacing = [1, 1, 1]
    grid = VectorGrid(v, spacing)

    seed = np.array([[0.5, 0.5, 0.5]])
    tracer = StreamTracer(10, 1.0)
    tracer.trace(seed, grid, direction=0)
    # Check that one step is out of bounds in the positive direction
    expected = np.roll(np.array([1.5, 0.5, 0.5]), dir)
    assert (tracer.xs[0][-1, :] == expected).all()
    # Check that one step is out of bounds in the negative direction
    expected = np.roll(np.array([-0.5, 0.5, 0.5]), dir)
    assert (tracer.xs[0][0, :] == expected).all()


def test_nthreads():
    # assert streamtracer.openmp_enabled
    streamtracer.set_num_threads(4)
    assert streamtracer.get_num_threads() == 4

    streamtracer.set_num_threads(2)
    assert streamtracer.get_num_threads() == 2
