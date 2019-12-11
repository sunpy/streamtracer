import numpy as np
import pytest

from streamtracer import StreamTracer, VectorGrid


@pytest.fixture
def tracer():
    return StreamTracer(1000, 0.1)


@pytest.fixture
def uniform_x_field():
    # A uniform field pointing in the x direction
    v = np.zeros((100, 100, 100, 3))
    # Make all vectors point in the x-direction
    v[:, :, :, 0] = 1
    spacing = [1, 1, 1]
    return VectorGrid(v, spacing)


@pytest.mark.parametrize('x0', [0, 50])
def test_different_seeds(tracer, uniform_x_field, x0):
    # Check that different seed points give sensible results
    seed = np.array([0, x0, x0])
    tracer.trace(seed, uniform_x_field)

    sline = tracer.xs[0]
    # Check that y, z coordinates are all zero
    np.testing.assert_equal(sline[:, 1], x0)
    np.testing.assert_equal(sline[:, 2], x0)


def test_uniform_field(tracer, uniform_x_field):
    # Check that tracing thought a uniform field gives sensible results
    seed = np.array([0, 0, 0])
    tracer.trace(seed, uniform_x_field)
    assert isinstance(tracer.xs[0], np.ndarray)
    assert len(tracer.xs) == len(tracer.ROT)

    sline = tracer.xs[0]
    # Check that y, z coordinates are all zero
    np.testing.assert_equal(sline[:, 1], 0)
    np.testing.assert_equal(sline[:, 2], 0)

    # Check that streamline always goes in a positive direction
    assert np.all(np.diff(sline[:, 0]) > 0)
    # Check that there are 100 * 0.1 = 1000 steps in the streamline
    assert sline.shape[0] == 1000


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
    seed = np.array([99.9, 50, 50])

    uniform_x_field.cyclic = [True, False, False]
    tracer.trace(seed, uniform_x_field)
    assert len(tracer.xs[0]) == (2 * maxsteps - 1)
    assert tracer.max_steps == maxsteps

    uniform_x_field.cyclic = [False, False, False]
    # Check that turning cyclic off interrupts the tracing
    tracer.trace(seed, uniform_x_field, direction=1)
    assert len(tracer.xs[0]) == 2


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
