import numpy as np
import pytest

from streamtracer import StreamTracer


@pytest.fixture
def tracer():
    return StreamTracer(1000, 0.1)


@pytest.fixture
def uniform_field():
    # Uniform field in the x direction
    v = np.zeros((100, 100, 100, 3))
    # Make all vectors point in the x-direction
    v[:, :, :, 0] = 1
    return v


def test_smoke(tracer):
    # Smoke test for StreamTracer creation
    pass


def test_uniform_field(tracer, uniform_field):
    seed = np.array([0, 0, 0])
    grid_spacing = [1, 1, 1]
    xc = [0, 0, 0]
    tracer.trace(seed, uniform_field, grid_spacing, xc)
    assert isinstance(tracer.xs[0], np.ndarray)
    assert len(tracer.xs) == len(tracer.ROT)
    assert len(tracer.xs) == len(tracer.ns)

    sline = tracer.xs[0]
    # Check that y, z coordinates are all zero
    np.testing.assert_equal(sline[:, 1], 0)
    np.testing.assert_equal(sline[:, 2], 0)

    # Check that streamline always goes in a positive direction
    assert np.all(np.diff(sline[:, 0]) > 0)


def test_trace_direction(tracer, uniform_field):
    seed = np.array([0, 0, 0])
    grid_spacing = [1, 1, 1]
    xc = [0, 0, 0]
    tracer.trace(seed, uniform_field, grid_spacing, xc, direction=1)
    sline = tracer.xs[0]
    assert np.all(sline[:, 0] >= 0)

    tracer.trace(seed, uniform_field, grid_spacing, xc, direction=-1)
    sline = tracer.xs[0]

    assert np.all(sline[:, 0] <= 0)


def test_bad_input(tracer, uniform_field):
    seed = np.array([0, 0, 0])
    grid_spacing = [1, 1, 1]
    xc = [0, 0, 0]
    with pytest.raises(ValueError, match='seeds must be a 2D array'):
        tracer.trace(np.array([[[1], [1]], [[1], [1]]]), uniform_field, grid_spacing, xc)

    with pytest.raises(ValueError, match='seeds must have shape'):
        tracer.trace(np.array([1, 1]), uniform_field, grid_spacing, xc)

    with pytest.raises(ValueError, match='field must be a 4D array'):
        tracer.trace(seed, np.array([1]), grid_spacing, xc)

    with pytest.raises(ValueError, match='field must have shape'):
        tracer.trace(seed, np.ones((10, 10, 10, 2)), grid_spacing, xc)

    with pytest.raises(ValueError, match='grid spacing must have shape'):
        tracer.trace(seed, uniform_field, 1, xc)
