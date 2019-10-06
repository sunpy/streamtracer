import numpy as np

from streamtracer import StreamTracer


def test_smoke():
    # Smoke test for StreamTracer creation
    StreamTracer(10000, 0.01)


def test_uniform_field():
    tracer = StreamTracer(10000, 0.1)

    seed = np.array([0, 0, 0])
    v = np.zeros((100, 100, 100, 3))
    # Make all vectors point in the x-direction
    v[:, :, :, 0] = 1
    grid_spacing = [1, 1, 1]
    xc = [0, 0, 0]
    tracer.trace(seed, v, grid_spacing, xc)
    assert isinstance(tracer.xs, np.ndarray)

    sline = tracer.xs[0]
    # Check that y, z coordinates are all zero
    np.testing.assert_equal(sline[:, 1], 0)
    np.testing.assert_equal(sline[:, 2], 0)

    # Check that streamline always goes in a positive direction
    assert np.all(np.diff(sline[:, 0]) > 0)
