# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
import numpy as np

from streamtracer import StreamTracer, VectorGrid


class TimeSuite:
    def setup(self):
        self.tracer = StreamTracer(1000, 0.1)
        # A uniform field pointing in the x direction
        v = np.zeros((100, 100, 100, 3))
        # Make all vectors point in the x-direction
        v[:, :, :, 0] = 1
        spacing = [1, 1, 1]
        self.grid = VectorGrid(v, spacing)

    def time_uniform_grid(self):
        # Check that tracing thought a uniform field gives sensible results
        seed = np.array([0, 0, 0])
        seed = np.tile(seed, (1000, 1))
        self.tracer.trace(seed, self.grid)


"""
class MemSuite:
    def mem_list(self):
        return [0] * 256
"""
