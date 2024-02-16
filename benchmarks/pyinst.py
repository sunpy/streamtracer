import numpy as np
from pyinstrument import Profiler

from streamtracer import StreamTracer, VectorGrid

profiler = Profiler()
profiler.start()

nsteps = 1000
step_size = 0.1
tracer = StreamTracer(nsteps, step_size)

field = np.random.rand(*(180, 360, 50, 3))
grid_spacing = [1, 2, 3]
grid = VectorGrid(field, grid_spacing)
seeds = np.repeat([[90, 180, 25]], 2**10, axis=0)
tracer.trace(seeds, grid)

profiler.stop()

profiler.print(show_all=True)
