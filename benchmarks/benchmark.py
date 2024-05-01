import time
import importlib.metadata

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from streamtracer import StreamTracer, VectorGrid

# Support old 1.x versions
__version__ = importlib.metadata.version("streamtracer")

nsteps = 1000
step_size = 0.1
tracer = StreamTracer(nsteps, step_size)

field = np.random.rand(*(180, 360, 50, 3))
grid_spacing = [1, 2, 3]
grid = VectorGrid(field, grid_spacing)

seedlist = 2 ** np.arange(12)
times = []
for nseeds in seedlist:
    t = time.time()
    # Trace from middle of field
    seeds = np.repeat([[90, 180, 25]], nseeds, axis=0)
    tracer.trace(seeds, grid)
    dt = time.time() - t
    assert len(tracer.xs) == nseeds
    times += [dt]
    print(nseeds, dt / nseeds, dt)


pd.DataFrame({"nseeds": seedlist, "time": times}).to_csv(
    f"v{__version__.replace('.', '')}.csv"
)

fig, ax = plt.subplots()

ax.plot(seedlist, times, marker=".")
ax.set_xlabel("n seeds")
ax.set_ylabel("time (s)")
ax.set_xscale("log")
ax.set_yscale("log")
plt.show()
