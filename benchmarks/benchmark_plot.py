from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

fig, ax = plt.subplots()

version_names = {
    "120": "v1.2 (FORTRAN)",
    "200": "v2.0 (Rust)",
    "210dev0": "v2.1 (Rust Parallel)",
}

files = Path(".").glob("v*.csv")

for file in files:
    label = version_names.get(file.stem[1:], file.stem)
    data = pd.read_csv(file)
    print(data)

    ax.plot(data["nseeds"], data["time"], label=label, marker="o")

ax.legend()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Number of seeds")
ax.set_ylabel("Time / seconds")
ax.set_title("Comparison of streamtracer performance for different versions")
plt.show()
