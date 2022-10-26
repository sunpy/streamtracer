import matplotlib.pyplot as plt
import pandas as pd

fig, ax = plt.subplots()

for v, label in zip(["120", "200"], ["v1.2 (FORTRAN)", "v2.0 (Rust)"]):
    data = pd.read_csv(f"v{v}.csv")
    print(data)

    ax.plot(data["nseeds"], data["time"], label=label, marker="o")

ax.legend()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Number of seeds")
ax.set_ylabel("Time / seconds")
ax.set_title("Comparison of streamtracer peformance for different versions")
plt.show()
