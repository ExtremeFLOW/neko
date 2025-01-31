import matplotlib.pyplot as plt
import numpy as np
from os.path import join

DATA_PATH = "."

n_probes = 200


# %% read data
data = np.genfromtxt(join(DATA_PATH, "./output.csv"), skip_header=n_probes+1, delimiter=",")


x = np.linspace(0, np.pi * 2, n_probes)
ic = np.sin(x)
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.plot(x, ic, '--k', label="ICs")
ax1.plot(x, data[200:, 1], label=r"Final solution")

ax1.set_xlabel(r"$x$")
ax1.legend()
ax1.set_xlim(0, np.pi * 2)



ax2.plot(x, ic - data[200:, 1], '--k', label="Absolute error")
ax2.set_xlabel(r"$x$")
ax2.set_ylabel(r"$absolute error$")

plt.tight_layout()
plt.savefig(join(DATA_PATH, "results.png"), dpi=300)
