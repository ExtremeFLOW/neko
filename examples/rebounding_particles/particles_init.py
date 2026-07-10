import numpy as np


x = np.array([0.0, 0.0, 0.0])
y = 0.0 * np.ones_like(x)
z = 0.0 * np.ones_like(x)

u = np.array([1.0, 1.5, 0.7])
v = 0.0 * np.ones_like(x)
w = 0.0 * np.ones_like(x)

diameter = np.array([1e-3, 1e-3, 1e-3])
rho = np.array([1e7, 1e7, 1e7])

data = np.stack([x, y, z, u, v, w, diameter, rho], axis=-1)

np.savetxt("particles.csv", data, delimiter=",")
