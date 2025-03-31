import numpy as np

x = np.linspace(-0.5, 1.75, 200)

p = np.column_stack((x, np.zeros(x.shape) + 1e-10, 0.05*np.ones(x.shape)))

np.savetxt("probes.csv", p, delimiter=',')
