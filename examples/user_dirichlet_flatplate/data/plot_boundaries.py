#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

o_d = np.genfromtxt("o.csv", delimiter = ",")
on_d = np.genfromtxt("on.csv", delimiter = ",")
v_d = np.genfromtxt("v.csv", delimiter = ",")

o = {"x": o_d[:,0],
     "y": o_d[:,1],
     "p": o_d[:,2],
     "u": o_d[:,3],
     "v": o_d[:,4]}
on = {"x": on_d[:,0],
     "y": on_d[:,1],
     "p": on_d[:,2],
     "u": on_d[:,3],
     "v": on_d[:,4]}

v = {"x": v_d[:,0],
     "y": v_d[:,1],
     "p": v_d[:,2],
     "u": v_d[:,3],
     "v": v_d[:,4]}


for data in [o, on, v]:
    plt.scatter(data["x"], data["y"], c=data["p"], marker = "+")

plt.colorbar()
plt.show()
plt.savefig("pressure.png")
