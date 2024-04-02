#!/bin/env python3

import matplotlib.pyplot as plt
from matplotlib import pylab
import numpy as np
from os.path import join

from HistPoint import HistPointList
from mathTools import autocorr, integrate

plt.close("all")

params = {'legend.fontsize': 13,
          'legend.loc': 'best',
          'figure.figsize': (10, 8),
          'lines.markerfacecolor': 'none',
          'axes.labelsize': 15,
          'axes.titlesize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15,
          'grid.alpha': 0.6}
pylab.rcParams.update(params)

if __name__ == "__main__":

    path = "/scratch/baconnet/simulations/neko/fst_sandbox/"

    hp = HistPointList(join(path, "output.csv"))  # Generate a list of history point objects

    dt = 2e-2

    x = -0.5
    y = -0.5
    z = -0.5
    deltax = 1
    deltay = 1
    deltaz = 1

    ds = hp.select_points(x-deltax, x+deltax, y-deltay, y+deltay, z-deltaz, z+deltaz)

    if len(ds) == 0:
        raise ValueError(f"No history points in ({x} +/- {deltax}, {y} +/- \
 {deltay}, {z} +/- {deltaz})")

    what = "u"

    fig, axs = plt.subplots(ncols=2)

    for d in ds:
        print(d)
        auto_correlation = autocorr(d[what])
        print("--> IL :",integrate(auto_correlation, dt))
        
        # Plot velocity
        d.plot(fig, axs[0], what = "mag")
        
        # Plot autocorrelation
        # axs[1].plot(d["t"], auto_correlation, label = d.coords_to_str())
        
        # Compute and plot fourier in time
        (freq, power) = d.power_density(nperseg=50)
        axs[1].loglog(freq*2.0*np.pi, power*2.0*np.pi, "+", label = d.coords_to_str())

    #axs[1].set_ylabel(f"Autocorrelation ({what})")
    #axs[1].set_xlabel("Lag [m]")
    axs[1].set_xlabel("Wavenumber [Hz]")
    axs[1].set_title(f"Power Spectral Density ({what})")
    axs[0].grid()
    axs[1].grid()
    #axs[2].grid()
    axs[0].legend()
    axs[1].legend()
    #axs[2].legend()

    plt.tight_layout()
    plt.savefig(f"pic.png")
    plt.show()
