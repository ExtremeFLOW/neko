import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pylab as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
import matplotlib.ticker as tick
#rcParams['text.usetex'] = True
rcParams['lines.linewidth'] = 1.5
rcParams['font.size'] = 16
import os
import sys
import csv

from matplotlib.animation import FuncAnimation

pts_filename = "output.csv"
n_probes = 100
n_fields = 2
line_counter = 0

x = np.zeros((n_probes,1))
y = np.zeros((n_probes,1))

t = []
fields_t = []
j = -1

if __name__ == "__main__":

    f = open(pts_filename)
    csv_reader = csv.reader(f)
    for line in csv_reader:

        i=(np.mod(csv_reader.line_num-1,n_probes))
        
        if line_counter < n_probes:
            
            x[i,0] = float(line[0])
            y[i,0] = float(line[1])
        
        else:

            if i == 0:
                t.append(float(line[0]))
                fld = np.zeros((n_probes,2))
                fields_t.append(fld)
                j += 1

            fields_t[j][i,0] = float(line[1])
            fields_t[j][i,1] = float(line[2])
            
        line_counter +=1


    # plot the mesh
    fig, ax = plt.subplots(1, 1,figsize=(7, 7))
    plt.plot(x.flatten(),y.flatten(),'.b')
    plt.gca().set_aspect('equal', adjustable='box')
    #plt.show()

    # Plot a field
    for isnap in range(0, len(t)):
        xx = x.flatten()
        yy = y.flatten()
        tt = isnap
        print(tt)
        i_field = 1
        zz = fields_t[tt][:,i_field].flatten()
        if i_field == 0:
            cmapp='RdBu_r'
            levels = 100
        else:
            cmapp='magma'
            levels = 100
        fig, ax = plt.subplots(1, 1,figsize=(7, 7))
        c1 = ax.tricontourf(xx,yy,zz,levels,cmap=cmapp)
        cbar=fig.colorbar(c1, ax=ax)
        plt.gca().set_aspect('equal', adjustable='box')
        ax.set_ylabel(r'x',labelpad=10)
        ax.set_xlabel(r'y',labelpad=10)
        path = "pics/"
        figname="xy_plane_"+repr(isnap)+".png"
        plt.savefig(path+figname, format="png", bbox_inches="tight")



