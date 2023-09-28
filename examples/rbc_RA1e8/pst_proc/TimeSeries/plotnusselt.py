import numpy as np

import matplotlib.pylab as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
import matplotlib.ticker as tick
rcParams['text.usetex'] = True
rcParams['lines.linewidth'] = 1.5
rcParams['font.size'] = 16


filename = ['nu_lx8.txt', 'nu_lx6.txt']
labels = [r'case 1', r'case 2']
colors = ['b','r']
ra = 1e11
plot_legend = False

if __name__ == "__main__":
    fig, ax = plt.subplots(1, 1,figsize=(5, 3), dpi=300)
    for filen in range(0, len(filename)):

        data = np.loadtxt(filename[filen])
        t = data[:,0]
        uzt = data[:,1]
        dtdz_top = abs(data[:,2])
        dtdz_bot = abs(data[:,3])
    

        Nu_v = 1 + np.sqrt(ra)*uzt
        #print(Nu_v)
        Nu_a = (dtdz_top + dtdz_bot)/2
        #print(Nu_a)


        #ax.plot(t, Nu_v,'--'+colors[filen],label = r'$Nu_V-'+labels[filen]+'$')
        ax.plot(t, Nu_a,'-'+colors[filen],label = r'$Nu_A-'+labels[filen]+'$')
        ax.set_xlabel(r'$t$')
        ax.set_ylabel(r'$Nu$')
        if plot_legend: plt.legend(loc='best')
    
    plt.tight_layout()
    plt.savefig("nu_timeseries.pdf", format="pdf", bbox_inches="tight")
    plt.show()

