import numpy as np

import matplotlib.pylab as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
import matplotlib.ticker as tick
#rcParams['text.usetex'] = True
rcParams['lines.linewidth'] = 1.5
rcParams['font.size'] = 16


filename = '../nusselt.txt'
ra = 1e8

if __name__ == "__main__":


    data = np.loadtxt(filename)
    t = data[:,0]
    uzt = data[:,1]
    dtdz_top = abs(data[:,2])
    dtdz_bot = abs(data[:,3])
    

    Nu_v = 1 + np.sqrt(ra)*uzt
    print(Nu_v)
    Nu_a = (dtdz_top + dtdz_bot)/2
    print(Nu_a)

    fig, ax = plt.subplots(1, 1,figsize=(3, 3), dpi=100)
    ax.plot(t, Nu_v,'-b',label = r'Nu_V')
    ax.plot(t, Nu_a,'-r',label = r'Nu_A')
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$Nu$')
    plt.legend(loc='best')
    plt.show()

