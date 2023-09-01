import numpy as np

import matplotlib.pylab as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
import matplotlib.ticker as tick
#rcParams['text.usetex'] = True
rcParams['lines.linewidth'] = 1.5
rcParams['font.size'] = 16


filename = '../../data/bc_heat_balance.txt'
ra = 1e11

if __name__ == "__main__":


    data = np.loadtxt(filename)
    t = data[:,0]
    dtdn_top = data[:,1]
    dtdn_bot = data[:,2]
    dtdn_sid = data[:,3]
    dive = data[:,5]
  
    print(dive)

    total = dtdn_top+dtdn_bot+dtdn_sid

    fig, ax = plt.subplots(1, 1,figsize=(3, 3), dpi=100)
    ax.plot(t, dtdn_bot,'-r',label = r'dtdn_bot')
    ax.plot(t, dtdn_top,'-b',label = r'dtdn_top')
    ax.plot(t, dtdn_sid,'-g',label = r'dtdn_side')
    ax.plot(t, total,'-k',label = r'total')
    ax.plot(t, dive,'--b',label = r'div(heat)')
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$dtdn$')
    plt.legend(loc='best')
    plt.show()

