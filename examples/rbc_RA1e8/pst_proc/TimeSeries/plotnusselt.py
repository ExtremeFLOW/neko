import numpy as np

import matplotlib.pylab as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
import matplotlib.ticker as tick
from statsmodels.tsa.stattools import adfuller
rcParams['text.usetex'] = True
rcParams['lines.linewidth'] = 0.5
rcParams['font.size'] = 16

def intitial_transient_index(Nu_a):
    #Augmented Dickey Fuller test
    m = 30 # Granulatity for iterations
    initial_transient = 0
    num_iter = int(np.ceil(len(Nu_a)/m))
    tolerance = 0.01

    for i in range(0,num_iter):
            
        pvalue = adfuller(Nu_a[:(i+1)*m])
        if pvalue[1] < tolerance:
            initial_transient = (i+1)*m
            break

    return initial_transient



filename = ['nu_lx8.txt', 'nu_lx6.txt']
labels = [r'case 1', r'case 2']
colors = ['b','r']
colors2 = ['k','g']
ra = 1e11
plot_legend = False

if __name__ == "__main__":
    fig, ax = plt.subplots(1, 2,figsize=(10, 3), dpi=300)
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

        #Get the initial transient with the Augmented Dickey fuller test
        it_nua = intitial_transient_index(Nu_a)
        it_nuv = intitial_transient_index(Nu_v)
        it_nuva = intitial_transient_index(Nu_v/Nu_a)


        #ax.plot(t, Nu_v,'-'+colors2[filen],label = r'$Nu_V-'+labels[filen]+'$')
        ax[0].plot(t, Nu_a,'-'+colors[filen],label = r'$Nu_A-'+labels[filen]+'$')
        ax[0].set_xlabel(r'$t$')
        ax[0].set_ylabel(r'$Nu_{\langle A \rangle}$')
        ax[0].axvline(x=t[it_nua], color=colors[filen], linestyle=':')

        ax[1].plot(t, Nu_v/Nu_a,'-'+colors[filen],label = r'$Nu_A-'+labels[filen]+'$')
        ax[1].set_xlabel(r'$t$')
        ax[1].set_ylabel(r'$Nu_{\langle V \rangle}/Nu_{\langle A \rangle}$')
        ax[1].set_ylim(0,4)
        ax[1].axvline(x=t[it_nuva], color=colors[filen], linestyle=':')
        

        if plot_legend: plt.legend(loc='best')
    

    plt.tight_layout()
    plt.savefig("nu_timeseries.pdf", format="pdf", bbox_inches="tight")
    plt.show()

