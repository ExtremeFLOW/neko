import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
import matplotlib.ticker as tick
rcParams['text.usetex'] = False
rcParams['lines.linewidth'] = 0.5
rcParams['font.size'] = 16

class time_series_plotter_c():
    
    def __init__(self):

        return

    def plot_time_series(self,ts):

        t = ts.t
        Nu_a = ts.nu_a
        Nu_v = ts.nu_v
        Nu_va = Nu_v/Nu_a
        
        fig, ax = plt.subplots(2, 2,figsize=(10, 6), dpi=300)
        
        ax[0,0].plot(t, Nu_a)
        ax[0,0].plot(t, Nu_v)
        ax[0,0].set_xlabel(r'$t$')
        ax[0,0].set_ylabel(r'$Nu_{\langle A \rangle}$')
        
        ax[0,1].plot(t, Nu_va,)
        ax[0,1].set_xlabel(r'$t$')
        ax[0,1].set_ylabel(r'$Nu_{\langle V \rangle}/Nu_{\langle A \rangle}$')

        ax[1,0].plot(t, ts.lx)
        ax[1,0].set_xlabel(r'$t$')
        ax[1,0].set_ylabel(r'$lx$')
        
        ax[1,1].plot(t, ts.Ra,)
        ax[1,1].set_xlabel(r'$t$')
        ax[1,1].set_ylabel(r'$Ra$')
        
        plt.tight_layout()
        plt.savefig("nu_timeseries.pdf", format="pdf", bbox_inches="tight")
        plt.show()

        return