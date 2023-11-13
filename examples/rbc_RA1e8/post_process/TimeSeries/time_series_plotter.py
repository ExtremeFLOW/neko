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


    def plot_ra_separate(self,ts_dict, keys, plot_transient_end):
        
        for which_ra in keys:
        
            print("------ Now writing for Ra = "+repr(which_ra)+ "---------------")
            
            fig, ax = plt.subplots(2, 2,figsize=(10, 6), dpi=300)
            i=-1
            for ts in ts_dict[which_ra]:
                
                i+=1
            
                ax[0,0].plot(ts.t, ts.nu_a)
                #ax[0,0].plot(ts.t, ts.nu_v)
                ax[0,0].set_xlabel(r'$t$')
                ax[0,0].set_ylabel(r'$Nu_{\langle A \rangle}$')
                if plot_transient_end:
                    ax[0,0].axvline(x=ts.t[ts.nu_a_transient_index], linestyle=':')
                    
                ax[0,1].plot(ts.t, ts.nu_v/ts.nu_a,)
                ax[0,1].set_xlabel(r'$t$')
                ax[0,1].set_ylabel(r'$Nu_{\langle V \rangle}/Nu_{\langle A \rangle}$')
                ax[0,1].set_ylim(0.5,1.5)
            
                ax[1,0].plot(ts.t, ts.lx)
                ax[1,0].set_xlabel(r'$t$')
                ax[1,0].set_ylabel(r'$lx$')
                    
                ax[1,1].plot(ts.t, ts.Ra,)
                ax[1,1].set_xlabel(r'$t$')
                ax[1,1].set_ylabel(r'$Ra$')
            
            plt.tight_layout()
            #plt.savefig("nu_timeseries.pdf", format="pdf", bbox_inches="tight")
            plt.show()
    
        return