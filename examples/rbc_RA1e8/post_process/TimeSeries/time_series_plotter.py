import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
import matplotlib.ticker as tick
rcParams['text.usetex'] = False
rcParams['lines.linewidth'] = 0.5
rcParams['font.size'] = 16

import numpy as np

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
        
            fig, ax = plt.subplots(4, 2,figsize=(10, 12), dpi=300)
            
            for ts in ts_dict[which_ra]:

                ra_string = '{:.0e}'.format(ts.Ra[-1])
                fig.suptitle(r'Time series for $Ra = $' + ra_string, fontsize=15)

                # Check if the time is always increasing, this could be broken if the run stopped suddenly and the time series was not cleaned.
                increasing = np.all(np.diff(ts.t) > 0)
                if increasing == False:
                    print("the time series is not strictly increasing, overwriting instances in which the time is larger than a subsequent simulation")
                    index_where_time_goes_back  = np.where(np.diff(ts.t) < 0)
                    print(index_where_time_goes_back[0][0])
                    index_that_should_be_deleted = np.where(ts.t[:index_where_time_goes_back[0][0]+1] >= ts.t[index_where_time_goes_back[0][0]+1])[0]
                    print("The following entries (rows) of the time series will be deleted to favor the newer simulation")
                    print(index_that_should_be_deleted)
                    #print(ts.t[index_that_should_be_deleted])
                    #print(ts.t[index_that_should_be_deleted[-1]:index_that_should_be_deleted[-1]+100])
                    ts.t = np.delete(ts.t, index_that_should_be_deleted, axis=0)
            
                # Check again if delete worked
                increasing = np.all(np.diff(ts.t) > 0)
                print(increasing)
            
                unique_lx = (np.unique(ts.lx))
                # Plot unique colors for different orders
                for lx_ in unique_lx:
                    ax[0,0].plot(ts.t[ts.lx == lx_] - ts.t[0], ts.nu_a[ts.lx == lx_])
                    #ax[0,0].plot(ts.t, ts.nu_v)
                    ax[0,0].set_xlabel(r'$t - t_{0}$')
                    ax[0,0].set_ylabel(r'$Nu_{\langle A \rangle}$')
                    #ax[0,0].set_xlim(50,600)
                    ax[0,0].set_ylim(min(ts.nu_eps_k)*0.9,max(ts.nu_eps_k)*1.1)
                            
                    #ax[0,1].plot(ts.t[ts.lx == lx_] - ts.t[0], ts.nu_v[ts.lx == lx_]/ts.nu_a[ts.lx == lx_],)
                    ax[0,1].plot(ts.t[ts.lx == lx_] - ts.t[0], ts.nu_v[ts.lx == lx_],)
                    ax[0,1].set_xlabel(r'$t - t_{0}$')
                    ax[0,1].set_ylabel(r'$Nu_{\langle V \rangle}$')
                    #ax[0,1].set_ylim(0.5,1.5)
                    ax[0,1].set_ylim(min(ts.nu_eps_k)*0.9,max(ts.nu_eps_k)*1.1)
            
                    ax[1,0].plot(ts.t[ts.lx == lx_] - ts.t[0], ts.nu_eps_t[ts.lx == lx_],)
                    ax[1,0].set_xlabel(r'$t - t_{0}$')
                    ax[1,0].set_ylabel(r'$eps_t$')
                    #ax[0,1].set_ylim(0.5,1.5)
                    ax[1,0].set_ylim(min(ts.nu_eps_k)*0.9,max(ts.nu_eps_k)*1.1)
            
                    ax[1,1].plot(ts.t[ts.lx == lx_] - ts.t[0], ts.nu_eps_k[ts.lx == lx_],)
                    ax[1,1].set_xlabel(r'$t - t_{0}$')
                    ax[1,1].set_ylabel(r'$eps_k$')
                    #ax[0,1].set_ylim(0.5,1.5)
                    ax[1,1].set_ylim(min(ts.nu_eps_k)*0.9,max(ts.nu_eps_k)*1.1)
                            
                    ax[2,0].plot(ts.t[ts.lx == lx_] - ts.t[0], ts.tke[ts.lx == lx_],)
                    ax[2,0].set_xlabel(r'$t - t_{0}$')
                    ax[2,0].set_ylabel(r'$tke$')
                    ax[2,0].set_ylim(min(ts.tke[ts.tke>0])*0.9,max(ts.tke[ts.tke>0])*1.1)
            
                    ax[2,1].axis('off')
                    
                    ax[3,0].plot(ts.t[ts.lx == lx_] - ts.t[0], ts.lx[ts.lx == lx_])
                    ax[3,0].set_xlabel(r'$t - t_{0}$')
                    ax[3,0].set_ylabel(r'$lx$')
                        
                    ax[3,1].plot(ts.t[ts.lx == lx_] - ts.t[0], ts.Ra[ts.lx == lx_],)
                    ax[3,1].set_xlabel(r'$t - t_{0}$')
                    ax[3,1].set_ylabel(r'$Ra$')

            if plot_transient_end==1:
                ax[0,0].axvline(x=ts.t[ts.nu_a_transient_index]-ts.t[0], linestyle=':', color= "red", linewidth = 1)
                ax[0,1].axvline(x=ts.t[ts.nu_v_transient_index]-ts.t[0], linestyle=':', color= "red", linewidth = 1)
                ax[1,0].axvline(x=ts.t[ts.nu_eps_t_transient_index]-ts.t[0], linestyle=':', color= "red", linewidth = 1)
                ax[1,1].axvline(x=ts.t[ts.nu_eps_k_transient_index]-ts.t[0], linestyle=':', color= "red", linewidth = 1)
                # Correct for tke index in case some files did not have tke on them
                tke_real_transient_index = ts.tke_transient_index
                if len(np.where(ts.tke == 0)[0]) > 0:
                    tke_real_transient_index += np.where(ts.tke == 0)[0][-1]
                ax[2,0].axvline(x=ts.t[tke_real_transient_index]-ts.t[0], linestyle=':', color= "red", linewidth = 1)
                    
                
            ax[0,0].grid(color = 'black', linestyle = '-', linewidth = 0.1)
            ax[0,1].grid(color = 'black', linestyle = '-', linewidth = 0.1)
            ax[1,0].grid(color = 'black', linestyle = '-', linewidth = 0.1)
            ax[1,1].grid(color = 'black', linestyle = '-', linewidth = 0.1)
            ax[2,0].grid(color = 'black', linestyle = '-', linewidth = 0.1)
            ax[2,1].grid(color = 'black', linestyle = '-', linewidth = 0.1)
            ax[3,0].grid(color = 'black', linestyle = '-', linewidth = 0.1)
            ax[3,1].grid(color = 'black', linestyle = '-', linewidth = 0.1)
            plt.tight_layout()
            #plt.savefig("nu_timeseries.pdf", format="pdf", bbox_inches="tight")
            plt.show()
    
        return