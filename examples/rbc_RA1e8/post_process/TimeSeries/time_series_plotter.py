class time_series_plotter_c():
    
    def __init__(self):

    def plot_time_series(ts):

        t = ts.t
        Nu_a = ts.nu_a
        Nu_v = ts.nu_v
        Nu_va = Nu_v/Nu_a
        
        fig, ax = plt.subplots(1, 2,figsize=(10, 3), dpi=300)
        ax[0].plot(t, Nu_a)
        ax[0].plot(t, Nu_v)
        ax[0].set_xlabel(r'$t$')
        ax[0].set_ylabel(r'$Nu_{\langle A \rangle}$')
        
        ax[1].plot(t, Nu_va,)
        ax[1].set_xlabel(r'$t$')
        ax[1].set_ylabel(r'$Nu_{\langle V \rangle}/Nu_{\langle A \rangle}$')
        
        plt.tight_layout()
        plt.savefig("nu_timeseries.pdf", format="pdf", bbox_inches="tight")
        plt.show()