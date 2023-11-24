import numpy as np
from stat_tools import adf_test
from stat_tools import get_lag1_cor
from stat_tools import get_SME_var


class time_series_data_c():
    
    def __init__(self, dat):

        # Assign the values from the data set
        self.t = dat[:,0]
        self.uzt = dat[:,1]
        self.dtdz_top = abs(dat[:,2])
        self.dtdz_bot = abs(dat[:,3])
        self.lx = dat[:,4]
        self.Ra = dat[:,5]
        self.nu_eps_t = dat[:,6]
        self.nu_eps_k = dat[:,7]
        self.tke = dat[:,8]

        ra = self.Ra[0]
        self.nu_v = 1 + np.sqrt(ra)*self.uzt
        self.nu_a = (self.dtdz_top + self.dtdz_bot)/2

        return 


    def get_transients(self, granularity_for_iterations):

        self.nu_v_transient_index = adf_test(self.nu_v, granularity_for_iterations)
        self.nu_a_transient_index = adf_test(self.nu_a, granularity_for_iterations)
        self.nu_va_transient_index = adf_test(self.nu_v/self.nu_a, granularity_for_iterations)
        self.nu_eps_t_transient_index = adf_test(self.nu_eps_t, granularity_for_iterations)
        self.nu_eps_k_transient_index = adf_test(self.nu_eps_k, granularity_for_iterations)
        self.tke_transient_index = adf_test(self.tke[np.where(self.tke > 0)], granularity_for_iterations)
        #self.tke_transient_index = self.tke_transient_index + np.where(self.tke[self.tke<=0])[0][-1]

        return


    def get_uncorrelated_batch(self, plot):

        self.nu_v_batch = get_lag1_cor(self.t[self.nu_v_transient_index:],self.nu_v[self.nu_v_transient_index:], len(self.nu_v[self.nu_v_transient_index:]), plot)

        self.nu_a_batch = get_lag1_cor(self.t[self.nu_a_transient_index:],self.nu_a[self.nu_a_transient_index:], len(self.nu_a[self.nu_a_transient_index:]), plot)

        self.nu_eps_t_batch = get_lag1_cor(self.t[self.nu_eps_t_transient_index:],self.nu_eps_t[self.nu_eps_t_transient_index:], len(self.nu_eps_t[self.nu_eps_t_transient_index:]), plot)

        self.nu_eps_k_batch = get_lag1_cor(self.t[self.nu_eps_k_transient_index:],self.nu_eps_k[self.nu_eps_k_transient_index:], len(self.nu_eps_k[self.nu_eps_k_transient_index:]), plot)

       
    def get_stats(self):
        
        mean, var, ci = get_SME_var(self.t[self.nu_v_transient_index:],self.nu_v[self.nu_v_transient_index:], self.nu_v_batch)
        self.stats_nu_v = [mean, var, ci]

        mean, var, ci = get_SME_var(self.t[self.nu_a_transient_index:],self.nu_a[self.nu_a_transient_index:], self.nu_a_batch)
        self.stats_nu_a = [mean, var, ci]

        mean, var, ci = get_SME_var(self.t[self.nu_eps_t_transient_index:],self.nu_eps_t[self.nu_eps_t_transient_index:], self.nu_eps_t_batch)
        self.stats_nu_eps_t = [mean, var, ci]

        mean, var, ci = get_SME_var(self.t[self.nu_eps_k_transient_index:],self.nu_eps_k[self.nu_eps_k_transient_index:], self.nu_eps_k_batch)
        self.stats_nu_eps_k = [mean, var, ci]
