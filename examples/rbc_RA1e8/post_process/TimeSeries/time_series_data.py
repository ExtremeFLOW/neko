class time_series_data_c():
    
    def __init__(self, dat):

        # Assign the values from the data set
        self.t = dat[:,0]
        self.uzt = dat[:,1]
        self.dtdz_top = abs(dat[:,2])
        self.dtdz_bot = abs(dat[:,3])
        self.lx = dat[:,4]
        self.Ra = dat[:,5]

        ra = self.Ra[0]
        self.nu_v = 1 + np.sqrt(ra)*self.uzt
        self.nu_a = (self.dtdz_top + self.dtdz_bot)/2

        return        