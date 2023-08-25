









import numpy as np
import sys
import matplotlib.pylab as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
import matplotlib.ticker as tick
rcParams['text.usetex'] = True
rcParams['lines.linewidth'] = 1.5
rcParams['font.size'] = 16

#####################################
# Case parameters
## Number of snapshots
nsnap=200
## Batch size
p=200 # Number of snapshots in a batch
## Pod of which quantity
qoi=2 #0=vx, 1=vy, 2=vz

# Data path
casename='rbc'
path='/scratch/adperez/200-streaming_paralell_codes/data/'+casename+'/'

# number of modes
k=200
maxk = 1000
# Options if autodecide number of modes
## Autoupdate?
ifautoupdate=False
## If autoupdate has been done, gather ALL modes?
ifget_all_modes=False
## Minimun percentage of new snapshots that need to be parallel
min_ratio=0.01
str_mr = str(min_ratio).replace('.','')
#str_mr='0.00'

# Update gbl or lcl?
ifgbl_update=True
num_recs = 10

#### Automatic, do not touch #####
if ifgbl_update==True: 
    outprefix='gbl'
else:
    outprefix='lcl'
#################################
# Overwrite the hyperparameters if command line arguments were given
if len(sys.argv) > 1:
    p=int(sys.argv[1])
    k=int(sys.argv[2])
    min_ratio=float(sys.argv[3])
    str_mr = str(min_ratio).replace('.','')
    ifgbl_update=bool(int(sys.argv[4]))
    ifautoupdate=bool(int(sys.argv[5]))

if ifgbl_update==True: 
    if ifautoupdate==True:
        outprefix='auto_gbl_p'+repr(p)+'_k'+repr(k)+'_mr'+str_mr
    else:
        outprefix='fixd_gbl_p'+repr(p)+'_k'+repr(k)+'_mr'+str_mr
else:
    if ifautoupdate==True:
        outprefix='auto_lcl_p'+repr(p)+'_k'+repr(k)+'_mr'+str_mr
    else:
        outprefix='fixd_lcl_p'+repr(p)+'_k'+repr(k)+'_mr'+str_mr

################################

if __name__ == "__main__":

    sigma_outprefix = path+'PODsig_'+outprefix+'_'+casename+'.npy'
    vtr_outprefix = path+'PODvtr_'+outprefix+'_'+casename+'.npy'
   
    #> Read and plot the singular values
    scaling = 1e3 #This scaling depens on the value that is in rayleigh.f90 for the mass matrix
    sigma = np.load(sigma_outprefix) / scaling
    ## Plot
    fig, ax = plt.subplots(1, 1,figsize=(5, 5), dpi=100)
    ax.plot(sigma,'-b',label = r'$\sigma$')
    ax.set_xlabel(r'$mod$')
    ax.set_ylabel(r'$\sigma$')
    plt.xscale('log')
    plt.yscale('log')
    #plt.legend(loc='best')
    plt.show()

    #> Cummulative sum to determine the energy percentage
    sigma_cumsum=np.cumsum(sigma)
    sigma_sum=np.sum(sigma)
    ## Plot
    fig, ax = plt.subplots(1, 1,figsize=(3, 3), dpi=100)
    ax.plot(sigma_cumsum/sigma_sum,'-b',label = r'$E$')
    ax.set_xlabel(r'$mod$')
    ax.set_ylabel(r'$E$')
    #plt.legend(loc='best')
    plt.show()

    #> Read right singular vectors
    vtr = np.load(vtr_outprefix)
    ##Get time coefficients
    T = np.diag(sigma)@vtr
    ## Plot
    fig, ax = plt.subplots(1, 1,figsize=(3, 3), dpi=100)
    for i in range(0,5):
        label = r'mod=' + repr(i)
        ax.plot(T[i,:],label=label)
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$a_i$')
    plt.legend(loc='best')
    plt.show()







