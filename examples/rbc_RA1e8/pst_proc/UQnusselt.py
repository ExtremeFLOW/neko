import numpy as np

import matplotlib.pylab as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
import matplotlib.ticker as tick
#rcParams['text.usetex'] = True
rcParams['lines.linewidth'] = 1.5
rcParams['font.size'] = 16
import os
import sys



#UQit related
UQit_ts_path=os.getcwd()+"/UQit_ts_module/UQit-ts/"  
sys.path.append(UQit_ts_path)
from SME import SME, SMEuncert
from NOBM import NOBM
from BMBC import BMBC
import math as mt
import plot_ts
import stats_ts

def get_lag1_cor(uqt,uqNu_v,max_batch_size):
    t = uqt
    fx = uqNu_v
    max_batch = max_batch_size
        
    batchsizes = np.arange(10, max_batch, 10)
    print(batchsizes)
    rho1_NOBM = np.zeros(len(batchsizes))
    var_NOBM = np.zeros(len(batchsizes))
    trueMean = np.mean(fx)
    sme_nobm = SME(t, fx, {'verbose':False,'conv':True})

    for j in range(len(batchsizes)):
        M = batchsizes[j]
        K1, tnobm, bmeans = NOBM.batchGen(t, fx, M)
        bmeans_mod = bmeans - trueMean
        s1 = np.sum(bmeans_mod[:-1] * bmeans_mod[1:])
        s0 = np.sum(bmeans_mod * bmeans_mod)
        rho1_NOBM[j] = np.divide(s1,s0)
        out_nobm = SMEuncert(sme_nobm,{'method':'NOBM','batchSize':M}).estim
        var_NOBM[j] = out_nobm['fSME_var'][-1]
    
    plt.figure(figsize=(8, 5))
    plt.plot(batchsizes, rho1_NOBM,'dodgerblue', linewidth=4, marker='o', markeredgecolor='black', markerfacecolor='lime')
    plt.axhline(0, color='k', linestyle=':', lw=4, label=r'$\rho_1$ = 0')
    plt.xlabel('M: the number of samples in each batch', fontsize=15)
    plt.ylabel('lag-1 autocorrelation for batch means series', fontsize=15)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc='best')
    plt.show()
    
    return

def get_SME_var(uqt,uqNu_v,Mb):

    t = uqt
    fy = uqNu_v 

    fy_SME = np.mean(fy)
    sme_y  = SME(t, fy, {'verbose':False,'conv':True})
    outvar = SMEuncert(sme_y,{'method':'NOBM','batchSize':Mb}).estim
    varnobm= outvar['fSME_var'][-1]  # last value is the final variance of all samples. This is the correct one
    cbi_SME=1.96*mt.sqrt(abs(varnobm))

    print(fy_SME)
    print(varnobm)
    print(cbi_SME)
    print("===")

    return

filename = '../nusselt.txt'
ra = 1e8
start_time = 0
if __name__ == "__main__":

    data = np.loadtxt(filename)
    t = data[:,0]
    uzt = data[:,1]
    dtdz_top = abs(data[:,2])
    dtdz_bot = abs(data[:,3])
    
    Nu_v = 1 + np.sqrt(ra)*uzt
    Nu_a = (dtdz_top + dtdz_bot)/2


    fig, ax = plt.subplots(1, 1,figsize=(8, 5), dpi=100)
    ax.plot(t, Nu_v,'-b',label = r'Nu_V')
    ax.plot(t, Nu_a,'-r',label = r'Nu_A')
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$Nu$')
    plt.legend(loc='best')
    plt.show()
    
    # Crop the time series to remove any transient
    uqt = t[t>=start_time]
    uqNu_v = Nu_v[t>=start_time]
    uqNu_a = Nu_a[t>=start_time]
    series_len = len(uqt)
    max_batch_size = int(series_len/1)

    # Get correlations to find best batch sizes for the SME
    get_lag1_cor(uqt,uqNu_v,max_batch_size)
    get_lag1_cor(uqt,uqNu_a,max_batch_size)
    
    # Get sample mean estimators and variance
    Mb = 5
    get_SME_var(uqt,uqNu_v,Mb)
    get_SME_var(uqt,uqNu_a,Mb)
