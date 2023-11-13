import numpy as np
from statsmodels.tsa.stattools import adfuller

# UQit time series module
import os
import sys
import time
UQit_ts_path=os.getcwd()+"/../UQ/UQit_ts_module/UQit-ts/"  
print("Using UQit in: " + UQit_ts_path)
sys.path.append(UQit_ts_path)
from SME import SME, SMEuncert
from NOBM import NOBM
from BMBC import BMBC
import math as mt
import plot_ts
import stats_ts

#Ploting
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
import matplotlib.ticker as tick
rcParams['text.usetex'] = False
rcParams['lines.linewidth'] = 0.5
rcParams['font.size'] = 16

from tqdm import tqdm
from scipy.signal import argrelextrema

def adf_test(variable, granularity_for_iterations):
        #Augmented Dickey Fuller test
        m = granularity_for_iterations # Granulatity for iterations
        initial_transient = 0
        num_iter = int(np.ceil(len(variable)/m))
        tolerance = 0.01

            
        pbar= tqdm(total=num_iter)
        for i in range(0,num_iter):
            pvalue = adfuller(variable[:(i+1)*m])

            pbar.update(1)            
            if pvalue[1] < tolerance: # Rejct the hypotyesis, which means that the series is stationary
                initial_transient = (i+1)*m
                pbar.close()
                break

        if initial_transient == 0:
            pbar.close()
    
        return initial_transient


def get_lag1_cor(t,fx,max_batch, plot):
    
    batchsizes = np.arange(5000, max_batch, 5000)
    rho1_NOBM = np.zeros(len(batchsizes))
    var_NOBM = np.zeros(len(batchsizes))
    trueMean = np.mean(fx)
    sme_nobm = SME(t, fx, {'verbose':False,'conv':True})

    pbar= tqdm(total=len(batchsizes))
    for j in range(len(batchsizes)):
        M = batchsizes[j]
        K1, tnobm, bmeans = NOBM.batchGen(t, fx, M)
        bmeans_mod = bmeans - trueMean
        s1 = np.sum(bmeans_mod[:-1] * bmeans_mod[1:])
        s0 = np.sum(bmeans_mod * bmeans_mod)
        rho1_NOBM[j] = np.divide(s1,s0)
        #out_nobm = SMEuncert(sme_nobm,{'method':'NOBM','batchSize':M}).estim
        #var_NOBM[j] = out_nobm['fSME_var'][-1]

        pbar.update(1)
    pbar.close()

    if plot==1:
        plt.figure(figsize=(8, 5))
        plt.plot(batchsizes, abs(rho1_NOBM),'dodgerblue', linewidth=4, marker='o', markeredgecolor='black', markerfacecolor='lime')
        plt.axhline(0, color='k', linestyle=':', lw=4, label=r'$\rho_1$ = 0')
        plt.xlabel('M: the number of samples in each batch', fontsize=15)
        plt.ylabel('lag-1 autocorrelation for batch means series', fontsize=15)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(loc='best')
        plt.show()


    # Find indices of local minima
    loc_min = (argrelextrema(abs(rho1_NOBM), np.less))[0]
    print(loc_min)
    if len(loc_min)==0:
        print("Warning!")
        print("Could not find local minima for autocorrelation")
        time.sleep(5)
        uncorrelated_batch_size = batchsizes[0]
        print("Continuing with first batch size (10)")
    else:    
        uncorrelated_batch_size = batchsizes[loc_min[0]]
      
    return uncorrelated_batch_size

def get_SME_var(t,fy,Mb):

    fy_SME = np.mean(fy)
    sme_y  = SME(t, fy, {'verbose':False,'conv':True})
    outvar = SMEuncert(sme_y,{'method':'NOBM','batchSize':Mb}).estim
    varnobm= outvar['fSME_var'][-1]  # last value is the final variance of all samples. This is the correct one
    cbi_SME=1.96*mt.sqrt(abs(varnobm))

    return fy_SME, varnobm, cbi_SME
