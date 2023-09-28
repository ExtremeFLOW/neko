import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pylab as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
import matplotlib.ticker as tick
from statsmodels.tsa.stattools import adfuller
rcParams['text.usetex'] = True
rcParams['lines.linewidth'] = 0.5
rcParams['font.size'] = 16
import os
import sys
import time

#From local Uqit folder
UQit_ts_path=os.getcwd()+"/UQit_ts_module/UQit-ts/"  
sys.path.append(UQit_ts_path)
from SME import SME, SMEuncert
from NOBM import NOBM
from BMBC import BMBC
import math as mt
import plot_ts
import stats_ts

# From Uqit pip package
from UQit.sobol import sobol
import UQit.analyticTestFuncs as analyticTestFuncs
from UQit.pce import pce, pceEval
import UQit.reshaper as reshaper
import UQit.sampling as sampling

## =========
## Functions
## =========

def get_lag1_cor(uqt,uqNu_v,max_batch_size):
    t = uqt
    fx = uqNu_v
    max_batch = max_batch_size
        
    batchsizes = np.arange(10, max_batch, 10)
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


    if 1==0:
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
    #print(loc_min[0])
    #print(batchsizes[loc_min[0]])
    if len(loc_min)==0:
        print("Warning!")
        print("Could not find local minima for autocorrelation")
        time.sleep(5)
        uncorrelated_batch_size = batchsizes[0]
        print("Continuing with first batch size (10)")
    else:    
        uncorrelated_batch_size = batchsizes[loc_min[0]]
    

    
    return uncorrelated_batch_size

def get_SME_var(uqt,uqNu_v,Mb):

    t = uqt
    fy = uqNu_v 


    fy_SME = np.mean(fy)
    sme_y  = SME(t, fy, {'verbose':False,'conv':True})
    outvar = SMEuncert(sme_y,{'method':'NOBM','batchSize':Mb}).estim
    varnobm= outvar['fSME_var'][-1]  # last value is the final variance of all samples. This is the correct one
    cbi_SME=1.96*mt.sqrt(abs(varnobm))

    return fy_SME, varnobm, cbi_SME


def intitial_transient_index(t,Nu_a):
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

    return initial_transient, t[initial_transient:], Nu_a[initial_transient:] 


def process_Nu_timeseries(ra, filename_list): 

    t_v_ts=[]
    t_a_ts=[]
    t_va_ts=[]
    Nu_v_ts=[]
    Nu_a_ts=[]
    Nu_va_ts=[]

    it_v=[]
    it_a=[]
    it_va=[]

    Nu_v_UQ = []
    Nu_a_UQ = []
    Nu_va_UQ = []

    for i in range(0, len(filename_list)):
 
        print("===================")
        print("Iteration: "+repr(i))

        filename = filename_list[i]
    
        data = np.loadtxt(filename)
        t = data[:,0]
        uzt = data[:,1]
        dtdz_top = abs(data[:,2])
        dtdz_bot = abs(data[:,3])
    
        Nu_v = 1 + np.sqrt(ra)*uzt
        Nu_a = (dtdz_top + dtdz_bot)/2
        Nu_va = Nu_v/Nu_a

        # Store the time series
        t_v_ts.append(t)
        t_a_ts.append(t)
        t_va_ts.append(t)
        Nu_v_ts.append(Nu_v)
        Nu_a_ts.append(Nu_a)
        Nu_va_ts.append(Nu_va)

        # Determine where the initial transient ends
        it_Nu_v, uqt_v, uqNu_v = intitial_transient_index(t,Nu_v)
        it_Nu_a, uqt_a, uqNu_a = intitial_transient_index(t,Nu_a)
        it_Nu_va, uqt_va, uqNu_va = intitial_transient_index(t,Nu_va)
        
        it_v.append(it_Nu_v)
        it_a.append(it_Nu_a)
        it_va.append(it_Nu_va)

        # Get correlations to find best batch sizes for the SME
        Nu_v_batch_size = get_lag1_cor(uqt_v,uqNu_v,len(uqt_v))
        Nu_a_batch_size = get_lag1_cor(uqt_a,uqNu_a,len(uqt_a))
        Nu_va_batch_size = get_lag1_cor(uqt_va,uqNu_va,len(uqt_va))
        
    
        # Get sample mean estimators and variance
        Nu_v_mean, Nu_v_var, Nu_v_ci = get_SME_var(uqt_v,uqNu_v,Nu_v_batch_size)
        Nu_a_mean, Nu_a_var, Nu_a_ci = get_SME_var(uqt_a,uqNu_a,Nu_a_batch_size)  
        Nu_va_mean, Nu_va_var, Nu_va_ci = get_SME_var(uqt_va,uqNu_va,Nu_va_batch_size)  
        

        print(Nu_v_mean)
        print(Nu_a_mean)
        print(Nu_va_mean)

        Nu_v_UQ.append((Nu_v_mean, Nu_v_var, Nu_v_ci))
        Nu_a_UQ.append((Nu_a_mean, Nu_a_var, Nu_a_ci))
        Nu_va_UQ.append((Nu_va_mean, Nu_va_var, Nu_va_ci))

    return t_v_ts, t_a_ts, t_va_ts, Nu_v_ts, Nu_a_ts, Nu_va_ts, it_v, it_a, it_va, Nu_v_UQ, Nu_a_UQ, Nu_va_UQ


## ======================================================
## Code Set-up
## ======================================================


filename_list = ['nu_lx8.txt', 'nu_lx6.txt']
labels = [r'case 1', r'case 2']
colors = ['b','r']
colors2 = ['k','g']
ra = 1e11
plot_legend = False

## ===========
## Main code
## ===========

if __name__ == "__main__":

    #> Time series analysis - Get Nusselt and CI
    t_v_ts, t_a_ts, t_va_ts, Nu_v_ts, Nu_a_ts, Nu_va_ts, it_v, it_a, it_va, Nu_v_UQ, Nu_a_UQ, Nu_va_UQ = process_Nu_timeseries(ra, filename_list)

    print(Nu_v_UQ) 
    print(Nu_a_UQ) 
    print(Nu_va_UQ) 
    
    fig, ax = plt.subplots(1, 2,figsize=(10, 3), dpi=300)
    for i in range(0, len(filename_list)):
        
        ax[0].plot(t_a_ts[i], Nu_a_ts[i],'-'+colors[i],label = r'$Nu_A-'+labels[i]+'$')
        ax[0].set_xlabel(r'$t$')
        ax[0].set_ylabel(r'$Nu_{\langle A \rangle}$')
        ax[0].axvline(x=t_a_ts[i][it_a[i]], color=colors[i], linestyle=':')

        ax[1].plot(t_va_ts[i], Nu_va_ts[i],'-'+colors[i],label = r'$Nu_A-'+labels[i]+'$')
        ax[1].set_xlabel(r'$t$')
        ax[1].set_ylabel(r'$Nu_{\langle V \rangle}/Nu_{\langle A \rangle}$')
        ax[1].set_ylim(0,4)
        ax[1].axvline(x=t_va_ts[i][it_va[i]], color=colors[i], linestyle=':')
        

        if plot_legend: plt.legend(loc='best')
    

    plt.tight_layout()
    plt.savefig("nu_timeseries.pdf", format="pdf", bbox_inches="tight")
    plt.show()

