import numpy as np
from scipy.signal import argrelextrema
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

    # Find indices of local minima
    loc_min = (argrelextrema(abs(rho1_NOBM), np.less))[0]
    uncorrelated_batch_size = batchsizes[loc_min[0]]
    
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


def process_Nu_timeseries(start_time, ra, filename_list, map_name_sizes, ifplot): 
    
    for i in range(0, len(filename_list)):
   
        filename = filename_list[i][:-1]

        dict_key = filename[3:6]
        print(dict_key)
        e_xy_l.append(map_name_sizes[dict_key][0])
        e_z_l.append(map_name_sizes[dict_key][1])
        lx_l.append(map_name_sizes[dict_key][2])
    
        data = np.loadtxt(filename)
        t = data[:,0]
        uzt = data[:,1]
        dtdz_top = abs(data[:,2])
        dtdz_bot = abs(data[:,3])
    
        Nu_v = 1 + np.sqrt(ra)*uzt
        Nu_a = (dtdz_top + dtdz_bot)/2

        if ifplot == True:
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
        Nu_v_batch_size = get_lag1_cor(uqt,uqNu_v,max_batch_size)
        Nu_a_batch_size = get_lag1_cor(uqt,uqNu_a,max_batch_size)
    
        # Get sample mean estimators and variance
        Nu_v_mean, Nu_v_var, Nu_v_ci = get_SME_var(uqt,uqNu_v,Nu_v_batch_size)
        Nu_a_mean, Nu_a_var, Nu_a_ci = get_SME_var(uqt,uqNu_a,Nu_a_batch_size)  


        print("elems xy = "+repr(e_xy_l[i]))
        print("elems z  = "+repr(e_z_l[i]))

        print("Batch size for Nu_v = " + repr(Nu_v_batch_size))
        print("Batch size for Nu_a = " + repr(Nu_a_batch_size))

        print('for i = ' + repr(i) + ' Nu_v_mean = ' + repr(Nu_v_mean) + ' diff = ' + repr(abs(Nu_v_mean - Nu_v_ref)) + ' CI =' + repr(Nu_v_ci) )
        print('for i = ' + repr(i) + ' Nu_a_mean = ' + repr(Nu_a_mean) + ' diff = ' + repr(abs(Nu_a_mean - Nu_a_ref)) + ' CI =' + repr(Nu_a_ci) )
        print("=================")


        Nu_v_mean_list.append(abs(Nu_v_mean-Nu_v_ref)/Nu_v_ref)
        Nu_a_mean_list.append(abs(Nu_a_mean-Nu_a_ref)/Nu_v_ref)

        
        Nu_v_dict[(e_xy_l[i],e_z_l[i], lx_l[i])] = Nu_v_mean 
        Nu_a_dict[(e_xy_l[i],e_z_l[i], lx_l[i])] = Nu_a_mean 

    if ifplot == True:
    
        x1 = np.array(e_xy_l)
        y1 = np.array(e_z_l)
        z1 = np.array(Nu_v_mean_list)
        cmapp='magma'
        fig, ax = plt.subplots(1, 1,figsize=(10, 5))
        ax2=ax
        c1 = ax.tricontourf(x1,y1,z1,50,cmap=cmapp)
        cbar=fig.colorbar(c1, ax=ax)
        cbar.set_label(r'$rel Nu_v error$',rotation=0,labelpad=40)
        ax.set_ylabel(r'e_z',labelpad=10)
        ax.set_xlabel(r'e_xy',labelpad=10)
        # Razterize
        for c in c1.collections:
            c.set_edgecolor("face")
            c.set_rasterized(True)
        plt.show()

        x1 = np.array(e_xy_l)
        y1 = np.array(e_z_l)
        z1 = np.array(Nu_a_mean_list)
        cmapp='magma'
        fig, ax = plt.subplots(1, 1,figsize=(10, 5))
        ax2=ax
        c1 = ax.tricontourf(x1,y1,z1,50,cmap=cmapp)
        cbar=fig.colorbar(c1, ax=ax)
        cbar.set_label(r'$Nu_a error$',rotation=0,labelpad=40)
        ax.set_ylabel(r'e_z',labelpad=10)
        ax.set_xlabel(r'e_xy',labelpad=10)
        # Razterize
        for c in c1.collections:
            c.set_edgecolor("face")
            c.set_rasterized(True)
        plt.show()

        fig, ax = plt.subplots(1, 1,figsize=(10, 5))
        ax.scatter(e_xy_l,e_z_l)   
        plt.show()


    return e_xy_l,e_z_l,lx_l,Nu_v_dict, Nu_a_dict


def GSA_3variables(exp_params1, exp_params2, exp_params3, Nu_dict, Nu_ref):

    #Experimental variables
    # The parameters need to be listed in ascending order and there should not be duplicates
    q1 = list(set(exp_params1))
    q2 = list(set(exp_params2))
    q3 = list(set(exp_params3)) 
    q3 = [6,8]
    q1.sort()
    q2.sort()
    q3.sort()
    params=[q1,q2,q3]
    experiments = []
    
    #Define the experiments
    experiments.append([(x,y,z) for x in q1 for y in q2 for z in q3])

    # Reasign lists
    exp = experiments[0] 
    exp_ind1 = 0
    exp_ind2 = 1
    exp_ind3 = 2
    q1 = params[exp_ind1]
    q2 = params[exp_ind2]
    q3 = params[exp_ind3]

    #allocate the vector
    n=[]
    n.append(len(q1))
    n.append(len(q2))
    n.append(len(q3))
    fEx = np.zeros(n)
    
    # Construct the quantity list
    q=[]
    q.append(np.array(q1))
    q.append(np.array(q2))
    q.append(np.array(q3))
 
    pdf=[]
    for i in range(len(q)):
        if n[i] != 1:
            pdf.append(np.ones(n[i])/(q[i][-1]-q[i][0]))   
        else:
            pdf.append(np.ones(n[i])/(1))    
    # Construct the function output
    for i in range(0,n[0]):
        for j in range(0,n[1]):
            for k in range(0,n[2]):
                fEx[i,j,k] = Nu_dict.get((q[0][i],q[1][j],q[2][k]), 1)
                #fEx[i,j,k] = abs(Nu_dict.get((q[0][i],q[1][j],q[2][k]), 2*Nu_ref) - Nu_ref)/Nu_ref

    print("=========")
    print("GSA") 
    print("TEMPORALLY COPYING THE SAME VALUES FOR LX 6 AND LX 8")
    print("=========")
    fEx[:,:,0] = fEx[:,:,1] ### This is a temporal step. Delete when you have more data points
    print('>q1 - ',q[0])
    print('>q2 - ',q[1])
    print('>q3 - ',q[2])
    print('> experiment Independent variables and Qoi value: ', Nu_dict)
    print('> Parameter matrix', fEx)
    print("=========")

    sobol_=sobol(q,fEx,pdf=pdf)

    return sobol_


## ======================================================
## Code Set-up
## ======================================================

master_filename = "nusselt_files.txt"
with open(master_filename) as f:
    filename_list = f.readlines()

# Code parameters
ra = 1e11
start_time = 200
Nu_v_ref = 222
Nu_a_ref = 229

# Mesh parameters
e_xy = 80
e_z  = 60
lx = 8
## Simulation code map to mesh size
map_name_sizes={}
map_name_sizes['001'] = np.array([e_xy   ,e_z   ,lx])
map_name_sizes['003'] = np.array([e_xy*2 ,e_z   ,lx])
map_name_sizes['004'] = np.array([e_xy*3 ,e_z   ,lx])
map_name_sizes['005'] = np.array([e_xy   ,e_z*2 ,lx])
map_name_sizes['006'] = np.array([e_xy*2 ,e_z*2 ,lx])
map_name_sizes['007'] = np.array([e_xy*3 ,e_z*2 ,lx])
map_name_sizes['008'] = np.array([e_xy   ,e_z*3 ,lx])
map_name_sizes['009'] = np.array([e_xy*2 ,e_z*3 ,lx])
map_name_sizes['010'] = np.array([e_xy*3 ,e_z*3 ,lx])

# Create list and directories to hold the data
Nu_v_mean_list = []
Nu_a_mean_list = []
e_xy_l = []
e_z_l  = []
lx_l  = []
Nu_v_dict={}
Nu_a_dict={}

## ===========
## Main code
## ===========

if __name__ == "__main__":

    #> Time series analysis - Get Nusselt and CI
    e_xy_l,e_z_l,lx_l,Nu_v_dict, Nu_a_dict = process_Nu_timeseries(start_time, ra, filename_list, map_name_sizes, False)
    
    #> GSA - Analize the effect of some parameters in the Nusselt calc
    sobol_nu_v_ = GSA_3variables(e_xy_l, e_z_l, [6,8], Nu_v_dict, Nu_v_ref) 
    sobol_nu_a_ = GSA_3variables(e_xy_l, e_z_l, [6,8], Nu_a_dict, Nu_a_ref) 

    print("=== GSA results ===")
    Si_v=sobol_nu_v_.Si
    Sij_v=sobol_nu_v_.Sij
    print('> nu_v - main indices Si by uqit:',Si_v)
    print('> nu_v main indices Sij by uqit:',Sij_v)
    print("=========")
        
    Si_a=sobol_nu_a_.Si
    Sij_a=sobol_nu_a_.Sij
    print('> nu_a - main indices Si by uqit:',Si_a)
    print('> nu_a main indices Sij by uqit:',Sij_a)
