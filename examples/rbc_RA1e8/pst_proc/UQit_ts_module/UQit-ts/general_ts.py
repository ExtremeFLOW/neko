#***************************************************
# General functions for time-series analysis
#***************************************************
# Saleh Rezaeiravesh, salehr@kth.se
#---------------------------------------------------
#
import os 
import sys
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
#UQit_ts_path=os.getenv("UQit_ts_path")
#sys.path.append(UQit_ts_path+'synDataGen_tUQ/')
import synDataGen_ts
#
#/////////////////////////////////////
def interp1D(t1,u1,t2,kind_='linear'):
    """ 
       1D interpolate (t1,u1) to (t2,u2) using method kind_
          Inputs: 
             t1: 1d numpy array of size n1, time
             u1: 1d numpy array of size n1, TS values corres. to t1
             t2: 1d numpy array of size n2, time samples (uniform)
             kind_: interpolation method
               'linear': linear interpolation
               ‘nearest’, ‘previous’, ‘next’
               ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’: 
                    0-th,1st,2nd,3rd- order Spline interpolation
          Outputs:
             u2: 1d numpy array of size n2, TS values corres. to t2
    """
    n1=len(t1);
    n2=len(t2);
    u2=np.zeros(n2)
    f=interpolate.interp1d(t1,u1,kind=kind_)
    u2=f(t2)    
    print('...... Interpolating %d samples to a series of %d samples.'%(n1,n2))    
    print('       using %s method.'%kind_)
    return u2
#
#////////////////////////
def interp1d_plot(t1,f1,t2,f2):
    """
       Plot original and mapped time-series
    """
    plt.figure(figsize=(10,4))
    plt.plot(t1,f1,'-b',lw=2,label='Original')
    plt.plot(t2,f2,'--r',label='Mapped')
    plt.legend(loc='best',fontsize=15)
    plt.xlabel('t',fontsize=17)
    plt.ylabel('f',fontsize=17)
    plt.tick_params(labelsize=16)
    plt.show()    
#
#############
# Ext Funcs
#############
#///////////////////////////
def interp1D_test():
    """ 
       Test interp1D()
       Interpolate a time-series with variable time-step size to a 
          series with fixed-size time steps.
       Number of samples in the original and interpolated series can be different.
    """
    #---- SETTINGS
    n1=1000          #no of original samples
    n2=600           #no of interpolated samples
    interpMethod='linear' #Interpolation method
    #-----------------------
    #(1) Generate original time-series with random time steps
    synOpts={'type':'synData3',
             'dt':(np.random.rand(n1-1)+0.1)*np.ones(n1-1),
             'noiseSdev':0.2}
    t1,f1=synDataGen_ts.syntheticDataGen(n1,synOpts)
    #(2) Interpolate the series to a fixed-size time-step series  
    #    with different number of samples
    t2=np.linspace(min(t1),max(t1),n2)
    f2=interp1D(t1,f1,t2,kind_=interpMethod)
    #(3) plot
    interp1d_plot(t1,f1,t2,f2)
#
#
