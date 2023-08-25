#*************************************************
# Tools for general statsitical analysis
#*************************************************
# Saleh Rezaeiravesh, salehr@kth.se
#-------------------------------------------------
# To Do: write and test autoCorrFunc_SME
#
import sys
import os
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf, acovf 
from statsmodels.tsa.stattools import pacf
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.graphics.tsaplots import plot_pacf
#UQit_ts_path=os.getenv("UQit_ts_path")
#sys.path.append(UQit_ts_path+'write_tUQ/')
#sys.path.append(UQit_ts_path+'synDataGen_tUQ/')
#sys.path.append(UQit_ts_path+'plot_tUQ/')
import synDataGen_ts
import write_ts
import plot_ts


#
#///////////////////////
def moments(t,f):
#def signalAnalyzer(t,f):
    """ 
        Statsitical analysis of a signal
    """
    print("Note: assumption is that the samples are equi-spaced over time!")
    #Statistical moments of u samples
    mean,std,skew,kurt,superSkew=statisticalMoments(f)
    write_ts.pw('    mean      std     skew      kurt    superkurt')
    write_ts.pw(write_ts.printRepeated('-', 50))
    write_ts.pw('%g  %g  %g  %g  %g' %(mean,std,skew,kurt,superSkew))
    out_={'mean':mean,'std':std,'skew':skew,'kurt':kurt,'superSkew':superSkew}
    return out_ 
#
#/////////////////////////
def statisticalMoments(f):
    """ 
       Compute statitical moments of a time-series using
           Numpy and scipy.stats methods
    """
    #f: time signal
    n=len(f)
    mean=np.mean(f)
    std=np.std(f)
    skew=st.skew(f)    #3rd central moment/std^3
    #4th central moment/std^4, if fisher=True => 3 is subtracted from calulated kurtosis 
    kurt=st.kurtosis(f,fisher=False)
    superSkew=st.moment(f,moment=6)/std**6.  #superskewness
    return mean,std,skew,kurt,superSkew
#
#//////////////////////////////////
def crossCorr(u1,u2,mode_='full'):
    """ 
       Compute cross correlations between signals u1 and u2 
       Inputs:
          u1, u2: two signals, each a numpy 1D array of size n.
                  The discrete times should be uniformly spaced
          mode_='full': discrete convolution
                'valid'
                'same'
       Outputs:
          r: cross correlation, numpy 1D array of size n or 1 (mode_=='valid')
       NOTE: for autocorrelation send in (u1,u1,mode='full')   
       NOTE: np.correlate only works for uniformly spaced samples (since works with the sample number)
    """
    #u1=u1.astype('Int64')
    #u2=u2.astype('Int64')
    r=np.correlate(u1,u2,mode=mode_) 
    if mode_ == 'full':
       r_= r[int(r.size/2):]
    else:
       r_=r
    return r_
#
#//////////////////////////////////////////////////
def autoCorrFunc(u,nlags_=1,fft_=False,plot=False):
    """
       Computes autocorrelation function, acf
       Inputs: 
           u    : time series data, numpy 1d array of size n
           nlag_: max number of lags to include in ACF   
       Ouputs:
           acf_ : ACF
       According to         
       "For very long time series it is recommended to use fft convolution instead. 
        When fft is False uses a simple, direct estimator of the autocovariances that only 
          computes the first nlag + 1 values. This can be much faster when the time series is 
          long and only a small number of autocovariances are needed.If unbiased is true, 
          the denominator for the autocovariance is adjusted but the autocorrelation is not 
          an unbiased estimator."
    """
    #write_ts.pw('...... Computing autocorrelations (ACF)')
    acf_=acf(u,fft=fft_,nlags=nlags_)#,alpha=0.05)
    if plot:
       plot_acf(u,lags=nlags_,alpha=0.05)
    return acf_
#
#//////////////////////////////////////
def autoCovFunc(u,nlags_=1,plot=False):
    """
       Computes autocorrelation function, acf
       Inputs: 
           u    : time series data, numpy 1d array of size n
           nlag_: max number of lags to includein ACF   
       Ouputs:
           acf_ : ACF
       According to         
       "For very long time series it is recommended to use fft convolution instead. 
        When fft is False uses a simple, direct estimator of the autocovariances that only 
          computes the first nlag + 1 values. This can be much faster when the time series is 
          long and only a small number of autocovariances are needed.If unbiased is true, 
          the denominator for the autocovariance is adjusted but the autocorrelation is not 
          an unbiased estimator."
    """
    #write_ts.pw('...... Computing autocorrelations (ACF)')
    acovf_=acovf(u,fft=False)#,alpha=0.05)
#    if plot:
#       plot_acf(u,lags=nlags_,alpha=0.05)
    return acovf_
#
#///////////////////////
def autoCorrFunc_SME(f):
    """
       Estimate ACF using Sample Mean Estimator
       Inputs: 
           f   : time series data, numpy 1d array of size n
           nLags: max lag 
       Output:
           acf_: ACF estimated by MSE
    """
    # See tests/main_tUQ.py, sample_ACF_analysis()
#
#///////////////////////
def partAutoCorrFunc(u,nlags_=1,method_='ywunbiased',plot=False):
    """
       Computes partial autocorrelation estimate, pacf
       Inputs: 
           u   : time series data, numpy 1d array of size n
           lag_: int, number of lags to consider  
           method_: methof to compute pacf: see
               https://www.statsmodels.org/stable/generated/statsmodels.tsa.stattools.pacf.html
       Outputs
           pacf=1d numpy array of PACFs
           confd=(1-alpha)% confidence for the estimated pacfs
    """
    write_ts.pw('...... Computing partial autocorrelations (PACF)')
    pacf_,ci_=pacf(u,method=method_,nlags=nlags_,alpha=0.05)   
    if plot:
       plot_pacf(u,lags=nlags_,alpha=0.05)
       plt.show()
    return pacf_,ci_
#
#
################
#External tests
################
def autoCorrFunc_test():
    """
       Test for autoCorrFunc
    """
    # ---- SETTINGS
    n=100000              #number of samples in the time-series
    dt_=1              #time step size in the time-series
    method_='synData2' #type of synthetic data  
    writeInFile_=False  #if time-series to be written in a file
    # ---------------------
    #(1) Generate a synthetic signal
    synDataOpts={'dt':dt_,   
                 'type':method_,
                 'noiseSdev':0.1,
                 'writeInFile':writeInFile_
                }
    t,f=synDataGen_ts.syntheticDataGen(n,synDataOpts)
    plot_ts.histoPlot_single(t,f,{})
    #(2) Compute the autocorrelation functions of the signal
    r=autoCorrFunc(f,int(n/5),plot=False)
    c=autoCovFunc(f,int(n/5),plot=False)
    #plot    
    plt.figure(figsize=(10,10))
    plt.subplot(2,1,1)
    plt.plot(r,'-ob')
    plt.xlabel('lag',fontsize=15)
    plt.ylabel('ACF',fontsize=14)
    plt.grid(alpha=0.3)
    plt.subplot(2,1,2)
    plt.plot(c,'-+r')
    plt.xlabel('lag',fontsize=15)
    plt.ylabel('ACovF',fontsize=14)
    plt.grid(alpha=0.3)
    print('AcovF[0]= ',c[0])
    print('var[f]= ',np.std(f)**2.)
    plt.show()
