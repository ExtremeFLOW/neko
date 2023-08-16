#*********************************************************************
#  Estimation of \sigma^2(\mu) where \mu is the sample mean estimator
#  Direct Autocovaraince Function (DACF) Method
#*********************************************************************
#  Saleh Rezaeiravesh, salehr@kth.se
#---------------------------------------------------------------------
#ToDo: move tests
#
#import os
#import sys
import numpy as np
import math as mt
#UQit_ts_path=os.getenv("UQit_ts_path")
##sys.path.append(UQit_ts_path+'stats_tUQ/')
#sys.path.append(UQit_ts_path+'synDataGen_tUQ/')
#sys.path.append(UQit_ts_path+'plot_tUQ/')
#import CLSC
#import NOBM
import SME
import stats_ts
import plot_ts
import synDataGen_ts


class DACF:
    """
    Direct Autocovariance Function to estimate uncerttainty in SME
    """
    def __init__(self,SMEuncert):
        self.t=SMEuncert.t
        self.f=SMEuncert.f
        self.n=SMEuncert.n
        self.opts=SMEuncert.opts
        self._uncertEstim()
    
    @classmethod
    def varEstim_dirAC(self,f_,acf_fft=False):
       """
       Estimator for \sigma^2(\mu) considering autocorrelations between the samples
       NOTE: The AC's are directly estimated. 
       NOTE: The direct estimators of ACF's are biased.
       Inputs:
          f_: 1d numpy array, samples of TS
          acf_fft: bol, if True: ACF are computed by FFT
       Outputs:
          var_mu: scalar, estimated varaince for \mu=E[f] using SME
       """
       n_=len(f_)
       var_f=np.std(f_)**2.   #estimation of var[f]
       acf=stats_ts.autoCorrFunc(f_,nlags_=n_,fft_=acf_fft,plot=False)
       k=np.arange(1,n_,1)
       var_mu=(1.+2*np.sum((1-k/n_)*acf[1:]))*var_f/float(n_)
       return var_mu

    def _uncertEstim(self):
    ##def estimator_tUQ(t,f,nChunk):
       """
       Direct Autocovariance Function (DACF) Estimator for \sigma^2(\mu)  
          Inputs:
             t: time samples, 1d numpy array of length n
             f: time-series values, 1d numpy array of length n
             nChunk: integer, number of chunks made out of the samples                      
                     convergence of estimators are plotted using the chunks
          Outputs:
             outList: a dictionary of the outputs         
       """
       methodName='Direct Autocovariance Function (DACF) Estimator'
       #Estimate the mean by sample-mean estimator
       tAvg,fMean,fVarDummy=SME.SME.sampEstim(self.t,self.f,self.opts['conv'])
       n=self.n
       nChunk=self.opts['nChunk']
       #Construct a set of data chunks (no averaging on the chunks)
       #NOTE: each chunk start at the first sample and contains more sample as the chunk number increases
       #  M = n/nChunk
       #  number of samples in chunks: {M,2M,3M,...,n}
       #  chunk1: 1,2,...,M
       #  chunk2: 1,2,...,M,M+1,...,2M
       M=int(n/nChunk)
       fMeanHat=[]
       fSMEVarHat=[]
       tAvgHat=[]
       t0=self.t[0]  #initial time
       for k in range(nChunk):  #chunkCounter
           nEndInChunk=(k+1)*M
           if k==(nChunk-1): #last chunk
              nEndInChunk=n
           #construct the chunks
           tCh=self.t[:nEndInChunk]
           fCh=self.f[:nEndInChunk]
           nCh=len(fCh)
           tAvgHat.append(tCh[-1]-t0)
           fMeanHat.append(np.mean(fCh))
           #Estimate variance of SME
           fSMEVarHat.append(self.varEstim_dirAC(fCh,acf_fft=False))   
       #    
       fMean_=fMeanHat[-1]
       fSMEVar_=fSMEVarHat[-1]
       #Make a database out of the outputs
       out={'t':self.t,'f':self.f,'tElapsed':tAvg,'tElapsedEstim':tAvgHat,'fSME':fMean,\
            'fSME_var':fSMEVarHat,'method':'DACF','figName':'DACF'+str(M),\
            'tDWN':tAvgHat+t0,'fDWN':fMeanHat};
       if self.opts['verbose']:
          print('... Estimating uncertainty in SME sing ')
          print('    '+methodName)
          print('...... Estimations using all samples:')
          print('       mean=%g, var[SME]=%g, std[SME]=%g' %(fMean_,fSMEVar_,mt.sqrt(fSMEVar_)))
          print('       N x var[SME] = %g' %(n*fSMEVar_))
       self.out=out   
#   
#
###############
# Ext Func
###############
#////////////////////////
def varEstim_dirAC_test():
    """
       Test varEstim_dirAC()
    """
    #------ SETTINGS
    n=10000   #number of samples in the TS
    nChunk=100 #number of chunks
    #-------------------------------------
    #(1) Geerate synthetic TS samples
    synDataOpts={'dt':1,
                 'type':'synData3',
                 'noiseSdev':0.1,
                 'writeInFile':False
                }
    t,f=synDataGen_ts.syntheticDataGen(n,synDataOpts)
    #(2) Estimate \sigma^2(\mu) using DACF
    outDACF=estimator_ts(t,f,nChunk)
    #(3) Plot the original TS and convergence of DACF estimator
    plot_ts.plotter_signal_estimates(outDACF,{})
#
