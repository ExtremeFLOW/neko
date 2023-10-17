#*************************************************************************
#  Estimation of uncertainty in the sample-mean estimator of a time-series
#  'BMBC': Batch Means and Batch Correlation Method 
#*************************************************************************
#  Saleh Rezaeiravesh, salehr@kth.se
#-------------------------------------------------------------------------
#
import numpy as np
import math as mt
import SME
import NOBM

class BMBC:
    """
    Batch Means and Batch Correlation (BMBC) Method to estimate uncertainty in a SME  

    Reference: S. Russo & P. Luchini, Journal of Computational Physics,
       Volume 347, 2017, Pages 328-340.

    Args:
       `SMEuncert`: Uncertainty estimator in a SME
       Required keys:
         * 'method':'BMBC',
         * 'batchSize': int

    Methods:
       `_uncertEstim`: Estimates the uncertainty in the SME 

    Attributes:   
       `out`: dict
          The results of the uncertainty estimation. The dict has the following keys:
            * `'t'`: numpy array of size n
              Sampled times in the original time-series
            * `'f'`: numpy array of size n
              Sampled values in the original time-series
            * `'tElapsed'`: numpy array of size n
              Elapsed times when computing cumulative averaging using (`t`,`f`)
            * `'fSME'`: numpy array of size n
              Values of the SME of `f` at `tElapsed` 
            * `'tDWN'`: numpy array of size K (number of batches)   
              Time samples associated to the non-overlapping batches
            * `'fDWN'`: numpy array of size K 
              Mean of the batches associated to `tDWN`
            * `'tElapsedEstim'`: numpy array of size K
              Elapsed times when computing cumulative averaging using (`tDWN`,`fDWN`)
            * `'fSME_var'`: numpy array of size K
              Variance of the SME of `f` at `tElapsedEstim` 
            * `'S1/S0'`: float
              Ratio of the lag-1 correlation to the variance of the batch means
            * `'method'`: string
              Method for estimating the uncertainty in the SME, here 'BMBC'
            * `'figName'`: string
              Figure name, when the data in `out` are plotted using the methods in 'plot_ts.py'

    Example:
      * To track the variation of the SME and its uncertainty with the sample size:
         `sme_=SME(t,f,{'verbose':False,'conv':True})`
         `out_=SMEuncert(sme_,{'method':'BMBC','batchSize':N}).estim`
      * To estimate the uncertainty considering the whole samples:
         `sme_=SME(t,f,{'verbose':False,'conv':False})`
         `out_=SMEuncert(sme_,{'method':'BMBC','batchSize':N}).estim`
    """
    def __init__(self,SMEuncert):
       self.t=SMEuncert.t
       self.f=SMEuncert.f
       self.opts=SMEuncert.opts
       self._uncertEstim()

    def _uncertEstim(self):
        methodName='Batch Means and Batch Correlations (BMBC) Estimator'
        tAvg,fMean,TMP=SME.SME.sampEstim(self.t,self.f,self.opts['conv'])
        #Compute Non-Overlapping Batch samples
        K,tNOB,fNOB=NOBM.NOBM.batchGen(self.t,self.f,self.opts['batchSize'])
        #Estimate mean and variance of f as a function of elapsed time
        tAvgBMBC,fMeanDummy,fVarDummy=SME.SME.sampEstim(tNOB,fNOB,self.opts['conv'])
        
        #Estimation of the variance of the SME as a function of the elapsed time
        fpNOB=fNOB-fMean[-1]
        fVarBMBC=np.zeros(K)
        S0=fpNOB[0]*fpNOB[0]
        S1=0.0
        S2=0.0
        fVarBMBC[0]=S0+2.*S1
        for k in range(1,K):
            S0+=fpNOB[k]**2.
            S1+=fpNOB[k-1]*fpNOB[k]
            if k>1:
               S2+=fpNOB[k-2]*fpNOB[k]
            km=k-1
            if km==0:
               km=k
            fVarBMBC[k]=(S0+2.*S1)/float(k*km)

        out={'t':self.t,'f':self.f,'tElapsed':tAvg,'tElapsedEstim':tAvgBMBC,'fSME':fMean,\
             'fSME_var':fVarBMBC,'method':'BMBC','figName':'BMBC'+str(self.opts['batchSize']),\
             'tDWN':tNOB,'fDWN':fNOB,'S1/S0':S1/S0,'S2/S0':S2/S0}

        if self.opts['verbose']:
           print('... Estimating uncertainty in SME[f] using ')
           print('    '+methodName)
           print('...... Estimations using all samples:')
           print('       SME=%g, var[SME]=%g, std[SME]=%g' %(fMean[-1],fVarBMBC[-1],mt.sqrt(fVarBMBC[-1])))
           print('       Value of S1/S0 = %g' %(S1/S0)) 
           print('       Value of S2/S0 = %g' %(S2/S0)) 
           print('       N x var[SME]   = %g' %(len(self.f)*fVarBMBC[-1]))
        self.out=out
