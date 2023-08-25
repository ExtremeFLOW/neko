#*************************************************************************
#  Estimation of uncertainty in the sample-mean estimator of a time-series
#  'NOBM': Non-Overlapping Batch Method 
#*************************************************************************
#  Saleh Rezaeiravesh, salehr@kth.se
#-------------------------------------------------------------------------
#
import numpy as np
import math as mt
import SME

class NOBM:
    """
    Non-overlapping Batch Mean (NOBM) method to estimate uncertainty in a SME

    Args:
       `SMEuncert`: Uncertainty estimator in a SME
       Required keys:
         * 'method':'NOBM',
         * 'batchSize': int

    Methods:
       `_uncertEstim`: Estimates the uncertainty in the SME 
       `batchGen`: Generates non-overlapping batches       

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
            * `'method'`: string
              Method for estimating the uncertainty in the SME, here 'NOBM'
            * `'figName'`: string
              Figure name, when the data in `out` are plotted using the methods in 'plot_ts.py'
    
    Example:
      * To track the variation of the SME and its uncertainty with the sample size:
         `sme_=SME(t,f,{'verbose':False,'conv':True})`
         `out_=SMEuncert(sme_,{'method':'NOBM','batchSize':N}).estim`
      * To estimate the uncertainty considering the whole samples:
         `sme_=SME(t,f,{'verbose':False,'conv':False})`
         `out_=SMEuncert(sme_,{'method':'NOBM','batchSize':N}).estim`
    """
    def __init__(self,SMEuncert):
        self.t=SMEuncert.t
        self.f=SMEuncert.f
        self.opts=SMEuncert.opts
        self._uncertEstim()

    @classmethod
    def batchGen(self,t,f,M):
        """ 
        Generates `K` non-overlapping batches of size `M` from imported time-series (`t`,`f`) 

        Args:
           `t`: numpy array of size n
              Time samples
           `f`: numpy array of size n
              Time-series values at times `t`
           `M`: int
              Batch size

        Returns:   
           `K`: int
              Number of batches
           `tNOB`: numpy array of size `K` 
              Time associated to the non-overlapping batches
           `fNOB`: numpy array of size `K` 
              Mean of the samples in the batches

        NOTE: Some of the original samples may remain unused since `K*M` is not necessarily equal to `n`. 
        """
        n_=len(t)        
        K=int(n_/M)      
        tNOB=np.zeros(K)
        fNOB=np.zeros(K)
        for k in range(K):  
            for i in range(M):   
                j=k*M+i  
                tNOB[k]+=t[j]
                fNOB[k]+=f[j]
        tNOB/=float(M)
        fNOB/=float(M)
        return K,tNOB,fNOB

    def _uncertEstim(self):
        methodName='Non-Overlapping Batch Mean (NOBM) Estimator'
        tAvg,fMean,TMP=SME.SME.sampEstim(self.t,self.f,self.opts['conv'])
        K,tNOB,fNOB=self.batchGen(self.t,self.f,self.opts['batchSize'])
        #Estimation of the variance of the SME as a function of the elapsed time
        tAvgNOB,fMeanNOB,fVarNOB=SME.SME.sampEstim(tNOB,fNOB,self.opts['conv'])

        out={'t':self.t,'f':self.f,'tElapsed':tAvg,'tElapsedEstim':tAvgNOB,'fSME':fMean,\
             'fSME_var':fVarNOB,'method':'NOBM','figName':'NOBM'+str(self.opts['batchSize']),\
             'tDWN':tNOB,'fDWN':fNOB};

        if self.opts['verbose']:
           print('... Estimating uncertainty in the SME[f] using ')
           print('    '+methodName)
           print('...... Estimations using all samples:')
           print('       SME=%g, var[SME]=%g, std[SME]=%g' %(fMean[-1],fVarNOB[-1],mt.sqrt(fVarNOB[-1])))
           print('       N x var(SME) = %g' %(len(self.f)*fVarNOB[-1]))
        self.out=out
