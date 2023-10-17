#*************************************************************************
#  Estimation of uncertainty in the sample-mean estimator of a time-series
#  'CLSC': Classical (Ensemble) Method
#*************************************************************************
#  Saleh Rezaeiravesh, salehr@kth.se
#-------------------------------------------------------------------------
#
import numpy as np
import math as mt
import SME

class CLSC:
    """
    Classical method to estimate uncertainty in a SME. 
    The samples are assumed to be i.i.d. and uncorrelated.

    Args:
       `SMEuncert`: Uncertainty estimator in a SME

    Methods:
       `_uncertEstim`: Estimates the uncertainty in the SME 
       Required keys:
          * 'method':'CLSC'

    Attributes:   
       `out`: dict
          The results of the uncertainty estimation. The dict has the following keys:
            * `'t'`: numpy array of size n
              Sampled times in the original time-series
            * `'f'`: numpy array of size n
              Sampled values of the original time-series
            * `'tElapsed'`: numpy array of size n
              Elapsed times when computing cumulative averaging using (`t`,`f`)
            * `'fSME'`: numpy array of size n
              Values of the SME of `f` at `tElapsed` 
            * `'fSME_var'`: numpy array of size n
              Values of the variance of the SME of `f` at `tElapsed` 
            * `'method'`: string
              Method for estimating the uncertainty in the SME, here 'CLSC'
            * `'figName'`: string
              Figure name, when the data in `out` are plotted using the methods in 'plot_ts.py'

    Example:
      * To track the variation of the SME and its uncertainty with the sample size:
         `sme_=SME(t,f,{'verbose':False,'conv':True})`
         `out_=SMEuncert(sme_,{'method':'CLSC'}).estim`
      * To estimate the uncertainty considering the whole samples:
         `sme_=SME(t,f,{'verbose':False,'conv':False})`
         `out_=SMEuncert(sme_,{'method':'CLSC'}).estim`
    """
    def __init__(self,SMEuncert):
        self.t=SMEuncert.t
        self.f=SMEuncert.f
        self.opts=SMEuncert.opts
        self._uncertEstim()

    def _uncertEstim(self):
        methodName='Classical Ensemble Estimator (CLSC)'
        tAvg,fMean,fVar=SME.SME.sampEstim(self.t,self.f,self.opts['conv'])   
        out={'t':self.t,'f':self.f,'tElapsed':tAvg,'fSME':fMean,'fSME_var':fVar,\
                 'method':'CLSC','figName':'CLSC'}

        if self.opts['verbose']:
           print('... Estimating uncertainty in SME[f] using ')
           print('    '+methodName)
           #Compute estimations of mean and variance as a function of the elapsed time
           print('...... Estimations using all samples:')
           print('       SME=%g, var[SME]=%g, std[SME]=%g' %(fMean[-1],fVar[-1],mt.sqrt(fVar[-1])))
           print('       N x var(SME) = %g' %(len(self.f)*fVar[-1]))
        self.out=out          
