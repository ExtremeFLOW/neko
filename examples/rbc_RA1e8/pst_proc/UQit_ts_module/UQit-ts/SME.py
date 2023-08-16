########################################################################
# Sample Mean Estimator (SME) and its Uncertainty Estimators (SMEuncert)
########################################################################
#-----------------------------------------------------------------------
# Saleh Rezaeiravesh, salehr@kth.se
#-----------------------------------------------------------------------
##ToDO:
# 1. change 'p' in ARM to 'order'
# 2. issue: outputs are out_
# 3. we may not need to call mean(), since uncert() retuns mean and tAvg as well. Later it might
#    be found that better to return values rather than out
# 4. DGX: the sampEstim would need to be retained in the old version. Error due to tElapsed exploding in the new version.
#    If only time is kept in the old version, the fSME_var array and tElapsed array are not of same length
# 5. DGX: I think a second division by n is needed for sampEstim if conv ='false'

import numpy as np
from CLSC import CLSC
from NOBM import NOBM
from BMBC import BMBC
from OBM import OBM
from ARM import ARM
from DACF import DACF

#
class SME:
    """
    SME constructor 

    Args:
      `t`: 1d numpy array of size n 
         Time samples in the time-series
      `f`: 1d numpy array of size n
         Samples of QoIs at samples `t`
      `opts`: dict (optional)
         Options for SME with the following keys:
           * `verbose`: bool (default: False)
                If True, actions are logged
           * `conv`: bool (default: True)
                If True, converegnce of the SME and its uncertainty over the elapsed time is considered.
                If False, SME and its uncertainty are computed considering the whole samples.

    Methods:
       `eval`: estimates the sample mean of `f`
       `uncert`: estimate SME and its uncertainty in the sample mean `val`

    Attributes:
       `tAvg`: 1d numpy array of size 1 (if `conv`:False) or (n-1) (if `conv`:True)
          Elapsed time when computing SME
       `mean`: 1d numpy array of size 1 (if `conv`:False) or (n-1) (if `conv`:True)
          Values of SME of f, <f>_n
       `uncert`: 1d numpy array of size 1 (if `conv`:False) or (n-1) (if `conv`:True)
          Uncertainty of SME of f
    """
    def __init__(self,t,f,opts={}):
        self.t=t
        self.f=f
        self.opts=opts
        self._info()
        self._check()

    def _info(self):
        self.verbose=False 
        if 'verbose' in self.opts.keys():
           self.verbose=self.opts['verbose']
        self.conv=True
        if 'conv' in self.opts.keys():
           self.conv=self.opts['conv']

    def _check(self):       
        if self.t.ndim>1:
           raise ValueError('t should be a 1D numpy array.') 
        if self.f.ndim>1:
           raise ValueError('f should be a 1D numpy array.') 
        self.n=len(self.t)
        if self.n!=len(self.f):
           raise ValueError('t and f should be of equal size.') 

    @classmethod
    def sampEstim(self,t_,f_,conv_):
        """
        Sample estimator of the mean and variance of the time-series
        """
        if conv_:
            tAvg=[0]   #elapsed time
            n = len(t_) #dgx changed time addition
            # fMean=[]  #estimator for mean of u: <f>_n
            # fVar=[]   #estimator for variance of f: <f^2-<f>_n^2>_n
            elapsedTim = 0.0
            # fSum=0.0
            # fSumSqr=0.0
            for i in range(1, n):#dgx changed time addition
                dt = t_[i]- t_[i-1]  # We need dt since the method is going to be also
                elapsedTim += dt
                tAvg.append(elapsedTim)   #To be used only for plotting
                # fSum+=f_[i]
                # fSumSqr+=f_[i]*f_[i]
                # fAvg_tmp=fSum/float(i)
                # fMean.append(fAvg_tmp)
                # fVar_tmp=fSumSqr/float(i)-fAvg_tmp**2.
                # fVar_tmp=fVar_tmp/(float(i)) ## second division by n
                # fVar.append(fVar_tmp)
           ##tAvg=np.cumsum(t_)
            # convert to numpy arrays
            tAvg = np.asarray(tAvg)
            # fMean = np.asarray(fMean)
            # fVar = np.asarray(fVar)
            fSum = np.cumsum(f_)
            fSumSqr = np.cumsum(f_**2.)
            n_ = np.arange(1, n+1)
            fMean = fSum/n_
            fVar = (fSumSqr/n_-fMean**2.)/n_
        else:
           tAvg=np.asarray([t_[-1]])    
           fMean=np.asarray([np.mean(f_)])
           fVar=np.asarray([np.std(f_)**2.])    

#        tAvg=np.asarray(tAvg)
#        fMean=np.asarray(fMean)
#        fVar =np.asarray(fVar)
        return tAvg,fMean,fVar

    def eval(self):
        """
        Estimates sample mean and variance of the time-series
        """
        tAvg,fMean,fVar=self.sampEstim(self.t,self.f,self.conv)
        self.tAvg=tAvg
        self.mean=fMean
        self.var=fVar
#
class SMEuncert:
    """    
    Class of estimators of the uncertainty in a SME
    """
    def __init__(self,SME,opts):
        self.t=SME.t
        self.f=SME.f
        self.n=SME.n
        self.SMEopts=SME.opts
        self.opts=opts
        self._check()
        self._estimate()

    def _messg(self,q_,L_,m_):
        if q_ not in L_:
           raise ValueError('%s is not provided for method %s.'%(q_,m_))
       
    def _check(self):
        if 'verbose' not in self.opts.keys():
           self.opts.update({'verbose':self.SMEopts['verbose']})
        if 'conv' not in self.opts.keys():
           self.opts.update({'conv':self.SMEopts['conv']})

        validMethods=['CLSC','NOBM','OBM','BMBC','ARM','DACF']
        if 'method' not in self.opts.keys():
            raise KeyError("'method' for estimating uncertainty is missing from SMEuncert.opts.")
        else:
            method_=self.opts['method']

        if method_ not in validMethods:
            print("Available methods for quantifying uncertainty in SME are: ",validMethods)
            raise ValueError("Provided 'method'=%s is invalid." %method_)

        if method_ in ['NOBM','OBM','BMBC']:
           self._messg('batchSize',self.opts.keys(),method_)
           if 'batchSize' in self.opts.keys():
              batchSize_=self.opts['batchSize']
              if batchSize_ >= self.n:
                 ValueError('Batch size should be smaller than the sample size.')               
        elif method_ in ['DACF','ARM']:
           self._messg('nChunk',self.opts.keys(),method_)
           if 'nChunk' in self.opts.keys():
              nChunk_=self.opts['nChunk']
              
           if method_=='ARM': #check for minimum required settings
              self._messg('nLagTrain',self.opts.keys(),method_)
              if isinstance(self.opts['nLagTrain'],int):
                 nLagTrain_=self.opts['nLagTrain']
              else:
                 raise ValueError("Value of 'nLagTrain' should be an integer!") 
              self._messg('AR-method',self.opts.keys(),method_)
              if 'p' not in self.opts.keys():
                 self._messg('maxlag',self.opts.keys(),method_)        

    def _estimate(self):
        """
        Estimator of the uncertainty in SME[f]
        """
        method_=self.opts['method']
        if method_=='CLSC':
           out_=CLSC(self)
        elif method_=='NOBM':
           out_=NOBM(self)
        elif method_=='OBM':
           out_=OBM(self)
        elif method_=='BMBC':
           out_=BMBC(self)
        elif method_=='ARM':
           out_=ARM(self)
        elif method_=='DACF':
           out_=DACF(self)
        self.estim=out_.out
