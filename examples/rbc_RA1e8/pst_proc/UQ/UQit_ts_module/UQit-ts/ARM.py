#*************************************************************************
#  Estimation of uncertainty in the sample-mean estimator of a time-series
#  'ARM': Auto-Regressive Model
#*************************************************************************
#  Saleh Rezaeiravesh, salehr@kth.se
#-------------------------------------------------------------------------
#
#ToDo: 1.change 'p' to 'order'
#      2. deal with what should be deprecated
#      3. fix the examples and move them to tests
#      4. use raise Errors 
#      5. check if options 'conv' and 'nChunk' are compatible
#
import os
import sys
import numpy as np
import math as mt
import subprocess
import shutil
import statsmodels.tsa.api as smt   #https://www.statsmodels.org/stable/api.html
import statsmodels.api as sm
import statsmodels
import time
import matplotlib.pyplot as plt

import SME
import synDataGen_ts
import stats_ts
import plot_ts

class ARM:
    """
    Autoregressive model (ARM) to estimate uncertainty in a SME

    Args:
       `SMEuncert`: Uncertainty estimator in a SME
       Required items:
         * 'method': 'ARM'
         * 'nChunk': int
              Number of chunks from the original samples
         * 'AR-method': string
              Method to fit an ARM to the samples: 'burg', 'yw', 'umle'
         * 'p': int
              Order of the ARM
         or
         * 'maxlag': int
              Maximum order of the ARM 
         * nLagTrain: int
              Number of lags at which sample-estimated autocorrelations are used to model the ACF.
       For other options for fitting the ARM, see ARM.fit.


    Methods:
       `_uncertEstim`: Estimates the uncertainty in the SME 
       
    

    """
    def __init__(self,SMEuncert):
       self.t=SMEuncert.t
       self.f=SMEuncert.f
       self.n=SMEuncert.n
       self.opts=SMEuncert.opts
       self._uncertEstim()

    @classmethod
#AR_fit -> fit
    def fit(self,f,opts):
        """
        Given n samples {x} of a TS, fit an Autoregressive model (ARM) to the demeaned samples.
        x:=f_demeaned=f-SME[f]
        SME[f] will be exported with the key 'orig_sampMean'
        ARM: x[i] = a0 + a1 x[i-1] + a2 x[i-2] + a3 x[i-3] + ... + ap x[i-p] + eps
             eps ~ N(0,sigma_d^2)
        Fitting includes estimating coefficients {(a0),a1,a2,...,ap}, the model order p, 
             and noise standard deviation sigma_d. a0 is the trend and is optional.
        Inputs:
           f: 1d numpy array of size n, samples of the TS
            NOTE: f is demeaned: i.e. f_new=f-SME[f]
           opts: settings with the description below; also see
            https://www.statsmodels.org/stable/generated/statsmodels.tsa.ar_model.AR.fit.html#statsmodels.tsa.ar_model.AR.fit
           https://www.statsmodels.org/stable/regression.html
           **Value or criterion to estimate p, the optimal lag        
             'p': integer, the value of p 
             if 'p' is not given, it needs to be estimated using a criterion and method
               'maxlag': integer, max possible value for p    
               'p-criter': (optional) criterion to select the optimal p
                  {'aic' [default],'bic','hqic','t-stat'}
              'p-method': (optional)
                 'cmle'- Conditional maximum likelihood using OLS                      
                 'mle' - [default] Unconditional (exact) maximum likelihood                  
              'p-trend': (optional) {'c','nc'} whether ('c') or not ('nc') [default] to 
                       include a0=const in ARM.
           **Method for estimating coefficients {a} of the AR with associated submethods
            'AR-method': 
               'umle': Unconditional MLE
               'yw' : Yule-Walker equation
               'burg': Burg's method
            'AR-submethod': (optional) submethod corresponding to a 'AR-method':
               if 'umle':   
                  'cmle': [default] Conditional maximum likelihood using OLS 
                  'mle' : Unconditional (exact) maximum likelihood 
               if 'yw': 
                  'mle': denominator of ACF estimator = n
                  'unbiased': [default] denominator of ACF estimator = (n-k)
               if 'burg':
                  empty (no 'AR-submethod')
         **Optional settings for a specific 'AR-method' (for a full set see the above links)        
           if 'AR-method'=='umle':
              'fit-trend': (optional) 'c' including trend (a0), 'nc' (no trend included) [default]
              'fit-optim': (optional) if 'AR-submethod'=='mle' optimzation method from scipy.optimize
                 {lbfgs'[default],'newton','nm','bfgs','powell','cg','ncg',...}
              'fit-maxiter': (optional) integer, max iterations in the optimization   
           if 'AR-method'=='yw':
              'demean': (optional) bool, True is [default]
           if 'AR-method'=='burg':
              'demean': (optional) bool, True is [default]
        Outputs:
           AR_coefs: coefficients {a} of the fitter AR(p) 
           AR_sig2d: sigma_d^2, variance of the noise in the AR model
           opts: the input opts which is updated with the optional settings
        """
        #(1) Assign the default values if opts is empty
        # (a) settings for p
        if 'p' not in opts.keys():
           if 'p-trend' not in opts.keys():
              opts.update({'p-trend':'nc'}) 
           if 'p-method' not in opts.keys():
              opts.update({'p-method':'mle'}) 
           if 'p-criter' not in opts.keys():
              opts.update({'p-criter':'aic'}) 
        # (b) settings for estimating the coefficients of AR(p)   
        if opts['AR-method']=='umle':
           if 'AR-submethod' not in opts.keys():
              opts.update({'AR-submethod':'cmle'})     
           if 'fit-trend' not in opts.keys():
              opts.update({'fit-trend':'nc'})     
           if 'fit-optim' not in opts.keys():
              opts.update({'fit-optim':'lbgfs'})#is used only for 'mle'     
           if 'fit-maxiter' not in opts.keys():
              opts.update({'fit-maxiter':100})     
        elif opts['AR-method']=='yw':
           if 'AR-submethod' not in opts.keys():
              opts.update({'AR-submethod':'unbiased'})     
           if 'demean' not in opts.keys():
              opts.update({'demean':True})     
        elif opts['AR-method']=='burg':
           if 'demean' not in opts.keys():
             opts.update({'demean':True})     
        elif opts['AR-method']!='ols':
            print('ERROR in ARM.fit(): Invalid "AR-method"') 
        # (c)    
        if 'verbose' in opts.keys() and opts['verbose']:
           print('... Fitting an AR to %d samples using %s method' %(len(f),opts['AR-method']))
        # 
        #(2) Assign or estimate p, the optimal lag
        if 'p' in opts.keys() and isinstance(opts['p'], (int)):
           p=opts['p'] 
           #print('    Assigned optimal lag p=',p)
        else:   
           p = smt.AR(f).select_order(maxlag=opts['maxlag'], ic=opts['p-criter'], \
                                  trend=opts['p-trend'],method=opts['p-method'])
           print('    Estimated optimal lag p=%d using method %s'%(p,opts['p-method']))
           opts.update({'p':p})       
        #
        #(3) Demean the input samples and save its sample-mean in opts
        f_mean=np.mean(f)
        opts.update({'orig_sampMean':f_mean})
        f=f-f_mean   #demean
        #(4) Estimate the AR(p) coefficients and noise varaince sigma_d^2
        if opts['AR-method']=='ols':
           AR_coefs,AR_sigd2=self.olsFit(f,p)     
        if opts['AR-method']=='umle':
           fit_=smt.AR(f).fit(maxlag=p,method=opts['AR-submethod'],maxiter=opts['fit-maxiter'],\
                          trend=opts['fit-trend'],solver=opts['fit-optim'])
           AR_coefs=fit_.params
           AR_sigd2=fit_.sigma2
        elif opts['AR-method']=='yw':   #Yule-Walker eqns
           AR_coefs,AR_sigd2=sm.regression.yule_walker(f,order=p,method=opts['AR-submethod'],\
                                                   demean=opts['demean'])               
           AR_sigd2=AR_sigd2**2.
        elif opts['AR-method']=='burg': #Burg's method
           AR_coefs,AR_sigd2=sm.regression.linear_model.burg(f,order=p,demean=opts['demean'])
        return AR_coefs,AR_sigd2,opts

    @classmethod    
#AR_sample -> sampGen
    def sampGen(self,fInit,n,AR_coefs,AR_noiseVar,opts):
        """
        Generating n samples from an AR(p) starting from initial samples fInit.
           ARM: x[i] = a1 x[i-1] + a2 x[i-2] + a3 x[i-3] + ... + ap x[i-p] + eps
               eps ~ N(0,sigma_d^2)
        Inputs:
           fInit: 1d numpy array of size nInit containing initial samples for AR(p)
              NOTE: nInit > p
           n: size of the samples to be taken from the AR(p)    
           AR_coefs: 1d numpy array of size p containing coeffcients of an AR model, 
                    AR_coefs:[a1,a2,a3,...,ap]
           NOTE: The AR is assumed to be demeaned.       
           AR_noiseVar: scalar, variance of the noise in the ARM    
           opts: options
        """
        nInit=len(fInit)
        f=np.zeros(n)
        nCoefs = len(AR_coefs)
        sig_=mt.sqrt(AR_noiseVar)
        eps=sig_*np.random.randn(n)
        for i in range(nInit):
            f[i]=fInit[i]
        for i in range(nInit,n):
            for j in range(nCoefs):
                k=(i-j-1)
                f[i]+=AR_coefs[j]*f[k]
            f[i]+=eps[i]
        return f

    @classmethod
#AR_ols_fit -> olsFit
    def olsFit(self,f,p):
        """
        OLS method to fit an AR(p) to samples f
           solve a linear system A a = R
           where A: pxp
                a:px1: coefficients of AR(p)
                R:px1
        Inputs: 
           f: 1d numpy array of size n, samples of a time-series
           p: order of the AR model
        Outputs: 
           a: coefficients of AR(p) estimated by OLS
           var_ols: variance of the sum of squares of residual =
                OLS estimate for varaince of the  iid noise in AR(p)
        """
        n=len(f)   #sample size
        R=np.zeros(p)     
        A=np.zeros((p,p)) 
        for k in range(p):
            sumR_=0.0
            k1=k+1
            for i in range(p,n):
                sumR_+=f[i-k1]*f[i]
            for l in range(p):
                l1=l+1
                sumA_=0.0
                for i in range(p,n):
                    sumA_+=f[i-k1]*f[i-l1]
                A[k,l]=sumA_
                R[k]=sumR_
        #solve the linear system to estimate AR coefficients       
        a=np.linalg.solve(A,R)
        #compute sum of sqaures of residuals =  SME of variance of AR noise
        resid=0
        for i in range(p,n):
            sum_=0
            for k in range(p):
                sum_+=a[k]*f[i-(k+1)]
            resid+=(f[i]-sum_)**2.
        var_ols=resid/float(n)    
        return a,var_ols       

    @classmethod
    def acf_AR_charactEq(self,a,p,verbose_=False):
        """
        Solving the characteristic equations for ACF's given the coefficients of an AR(p).
        Ansatz: \rho_k=ACF(k)=\lambda^k = \sum_{i=1}^p a_i \lambda^(k-i)    
        AR(p): x[i] = a0 + a1 x[i-1] + a2 x[i-2] + a3 x[i-3] + ... + ap x[i-p] + eps
              eps ~ N(0,sigma_d^2)
              a0: trend, optional
        Charactersitic equation:
          (a0-1)\lambda^p+a1\lambda^(p-1)+a2\lambda^(p-2)+...+ap = 0
          solve this p-order polynomial eqn for p lambda's
        Inputs:
           a: 1d numpy array, coeffcients of AR(p): (a0),a1,a2,a3,...,ap
              (It may or may not include a0, it will be figured out using p)
           p: order of the AR   
         Outputs:
           lambda: 1d numpy array of size p containing the lambda's in the charactersitic equation
        """
        if verbose_:
           print('...... Estimating solutions for charactersitic equations of ACF\'s')
        #(1) Update array {a} acc. to the charactersitic eqn: a=[(a0-1),a1,a2,...,ap]
        # (a) check if a0 (AR's trend) is in input {a}
        if (len(a)==(p+1)): #if a0 already exsists in {a}
           a[0]-=1.0 
        elif (len(a)==p):
           a=np.insert(a,0,-1.0)
        #(2) solve the characteristic equation  
        lambda_=np.roots(a)
        return lambda_

    @classmethod
    def acf_powModel_fit(self,a,p,acf_samp,verbose_=False):
        """
        Making a power-law function for ACF in terms of lag, given coeffcients of an
           AR(p) fitted to TS samples along with a set of ACF estimated for (K-1) initial lags.
           Model ACF: \rho_k=\sum_{i=1}^p C_i \lambda_i^k
        Tuning the Model: 
           1. Given {a}'s, estimate \lambda's by solving a characteristic eqn.
           2. Given K \rho, estimate p coeffcients {C}'s. These K training {\rho}'s
              are given for instance, by sample estimator.
           3. The tuned ACF model can be used for prediction of ACF's up to any lag.    
        Inputs:
           a: 1d numpy array containing the coeffcients of AR(p)
              a:{(a0),a1,a2,...,ap}, a0 (trend) is optional
           p: integer, order of the AR fitted to the samples   
           acf_samp: 1d numpy array of size K, containing a set of training ACF's estimated by 
              sampling, for instance, with (K-1) number of lags.  
              acf_samp={acf[0]=1,acf[1],...,acf[K-1]}
        Outputs:
           C: 1d numpy array of size p, containing C's: C1,C2,...,Cp
           lam: 1d numpy array of size p, containing \lambda's
        """
        if verbose_:
           print('... Fitting a power-law model for ACF vs. lag')
        #(1) Compute p lambda's from a characteristic eqn
        lam=self.acf_AR_charactEq(a,p,verbose_)
        #(2) Construc the linear-set of equations to compute C's
        K=len(acf_samp)   #K-1: max-lag included in the training ACF samples
        A=np.zeros((K,p),dtype=np.complex_)
        for i in range(K):
            A[i,:]=lam**i
        #Make the eqn normal and solve for {C}'s
        M=np.dot(A.T,A)
        R=np.dot(A.T,acf_samp)
        #(3) Solve the normal linear system for {C}'s
        C=np.linalg.solve(M,R)
        return C,lam

    @classmethod
    def acf_powModel_pred(self,C,lambda_,K):
        """
        Predict ACF using a fitted power-law model
        Model ACF: \rho_k=\sum_{i=1}^p C_i \lambda_i^k
                 where, k=0,1,2,...,K
        Inputs: 
           C: 1d numpy array of size p, coeffcients {C} in the model
           lambda_: 1d numpy array of size p, values of lambda in the model
           K: integer, maximum lag at which ACFs are predicted
        Outputs:
           acf_model: 1d numpy arrau of size K, predicted ACFs by the model
        """
        acf_model=np.zeros(K) 
        #NOTE: The imaginary parts are discarded below since they are too small
        for k in range(K):
           tmp=lambda_**k
           acf_model[k]=np.sum(np.real(C*tmp))
        return acf_model    

    @classmethod
    def varEstim(self,f,nLagTrain,AR_opts):
        """
        Estimator for \sigma^2(\mu) using ARM 
        NOTE: The ACF's are predicted by a power-law model which is
             tuned-up in advance using some training ACF's 
        Inputs:
           f: 1d numpy array, samples of TS
           nLagTrain: integer, maximum lag for training ACFs in the fitted model
           Note: it can be overriden by the lag corresponding to `acfThresh` that can be provided in `AR_opts`
           AR_opts: dictionary containing options for the AR to be fitted
               For main associated keys, see ARM.fit()
               Extra keys:
               `acfThresh`: float btween 0 and 1
                   Threshhold for autocorrelation; if provided sample-estimated ACF up to that are used to fit the ACF model
               `acfPLot`: bool,
                   If True, the sample-estimated ACF and modeled ACF are plotted.
                   
        Outputs:
           var_mu: scalar, estimated variance for \mu=E[f] where E is SME
        """
        if 'verbose' in AR_opts.keys() and AR_opts['verbose']:
           print('... Estimating uncertainty in Sample-mean estimator using AR model')
        n=len(f)   #number of TS samples 
        #(1) Fit AR to the samples
        coefs_AR,sig2d_AR,AR_opts=self.fit(f,AR_opts)
#        print('AR',coefs_AR)
        #(2) Fit a power-law model to the ACF

        # (a) Estimate (sample-mean) training ACF's 
        if 'acfThresh' not in AR_opts.keys():
           # (a) Estimate (sample-mean) training ACF's with max-lag=(nLagTrain-1)
           acf_train=stats_ts.autoCorrFunc(f,nlags_=nLagTrain,fft_=False,plot=False)            
        else:  #override nLagTrain by the lag corresponding to `acfThresh`   
           acfThresh=AR_opts['acfThresh'] 
           acf_train=stats_ts.autoCorrFunc(f,nlags_=int(n/2),fft_=False,plot=False)
           iACF_thresh=np.argmin(abs(acf_train-acfThresh))  
#           if iACF_thresh < nLagTrain:
           print('Overriding nLagTrain=%d by lag=%d associated to acfThresh=%g' %(nLagTrain,iACF_thresh,acfThresh))
#              nLagTrain=min(nLagTrain,iACF_thresh)
#              acf_train=acf_train[:nLagTrain]           
           nLagTrain=int(iACF_thresh)
           acf_train=acf_train[:nLagTrain]           

        #(4) Fit a power-law model to ACF'
        verbose_=False 
        if 'verbose' in AR_opts.keys():
           verbose_=AR_opts['verbose']              
        C,lambda_=self.acf_powModel_fit(coefs_AR,AR_opts['p'],acf_train,verbose_)
        #(5) Use the fitted model to estimate ACF's at all n-samples
        acf_model=self.acf_powModel_pred(C,lambda_,n)

        if 'acfPlot' in AR_opts.keys() and AR_opts['acfPlot']: #plot modeled and sampled ACF in ARM estimator
           acf_full=stats_ts.autoCorrFunc(f,nlags_=n,fft_=False,plot=False)
           k_test=np.arange(n)
           k_train=np.arange(len(acf_train))

           plt.figure(figsize=(10,4))
           plt.plot(k_test,acf_model,'--r',lw=2,label='Modelled')
           plt.plot(k_test,acf_full,'-b',mfc='none',ms=5,label='Sample-estimate')
           plt.semilogx(k_train,acf_train,'sk',mfc='y',ms=5,label='Training for modeling')
           #plt.plot(k_test,MACF_fun(k_test,q[0],q[1],q[2]),':k',label='curve-fit, Newton')
           plt.legend(loc='best',fontsize=14)
           plt.xlabel('Lag',fontsize=14)
           plt.ylabel('ACF',fontsize=14)
           plt.grid(alpha=0.3)
           plt.show()

        #(6) Estimate \sigma^2(\mu), the varaince of the 
        var_f=np.std(f)**2.   #estimation of var[f]
        k=np.arange(1,n,1)
        var_mu=(1.+2*np.sum((1-k/n)*acf_model[1:]))*var_f/float(n)
        return var_mu

    def _uncertEstim(self):
##def estimator_tUQ(t,f,nChunk,nLagTrain,AR_opts):
        """
        Estimator for \sigma^2(\mu) using ARM
        Inputs:
           t: time samples, 1d numpy array of length n
           f: time-series values, 1d numpy array of length n
           nChunk: integer, number of chunks made out of the samples                      
                 convergence of estimators are plotted using the chunks       
           nLagTrain: integer, maximum lag for training ACFs in the fitted model
           acfThresh: float (optional), between 0 and 1
                  Maximum sample ACF which is considered in modeling the ACF. When using this, choose `nLagTrain` large enough, so that
                  `acfThresh` is <= sampled-ACF[`nLagTrain`] 
           AR_opts: dictionary containing options for the AR to be fitted
                  For associated keys, see ARM.fit()
        Outputs:
           outList: a dictionary of the outputs
        """
        methodName='Autoregressive Model (ARM) Estimator'
        #Estimate the mean by sample-mean estimator
        tAvg,fMean,fVarDummy=SME.SME.sampEstim(self.t,self.f,self.opts['conv'])
        n=self.n
        nChunk=self.opts['nChunk']
        #Construct a set of data chunks (no averaging on the chunks)
        #NOTE: each chunk start at the first sample and contains more sample as the chunk number increases
        #  M = N/nChunk
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
               nEndInChunk=n;
            #construct the chunks
            tCh=self.t[:nEndInChunk]
            fCh=self.f[:nEndInChunk]
            nCh=len(fCh)
            tAvgHat.append(tCh[-1]-t0)
            fMeanHat.append(np.mean(fCh))
            #Estimate variance of SME
            fSMEVarHat_=self.varEstim(fCh,self.opts['nLagTrain'],self.opts)
            fSMEVarHat.append(fSMEVarHat_)
        #
        fMean_=fMeanHat[-1]
        fSMEVar_=fSMEVarHat[-1]
        if self.opts['verbose']:
           print('... Estimating uncertainty in SME[f] using ')
           print('    '+methodName)
           print('...... Estimations using all samples:')
           print('       SME=%g, var[SME]=%g, std[SME]=%g' %(fMean_,fSMEVar_,mt.sqrt(fSMEVar_)))
           print('       N x var[SME] = %g' %(n*fSMEVar_))
        #Make a database out of the outputs
        out={'t':self.t,'f':self.f,'tElapsed':tAvg,'tElapsedEstim':tAvgHat,'fSME':fMean,\
             'fSME_var':fSMEVarHat,'method':'ARM','figName':'ARM'+str(M),\
             'tDWN':tAvgHat+t0,'fDWN':fMeanHat};
        self.out=out

###################
# SHOULD be deprecated
###################
#
#////////////////////////////////////////////////////
#deprec ????
def AR_estimator(f,maxlag_=40,ic_='bic',trend_='nc'):
    """
       Autoregressive Estimator
       See http://www.blackarbs.com/blog/time-series-analysis-in-python-linear-models-to-garch/11/1/2016
       ARM: x[i] = a1 x[i-1] + a2 x[i-2] + a3 x[i-3] + ... + ap  x[i-p] +eps
            eps ~ N(0,sig^2)
    """
    #(1) construct a AR Model (estimate coeffcients in ARM)
    arm=smt.AR(f).fit(maxlag=maxlag_,ic=ic_,trend=trend_,missing='raise')
    coefs=arm.params[:]  #[a1,a2,...,ap] 
    #(2) estimate best lag order = p
    best_lag = smt.AR(f).select_order(maxlag=maxlag_, ic=ic_, trend=trend_)
    #(3) Yule-Walker Estimator for autocorrelations and noise sdev in AR model
    coefs_yw,sigma_yw=sm.regression.yule_walker(f,order=best_lag,method='mle')
    #(4) Burg's Estimator for autocorrelations and noise sdev in AR model
    coefs_bg,sigma_bg=sm.regression.linear_model.burg(f,order=best_lag)
    print('est_lag',best_lag)
    print('ARM coefs:',coefs)
    print('rho,sigma (Yule-walker):', coefs_yw, sigma_yw)
    print('rho,sigma (Burg):', coefs_bg, sigma_bg)
    #(4) PACF
    pacf_,ci_=stats_ts.partAutoCorrFunc(f,nlags_=100,method_='ywunbiased',plot=True)

    #(5) Estimate autocovariances
    #(a) sample-estimators

    #(b) using
    print('PACFs=',pacf_)
    return arm,coefs,best_lag
#
#/////////////////////////////////////
#incomplete????
def estimator_AR_ext(t,f,chunkLength,ARpkgPath):
    """ 
       Autoregressive Model to estimate correlations and variance of SME
       This method uses the C++ library by RhysU: https://github.com/RhysU/ar/
    """
    methodName='Autoregressive Model Estimator'
#    figName='AR_estimator'
#    print( printRepeated('-', 60))
#    print('>>>> '+methodName+' <<<<')
    print('... '+methodName)
    n=len(t)

    #SETTINGS
    #nChunk=25;   #no of chunks of input data for which we evaluate estimators    
    nChunk=int(n/chunkLength)+1
    print('...... nChunk=%d' %nChunk)
    #NOTE: nChunk should be chosen so that enough no of samples are each in chunk. otherwise the estimated var can be NaN
    ARpkgPath="/home/salre674/Documents/installedSoftware/AR/ar" #where the AR library is complied
    ARcode='arsel'

    #Assignments and Settings for the AR library
    #create a temp directory in AR-library working directory
    tempDir=ARpkgPath+'/saleh/tempFiles/'
    if os.path.exists(tempDir):
       shutil.rmtree(tempDir)    #delete tempDir if it exists
    os.makedirs(tempDir)   #create tempDir
      
    inDataFileName='timeSeriesData.in'  #argument 2 for bash
    outAR='ARModelResults.out'

    #Estimate the mean by sample-Mean Estimator
    [timAvg,fMean,fVarDummy]=CLSC.sampleMeanEstimator(t,f)

    #construct a set of data chunks (no averaging on the chunks)
    #NOTE: each chunk start at the first sample and contains more sample as the chunk number increases
    #  M = N/nChunk
    #  number of samples in chunks: {M,2M,3M,...,n}
    #  chunk1: 1,2,...,M
    #  chunk2: 1,2,...,M,M+1,...,2M
    M=int(n/nChunk)
    fMeanHat=[];
    fVarAR=[];
    timAvgAR=[];
    fMeanHat.append(f[0]);  #these assignments are required to make CI95 plots statrting from tElapsed=0
    timAvgAR.append(timAvg[0]);
    fVarAR.append(fVarDummy[1])

    for k in range(nChunk):  #chunkCounter
        nEndInChunk=(k+1)*M
        if k==(nChunk-1): #last chunk
           nEndInChunk=n;
        #construct the chunk
        [tCh,fCh]=chunkMaker(t,f,nEndInChunk)
        nCh=len(fCh)
        #write the time-series data in tempDir to be fed into the AR-library    
        inDataFile=open(tempDir+inDataFileName,'w')
        for i in range(nCh):
            inDataFile.write(" %g \n" %(fCh[i]))
        #run the AR model through a bash        
        subprocess.call("bash ./driver/ARModelDriver.sh "+ARpkgPath+" "+inDataFileName+" "+outAR, shell=True)

        #grab the required infor fro mthe results of the AR which are written in outAR
        AROutFile=open(tempDir+outAR,'r')
        ain=AROutFile.readlines();   #read all lines
        ain_sep=[];
        for i in range(len(ain)):
            ain_sep.append(ain[i].split())

        #mu: sample mean estimator
        mu=float(ain_sep[7][2]);  #line 8, column 3
        #mu_sigma: sdev in sample-mean estimator
        mu_sigma=float(ain_sep[8][2]);  #line 8, column 3
        #sigma2x: estimation of the variance
        sigma2x=float(ain_sep[12][2])

        #AR estimators
        fMeanHat.append(mu);
        fVarAR.append(mu_sigma**2.)
        timAvgAR.append(tCh[nCh-1])

        inDataFile.close()
        AROutFile.close()

    print('...... Estimated (at time(n)) mean=%g , var=%g, sdev=%g' %(fMean[len(fMean)-1],fVarAR[len(fVarAR)-1],mt.sqrt(fVarAR[len(fVarAR)-1])))
    print('       N x var(\mu) = %g' %(len(f)*fVarAR[len(fVarAR)-1]))

    #make a data list to plot
    outList={'t':t,'f':f,'tElapsed':timAvg,'tElapsedEstim':timAvgAR,'fMean':fMean,'fVar':[],'fVarEstim':fVarAR,'method':'ARM','figName':'AR'+str(chunkLength),'tDWN':timAvgAR,'fDWN':fMeanHat}; #the last two are dummy
    #print(outList)
#    print(outList['tDWN'])
#    if plot:
#       #>>>> Plot
#       plotter_signal_estimations(outList)
    return outList
#
#//////////////////////////////
#deprec???
def chunkMaker(t,u,nEnd):
    n=len(t);  #total number of samples
    tCh=[];
    uCh=[];
    for i in range(nEnd):
        tCh.append(t[i]);
        uCh.append(u[i]);
    return tCh,uCh
#
#
