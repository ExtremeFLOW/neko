#*****************************************************************************
# Analysis of time-series taken from turbulent flow simulations
#*****************************************************************************
# Saleh Rezaeiravesh, salehr@kth.se
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
import time
import sys
import os
import math as mt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 15,
          'legend.loc':'best',
#          'figure.figsize': (14,5),
          'lines.markerfacecolor':'none',
         'axes.labelsize': 15,
         'axes.titlesize': 15,
         'xtick.labelsize':15,
         'ytick.labelsize':15,
         'grid.alpha':0.3}
pylab.rcParams.update(params)


from UQit.stats import pdfFit_uniVar
#
#
sys.path.append(os.getenv("UQit_ts_path"))
sys.path.append(os.getenv("UQit_ts_app_path")+"/extDataHandler_tUQ/")
#
#import ARM
#import DACF
#import SME_uncert
from SME import SME,SMEuncert
import dbMakerChan_tUQ
from ARM import ARM
#import synDataGen_tUQ
#import write_tUQ
#import stats_tUQ
#import plot_tUQ
#import general_tUQ
#
#####################
# MAIN
#####################
os.system('clear')
#fLog=open('tUQ.log','w')  #logfile: not working now
#----- SETTINGS -----------------------
workDir='/home/saleh/Desktop/temp_noBackUp/uqTimeErrorBars/outputs/'  #working directory to dump figs and data
#--------------------------------------

ARpkgPath='/home/salre674/Documents/installedSoftware/AR/ar/'  #where the AR package is located
#
#
################
# MAIN
################
#default set of TS data (they can be overwritten in an external function)
#-------- SETTINGS ----------------
#Directory contains the data (both raw and processed databases)
dataDir="/home/saleh/Desktop/DATA_mustBackUp/uqTimeErrorBars/nekRuns/chan300/C8gp3_trans"
#raw database name
rawData="uq_timeseries.dat"
#caseName (when dumping figs/outdata)
caseName='C8gp3_trans'
doPickle=True   #True: if you want to create a pickle from raw channel time series
#----------------------------------
#(1) create pickle databases on the disk from rawdata achieved by running nek5000
if doPickle:
   dbMakerChan_tUQ.createPickleDB(dataDir,rawData)
#(2) read in existing databases
tim,y,uTau,vel=dbMakerChan_tUQ.readPickleDB(dataDir+'/database')
print(vel[0].keys())
print('a1',tim.shape)
ypls=y*298
ny=len(y)
n=len(tim)

#test plot (plot uTau-t, Uc-t)
if 0==0:
   t=tim 
   plt.figure(figsize=(12,5))
   plt.subplot(2,1,1)
   plt.plot(t,uTau,'-k')
   plt.plot(t[60000:],uTau[60000:],'-r')
   plt.xlabel(r'$tU_b/\delta$',fontsize=15)
   plt.ylabel(r'$u_\tau/U_b$',fontsize=15)
   plt.ylim([0.0525,0.065])
   plt.vlines(300,0.0525,0.065,colors='b',linestyle='dashed')
   plt.grid()
   plt.subplot(2,1,2)
   plt.plot(t,vel[-1]["u"],'-k')
   plt.xlabel(r'$tU_b/\delta$',fontsize=15)
   plt.ylabel(r'$u_c/U_b$',fontsize=15)
   plt.ylim([1.1,1.2])
   plt.vlines(300,1.1,1.2,colors='b',linestyle='dashed')
   plt.grid()
   plt.show()


print('ss',t.shape)
#
#
###################
# External Funcs
###################
def chan_ARfit_diffMethods():
    """
    Test different methods for fitting AR to channel flow time-series samples
    """
    #----- SETTINGS
    maxLag=100
    #--------------------------
    #(0) Training data
    f=vel[40]["u"]/uTau
    print('y+=',ypls[40])
    #use different methods to fit AR
    #p_opts={'maxlag':40}
    p_opts={'p':maxLag}
    # Unconditional MLE
    fit_opts=p_opts.copy()
    fit_opts.update({'AR-method':'umle'})
    fit_coefs,fit_sig2d,fit_opts=ARM.AR_fit(f,fit_opts)
    # OLS method
    ols_opts={'p':fit_opts['p']}
    ols_opts.update({'AR-method':'ols'})
    ols_coefs,ols_sig2d,ols_opts=ARM.AR_fit(f,ols_opts)
    # Burg's method
    bg_opts={'p':fit_opts['p']}
    bg_opts.update({'AR-method':'burg'})
    bg_coefs,bg_sig2d,bg_opts=ARM.AR_fit(f,bg_opts)
    #YW method
    yw_opts={'p':fit_opts['p']}
    yw_opts.update({'AR-method':'yw'})
    yw_coefs,yw_sig2d,yw_opts=ARM.AR_fit(f,yw_opts)
    #(3) plot
    plt.figure(figsize=(10,4))
    plt.plot(fit_coefs,'+m',ms=9,mfc='none',label='OLS')
    plt.plot(fit_coefs,'ob',ms=9,mfc='none',label='Unconditional MLE')
    plt.plot(yw_coefs,'sg',ms=9,mfc='none',label='Yule-Walker eqns')
    plt.plot(bg_coefs,'xr',ms=9,label='Burg\'s method')
    plt.ylabel('AR coefficients',fontsize=17)
    plt.xlabel('Lag-1',fontsize=15)
    plt.tick_params(labelsize=16)
    plt.legend(loc='best',fontsize=15)
    plt.grid(alpha=0.3)
    plt.show()
#
#//////////////////////////
def chan_ARfit_diff_y():
    """
    Compare the coeffcients of the AR fitted to TS of flow quantities at different y's    
    """
    #(3) select the time-series to analyze
    fig=plt.figure(figsize=(15,7))
    for i in range(1,ny,1):
        #specify quantitiy and its TS samples
        f=vel[i]["u"]/uTau
        #Fit an ARM using different methods
        p_opts={'p':100}   #max lag
        bg_opts={'p':p_opts['p']}
        bg_opts.update({'AR-method':'burg'})
        bg_coefs,bg_sig2d,bg_opts=ARM.AR_fit(f,bg_opts)
        #
        h_=1+np.arange(len(bg_coefs))
        ltyp_='-'
        if i > 36:
           ltyp_='--' 
        ypls_=f"{ypls[i]:.1f}"
        plt.semilogx(h_,bg_coefs,ltyp_,ms=9,label=r'$y^+$='+ypls_)
    plt.grid(alpha=0.4)
    plt.legend(loc='best',ncol=3,fontsize=14)
    plt.xlabel('Lag',fontsize=16)
    plt.ylabel('ARM Coefficients',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.show()
#
#
def chan_ARfit_convergence():
    """
       Convergence of ARM coeffcients with the size of the training samples
    """
    #----- SETTINGS
    nTrain=[2000,5000,7500,10000,20000,30000,40000,50000,70000,80000,n]     #size of training samples, should be less than n
    maxLag=100 #=p, order of ARM
    #----------------------------
    #(0) Training data
    f=vel[40]["u"]/uTau
    #
    ar_coefs=[]
    ar_noiseVar=[]
    for n_ in nTrain:
        # select training samples
        f_=f[:n_]
        #Construct AR based on traiing data
        p_opts={'p':maxLag}   #max lag
        bg_opts={'p':p_opts['p']}
        bg_opts.update({'AR-method':'burg'})
        bg_coefs,bg_sig2d,bg_opts=ARM.AR_fit(f_,bg_opts)
        ar_coefs.append(bg_coefs)
        ar_noiseVar.append(bg_sig2d)
    ar_coefs=np.asarray(ar_coefs)    
    plt.figure(figsize=(10,4))
    plt.semilogx(nTrain,ar_coefs[:,0:10])   
    plt.ylabel('ARM Coefficients',fontsize=16)
    plt.xlabel(r'Sample Size $n$',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)    
    plt.grid(alpha=0.3)
    plt.show()
#
#
#///////////////////////////
import statsmodels.api as sm
def chan_ARfit_pred_uncert():
    """
    Compare the uncertainty in MSE estimated using the original samples and 
       the synthetic samples generated from the AR fitted to the original samples.
    """
    #----- SETTINGS
    nTrain=30000     #size of training samples, should be less than n
    nTest=300000     #size of test samples
    nInit=200        #number of initial samples from original samples when using ARM for prediction
    maxLag=100 #=p, order of ARM
    #----------------------------
    #(0) Training data
    f_=vel[40]["u"]/uTau
    f=f_[:nTrain]
    #(1) Construct AR based on traiing data
    p_opts={'p':maxLag}   #max lag
    bg_opts={'p':p_opts['p']}
    bg_opts.update({'AR-method':'burg'})
    bg_coefs,bg_sig2d,bg_opts=ARM.AR_fit(f,bg_opts)
    #(2) Draw nTest samples from a fitted ARM
    f_mean=bg_opts['orig_sampMean']   #mean of the original samples
    gInit=f[:nInit]-f_mean
    g=ARM.AR_sample(gInit,nTest,bg_coefs,bg_sig2d,{})    
    #(3) plot original TS and ARM test samples
    plt.figure(figsize=(15,4))
    plt.plot(f,'-b',lw=3,label='Training Sample')
    plt.plot(f_mean*np.ones(nTrain),'--b')
    plt.plot(f_,'-k',lw=1,label='Flow Simulation')
    plt.plot((np.mean(f_))*np.ones(n),'--k')
    plt.plot(g+f_mean,'-r',label='AR Emulation')
    plt.plot((np.mean(g)+f_mean)*np.ones(nTest),'--r')
    plt.xlabel('Lag',fontsize=16)
    plt.ylabel('ARM Coefficients',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(alpha=0.3)
    plt.legend(loc='best',fontsize=17)
    plt.show()
    #(3b) plot histogram of the original and synthetic samples
    #training
    plt.figure(figsize=(8,5))
    plt.hist(f_,bins=50,density=True,color='steelblue',alpha=0.4,edgecolor='b',label='Training')
    kde = sm.nonparametric.KDEUnivariate(f)
    kde.fit()
    plt.plot(kde.support,kde.density,'-b',lw=2)
    #flow simulation
    plt.hist(f_,bins=50,density=True,color='grey',alpha=0.3,edgecolor='k',label='Flow Simulation')
    kde = sm.nonparametric.KDEUnivariate(f_)
    kde.fit()
    plt.plot(kde.support,kde.density,'-k',lw=2)
    #AR test emulation
    plt.hist(g+f_mean,bins=50,density=True,color='salmon',alpha=0.4,edgecolor='r',label='AR Emulation')
    kde = sm.nonparametric.KDEUnivariate(g+f_mean)
    kde.fit()
    plt.plot(kde.support,kde.density,'-r',lw=2)
    plt.legend(loc='best',fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.show()

    #(4) Compare uncerainties in SME
    #NOBM Estimator
    out_nobm=SME_uncert.uncertEstimator(t[:nTrain],f,{'method':'NOBM','batchSize':1000})
    plot_tUQ.plotter_signal_estimates(out_nobm,{})
    #ARM Estimator
    nChunk_=20
    nLagTrain=200
    outAR=ARM.estimator_tUQ(t[:nTrain],f,nChunk_,nLagTrain,bg_opts)
    plot_tUQ.plotter_signal_estimates(outAR,{})
    #DACF Estimator
    nChunk_=20
    out_DACF=DACF.estimator_tUQ(t[:nTrain],f,nChunk_)
    plot_tUQ.plotter_signal_estimates(out_DACF,{})   
#
#
#///////////////////////////
def chan_AR_modelACF():
    """
    Model ACF for channel flow TS data
    """
    #(0) Training data
    f=vel[56]["u"]/uTau
    #(1) fit AR
    p_opts={'p':100}
    bg_opts={'p':p_opts['p']}
    bg_opts.update({'AR-method':'burg'})
    bg_coefs,bg_sig2d,bg_opts=ARM.AR_fit(f,bg_opts)
    #(2) model ACF
    nLagTrain=200  # max-lag of the training ACFs when fitting a power-law model
    acf_samp=stats_tUQ.autoCorrFunc(f,nlags_=nLagTrain,fft_=False,plot=False)
    C,lambda_=ARM.acf_powModel_fit(bg_coefs,bg_opts['p'],acf_samp)
    acf_model=ARM.acf_powModel_pred(C,lambda_,n)
    #(3) plot
    plt.figure(figsize=(15,5))
    acf_samp_full=stats_tUQ.autoCorrFunc(f,nlags_=int(n),fft_=False,plot=False)
    plt.plot(acf_samp,'ob',mfc='none',ms=6,label='Training sample with max-lag, K=%d'%(nLagTrain))
    plt.plot(acf_samp_full,'-',color='skyblue',lw=2,label='Full-lag sample-estim. (reference)')
    plt.semilogx(np.arange(1,n,1),acf_model[1:],'-r',lw=2,label='Modelled')
    plt.plot(acf_samp_full[:nLagTrain],'o',color='skyblue',ms=2)
    plt.ylabel('ACF',fontsize=16)
    plt.xlabel('Lag',fontsize=16)
    plt.tick_params(labelsize=16)
    plt.grid(alpha=0.3)
    plt.legend(loc='best',fontsize=15)
    plt.show()
#
#    
#///////////////////////////
def chan_Profile_ARM_uncert():
    """
       Using ARM to predict uncertainty at all points along a profile
    """
    #----- SETTINGS
    qoiName="u"
    maxLag=100 #=p, order of ARM
    nLagTrain=200  #number of lags to train modeled ACF (200 seems to be optimal)
    #----------------------------
    f_SME=[]   #sample mean of qoi
    ci_SME=[] #95% CI of the sample mean
    ypls_=[]   #ypls

    t=tim

    for i in range(1,ny):
        #(0) Training data
        f_=vel[i][qoiName]#/uTau
        f=f_
        #a subset of time-series
        f=f_[60000:]
        t=tim[60000:]

        f_SME.append(np.mean(f))
        ypls_.append(ypls[i])

        print(i)

        sme_=SME(t,f,{'verbose':False,'conv':False})  #new
        out_=SMEuncert(sme_,{'method':'ARM','nChunk':1,'nLagTrain':nLagTrain,'p':maxLag,'AR-method':'burg'}).estim
#        out_=SMEuncert(sme_,{'method':'NOBM','batchSize':40000}).estim

        var_=out_['fSME_var']
        print(i,var_)


        #(1) Set the options of AR(maxLag)
#        p_opts={'p':maxLag}   #max lag
#        AR_opts={'p':p_opts['p']}
#        AR_opts.update({'AR-method':'burg'})
#        #(2) estimate variance of uncertainty in the SME fitting an AR
#        var_=ARM.varEstim(f,nLagTrain,AR_opts)


        ci_SME.append(1.96*mt.sqrt(abs(var_[0])) * 20)

    ci_SME=np.asarray(ci_SME)    
    print(ci_SME)
    f_SME=np.asarray(f_SME)    
    ypls_=np.asarray(ypls_)
    #(3) plot   
    plt.figure(figsize=(10,3))
    plt.semilogx(ypls_,f_SME,'-b',label=r'$\hat{\mathbb{E}}[\cdot]$')    
    ax=plt.gca()
    ax.fill_between(ypls_,f_SME-ci_SME,f_SME+ci_SME,color='salmon',alpha=0.4,label=r'$20\times 95\% CI$')
    plt.legend(loc='best',fontsize=20)
    plt.xlabel(r'$y^+$',fontsize=20)
    plt.ylabel(r'$\langle u\rangle$',fontsize=20)
    #plt.ylabel(r'$\langle u^\prime v^\prime \rangle^+$',fontsize=20)
    #plt.ylabel(r'$v^{\prime^+}_{rms}$',fontsize=20)
    plt.grid(alpha=0.3)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.show()
#
#
#///////////////////////////
#import NOBM
def chanDataAnalysis_errorEstim():
    """
        Estimate uncertainties in the time-averaged data of channel flow
    """    
    print('>>> chanDataAnalysis()')
    #-------- SETTINGS ----------------
    #Directory contains the data (both raw and processed databases)
    dataDir="/home/saleh/Desktop/temp_noBackUp/uqTimeErrorBars/nekRuns/chan300/C8_25"
    #raw database name
    rawData="uq_timeseries.dat"
    #caseName (when dumping figs/outdata)
    caseName='C8_25'
    doPickle=False   #True: if you want to create a pickle from raw channel time series
    #----------------------------------
    #(1) create pickle databases on the disk from rawdata achieved by running nek5000
    if doPickle:
       dbMakerChan_tUQ.createPickleDB(dataDir,rawData)

    #(2) read in existing databases
    t,y,uTau,vel=dbMakerChan_tUQ.readPickleDB(dataDir+'/database')
    ny=len(y)
    nt=len(t)
    #
#    plt.plot(t,uTau)
#    plt.show()
    f_=uTau
#    print(vel[50]["u'"].shape)
    #(2) Preprocessing
    # (a) compute autocorrelations
    ##r=stats_tUQ.crossCorr(f_,f_,'full')
    #print('r',r)
    #print(r.size)
    r=stats_tUQ.autoCorrFunc(f_,nt-2,plot=False)
    #plt.plot(t,r,'ob')
    #plt.show()
    print(r)
    time.sleep(5)

    ### plot signal and NOBM reduced version 
    #PLAN: use this part to compute effect of autocorrelation 
    Kbm_,tbm_,fbm_=NOBM.NOBM_batchGen(t,f_,1000)
    plt.figure(figsize=(10,4))
    plt.plot(t,f_,'-b',lw=2,label='original signal')
    plt.plot(tbm_,fbm_,'--or',label='NOBM signal')
    plt.legend(loc='best')
    plt.xlabel('Time')
    plt.ylabel('Signal')
    plt.grid(alpha=0.3)
    plt.show()
    rbm=stats_tUQ.autoCorrFunc(fbm_,Kbm_-1,plot=False)

    print('size of the orignal signal:',len(t),len(r))
    print('size of the NOBM signal:',len(tbm_),len(rbm))
    plt.figure(figsize=(10,4))
    plt.plot(t[1:],r,'-ob',label='Original signal')
    plt.plot(tbm_,rbm,'-xr',label='NOBM signal')
    plt.legend(loc='best')
    plt.xlabel('Time Lag')
    plt.ylabel('ACF')
    plt.grid(alpha=0.3)
    plt.show()

    for it in range(30,30):  #time
        u_=np.array([vel[i]['u'][it] for i in range(Ny)])
        for i in range(ny):
            plt.semilogx(y[i],u_/uTau[it],'.')
    plt.show()

    #(3) run the UQ estimator for all points on a channel profile
    batchSize=1000
    qoiName='u'
    ny=len(y)
    sdev_list=[]
    mean_list=[]
    pltOpts={'figDir':workDir+caseName,   #plot options
             'figName':caseName+'_nobm'}
    iPltList=[5,25,40,50,60]   #y-locations at which we plot signal processing
    for i in range(ny):
        #Estimate uncertainties
        f_=vel[i][qoiName]/uTau
        out_nobm=SME_uncert.uncertEstimator(t,f_,{'method':'NOBM','batchSize':batchSize})
        #list of estmated mean and sdev of \mu
        mean_list.append(out_nobm['fMean'][-1])
        sdev_list.append(mt.sqrt(out_nobm['fVarEstim'][-1]))
        if i in iPltList:
           print(i,y[i]*298)
           pltOpts.update({'figName':pltOpts['figName']+'_'+str(i)}) 
           plot_tUQ.plotter_signal_estimates(out_nobm,pltOpts)
    mean_list=np.asarray(mean_list)
    sdev_list=np.asarray(sdev_list)
    #
    plt.figure(figsize=(10,6));
    ax=plt.gca();
    scl_=20   #magnify uncertainty for visualization purpose
    ax.fill_between(y*298,mean_list-1.96*sdev_list*scl_,mean_list+1.96*sdev_list*scl_,color='skyblue',alpha=0.4,label=r'$95\%CI \times 20$')    
    plt.semilogx(y*298,mean_list,'-b',label='Estimated Mean')    
    plt.legend(loc='best',fontsize=17)
    plt.xlabel(r'$y^+$',fontsize=18)
    plt.ylabel(r'$\langle u\rangle^+$',fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(alpha=0.3)
    plt.show()
    #
    for i in range(ny):
        print(y[i])
    print('-----------------')
    for i in range(ny):
        print(sdev_list[i]/(mean_list[i]+1e-6))
    plt.figure(figsize=(10,5))
    plt.semilogx(300*y,sdev_list/(mean_list+1e-6),'-r')
    plt.xlabel(r'$y$',fontsize=16)
    plt.ylabel(r'$\hat{\sigma}[\langle u \rangle]/\hat{\mathbb{E}}[\langle u \rangle]$',fontsize=16)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(alpha=0.4)
    plt.show()

    #(4) plot estimates for the signal (last on the above loop)
#    pltOpts={'figDir':workDir+caseName,   #plot options
#             'figName':caseName+'_nobm'}
#    plot_tUQ.plotter_signal_estimates(out_nobm,pltOpts)
#
#///////////////////////
def plot_ACF_ARM_y():
    """
    Plot ACF (sampled and modeled) and ARM coeffcients for channel flow data
    at different wall-distances    
    """
    #----SETTINGS
    qoi="u"   #QoI
    pMax=200  #order of ARM
    nLagTrain=200  # max-lag of the training ACFs when fitting a power-law model
    #-------------------
    #
    coef1_list=[] #list of first coeffcients of the ARM at any wall-distance
    yp_list=[]
    lagACFdev_list=[] #list of first big deviation of ACF from unity
    for i in range(1,ny,1):
        plt.figure(figsize=(12,10))
        #specify quantitiy and its TS samples
        f=vel[i][qoi]#/uTau   #QoI
        #1. Fit an ARM using different methods
        p_opts={'p':pMax}   #max lag
        bg_opts={'p':p_opts['p']}
        bg_opts.update({'AR-method':'burg'})
        bg_coefs,bg_sig2d,bg_opts=ARM.AR_fit(f,bg_opts)
        coef1_list.append(bg_coefs[0])
        yp_list.append(y[i]*np.mean(uTau)/1.9922E-4)
        #2. Sample ACF
        acf_samp_full=stats_tUQ.autoCorrFunc(f,nlags_=f.size,fft_=False,plot=False)
        # find first lag at which the ACF becomes deviated from 1.0 more than eps
        for j in range(200):
            if abs(1.0-acf_samp_full[j])>0.005:
               iACF_dev=j               
               break     
        lagACFdev_list.append(j)   
        #3. Fit ACF using ARM
        acf_samp=stats_tUQ.autoCorrFunc(f,nlags_=nLagTrain,fft_=False,plot=False)
        C,lambda_=ARM.acf_powModel_fit(bg_coefs,bg_opts['p'],acf_samp)
        acf_model=ARM.acf_powModel_pred(C,lambda_,n)
        #
        h_=1+np.arange(len(bg_coefs))
        ltyp_='-b'
        #if i > 36:
        #   ltyp_='--' 
        ypls_=f"{ypls[i]:.1f}"
        #Plots
        #
        plt.subplot(2,1,1)
        plt.plot(acf_samp_full,'-k',label='Sample-estimated')
        plt.semilogx(acf_model,'--r',mfc='none',lw=2,label='Modeled')
        print('A',iACF_dev)
        plt.axvline(x=iACF_dev,color='m',ls='-.')
        plt.legend(loc='best',fontsize=13)
        plt.grid(alpha=0.4)
        plt.xlabel('Lag',fontsize=18)
        plt.ylabel('ACF',fontsize=18)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.ylim([-0.3,1.1])
        #
        plt.subplot(2,1,2)
        plt.semilogx(h_,bg_coefs,ltyp_,ms=9,lw=2,label=r'$y^+$='+ypls_)
        plt.axvline(x=iACF_dev,color='m',ls='-.')
        plt.legend(loc='upper right',fontsize=18)
        plt.grid(alpha=0.4)
        plt.xlabel('Lag',fontsize=18)
        plt.ylabel('ARM Coefficients',fontsize=18)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.ylim([-0.21,0.81])
        #save the fig
        if 0==0:
           acf_samp_full=stats_tUQ.autoCorrFunc(f,nlags_=f.size,fft_=False,plot=False)
           fig = plt.gcf()
           DPI = fig.get_dpi()
           fig.set_size_inches(800/float(DPI),900/float(DPI))
           figDir='../testFigs/'+caseName+'/ARM_ACF_'+qoi+'_y/'
           if not os.path.exists(figDir):
              os.makedirs(figDir)
           plt.savefig(figDir+caseName+'_'+qoi+'_'+str(i)+'.png',bbox_inches='tight')           
        else:  
           plt.show()
    #       
    # plot2: 1st ARM coef vs. y+
    plt.figure(figsize=(12,10))
    plt.subplot(2,1,1)
    plt.semilogx(yp_list,coef1_list,'-ob')
    plt.xlabel(r'$y^+$',fontsize=18)
    plt.ylabel('First ARM Coefficient',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(alpha=0.4)
    # plot3: 1st ARM coef vs. y+
    plt.subplot(2,1,2)
    plt.semilogx(yp_list,lagACFdev_list,'-.sm')
    plt.xlabel(r'$y^+$',fontsize=18)
    plt.ylabel('Max p (by ACF)',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(alpha=0.4)
    #save the fig
    if 0==0:
       acf_samp_full=stats_tUQ.autoCorrFunc(f,nlags_=f.size,fft_=False,plot=False)
       fig = plt.gcf()
       DPI = fig.get_dpi()
       fig.set_size_inches(800/float(DPI),1000/float(DPI))
       figDir='../testFigs/'+caseName+'/ARM_ACF_'+qoi+'_y/'
       if not os.path.exists(figDir):
          os.makedirs(figDir)
       plt.savefig(figDir+caseName+'_'+qoi+'_'+str(i)+'pMAX'+'.png',bbox_inches='tight')           
    else:  
       plt.show()
#
#///////////////////////////
def chan_profs_2ndOMoments_Uncert():
    """
       Accurate estimation of uncertainty in 2nd-order statsitical moments and rms of channel flow QoIs
       E[u]   = <u>   + e_u,     where e_u~N(0,sig_u)
       E[u^2] = <u^2> + e_{u^2}, where e_{u^2}~N(0,sig_{u^2}) 
       => 
       E[u'^2] = E[u^2]-E[u]^2 + e_{u^2} - e_u^2 - 2 E[u]*e_u
    """
    def cumSME(X):
        """
        cumulative SME:
        <X>_k=k^{-1} \sum_{1}^k X_k
        """
        N=len(X)
        Xn=np.zeros(N)
        sum_=X[0]
        Xn[0]=sum_
        for i in range(1,N):
            sum_+=X[i]
            Xn[i]=sum_/float(i+1)
        return Xn    

    #----- SETTINGS
    qoiName="u"
    maxLag=100 #=p, order of ARM
    nLagTrain=200  #number of lags to train modeled ACF (200 seems to be optimal)
    #----------------------------
    f_mean=[]   #sample mean of qoi
    f_ci=[]  #95% CI of the sample mean
    f2_mean=[]   #sample mean of qoi
    f2_ci=[]  #95% CI of the sample mean
    f2prime_mean=[]
    f2prime_ci=[]
    frms_mean=[]
    frms_ci=[]
    ypls_=[]   #ypls


    for i in range(1,ny,1):
#    for i in range(40,46):
        #(0) Training data
        f=vel[i][qoiName]#/uTau   
        f2=vel[i][qoiName+qoiName]#/uTau**2.       

        #compute correlation matrix of cumlative SME of f and f^2
        Fn=cumSME(f)
        F2n=cumSME(f2)
        corrMat=np.corrcoef(Fn,F2n) 
        #print(corrMat)

        #compute SME over all samples
        f_mean.append(np.mean(f))
        f2_mean.append(np.mean(f2))
        #print('A',np.mean(f2)-np.mean(f)**2.,np.mean(f2-np.mean(f)**2.),mt.sqrt(np.mean(f2-np.mean(f)**2.)))
        ypls_.append(ypls[i])

        #(1) Set the options of AR(maxLag)
        p_opts={'p':maxLag}   #max lag
        AR_opts={'p':p_opts['p']}
        AR_opts.update({'AR-method':'burg'})

        #(2) estimate variance of uncertainty in the SME of f and f^2 
        var_f_=ARM.varEstim_ARM(f,nLagTrain,AR_opts)
        var_f2_=ARM.varEstim_ARM(f2,nLagTrain,AR_opts)
        f_ci.append(1.96*mt.sqrt(abs(var_f_)))
        f2_ci.append(1.96*mt.sqrt(abs(var_f2_)))

        #(3) estimate mean and uncertainty in f'^2
        f2prm_mean=f2_mean[-1]-f_mean[-1]**2.
        #we generate MC samples to compute uncertainty in the SME of f'^2
        #Note that the samples of SME of f and f^2 are correlated
        nMC=100000
        corrMat[0,0]=corrMat[0,0]*var_f_
        corrMat[1,1]=corrMat[1,1]*var_f2_
        corrMat[0,1]=corrMat[0,1]*mt.sqrt(var_f_)*mt.sqrt(var_f2_)
        corrMat[1,0]=corrMat[1,0]*mt.sqrt(var_f_)*mt.sqrt(var_f2_)
        mean_=np.asarray([0,0])
        e_f,e_f2=np.random.multivariate_normal(mean_, corrMat, nMC).T
        #plt.plot(e_f,e_f2,'.',alpha=0.2)
        #plt.show()

        e_f2prm=e_f2-e_f**2.-2*f_mean[-1]*e_f
        #kde=pdfFit_uniVar(e_f2prm,True,{})
        f2prime_mean.append(f2prm_mean+np.mean(e_f2prm))
        #f2prime_mean.append(f2prm_mean)
        f2prime_ci.append(1.96*np.std(e_f2prm))
        #rms mean and associated ci
        #e_rms=np.sqrt(abs(f2prm_mean+e_f2prm))
        #e_rms=np.sqrt(abs(f2prm_mean))
        e_rms=f2_mean[-1]-f_mean[-1]**2.
        frms_mean.append(np.mean(e_rms)) 
        frms_ci.append(1.96*np.std(e_rms)) 
        print('var_f, var_f2: ',mt.sqrt(var_f_),mt.sqrt(var_f2_))
        print("f'^2 mean, e_f'^2 mean, e_f'^2_std",f2prm_mean,np.mean(e_f2prm),np.std(e_f2prm))
#        e_rms_std=np.std(e_f2prm)
#        print(mt.sqrt(f2prm_mean),0.5*e_rms_mean/mt.sqrt(f2prm_mean),0.5*e_rms_std/mt.sqrt(f2prm_mean))
        print('mean_erms, std_erms: ',np.mean(e_rms),np.std(e_rms))
#
#   
    f_mean=np.asarray(f_mean)    
    f2_mean=np.asarray(f2_mean)    
    f_ci=np.asarray(f_ci)    
    f2_ci=np.asarray(f2_ci)    
    f2prime_mean=np.asarray(f2prime_mean)
    f2prime_ci=np.asarray(f2prime_ci)
    frms_mean=np.asarray(frms_mean)
    frms_ci=np.asarray(frms_ci)

    ypls_=np.asarray(ypls_)
    #(3) plot   
    plt.figure(figsize=(10,8))
    plt.subplot(2,1,1)
    ax=plt.gca()
    plt.semilogx(ypls_,f2prime_mean,'-b',label=r'$\hat{\mathbb{E}}[\cdot]$')    
    ax.fill_between(ypls_,f2prime_mean-f2prime_ci,f2prime_mean+f2prime_ci,color='skyblue',alpha=0.4,label=r'$20\times 95\% CI$')
    plt.ylabel(r'$\langle u{^\prime^2}\rangle^+$',fontsize=20)
    plt.subplot(2,1,2)
    ax=plt.gca()
    plt.semilogx(ypls_,np.sqrt(frms_mean),'-b',label=r'$\hat{\mathbb{E}}[\cdot]$')    
    ax.fill_between(ypls_,frms_mean-frms_ci,frms_mean+frms_ci,color='skyblue',alpha=0.4,label=r'$20\times 95\% CI$')

    plt.legend(loc='best',fontsize=20)
    plt.xlabel(r'$y^+$',fontsize=20)
#    #plt.ylabel(r'$\langle u^\prime v^\prime \rangle^+$',fontsize=20)
    plt.ylabel(r'$u^{\prime^+}_{rms}$',fontsize=20)
    plt.grid(alpha=0.3)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.show()
