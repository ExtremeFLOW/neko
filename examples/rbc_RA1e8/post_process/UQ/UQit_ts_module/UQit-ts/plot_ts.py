#*************************************************************
#    Plotters for UQit time-averge
#*************************************************************
#  Saleh Rezaeiravesh, salehr@kth.se
#-------------------------------------------------------------
#
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
#UQit_ts_path=os.getenv("UQit_ts_path")
import general_ts
#
#//////////////////////////////////////////////
def plotter_signal_estimates(dataList,pltOpts):
    """ 
       Plots a time-series and the estimations for its mean and variance 
    """
    print('...... Plotting UQ estimates for the signal')
    t=dataList['t']  #original signal: (t,f)
    f=dataList['f']
    timAvg=dataList['tElapsed']   #elapsed time when computing the mean (the same size of the original signal)
    if 'tElapsedEstim' in dataList.keys():
        timAvgEstim=dataList['tElapsedEstim']   #elapsed time when estimating the variance
    else:
        timAvgEstim=timAvg;
    fMean=dataList['fSME']       #sequence of estimated mean of f
    fVar=dataList['fSME_var']         #sequence of estimated variance of f (mapped to have the size of uMean)
    if 'fVarEstim' in dataList.keys():
        fVarEstim=dataList['fVarEstim']         #sequence of original estimated variance of u
    else:
        fVarEstim=fVar
    method=dataList['method']
    figName=dataList['figName']
    dwnData=False
    if 'tDWN' in dataList.keys():   #Only for NOBM, OBM, BMBC
       dwnData=True
       tDWN=dataList['tDWN'];
       fDWN=dataList['fDWN'];
       #linear Interpolate fVar and fMean to a vector of size t (to be used only for plotting the CI95%)
       timAvgMap=np.linspace(min(timAvgEstim),max(timAvgEstim),num=len(t))
       fVarMap=general_ts.interp1D(timAvgEstim,fVarEstim,timAvgMap,'linear')
       fMeanMap=general_ts.interp1D(timAvg,fMean,timAvgMap,'linear')
    else:
       timAvgMap=timAvg
       fVarMap=fVar;
       fMeanMap=fMean;
    #>>>> Plot 
    plt.figure(figsize=(10,12));
    #(a) original signal (+ batch signal)
    plt.subplot(311)
    ax=plt.gca();
    plt.plot(t,f,'-k',linewidth=1,label='Original')
    if dwnData:
       if method=='ARM' or method=='DACF':
          sym_='or'
          lab_='Chunk'
          f_=np.mean(f)*np.ones(len(tDWN))
       else:   
          sym_='-r'
          lab_='Batch Mean'
          f_=fDWN
       plt.plot(tDWN,f_,sym_,linewidth=1,markersize='3',label=lab_)
       plt.legend(loc='best')
    plt.xlabel('$t$',fontsize=16)
    plt.ylabel(r'$ f $',fontsize=24)
    plt.tick_params(labelsize=16)
    plt.grid(alpha=0.3)
    # plt.xlim((1e-05, max(timAvg)))
    #(b) estimation of the mean+95% confidence interval
    ci95=1.96*np.sqrt(abs(fVarMap))
    plt.subplot(312)
    ax=plt.gca()
    ax.fill_between(timAvgMap,(fMeanMap + ci95),(fMeanMap -ci95), color='powderblue')
    plt.semilogx(timAvgMap,fMeanMap,'-b',linewidth=2)
    plt.xlabel(r'Elapsed $t$',fontsize=16)
    plt.ylabel(r'$\hat{\mu}\pm{CI}_{95\%}$',fontsize=24)
    plt.tick_params(labelsize=16)
    plt.grid(alpha=0.3)
    # plt.xlim((min(timAvgEstim), max(timAvgEstim)))
    plt.xlim((24, max(timAvgMap)))
    ##plt.xlim((min(timAvgMap),max(timAvgMap)))
    #(c) estimation of the variance
    plt.subplot(313)
    plt.loglog(timAvgEstim[1:],fVarEstim[1:],'-r',linewidth=2)
    plt.xlabel(r'Elapsed $t$',fontsize=16)
    plt.ylabel(r'$\hat{\sigma}^2(\hat{\mu})$',fontsize=24)
    plt.tick_params(labelsize=16)
    plt.grid(alpha=0.3)
    plt.xlim((min(timAvgEstim),max(timAvgEstim)))

    if 'title' in pltOpts.keys():
        plt.suptitle(pltOpts['title'],fontsize=16)
    # 
    plt.subplots_adjust(hspace=0.3)
    # save the figure
    if 'saveFig' in pltOpts.keys() and pltOpts['saveFig']:
       fig = plt.gcf()
       DPI = fig.get_dpi()
       fig.set_size_inches(1000.0/float(DPI),+900.0/float(DPI))
       figDir_=pltOpts['figDir']   #output path
       if not os.path.exists(figDir_):
          os.makedirs(figDir_)
       figOut=figDir_+'/'+pltOpts['figName']+'.pdf'
       plt.savefig(figOut,bbox_inches='tight')
    plt.show()
    return timAvgMap,fMeanMap,ci95,timAvgEstim,fVarEstim  #in case comparison is needed
#
#/////////////////////////////////
def histoPlot_single(t,f,pltOpts):
    """ 
       plot time-series along with associated histogram
         Inputs:
           t: 1d numpy array of size n containing time samples
           f: 1d numpy array of size n containing sample values
           pltOpts: plot options containing the following keys:
             
         Output:           
    """
    plt.figure(figsize=(10,7))
    #(a) plot f vs t
    plt.subplot(2,1,1)
    plt.plot(t,f,'-',color='steelblue',label='Samples')
    X=np.linspace(min(t),max(t),10)
    plt.plot(X,np.mean(f)*np.ones(10),'--k',linewidth=2,label=r'$\hat{E}[f]$')
    plt.xlabel(r'$t$',fontsize=17)
    plt.ylabel(r'$f$',fontsize=17)
    plt.tick_params(labelsize=16)
    plt.legend(loc='best',fontsize=15)    
    #(b) plot histogram of f
    plt.subplot(2,1,2)
    plt.hist(f,bins=50,density=True,color='steelblue',alpha=0.4,edgecolor='b',label='histogram')
    kde = sm.nonparametric.KDEUnivariate(f)
    kde.fit()
    plt.plot(kde.support,kde.density,'-r',lw=2,label='kde')
    plt.axvline(np.mean(f),color='k', linestyle='dashed', linewidth=2,label=r'$\hat{E}[f]$')
    plt.xlabel(r'$f$',fontsize=17)
    plt.ylabel('PDF',fontsize=17)
    plt.tick_params(labelsize=16)
    plt.legend(loc='best',fontsize=15)   
    plt.subplots_adjust(hspace=0.3)
    #save fig
    if 'saveFig' in pltOpts.keys() and pltOpts['saveFig']:
       fig = plt.gcf()
       DPI = fig.get_dpi()
       fig.set_size_inches(400/float(DPI),800/float(DPI))
       figDir=pltOpts['figDir']
       if not os.path.exists(figDir):
          os.makedirs(figDir)
       plt.savefig(figDir+pltOpts['figName']+'.pdf',bbox_inches='tight')
    plt.show()
#       
