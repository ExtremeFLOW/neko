#****************************************************
# Generating synthetic time-series (TS) samples
#****************************************************
# Saleh Rezaeiravesh, salehr@kth.se
#----------------------------------------------------
#
import sys
import os
import numpy as np
import math as mt
import matplotlib.pyplot as plt
#UQit_ts_path=os.getenv("UQit_ts_path")
#sys.path.append(UQit_ts_path+'/write_tUQ/')
#ys.path.append('./')
#import __init__
import write_ts



class syntheticData:
    """
    Generating a synthetic time-series of size n
    Inputs:
       n: int, size of the time-series
       opts: options 
          {'type':[],         #method of geneating the TS
           'dt':scalar or numpy array,  #time step size of the TS samples
                scalar: dt is fixed for all samples
                1d numpy array: variable dt for samples, size = n-1
           'noiseSdev':       #float, noise standard deviation
           'writeInFile':[],  #bool, if True, TS data are written in file
           'writeFile':[]}    #File to write the TS data
    Outputs:               
       t2: 1d numpy array of size n, time 
       u2: 1d numpy array of size n, values of the TS 
    """
    def __init__(self,n,opts):
        self.n=n
        self.opts=opts

    def synData1(self,N,noiseSdev_,dt):
        t=np.zeros(N)
        u=np.zeros(N)
        t[0]=0.0
        t[1]=t[0]+dt[0]
        u[0]=0.0
        u[1]=0.1
        r=np.random.normal(size=N)   #samples with ~ N(0,1)
        for i in range(2,N):
            t[i]=t[i-1]+dt[i-1]
            u[i]=0.2*u[i-2]+0.8*u[i-1]+noiseSdev_*r[i]
        return t,u

    def synData2(self,N,noiseSdev_,dt):
        t=np.zeros(N)
        u=np.zeros(N)
        t[0]=0.0
        #u[0]=0.5
        u[0]=np.random.rand()
        r=np.random.uniform(0.,1.,N)   #samples with U[0,1]
        for i in range(1,N):
            t[i]=t[i-1]+dt[i-1]
            #Eq.(33) in Russo-Luchini 2017
            u[i]=0.9*u[i-1]+noiseSdev_*r[i]
        return t,u    
        
    def synData3(self,N,noiseSdev_,dt,opts):
        t=np.zeros(N)
        u=np.zeros(N)
        if 'coefs' in opts.keys():
           coefs_=opts['coefs']     
        else:   
           #default, from Example 4.1 in https://arxiv.org/pdf/1802.01056.pdf 
           coefs_=[3.1378,-3.9789,2.6788,-1.0401,0.2139,-0.0133]
        nor=noiseSdev_*np.random.normal(0,1,N)
        t[0]=0.0
        u[0]=0.0
        nCoefs=len(coefs_)
        for i in range(1,nCoefs):
            t[i]=t[i-1]+dt[0]
            u[i]=-0.1+.2*np.random.rand()
        for i in range(nCoefs,N):
            t[i]=t[i-1]+dt[i-1]
            tmp_=0.0
            for j in range(nCoefs):
                tmp_+=coefs_[j]*u[i-(j+1)]
            u[i]=tmp_+nor[i]  
        return t,u    
    
    def gen(self):
##def syntheticDataGen(n,opts)
       n=self.n
       #(1) Assignments
       nExtra=int(0.7*n)   #extra number of samples for the initial burn-in period
       N=n+nExtra          #total number of samples to generate
       dt_=self.opts['dt']       #time step size
       #check if dt is fixed or varaible for the samples
       if (isinstance(dt_,int) or isinstance(dt_,float)): #fixed dt
          dt=dt_*np.ones(N-1) 
       elif (isinstance(dt_,np.ndarray)): #variable dt
           #check the length of the dt array
           ndt_=len(dt_)
           if (n-len(dt_)) != 1:
              raise ValueError('Length of the dt array should be one less than the number of samples.')
           #extend the dt size by adding nExtra dt[0] at the beginning of the array        
           dt=np.zeros(N)
           for i in range(N):
               if i<nExtra:
                  dt[i]=dt_[0]
               else:
                  dt[i]=dt_[i-(nExtra+1)] 

       type_=self.opts['type']           
       if 'noiseSdev' in self.opts.keys():
           noiseSdev_=self.opts['noiseSdev']
       else:
           noiseSdev_=0.1 #default noise
       write_ts.pw('... Generating n=%d synthetic time-series samples.' %n)
       write_ts.pw('    using method %s' %type_)
       write_ts.pw('    with noise sdev = %g' %noiseSdev_)

       #(2) Generating Samples
       if type_=='synData1':   
          t_,f_=self.synData1(N,noiseSdev_,dt)
       elif type_=='synData2':
          t_,f_=self.synData2(N,noiseSdev_,dt)
       elif type_=='synData3': 
          t_,f_=self.synData3(N,noiseSdev_,dt,self.opts)
       else:
          print("Available 'type's: ", ['synData1','synData2','synData3']) 
          raise ValueError("Invalid 'type'.")
       #(3) Remove the initial burn-in samples in the series
       t=t_[(nExtra-1):(n+nExtra-1)]-t_[nExtra-1]
       f=f_[(nExtra-1):(n+nExtra-1)]
       n=len(f)
       #(4) Write the synthetic samples in file
       if 'writeInFile' in self.opts.keys():
          if self.opts['writeInFile']:
             fSynData=open(self.opts['writeFile'],'w')
             fSynData.write('# Time-series samples used as input in the code\n')
             fSynData.write('# no. of smaples= %d \n' %n)
             for i in range(n):
                 fSynData.write("%E \t %E \n" %(t[i],f[i]))
       return t,f          
#
#
###############
# Ext Funcs
###############
def syntheticDataGen_test():
    """
       Test syntheticDataGen()
       Generate synthetic data for time-series using fixed and varying time-steps
    """
    #---- SETTINGS
    n=1000          #no of samples
    noiseSdev_=0.15  #noise standard deviation
    #-----------------------
    #time-series with fixed dt
    synOpts={'type':'synData3',
             'dt':1.0,
             'noiseSdev':noiseSdev_}
    t1,f1=syntheticData(n,synOpts).gen()
    #time-series with random dt
    synOpts={'type':'synData3',
             'dt':(np.random.rand(n-1)+0.1)*np.ones(n-1),
             'noiseSdev':noiseSdev_}
    t2,f2=syntheticData(n,synOpts).gen()
    #plot
    plt.figure(figsize=(10,6))
    plt.subplot(211)
    plt.plot(t1,f1,'-b',label=r'Fixed $\Delta t$')
    plt.subplot(212)
    plt.plot(t2,f2,'-r',label=r'Variable $\Delta t$')
    plt.legend(loc='best')
    plt.show()
#
#
