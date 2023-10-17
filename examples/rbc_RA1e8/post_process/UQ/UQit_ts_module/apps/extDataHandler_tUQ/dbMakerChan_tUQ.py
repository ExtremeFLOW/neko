################################################################################
#   Interface with Nek5000
#      to analyze the time-series of turbulent channel flow data 
#
#   * The time-series data dumped during the simulations of turbulent flows 
#     are read in.
#   * The read-in data are written in pickle format to speed up their future use
################################################################################
#  Saleh Rezaeiravesh, salehr@kth.se
#---------------------------------------------------------------------------------
#
import sys
import os
import numpy as np
import math as mt
import pickle
#
#////////////////////////////////////
def readNekRawData(dirName,fileName):
    """ 
       Read in the raw data of Nek5000 saved at each simulation time step.
       There are 'ny' points in the channel half height at each the time-series of
           velocity, ... are dumped on the course of flow simulation. 
       NOTE: The function reads the input file line-by-line to avoid memory issues. 
       Inputs:
          dirName : Path at which Nek5000 data are located
          fileName: Name of the Nek5000 data file
       Outputs:
          y: 1D numpy array of size ny holding wall-normal distance of a profile
          time: 1D numpy array of size nTime holding time associated to each sample
          uTau: 1D numpy array of size nTime holding <uTau> at each time step
          vel_db: list (of length ny) of 1D numpy arrays of length nTime holding samples of 
                'u', 'v', 'w' (velocity components)
                'uu','vv','ww','uv' (multiplication of the velocity components)
                  Ex. vel_db[1]['u']: u-TS at y[1] (NOTE: y[0]=0 specifies the wall)
                  Ex. vel_db[3]['v']: v-TS at y[3]                  
    """ 
    db={'time':[],'uTau':[]};  #initialize the database for time and uTau
    vel_db=[];   #list of dictionaries (databases) to keep velocity time series at \
                 #  each wall-normal location. 
                 #NOTE: each line is devoted to one wall-normal location
    iHit=False     
    ic=0                 
    with open(dirName+'/'+fileName) as infile:
        iLine=0
        for line in infile:
            ain_sep=line.split() 
            #Find ny: The num of points in the channel half height and associated wall-distance y
            if iLine==0:
               ny=int(ain_sep[len(ain_sep)-1])
               y=[]; #wall-normal coordinates of GLL pts
               for i in range(ny):
                   db_={'u':[],'v':[],'w':[],'uu':[],'vv':[],'ww':[],'uv':[]} 
                   vel_db.append(db_)              
            #Wall-normal coordinates of the GLL point over channel half-height    
            if iLine>0 and iLine<=ny:
               y.append(float(ain_sep[0]))
               if iLine==ny:
                  db.update({'y':y})
            #Update time and uTau list in db
            if ain_sep[0]=='time,':
               iHit= True
               if 'time' in db.keys():
                  time_list=db['time'] #read the current list
                  uTau_list=db['uTau'] 
                  time_list.append(float(ain_sep[2]))
                  uTau_list.append(float(ain_sep[3]))
               db.update({'time':time_list})  #update the dict
               db.update({'uTau':uTau_list})  
            #Read velocity at the current time and at all y locations            
            if iHit and ain_sep[0]!='time,':
               #read new velocity at y[ic]
               # and update the ic-th row of the lis_db          
               for iq_ in range(len(db_.keys())):
                   dbList_update(vel_db,ic,list(db_.keys())[iq_],float(ain_sep[iq_]))
               if ic==ny-1:
                  iHit=False
                  ic=0
               else:
                  ic+=1                                                      
            #update no of lines 
            iLine+=1
    #convert lists to numpy arrays        
    print('...... Converting lists to numpy arrays in readNekRawData()')
    y   =np.asarray(db['y'])
    t   =np.asarray(db['time'])    
    uTau=np.asarray(db['uTau'])
    keyList=['u','v','w','uu','vv','ww','uv']
    for i in range(len(y)):
        for key_ in keyList:
            vel_=np.asarray(vel_db[i][key_])
            vel_db[i].update({key_:vel_})
    return vel_db,t,y,uTau
#
#
#/////////////////////////////////////
def dbList_update(vel_db,ic,key,val):
    """
       Update the ic-th row of a list of dictionaries (vel_db) 
           with a new value 'val', corresponding to the key
    """
    db_=vel_db[ic]     #get the dict at row ic 
    f_list=db_[key]    #get the dict values corresponding to key
    f_list.append(val)  #update the list
    db_.update({key:f_list})   #update the dictionary
    vel_db[ic]=db_    #update row ic of the list of dictionaries    
    return vel_db
#
#
#///////////////////////////////////////
def createPickleDB(dataDir,rawDataName):
    """
       Create pickle database from raw time series of turbulent channel flow. 
       Inputs:
           dataDir: Path to the folder containing the data (both raw and processed)
           rawdataName: Name of the file which contains the time-series data
           #>> Guide for data:
           # t: list of length nTime containing the times at which data have been saved in NEK
           # uTau: list of length nTime containing uTau (avg on both walls) samples
           # y: list of length ny(=no of GLL pts normal to the wall across the channel half height)
           # vel: a list of dictionaries containing time-series (1D numpy arrays) of 
                    velocity components u,v,w at each of the ny GLL pts. 
                  Each row of vel belongs to one y:
           #      ex. vel[1]['u']: u-signal at y[1] (y[0]=0 specifies the wall)
           #      ex. vel[3]['v']: v-signal at y[3]
       Outputs:
           written databases in pickle format
    """
    print('... Creating a pickle database from raw time-series data')
    # Read in Nek5000 data
    vel,t,y,uTau=readNekRawData(dataDir+'/rawData',rawDataName)
    print('...... Data at %d time instants and %d wall-normal locations are read!' %(len(t),len(y)))
    # Create the ouput folder if does not exist
    dataDir_=dataDir+'/database/'   #output path
    if not os.path.exists(dataDir_):
       os.makedirs(dataDir_)
    # Save the read-in data as pickle
    print('...... Dumping data (using pickle) at ')
    print('       %s'%dataDir_)
    with open(dataDir_+'t','wb') as f:
        pickle.dump(t,f);
    with open(dataDir_+'y','wb') as f:
        pickle.dump(y,f);
    with open(dataDir_+'uTau','wb') as f:
        pickle.dump(uTau,f);
    with open(dataDir_+'vel','wb') as f:
        pickle.dump(vel,f);
#
#
#//////////////////////////
def readPickleDB(dirName):
    """ 
       Read in pickle databases created by createPickleDB()
    """
    print('... Reading pickle databases located at:')
    print('    %s' %dirName)
    with open(dirName+'/'+'t', 'rb') as f:
         t= pickle.load(f)
    print('... Size of time series = %d' %(t.size))
    with open(dirName+'/'+'uTau', 'rb') as f:
         uTau= pickle.load(f)
    with open(dirName+'/'+'vel', 'rb') as f:
         vel= pickle.load(f)
    with open(dirName+'/'+'y', 'rb') as f:
         y= pickle.load(f)
    return t,y,uTau,vel
#
#
