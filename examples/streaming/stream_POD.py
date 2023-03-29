import sys
import os
import copy
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6

from mpi4py import MPI #equivalent to the use of MPI_init() in C
import numpy as np
import scipy.optimize
import time
import matplotlib.pyplot as plt
import pymech as pm
from pymech.neksuite import readnek
from pymech.neksuite import writenek
from tqdm import tqdm
import mpi_spSVD as spSVD
import pandas as pd

#####################################
# Case parameters
## Number of snapshots
nsnap=500
## Batch size
p=10 # Number of snapshots in a batch
## Pod of which quantity
qoi=0 #0=vx, 1=vy, 2=vz

# Data path
casename='mixlay'
path='/scratch/adperez/TOOLS_iostream_POD/data/'+casename+'/'

# number of modes
k=500
maxk = 1000
# Options if autodecide number of modes
## Autoupdate?
ifautoupdate=False
## If autoupdate has been done, gather ALL modes?
ifget_all_modes=False
## Minimun percentage of new snapshots that need to be parallel
min_ratio=0.01

# Update gbl or lcl?
ifgbl_update=True
num_recs = 10

#### Automatic, do not touch #####
if ifgbl_update==True: 
    outprefix='gbl'
else:
    outprefix='lcl'
#################################
# Overwrite the hyperparameters if command line arguments were given
if len(sys.argv) > 1:
    p=int(sys.argv[1])
    k=int(sys.argv[2])
    min_ratio=float(sys.argv[3])
    str_mr = str(min_ratio).replace('.','')
    ifgbl_update=bool(int(sys.argv[4]))
    ifautoupdate=bool(int(sys.argv[5]))

    if ifgbl_update==True: 
        if ifautoupdate==True:
            outprefix='auto_gbl_p'+repr(p)+'_k'+repr(k)+'_mr'+str_mr
        else:
            outprefix='fixd_gbl_p'+repr(p)+'_k'+repr(k)+'_mr'+str_mr
    else:
        if ifautoupdate==True:
            outprefix='auto_lcl_p'+repr(p)+'_k'+repr(k)+'_mr'+str_mr
        else:
            outprefix='fixd_lcl_p'+repr(p)+'_k'+repr(k)+'_mr'+str_mr

################################

if __name__ == "__main__":
    
    # MPI - MPMD
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    if rank ==0: print(outprefix)
    
    #define the class
    spPOD = spSVD.spSVD(k=k, p=p, nsnap=nsnap, qoi= qoi, path =path, casename=casename, ifautoupdate=ifautoupdate, min_ratio=min_ratio, ifget_all_modes=ifget_all_modes,ifgbl_update=ifgbl_update)
    #retreive case information
    spPOD.get_case_info(comm)
    nel = spPOD.nel
    lx1 = spPOD.lx1
    ly1 = spPOD.ly1
    lz1 = spPOD.lz1
    nxyz= spPOD.nxyz
    nxyze= spPOD.nxyze

    #Create buffer for reading data and redistributing it
    if (rank==0): 
        v=np.empty((nxyze,1))
        bm1sqrt=np.empty((nxyze,1))
    else:
        v=None
        bm1sqrt=None

    #create local recieve buffer in each MPI rank
    Xi = np.empty((int(nxyze/size),1))
    buff=np.empty((int(nxyze/size),p))
    bm1sqrti = np.empty((int(nxyze/size),1))

    # Create dummies for the streaming POD
    U_1t = None
    D_1t = None
    Vt_1t = None

    # Change k to a low value if autoupdate is required
    if ifautoupdate==True: spPOD.k=p

    running_ra=[]

    #This is the streaming process
    pbar= tqdm(total=nsnap)
    for isnap in range(0,nsnap):

        #read mass matrix
        if (isnap==0):
            # Read mass matrix
            if (isnap==0):
                if (rank==0):
                    bm1sqrt=spPOD.read_mass_data(bm1sqrt,isnap)
		        #Scatter the data red from rank0
                comm.Scatter(bm1sqrt, bm1sqrti, root=0)

        # Read vel matrix
        if (rank==0):
            v = spPOD.read_vel_data(v,isnap)
        #Scatter the data red from rank0
        comm.Scatter(v, Xi, root=0)
		
        #Scale the data with the mass matrix
        Xi=spPOD.scale_data(Xi,bm1sqrti,int(nxyze/size),1,'mult')

        # Update the number of snapshots processed up to now
        m=isnap+1

        # fill in the buffer
        buff[:,np.mod(isnap,spPOD.p)] = Xi[:,0]
        cbf = np.mod(isnap,spPOD.p)


        # Calculate the residual and check if basis needs to be expanded 
        if isnap>=spPOD.p: 
            if ifautoupdate==True:
                if ifgbl_update == False:
                    ra=spPOD.get_perp_ratio(U_1t,Xi.reshape((-1,1)))
                    running_ra.append(ra)
                else:
                    ra=spPOD.mpi_get_perp_ratio(U_1t,Xi.reshape((-1,1)),comm)
                    running_ra.append(ra)
            else:
                ra=0
                running_ra.append(ra)
            if spPOD.ifautoupdate==True and ra>=spPOD.min_ratio and spPOD.k<maxk: spPOD.k+=1

        # Update
        if np.mod(isnap+1,spPOD.p) == 0 or (isnap+1)==nsnap:
            # Print the local state
            print("rank "+repr(rank)+" has k="+repr(spPOD.k))
            # Update the svd with each new snapshot 
            if ifgbl_update==True:
                U_1t,D_1t,Vt_1t = spPOD.gblSVD_update_fromBatch(U_1t,D_1t,Vt_1t,buff[:,:(cbf+1)],isnap,comm)
            else:
                U_1t,D_1t,Vt_1t = spPOD.lclSVD_update_fromBatch(U_1t,D_1t,Vt_1t,buff[:,:(cbf+1)],isnap)

        pbar.update(1)
    pbar.close() 

    #Do this in case you want to keep more data than you actually collected
    if spPOD.setk>=m: spPOD.setk=m

    ##Scale the modes back before gathering them
    U_1t=spPOD.scale_data(U_1t,bm1sqrti,int(nxyze/size),spPOD.k,'div')

    ## If local update, reconstruct for later comparison
    if ifgbl_update == False: Si = U_1t@np.diag(D_1t)@Vt_1t[:,0:num_recs]

    # Obtain a global svd from a local one only if the svd has been updated locally
    if ifgbl_update == False: U_1t,D_1t,Vt_1t = spPOD.lclSVD_to_gblSVD(U_1t,D_1t,Vt_1t,comm)
    
    if ifgbl_update == False: 
        ortho = comm.gather(running_ra,root=0)
    else:
        ortho = running_ra

    print(U_1t.shape)

    #Once the streaming is done. Gather the modes at one rank
    if ifget_all_modes==True: kk=m
    if ifget_all_modes==False:
        if ifgbl_update == True:
            if spPOD.setk <= spPOD.k: # This line here is to protect the code from truncating when self updating, other wise they are the same
                kk=spPOD.setk
            else:
                kk=spPOD.k
        else:
            kk=spPOD.setk

        U_1t=np.copy(U_1t[:,0:kk])
        D_1t=np.copy(D_1t[0:kk])
        Vt_1t=np.copy(Vt_1t[0:kk,:])

    U, bm1sqrt2 = spPOD.gatherModesandMass_atRank0(U_1t,bm1sqrti,spPOD.nxyze,kk,comm)
 
    if ifgbl_update == False: S, _ = spPOD.gatherModesandMass_atRank0(Si,bm1sqrti,spPOD.nxyze,num_recs,comm)

    #Write the modes
    if rank==0:
        
        # Save the running orthogonality and others
        np.save(spPOD.path+'PODrra'+'_'+outprefix+'_'+spPOD.casename,np.array(ortho))
    
        print(U.shape)
        print(D_1t.shape)
        print(Vt_1t.shape)

        ## Do this because when globaly autoupdating, you are constrained
        #if ifgbl_update== True and ifautoupdate== True:
        #    kkk=spPOD.k
        #else:
        #    kkk=spPOD.setk

        spPOD.write_modes(U,0,kk,outprefix)
        spPOD.write_timecoeff(D_1t,Vt_1t,outprefix)
        
        if ifgbl_update == False:
            spPOD.write_rec(S,0,num_recs,outprefix)
