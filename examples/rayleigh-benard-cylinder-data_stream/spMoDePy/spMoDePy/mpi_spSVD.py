import sys
import os
import copy
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6

from mpi4py import MPI
import numpy as np
import scipy.optimize
NoneType = type(None)

class spSVD_c():
    def __init__(self, logger):
        
        self.log = logger
        self.ifget_all_modes = False
        logger.write("warning","ifget_all_modes is hard coded to False. This parameter applies to lcl updates. It controls if one gets all modes in the global rotation, despite keeping less modes locally. I do not see a use for this in production runs. Thus it is set to false. If needed, activate in mpi_spSVD.py module")
        #logger.write("debug","This parameter applies to lcl updates. It controls if one gets all modes in the global rotation, despite keeping lees modes locally. I do not see a use for this in production runs. Thus it is set to false. If needed, activate in mpi_spSVD.py module")

    def lclSVD_to_gblSVD(self,Uii,Dii,Vtii,k_set,comm):


        # Get information from the communicator
        rank = comm.Get_rank()
        size = comm.Get_size()

        n=Uii.shape[0]
        Yii=np.diag(Dii)@Vtii # This is size KxM
        k=Yii.shape[0]
        k_set=k_set # The number of modes that want to be kept. This is different than k because in local updates, each rank can have different k that will provide with ALL global modes.
        m=Yii.shape[1]

        #print("rank "+repr(rank)+" has n="+repr(n)+" m="+repr(m)+" k="+repr(k))

        #Yii could have different k accross ranks, so pad them to standard size
        Yi=np.zeros((m,m),dtype=Yii.dtype)
        if self.ifget_all_modes==True:
            Ui=np.zeros((n,m),dtype=Uii.dtype)  #Pad the modes to be able to obtain all global rotations
        else:
            Ui=np.zeros((n,k),dtype=Uii.dtype)  #Do not pad the modes, since it requires much memory
        Ui[:,:k]=Uii[:,:]
        Yi[:k,:]=Yii[:,:]

        #Gather Yi into Y in rank 0
        #prepare the buffer for recieving
        Y = None
        if rank == 0:
            #Generate the buffer to gather in rank 0
            Y = np.empty((m*size,m),dtype=Yii.dtype)
        comm.Gather(Yi, Y, root=0)

        if rank == 0:
            #If tank is zero, calculate the svd of the combined eigen matrix
            #Perform the svd of the combined eigen matrix
            #tic_in = time.perf_counter()
            Uy,Dy,Vty=np.linalg.svd(Y, full_matrices=False)
            #toc_in = time.perf_counter()
            #print(f"Time for SVD of Y in rank {rank}: {toc_in - tic_in:0.4f} seconds")
        else:
            #If the rank is not zero, simply create a buffer to recieve the Uy Dy and Vty
            Uy  = np.empty((m*size,m),dtype=Uii.dtype)
            Dy  = np.empty((m),dtype=Dii.dtype)
            Vty = np.empty((m,m),dtype=Vtii.dtype)
        comm.Bcast(Uy, root=0)
        comm.Bcast(Dy, root=0)
        comm.Bcast(Vty, root=0)
    
        #Perform the rotation of the local bases to obtain the global one
        if self.ifget_all_modes==True:
            Ui_global=Ui@Uy[rank*m:(rank+1)*m,:]
        else:
            Ui_global=Ui@Uy[rank*m:rank*m+k,:k_set]

        string = 'Ui is of shape[%d,%d]' % (Ui.shape[0], Ui.shape[1])
        self.log.write("debug",string)

        string = 'Ui_global is of shape[%d,%d]' % (Ui_global.shape[0], Ui_global.shape[1])
        self.log.write("debug",string)
         
        if self.ifget_all_modes!=True:
            Ui_global = np.ascontiguousarray(Ui_global[:,:k_set])
            Dy = np.ascontiguousarray(Dy[:k_set])
            Vty = np.ascontiguousarray(Vty[:k_set,:])

        return Ui_global,Dy,Vty


    def lclSVD_update_fromBatch(self,U_1t,D_1t,Vt_1t,Xi,k): 
        '''
        Xi: new data batch

        '''

        if type(U_1t)==NoneType:
            #Perform the distributed SVD and don't accumulate
            U_1t,D_1t,Vt_1t=np.linalg.svd(Xi,full_matrices=False)
        else:
            #Find the svd of the new snapshot
            U_tp1,D_tp1,Vt_tp1=np.linalg.svd(Xi,full_matrices=False)
            #2 contruct matrices to Do the updating
            V_tilde=scipy.linalg.block_diag(Vt_1t.conj().T,Vt_tp1.conj().T)
            W=np.append(U_1t@np.diag(D_1t),U_tp1@np.diag(D_tp1),axis=1)
            Uw,Dw,Vtw=np.linalg.svd(W,full_matrices=False)
            #3 Update
            U_1t=Uw
            D_1t=Dw
            Vt_1t=(V_tilde@Vtw.conj().T).conj().T

        #Truncate the matrices if needed.
        #Eliminate the lower energy mode, which should be the last ones
        if U_1t.shape[1]>k:
            U_1t=np.copy(U_1t[:,0:k])
            D_1t=np.copy(D_1t[0:k])
            Vt_1t=np.copy(Vt_1t[0:k,:])

        return U_1t,D_1t,Vt_1t


    def gblSVD(self,Xi,comm):

        # Get information from the communicator
        rank = comm.Get_rank()
        size = comm.Get_size()
        # Get some set up data
        m=Xi.shape[1]

        #Perfrom Svd in all ranks
        #tic_in = time.perf_counter()
        Ui,Di,Vti=np.linalg.svd(Xi, full_matrices=False)
        #toc_in = time.perf_counter()
        Yi=np.diag(Di)@Vti
        #print(f"Time for SVD of Xi in rank {rank}: {toc_in - tic_in:0.4f} seconds")

        #Gather Yi into Y in rank 0
        #prepare the buffer for recieving
        Y = None
        if rank == 0:
            #Generate the buffer to gather in rank 0
            Y = np.empty((m*size,m),dtype=Xi.dtype)
        comm.Gather(Yi, Y, root=0)

        if rank == 0:
            #If tank is zero, calculate the svd of the combined eigen matrix
            #Perform the svd of the combined eigen matrix
            #tic_in = time.perf_counter()
            Uy,Dy,Vty=np.linalg.svd(Y, full_matrices=False)
            #toc_in = time.perf_counter()
            #print(f"Time for SVD of Y in rank {rank}: {toc_in - tic_in:0.4f} seconds")
        else:
            #If the rank is not zero, simply create a buffer to recieve the Uy Dy and Vty
            Uy  = np.empty((m*size,m),dtype=Xi.dtype)
            Dy  = np.empty((m),dtype=Xi.dtype)
            Vty = np.empty((m,m),dtype=Xi.dtype)
        comm.Bcast(Uy, root=0)
        comm.Bcast(Dy, root=0)
        comm.Bcast(Vty, root=0)
        #Now matrix multiply each Ui by the corresponding entries in Uy
        U_local=Ui@Uy[rank*m:(rank+1)*m,:]

        return U_local, Dy, Vty


    def gblSVD_update_fromBatch(self,U_1t,D_1t,Vt_1t,Xi,k,comm): 
        '''

        '''

        if type(U_1t)==NoneType:
            self.log.write("info","Initialized the updating arrays")
            #Perform the distributed SVD and don't accumulate
            U_1t,D_1t,Vt_1t=self.gblSVD(Xi,comm)
        else:
            #Find the svd of the new snapshot
            U_tp1,D_tp1,Vt_tp1=self.gblSVD(Xi,comm)
            #2 contruct matrices to Do the updating
            V_tilde=scipy.linalg.block_diag(Vt_1t.conj().T,Vt_tp1.conj().T)
            W=np.append(U_1t@np.diag(D_1t),U_tp1@np.diag(D_tp1),axis=1)
            Uw,Dw,Vtw=self.gblSVD(W,comm)
            #3 Update
            U_1t=Uw
            D_1t=Dw
            Vt_1t=(V_tilde@Vtw.conj().T).conj().T

        #Truncate the matrices if needed.
        #Eliminate the lower energy mode, which should be the last ones
        if U_1t.shape[1]>=k:
            U_1t=np.copy(U_1t[:,0:k])
            D_1t=np.copy(D_1t[0:k])
            Vt_1t=np.copy(Vt_1t[0:k,:])

        return U_1t,D_1t,Vt_1t



