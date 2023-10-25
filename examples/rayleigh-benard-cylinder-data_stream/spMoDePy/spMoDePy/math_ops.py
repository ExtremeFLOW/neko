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

class math_ops_c():
    def __init__(self):
        a=1

    def scale_data(self, X,bm,rows,columns,scale):
        
        #Scale the data with the mass matrix
        if scale=='mult':
            for i in range(0,columns):
                for j in range(0,rows):
                    X[j,i]= X[j,i]*bm[j,0]
        elif scale=='div':
            for i in range(0,columns):
                for j in range(0,rows):
                    X[j,i]= X[j,i]/bm[j,0]

    
    def get_perp_ratio(self,U,v):
        m=U.conj().T@v
        v_parallel=U@m
        v_orthogonal=v-v_parallel
        nrm_o=np.linalg.norm(v_orthogonal)
        nrm=np.linalg.norm(v)    
        return nrm_o/nrm

    def mpi_get_perp_ratio(self,U,v,comm):
        # Get information from the communicator
        rank = comm.Get_rank()
        size = comm.Get_size()
        
        # Get the local partial step
        mi=U.conj().T@v

        # Use MPI_SUM to get the global one by aggregating local
        m = np.zeros_like(mi,dtype=mi.dtype)
        comm.Allreduce(mi, m, op=MPI.SUM)
        
        # Get local parallel component
        v_parallel=U@m
        # Get local orthogonal component
        v_orthogonal=v-v_parallel

        #  Do a local sum of squares
        nrmi=np.zeros((2))
        nrmi[0]=np.sum(v_orthogonal**2)
        nrmi[1]=np.sum(v**2)    

        # Then do an all reduce
        nrm = np.zeros_like(nrmi,dtype=nrmi.dtype)
        comm.Allreduce(nrmi, nrm, op=MPI.SUM)
 
        # Get the actual norm
        nrm_o=np.sqrt(nrm[0])
        nrm=np.sqrt(nrm[1])

        return nrm_o/nrm
    
    def gatherModesandMass_atRank0(self,U_1t,BM1, n,comm):
        
        # Get information from the communicator
        rank = comm.Get_rank()
        size = comm.Get_size()
        m = U_1t.shape[1]

        U = None #prepare the buffer for recieving
        bm1sqrt = None #prepare the buffer for recieving
        if rank == 0:
           #Generate the buffer to gather in rank 0
            U = np.empty((n,m),dtype=U_1t.dtype)
            bm1sqrt = np.empty((n,1))
        comm.Gather(U_1t, U, root=0)
        comm.Gather(BM1, bm1sqrt, root=0)
        return U, bm1sqrt
