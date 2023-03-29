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
import time
import matplotlib.pyplot as plt
import pymech as pm
from pymech.neksuite import readnek
from pymech.neksuite import writenek
from tqdm import tqdm


class spSVD():
    def __init__(self,
            k: int=0,
            p: int=0,
            nsnap: int=0,
            qoi: int=0,
            path: str='',
            casename: str='',
            ifautoupdate: bool=False,
            min_ratio: float=0.0,
            ifget_all_modes: bool=False,
            ifgbl_update: bool=False):

        super().__init__()

        self.setk = k
        self.k = k
        self.p = p
        self.nsnap=nsnap
        self.qoi=qoi
        self.path=path
        self.casename=casename
        self.ifautoupdate=ifautoupdate
        self.min_ratio=min_ratio
        self.ifget_all_modes=ifget_all_modes
        self.ifgbl_update=ifgbl_update

    
    def get_case_info(self,comm):
        
        # Get information from the communicator
        rank = comm.Get_rank()
        size = comm.Get_size()
        path=self.path
        casename=self.casename

        if rank==0:
            buff  = np.zeros((10))
            # Read information about the case (number of elements, etc)
            filename=path+casename+'0.f'+str(0+1).zfill(5)
            data=readnek(filename)
            # simulation parameters
            buff[0]=data.nel #nel
            buff[1]=data.lr1[0] #lx1
            buff[2]=data.lr1[1] #ly1
            buff[3]=data.lr1[2] #lz1
            buff[4]=np.prod(data.lr1) #nxyz
            buff[5]= data.nel*data.lr1[0]*data.lr1[1]*data.lr1[2] #nxyze
            if data.lr1[2]==1:
                buff[6]=2
            else:
                buff[6]=3
        else:
            buff  = np.empty((10))

        comm.Bcast(buff, root=0)

        self.nel=int(buff[0])
        self.lx1=int(buff[1])
        self.ly1=int(buff[2])
        self.lz1=int(buff[3])
        self.nxyz=int(buff[4])
        self.nxyze=int(buff[5])
        self.dim=int(buff[6])

        return


    def lclSVD_to_gblSVD(self,Uii,Dii,Vtii,comm):

        # Get information from the communicator
        rank = comm.Get_rank()
        size = comm.Get_size()

        n=Uii.shape[0]
        Yii=np.diag(Dii)@Vtii # This is size KxM
        k=Yii.shape[0]
        k_set=self.setk
        m=Yii.shape[1]

        #print("rank "+repr(rank)+" has n="+repr(n)+" m="+repr(m)+" k="+repr(k))

        #Yii could have different k accross ranks, so pad them to standard size
        Yi=np.zeros((m,m))
        if self.ifget_all_modes==True:
            Ui=np.zeros((n,m))  #Pad the modes to be able to obtain all global rotations
        else:
            Ui=np.zeros((n,k))  #Do not pad the modes, since it requires much memory
        Ui[:,:k]=Uii[:,:]
        Yi[:k,:]=Yii[:,:]

        #Gather Yi into Y in rank 0
        #prepare the buffer for recieving
        Y = None
        if rank == 0:
            #Generate the buffer to gather in rank 0
            Y = np.empty((m*size,m))
        comm.Gather(Yi, Y, root=0)

        if rank == 0:
            #If tank is zero, calculate the svd of the combined eigen matrix
            #Perform the svd of the combined eigen matrix
            tic_in = time.perf_counter()
            Uy,Dy,Vty=np.linalg.svd(Y, full_matrices=False)
            toc_in = time.perf_counter()
            #print(f"Time for SVD of Y in rank {rank}: {toc_in - tic_in:0.4f} seconds")
        else:
            #If the rank is not zero, simply create a buffer to recieve the Uy Dy and Vty
            Uy  = np.empty((m*size,m))
            Dy  = np.empty((m))
            Vty = np.empty((m,m))
        comm.Bcast(Uy, root=0)
        comm.Bcast(Dy, root=0)
        comm.Bcast(Vty, root=0)
    
        #Perform the rotation of the local bases to obtain the global one
        if self.ifget_all_modes==True:
            Ui_global=Ui@Uy[rank*m:(rank+1)*m,:]
        else:
            Ui_global=Ui@Uy[rank*m:rank*m+k,:k_set]

        print('shape of Ui:')
        print(Ui.shape)
        print('shape of Ui_global:')
        print(Ui_global.shape)
    
        if self.ifget_all_modes!=True:
            Ui_global = np.ascontiguousarray(Ui_global[:,:k_set])
            Dy = np.ascontiguousarray(Dy[:k_set])
            Vty = np.ascontiguousarray(Vty[:k_set,:])

        return Ui_global,Dy,Vty


    def lclSVD_update_fromBatch(self,U_1t,D_1t,Vt_1t,Xi,j): 
        '''
        Xi: new data batch
        j:  Current snapshot

        '''
        k=self.k
        p=self.p

        if j==(p-1):
            #Perform the distributed SVD and don't accumulate
            U_1t,D_1t,Vt_1t=np.linalg.svd(Xi,full_matrices=False)
        else:
            #Find the svd of the new snapshot
            U_tp1,D_tp1,Vt_tp1=np.linalg.svd(Xi,full_matrices=False)
            #2 contruct matrices to Do the updating
            V_tilde=scipy.linalg.block_diag(Vt_1t.T,Vt_tp1.T)
            W=np.append(U_1t@np.diag(D_1t),U_tp1@np.diag(D_tp1),axis=1)
            Uw,Dw,Vtw=np.linalg.svd(W,full_matrices=False)
            #3 Update
            U_1t=Uw
            D_1t=Dw
            Vt_1t=(V_tilde@Vtw.T).T

            #Truncate the matrices if needed.
            #Eliminate the lower energy mode, which should be the last ones
            if (j+1)>=k:
                U_1t=np.copy(U_1t[:,0:k])
                D_1t=np.copy(D_1t[0:k])
                Vt_1t=np.copy(Vt_1t[0:k,:])

        return U_1t,D_1t,Vt_1t


    def gblSVD(self,Xi,comm):

        # Get information from the communicator
        rank = comm.Get_rank()
        size = comm.Get_size()
        # Get some set up data
        n=self.nxyze
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
            Y = np.empty((m*size,m))
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
            Uy  = np.empty((m*size,m))
            Dy  = np.empty((m))
            Vty = np.empty((m,m))
        comm.Bcast(Uy, root=0)
        comm.Bcast(Dy, root=0)
        comm.Bcast(Vty, root=0)
        #Now matrix multiply each Ui by the corresponding entries in Uy
        U_local=Ui@Uy[rank*m:(rank+1)*m,:]

        return U_local, Dy, Vty


    def gblSVD_update_fromBatch(self,U_1t,D_1t,Vt_1t,Xi,j,comm): 
        '''
        Xi: new data batch
        j:  Current snapshot

        '''
        k=self.k
        kset=self.setk
        p=self.p
        nsnap=self.nsnap
        ifgbl_update=self.ifgbl_update

        if j==(p-1):
            #Perform the distributed SVD and don't accumulate
            U_1t,D_1t,Vt_1t=self.gblSVD(Xi,comm)
            
            if p == nsnap and ifgbl_update == True:
                U_1t=np.copy(U_1t[:,0:k])
                D_1t=np.copy(D_1t[0:k])
                Vt_1t=np.copy(Vt_1t[0:k,:])

        else:
            #Find the svd of the new snapshot
            U_tp1,D_tp1,Vt_tp1=self.gblSVD(Xi,comm)
            #2 contruct matrices to Do the updating
            V_tilde=scipy.linalg.block_diag(Vt_1t.T,Vt_tp1.T)
            W=np.append(U_1t@np.diag(D_1t),U_tp1@np.diag(D_tp1),axis=1)
            Uw,Dw,Vtw=self.gblSVD(W,comm)
            #3 Update
            U_1t=Uw
            D_1t=Dw
            Vt_1t=(V_tilde@Vtw.T).T

            #Truncate the matrices if needed.
            #Eliminate the lower energy mode, which should be the last ones
            if (j+1)>=k:
                U_1t=np.copy(U_1t[:,0:k])
                D_1t=np.copy(D_1t[0:k])
                Vt_1t=np.copy(Vt_1t[0:k,:])

        return U_1t,D_1t,Vt_1t



##################### From here on can be in another file or class #######################

    
    def gatherModes_atRank0(self,U_1t, n, m,comm):
        
        # Get information from the communicator
        rank = comm.Get_rank()
        size = comm.Get_size()
        
        U = None #prepare the buffer for recieving
        if rank == 0:
            #Generate the buffer to gather in rank 0
            U = np.empty((n,m))
        comm.Gather(U_1t, U, root=0)
        return U

    def gatherModesandMass_atRank0(self,U_1t,BM1, n, m,comm):
        
        # Get information from the communicator
        rank = comm.Get_rank()
        size = comm.Get_size()
        
        U = None #prepare the buffer for recieving
        bm1sqrt = None #prepare the buffer for recieving
        if rank == 0:
           #Generate the buffer to gather in rank 0
            U = np.empty((n,m))
            bm1sqrt = np.empty((n,1))
        comm.Gather(U_1t, U, root=0)
        comm.Gather(BM1, bm1sqrt, root=0)
        return U, bm1sqrt


    
    def write_modes(self,U, first_mode, number_of_modes,prefix):

        path=self.path
        casename=self.casename
        filename=path+casename+'0.f'+str(0+1).zfill(5)
        data=readnek(filename)

        # Copy the data to overwrite
        dataout=copy.deepcopy(data)

        # Retrieve the information
        path=self.path
        casename=self.casename
        qoi=self.qoi
        nxyz=self.nxyz
        nel=self.nel

        filename=path+'PODmod'+'_'+prefix+'_'+casename+'0.f'+str(0+1).zfill(5)
        print('Writing POD modes in: '+filename)

        pbar= tqdm(total=number_of_modes)
        for j in range(first_mode,first_mode+number_of_modes):
            for e in range(0,nel):
                if qoi != 3:
                    #dataout.elem[e].vel[qoi,:,:,:]=np.reshape(U[e*nxyz:e*nxyz+nxyz,j],(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
                    dataout.elem[e].vel[qoi,:,:,:]=np.reshape(U[e*nxyz:e*nxyz+nxyz,j],(data.lr1[2],data.lr1[1],data.lr1[0]))
                elif qoi ==3:
                    dataout.elem[e].pres[0,:,:,:]=np.reshape(U[e*nxyz:e*nxyz+nxyz,j],(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')

            filename=path+'PODmod'+'_'+prefix+'_'+casename+'0.f'+str(j+1).zfill(5)
            writenek(filename,dataout)
            pbar.update(1)
        pbar.close()

    def write_rec(self,U, first_mode, number_of_modes,prefix):

        path=self.path
        casename=self.casename
        filename=path+casename+'0.f'+str(0+1).zfill(5)
        data=readnek(filename)

        # Copy the data to overwrite
        dataout=copy.deepcopy(data)

        # Retrieve the information
        path=self.path
        casename=self.casename
        qoi=self.qoi
        nxyz=self.nxyz
        nel=self.nel

        filename=path+'PODrec'+'_'+prefix+'_'+casename+'0.f'+str(0+1).zfill(5)
        print('Writing POD rec in: '+filename)


        pbar= tqdm(total=number_of_modes)
        for j in range(first_mode,first_mode+number_of_modes):
            for e in range(0,nel):
                if qoi != 3:
                    dataout.elem[e].vel[qoi,:,:,:]=np.reshape(U[e*nxyz:e*nxyz+nxyz,j],(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
                elif qoi ==3:
                    dataout.elem[e].pres[0,:,:,:]=np.reshape(U[e*nxyz:e*nxyz+nxyz,j],(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
            filename=path+'PODrec'+'_'+prefix+'_'+casename+'0.f'+str(j+1).zfill(5)
            writenek(filename,dataout)
            pbar.update(1)
        pbar.close()

    def write_timecoeff(self,D,Vt,prefix):

        # Retrieve the information
        path=self.path
        casename=self.casename

        filename1=path+'PODsig'+'_'+prefix+'_'+casename
        filename2=path+'PODvtr'+'_'+prefix+'_'+casename
    
        np.save(filename1,D) 
        np.save(filename2,Vt) 


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
        return X

    def get_perp_ratio(self,U,v):
        m=U.T@v
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
        mi=U.T@v

        # Use MPI_SUM to get the global one by aggregating local
        m = np.zeros_like(mi)
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
        nrm = np.zeros_like(nrmi)
        comm.Allreduce(nrmi, nrm, op=MPI.SUM)
 
        # Get the actual norm
        nrm_o=np.sqrt(nrm[0])
        nrm=np.sqrt(nrm[1])

        return nrm_o/nrm

    def read_mass_data(self,bm1sqrt,isnap):
        path=self.path
        casename=self.casename
        qoi=self.qoi
        nxyz=self.nxyz    
        nel=self.nel
    
        filename=path+'bm1'+casename+'0.f'+str(isnap+1).zfill(5)
        bm1data=readnek(filename)
        for e in range(0,nel):
            #Rearange 
            if qoi != 3:
                bm1e=bm1data.elem[e].vel[qoi,:,:,:].reshape((nxyz,1),order='F')
            elif qoi==3:
                bm1e=bm1data.elem[e].vel[0,:,:,:].reshape((nxyz,1),order='F')
            #Copy into a bigger vector to scatter
            bm1sqrt[e*nxyz:e*nxyz+nxyz,0]=np.copy(np.sqrt(bm1e[:,0]))
        return bm1sqrt

    def read_vel_data(self,v,isnap):
        path=self.path
        casename=self.casename
        qoi=self.qoi
        nxyz=self.nxyz
        nel=self.nel

        filename=path+casename+'0.f'+str(isnap+1).zfill(5)
        data=readnek(filename)
        for e in range(0,nel):
            if qoi != 3:
                #Rearange 
                ve=data.elem[e].vel[qoi,:,:,:].reshape((nxyz,1),order='F')
                #Copy into a bigger vector to scatter
                v[e*nxyz:e*nxyz+nxyz,0]=np.copy((ve[:,0]))
            elif qoi ==3:
                #Rearange 
                ve=data.elem[e].pres[0,:,:,:].reshape((nxyz,1),order='F')
                #Copy into a bigger vector to scatter
                v[e*nxyz:e*nxyz+nxyz,0]=np.copy((ve[:,0]))



        return v
