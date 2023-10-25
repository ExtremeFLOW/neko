import sys
import os
import copy

import numpy as np
from pymech.neksuite import readnek
from pymech.neksuite import writenek
from tqdm import tqdm
import h5py
NoneType = type(None)
        

class io_controller_c():
    def __init__(self,
                params_c,
                case_c,
                logger,
                comm):

        self.params = params_c
        self.case = case_c
        self.comm = comm
        self.log  = logger
        global_alloc(self,params_c,case_c,comm)

    # first step
    def first_step(self,params,case,iterator,md,comm):

        if params.ifrestart:
            io_read_restart(self,params,case,iterator,md,comm)

        io_allocate(self,params,case,md,comm)
        io_read(self,params,case,iterator,comm)

    def step(self,params,case,iterator,comm):
        
        io_read(self,params,case,iterator,comm)

    def close(self):
        self.log.write("warning", "Stream of data is over")

    def write(self,params,case,iterator,md,comm):
    
        rank     = comm.Get_rank()
        size     = comm.Get_size()

        if rank == 0:

            # Write the modes
            write_modes(self,params, case,md,params.update_type)

            # Write the coefficients
            write_timecoeff(self,params,md,params.update_type)

            # Write restart file
            write_restart_file(self,params,iterator,md)

# Support functions

def global_alloc(self,params,case,comm):

    # Get information from the communicator
    rank     = comm.Get_rank()
    size     = comm.Get_size()

    nxyzef = case.nxyzef
    nel = case.nel


    if (rank==0): 
        self.v=np.empty((nxyzef,1))
        self.bm1sqrt=np.empty((nxyzef,1))
        self.lglel=np.empty((nel,1),dtype=np.int32)
    else:
        self.v=None
        self.bm1sqrt=None
        self.lglel=None


def io_allocate(self,params,case,md,comm):
       
    #Associate
    rank = comm.Get_rank()
    size = comm.Get_size()
    nel = case.nel
    nxyzef= case.nxyzef
    p = params.batch_size
    
    #create local recieve buffer in each MPI rank
    self.Xi = np.zeros((int(nxyzef/size),1))
    md.buffer_allocator(params,case,int(nxyzef/size))
    self.bm1sqrti = np.zeros((int(nxyzef/size),1))
    self.lgleli = np.zeros((int(nel/size),1), dtype=np.int32)
        
    case.mynxyzef = int(nxyzef/size)

def io_read(self,params, case, iterator, comm):
        
    rank = comm.Get_rank()
    size = comm.Get_size()
        
    #read mass matrix
    if (iterator.first_iteration):
        if (rank==0):
            read_mass_data(self,params,case,0)
	#Scatter the data red from rank0
        comm.Scatter(self.bm1sqrt, self.bm1sqrti, root=0)

    # Read vel matrix
    if (rank==0):
        if params.decomposition == "POD": isnap = iterator.iteration
        if params.decomposition == "SPOD": isnap = iterator.iteration
        if params.decomposition == "DMD": isnap = iterator.iteration
        isnap = iterator.iteration
        read_vel_data(self,params,case,isnap)
    #Scatter the data red from rank0
    comm.Scatter(self.v, self.Xi, root=0)

    # Read global element number, which is just looping when reading from file
    if (rank==0):
        for e in range(0, case.nel):
            self.lglel[e] = e + 1
    #Scatter the data red from rank0
    comm.Scatter(self.lglel, self.lgleli, root=0)

    
def read_mass_data(self,params,case,isnap):
    path=params.dataPath
    casename=params.casename
    qoi=params.fields
    num_fields=params.num_fields
    nxyz=case.nxyz    
    nel=case.nel
    nxyze = case.nxyze
    
    filename=path+'bm1'+casename+'0.f'+str(isnap+1).zfill(5)
    bm1data=readnek(filename)
    for e in range(0,nel):
        for field in range(0,num_fields):
            qo = qoi[field]
            if qo != 3:
                bm1e=bm1data.elem[e].vel[qo,:,:,:].reshape((nxyz,1),order='F')
            elif qo==3:
                bm1e=bm1data.elem[e].vel[0,:,:,:].reshape((nxyz,1),order='F')                
            #Copy into a bigger vector to scatter
            self.bm1sqrt[(e*nxyz+field*nxyze):(e*nxyz+field*nxyze)+nxyz,0]=np.copy(np.sqrt(bm1e[:,0]))


def read_vel_data(self,params,case,isnap):
    path=params.dataPath
    casename=params.casename
    qoi=params.fields
    num_fields=params.num_fields
    nxyz=case.nxyz
    nel =case.nel
    nxyze = case.nxyze

    filename=path+casename+'0.f'+str(isnap+1).zfill(5)
    self.log.write("info", "Loading velocity data from: "+filename)
    data=readnek(filename)
    for e in range(0,nel):
        for field in range(0,num_fields):
            qo = qoi[field]
            if qo != 3:
                ve=data.elem[e].vel[qo,:,:,:].reshape((nxyz,1),order='F')
            elif qo==3:
                ve=data.elem[e].pres[0,:,:,:].reshape((nxyz,1),order='F')
            #Copy into a bigger vector to scatter
            self.v[(e*nxyz+field*nxyze):(e*nxyz+field*nxyze)+nxyz,0]=np.copy((ve[:,0]))
    

def write_modes(self,params, case,md,prefix):

    path=params.dataPath
    casename=params.casename
    filename=path+casename+'0.f'+str(0+1).zfill(5)
    data=readnek(filename)
    first_mode = 0
    number_of_modes = md.k

    if number_of_modes > md.U_1t[0].shape[1]:
        self.log.write("warning", "Number of actual modes is less than the specified K value (This can happen in DMD). Setting writing only the available modes")
        number_of_modes = md.U_1t[0].shape[1]

    # Copy the data to overwrite
    dataout=copy.deepcopy(data)

    # Retrieve the information
    qoi=params.fields
    nxyz=case.nxyz
    nxyze=case.nxyze
    nel=case.nel
    num_fields=params.num_fields

    for ii in range(0,md.number_of_decompositions):
        filename=path+'PODmod'+'_'+prefix+'_md_'+repr(ii)+'_'+casename+'0.f'+str(0+1).zfill(5)
        self.log.write("info", "Writing POD modes to: "+filename)
        pbar= tqdm(total=number_of_modes)
    
        for j in range(first_mode,first_mode+number_of_modes):
            for e in range(0,nel):
                gl_e = int(self.lglel[e]-1)
                for field in range(0,num_fields):
                    qo = qoi[field]
                    if qo != 3:
                        dataout.elem[gl_e].vel[qo,:,:,:]=np.reshape(md.U[ii][(e*nxyz+field*nxyze):(e*nxyz+field*nxyze)+nxyz,j],(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
                    elif qo==3:
                        dataout.elem[gl_e].pres[0,:,:,:]=np.reshape(md.U[ii][(e*nxyz+field*nxyze):(e*nxyz+field*nxyze)+nxyz,j],(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')

            filename=path+'PODmod'+'_'+prefix+'_md_'+repr(ii)+'_'+casename+'0.f'+str(j+1).zfill(5)
            writenek(filename,dataout)
            pbar.update(1)
        pbar.close()


def write_timecoeff(self,params,md,prefix):

    self.log.write("info", "Writing singular values and right singular vectors")

    # Retrieve the information
    path=params.dataPath
    casename=params.casename

    for ii in range(0,md.number_of_decompositions):
    
        filename1=path+'PODsig'+'_'+prefix+'_md_'+repr(ii)+'_'+casename
        filename2=path+'PODvtr'+'_'+prefix+'_md_'+repr(ii)+'_'+casename
    
        np.save(filename1,md.D_1t[ii]) 
        np.save(filename2,md.Vt_1t[ii]) 

        try:
            getattr(md, "Lambda")
            self.log.write("info","Writing eigen values of At")
            filename3=path+'PODLambda'+'_'+prefix+'_md_'+repr(ii)+'_'+casename
            np.save(filename3,md.Lambda[ii])
        except AttributeError:
            ''' Do nothing '''

    
def write_restart_file(self,params,iterator,md):

    self.log.write("info", "Writing restart file")

    # Retrieve the information
    path=params.dataPath
    casename=params.casename

    # Set up the writer
    filename=path+'restartPOD'+'_'+casename+'.hdf5'
    fw = h5py.File(filename,'w')

    # Write some attributes
    fw.attrs['isnap'] = iterator.iteration - 1
    fw.attrs['updates_to_b1'] = md.number_of_updates_buffer1
    fw.attrs['updates_to_b2'] = md.number_of_updates_buffer2
    fw.attrs['U_shape'] = md.U[0].shape
    fw.attrs['D_shape'] = md.D_1t[0].shape
    fw.attrs['Vt_shape'] = md.Vt_1t[0].shape

    # Create the data sets
    dset1 = fw.create_dataset('U',data=md.U[0], compression="gzip", compression_opts=4)
    dset2 = fw.create_dataset('D',data=md.D_1t[0], compression="gzip", compression_opts=4)
    dset3 = fw.create_dataset('Vt',data=md.Vt_1t[0], compression="gzip", compression_opts=4)

    fw.close()

def io_read_restart(self,params,case,iterator,md,comm):
        
    rank = comm.Get_rank()
    size = comm.Get_size()

    filename=params.restart_file_path
    f = h5py.File(filename, "r")
        
    #Update isnap
    iterator.iteration = f.attrs['isnap']
    md.number_of_updates_buffer1 = f.attrs['updates_to_b1']
    md.number_of_updates_buffer2 = f.attrs['updates_to_b2']
    U_shape = f.attrs['U_shape']
    D_shape = f.attrs['D_shape']
    Vt_shape = f.attrs['Vt_shape']

    #Create buffers
    U = None
    md.U_1t[0]  = np.empty((int(U_shape[0]/size),U_shape[1]))
    md.D_1t[0]  = np.empty((D_shape))
    md.Vt_1t[0] = np.empty((Vt_shape[0],Vt_shape[1]))

    if rank == 0:
        #Read the modes
        U = np.array(f['U'][:])
        md.D_1t[0] = np.array(f['D'][:])
        md.Vt_1t[0] = np.array(f['Vt'][:])

    #Communicate from rank 0
    comm.Scatter(U, md.U_1t[0], root=0)
    comm.Bcast(md.D_1t[0], root=0)
    comm.Bcast(md.Vt_1t[0], root=0)

    f.close()


    self.log.write("warning","restart file found")
    string="Last finished iteration: "+repr(iterator.iteration)
    self.log.write("warning",string)
    iterator.iteration += 1
    string="restarting from iteration: "+repr(iterator.iteration)
    self.log.write("warning",string)
    string = 'The shape of of the saved modes is U[%d,%d]' % (md.U_1t[0].shape[0], md.U_1t[0].shape[1])
    self.log.write("warning",string)

