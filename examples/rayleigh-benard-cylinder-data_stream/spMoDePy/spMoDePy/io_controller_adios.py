import sys
import os
import copy

import numpy as np
from pymech.neksuite import readnek
from pymech.neksuite import writenek
from tqdm import tqdm
import h5py
NoneType = type(None)
        
import adios2

class io_controller_adios_c():
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

        #Adios status
        self.okstep = adios2.StepStatus.OK
        self.endstream = adios2.StepStatus.EndOfStream

        # ADIOS2 instance
        self.adios = adios2.ADIOS(comm)
        # ADIOS IO - Engine
        self.ioRead = self.adios.DeclareIO("ioReader")
        self.ioRead.SetEngine('SST')
        
        # Open the stream
        self.ibpStream = self.ioRead.Open("globalArray", adios2.Mode.Read, comm)

    # first step
    def first_step(self,params,case,iterator,md,comm):

        if params.ifrestart:
            io_read_restart(self,params,case,iterator,md,comm)

        self.log.write("warning", "Remember that ADIOS2 reads the data at the end of the step")

        self.stepStatus = self.ibpStream.BeginStep()

        if self.stepStatus == adios2.StepStatus.OK:

            self.var_inVX = self.ioRead.InquireVariable("VX")
            self.var_inVY = self.ioRead.InquireVariable("VY")
            self.var_inVZ = self.ioRead.InquireVariable("VZ")
            self.var_inBM1 = self.ioRead.InquireVariable("BM1")
            self.var_inLGLEL = self.ioRead.InquireVariable("LGLEL")

            io_allocate(self,params,case,md,comm)
            io_read(self,params,case,iterator,comm)
        
            self.ibpStream.EndStep() # The data is actually read when the step ends!
        
            for field in range(0, params.num_fields):
                self.Xi[(field*case.my_count):(field*case.my_count)+case.my_count]=np.copy(self.Xtemp[field])
                self.bm1sqrti[(field*case.my_count):(field*case.my_count)+case.my_count]=np.copy(np.sqrt(self.bm1i))


        elif self.stepStatus == adios2.StepStatus.EndOfStream:
            self.log.write("warning", "Signal to end the stream given from Adios. Skipping the update loop")
            iterator.run_update_loop = False

    def step(self,params,case,iterator,comm):
        
        self.stepStatus = self.ibpStream.BeginStep()

        if self.stepStatus == adios2.StepStatus.OK:

            self.var_inVX = self.ioRead.InquireVariable("VX")
            self.var_inVY = self.ioRead.InquireVariable("VY")
            self.var_inVZ = self.ioRead.InquireVariable("VZ")
            self.var_inBM1 = self.ioRead.InquireVariable("BM1")
            self.var_inLGLEL = self.ioRead.InquireVariable("LGLEL")

            io_read(self,params,case,iterator,comm)
        
            self.ibpStream.EndStep() # The data is actually read when the step ends!
        
            for field in range(0, params.num_fields):
                self.Xi[(field*case.my_count):(field*case.my_count)+case.my_count]=np.copy(self.Xtemp[field])
                self.bm1sqrti[(field*case.my_count):(field*case.my_count)+case.my_count]=np.copy(np.sqrt(self.bm1i))


        elif self.stepStatus == adios2.StepStatus.EndOfStream:
            iterator.run_update_loop = False
            self.log.write("warning", "Signal to end the stream given from Adios. Skipping the update loop")

    def close(self):
        self.log.write("warning", "Stream of data is over")
        self.ibpStream.Close()

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

def io_allocate(self, params,case,md, comm):
    
    # Get information from the communicator
    rank = comm.Get_rank()
    size = comm.Get_size()
    num_fields = params.num_fields
    p = params.batch_size

    if (self.var_inVX) is not None:
        total_size=self.var_inVX.Shape()[0]
        my_count= int(total_size/size)
        my_start= int(my_count*rank)
        if total_size % size != 0:
            if rank < (total_size % size):
                my_count += 1
                my_start += rank
            else:
                my_start += (total_size%size)    

        #Allocate arrays
        Xtemp1 = np.zeros((my_count,1), dtype=np.double)
        self.Xtemp = []
        for i in range(0,num_fields):
            self.Xtemp.append(Xtemp1.copy())

        self.Xi = np.zeros((my_count*num_fields,1), dtype=np.double) 
        md.buffer_allocator(params,case,int(my_count*num_fields))
        
        #Allocate arrays
        self.bm1i = np.zeros((my_count,1), dtype=np.double)
        self.bm1sqrti = np.zeros((my_count*num_fields,1), dtype=np.double)
            
    if self.var_inLGLEL is not None:
        total_size2=self.var_inLGLEL.Shape()[0]
        my_count2= int(total_size2/size)
        my_start2= int(my_count2*rank)
        if total_size2 % size != 0:
            if rank < (total_size2 % size):
                my_count2 += 1
                my_start2 += rank
            else:
                my_start2 += (total_size2%size)
                        
        #Allocate arrays
        self.lgleli = np.zeros((my_count2,1), dtype=np.int32)

        case.my_start  = my_start
        case.my_count  = my_count
        case.my_start2 = my_start2
        case.my_count2 = my_count2
        
        case.mynxyzef = int(my_count*num_fields)


def io_read(self, params,case,iterator, comm):
        
        # Associate
        qoi=params.fields
        num_fields=params.num_fields
        nxyz=case.nxyz    
        nel=case.nel
        mycount = case.my_count
        mynxyzef = case.mynxyzef

        #Select the data in the global array that belongs to me
        self.var_inVX.SetSelection([[case.my_start], [case.my_count]])
        self.var_inVY.SetSelection([[case.my_start], [case.my_count]])
        self.var_inVZ.SetSelection([[case.my_start], [case.my_count]])
        self.var_inBM1.SetSelection([[case.my_start], [case.my_count]])
        self.var_inLGLEL.SetSelection([[case.my_start2], [case.my_count2]])
    
        possible_fields = {0: self.var_inVX,
                           1: self.var_inVY,
                           2: self.var_inVZ
                          }

        #Read the variable into array
        for field in range(0, num_fields):
            qo = qoi[field]
            self.ibpStream.Get(possible_fields[qo], self.Xtemp[field])
            #self.Xi[(field*mycount):(field*mycount)+mycount]=np.copy(self.Xtemp[field])

        #Read and process mass
        self.ibpStream.Get(self.var_inBM1, self.bm1i)
        #for field in range(0, num_fields):
            #self.bm1sqrti[(field*mycount):(field*mycount)+mycount]=np.copy(np.sqrt(self.bm1i))

        #Read global element numbers
        self.ibpStream.Get(self.var_inLGLEL, self.lgleli)

def write_modes(self,params, case,md,prefix):

    path=params.dataPath
    casename=params.casename
    filename=path+casename+'0.f'+str(0).zfill(5)
    data=readnek(filename)
    first_mode = 0
    number_of_modes = md.k
    
    if number_of_modes > md.U_1t[0].shape[1]:
        self.log.write("warning", "Number of actual modes is less than the specified K value (This can happen in DMD). Setting writing only the available modes")

    # Copy the data to overwrite
    dataout=copy.deepcopy(data)

    # Retrieve the information
    qoi=params.fields
    nxyz=case.nxyz
    nxyze=case.nxyze
    nel=case.nel
    num_fields=params.num_fields

    for ii in range(0,md.number_of_decompositions):
        
        filename=path+'PODmod'+'_'+prefix+'_md_'+repr(ii)+'_'+casename+'0.f'+str(0).zfill(5)
        self.log.write("info", "Writing POD modes to: "+filename)
        pbar= tqdm(total=number_of_modes)
    
        my_count = case.my_count
        my_count2 = case.my_count2
            
        for j in range(first_mode,first_mode+number_of_modes):
            for field in range(0,num_fields):
                rank_id = -1
                qo = qoi[field]
                
                for e in range(0,nel):
                    gl_e = int(self.lglel[e]-1)
                    module = (np.mod(e,case.my_count2))
                    if module==0:
                        rank_id = rank_id + 1
                    
                    low_lim = (e*nxyz+rank_id*my_count*(num_fields-1)) + field * my_count2 * nxyz
                    hig_lim = (e*nxyz+rank_id*my_count*(num_fields-1)) + nxyz + field * my_count2 * nxyz
                    #print(repr(rank_id)+","+repr(gl_e)+","+repr(field)+","+repr(low_lim)+","+repr(hig_lim))
                    dataout.elem[gl_e].vel[qo,:,:,:]=np.reshape(md.U[ii][low_lim:hig_lim,j],(data.lr1[2],data.lr1[1],data.lr1[0]))

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

