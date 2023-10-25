
import numpy as np
from mpi4py import MPI #equivalent to the use of MPI_init() in C
from spMoDePy.ModalDecompositionInterface import ModalDecompositionInterface
from spMoDePy.mpi_spSVD import spSVD_c
from spMoDePy.math_ops import math_ops_c
import sys

class SPOD_c(ModalDecompositionInterface):

    def __init__(self,
                params_c, logger, comm,
                number_of_decompositions: int=1):

        rank     = comm.Get_rank()
        size     = comm.Get_size()

        self.log = logger
        
        # Parse from params file
        self.k         = params_c.keep_modes
        self.setk      = params_c.keep_modes
        self.p      = params_c.batch_size
        self.number_of_decompositions = params_c.number_of_decompositions
        self.window_size = params_c.window_size
        self.frequency_locations = params_c.frequency_locations

        # Change k to a low value if autoupdate is required
        if params_c.ifautoupdate==True: self.k=self.p 
        self.running_ra=[]

        # Define some buffers
        U_1t = [] 
        D_1t = []
        Vt_1t = []
        for i in range(0,self.number_of_decompositions):
            U_1t.append(None) 
            D_1t.append(None)
            Vt_1t.append(None)

        self.U_1t = U_1t
        self.D_1t = D_1t
        self.Vt_1t = Vt_1t
         
        self.svd  = spSVD_c(logger)
        self.math = math_ops_c()

        self.buffer1_index = 0
        self.buffer1_max_index = self.window_size - 1
        
        self.buffer2_index = 0
        self.buffer2_max_index = self.p - 1
        
        self.number_of_updates_buffer1 = 0 
        self.number_of_updates_buffer2 = 0 

    def buffer_allocator(self,params,case,rows):
        
        self.buff = [] 
        buff_temp=np.zeros((rows,self.window_size))
        for i in range(0,1):
            self.buff.append(buff_temp.copy())
        
        self.buff2 = [] 
        buff_temp=np.zeros((rows,self.p),dtype=np.complex_)
        for i in range(0,(self.number_of_decompositions)):
            self.buff2.append(buff_temp.copy())
    

    def firststep(self, params,case,iterator, io_controller,comm):
        if params.ifrestart:
            self.log.write("critical","Restart not properly implemented for SPOD. Finalizing run")
            sys.exit(1)


    def load_buffer1(self, params, case, iterator, io_controller, comm):
        
        if self.buffer1_index > self.buffer1_max_index:
            self.buffer1_index = 0

        #Scale the data with the mass matrix
        self.math.scale_data(io_controller.Xi,io_controller.bm1sqrti,case.mynxyzef,1,'mult')

        #Fill the buffer
        self.buff[0][:,self.buffer1_index] = np.copy(io_controller.Xi[:,0]) 
        
        if self.buffer1_index == self.buffer1_max_index:
            self.log.write("info","Loaded snapshot in buffer 1 in pos: "+repr(self.buffer1_index))
            self.log.write("info", "Buffer one is full, proceed to update")
            self.buffer1_index += 1
            iterator.update_from_buffer1 = True
        else: 
            self.log.write("info","Loaded snapshot in buffer 1 in pos: "+repr(self.buffer1_index))
            self.buffer1_index += 1
            iterator.update_from_buffer1 = False

    def update_from_buffer1(self, params, case, iterator, io_controller, comm): 
        # Get rank info
        rank     = comm.Get_rank()
        size     = comm.Get_size()
                
        self.fftbuff = np.fft.fft(self.buff[0], n=None, axis=-1, norm=None)/self.window_size    

        self.log.write("info","FFT performed on buffer 1")

        self.number_of_updates_buffer1 += 1
             
        self.log.write("info","proceding to load frquency data to buffer2")
        iterator.load_buffer2 = True
    
    def load_buffer2(self, params, case, iterator, io_controller, comm):
        
        if self.buffer2_index > self.buffer2_max_index:
            self.buffer2_index = 0

        #Fill the buffer
        for i in range(0,(self.number_of_decompositions)):
            self.buff2[i][:,self.buffer2_index] = np.copy(self.fftbuff[:,self.frequency_locations[i]]) 
        
        # Reset the iterator loading sequence
        iterator.load_buffer2 = False

        if self.buffer2_index == self.buffer2_max_index:
            self.log.write("info","Loaded snapshot in buffer 2 in pos: "+repr(self.buffer2_index))
            self.log.write("info", "Buffer two is full, proceed to update")
            self.buffer2_index += 1
            iterator.update_from_buffer2 = True
        else: 
            self.log.write("info","Loaded snapshot in buffer 2 in pos: "+repr(self.buffer2_index))
            self.buffer2_index += 1
            iterator.update_from_buffer2 = False
        
        # Calculate the residual and check if basis needs to be expanded 
        if self.number_of_updates_buffer1 >= 1:
            if params.ifautoupdate==True:
                self.log.write("critical","Restart not properly implemented for SPOD. Finalizing run")
                sys.exit(1)
    
    def update_from_buffer2(self, params, case, iterator, io_controller, comm): 
        
        # Get rank info
        rank     = comm.Get_rank()
        size     = comm.Get_size()
                 
        # Perform the update
        if params.ifgbl_update==True:
            for i in range(0,(self.number_of_decompositions)):
                self.U_1t[i],self.D_1t[i],self.Vt_1t[i] = self.svd.gblSVD_update_fromBatch(self.U_1t[i],self.D_1t[i],self.Vt_1t[i],self.buff2[i][:,:(self.buffer2_index)],self.k, comm)
        else:
            for i in range(0,(self.number_of_decompositions)):
                self.U_1t[i],self.D_1t[i],self.Vt_1t[i] = self.svd.lclSVD_update_fromBatch(self.U_1t[i],self.D_1t[i],self.Vt_1t[i],self.buff2[i][:,:(self.buffer2_index)],self.k)

        string = 'The shape of the modes after this update is U[%d,%d]' % (self.U_1t[0].shape[0], self.U_1t[0].shape[1])
        self.log.write("info",string)

        # Reset iterator loading sequence
        iterator.update_from_buffer2 = False

        self.number_of_updates_buffer2 += 1

    def poststream(self, params, case, iterator, io_controller, comm):
        
        # Check if there is information in the buffer that should be taken in case the loop exit without flushing
        if self.buffer1_index > self.buffer1_max_index:
            self.log.write("info","All snapshots in buffer 1 were properly included in the updates")
        else: 
            self.log.write("critical","Last loaded snapshot to buffer was: "+repr(self.buffer1_index-1))
            self.log.write("critical","The buffer updates when it is full to position: "+repr(self.buffer1_max_index))
            self.log.write("critical","Can not perform FFT in incomplete realization. The data in the buffer will be lost")
        
        if self.buffer2_index > self.buffer2_max_index:
            self.log.write("info","All snapshots in buffer 2 were properly included in the updates")
        else: 
            self.log.write("warning","Last loaded snapshot to buffer was: "+repr(self.buffer2_index-1))
            self.log.write("warning","The buffer updates when it is full to position: "+repr(self.buffer2_max_index))
            self.log.write("warning","Data must be updated now to not lose anything,  Performing an update with data in buffer 2")
            self.update_from_buffer2(params, case, iterator, io_controller, comm)

        self.log.write("info","Rescaling the obtained modes")
        # Scale the modes back before gathering them
        for i in range(0,(self.number_of_decompositions)):
            self.math.scale_data(self.U_1t[i],io_controller.bm1sqrti,case.mynxyzef,self.k,'div')


        # If local updates where made
        if params.ifgbl_update == False: 
            self.log.write("info","Obtaining global modes from local ones")
            ## Obtain global modes
            for i in range(0,(self.number_of_decompositions)):
                self.U_1t[i],self.D_1t[i],self.Vt_1t[i] = self.svd.lclSVD_to_gblSVD(self.U_1t[i],self.D_1t[i],self.Vt_1t[i],params.keep_modes,comm)
            ## Gather the orthogonality record
            ortho = comm.gather(self.running_ra,root=0)
        else:
            ortho = self.running_ra


        # Gather the modes in rank 0
        self.U= []
        for i in range(0,(self.number_of_decompositions)):
            U_, io_controller.bm1sqrt2 = self.math.gatherModesandMass_atRank0(self.U_1t[i],io_controller.bm1sqrti,case.nxyzef,comm)
            if comm.Get_rank() == 0: self.U.append(U_.copy())

        # Gather the global element numbers in rank 0 
        comm.Gather([io_controller.lgleli,MPI.INT], [io_controller.lglel, MPI.INT], root=0)
    
