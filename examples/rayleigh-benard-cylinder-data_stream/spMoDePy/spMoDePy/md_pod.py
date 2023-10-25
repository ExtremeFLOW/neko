
import numpy as np
from mpi4py import MPI #equivalent to the use of MPI_init() in C
from spMoDePy.ModalDecompositionInterface import ModalDecompositionInterface
from spMoDePy.mpi_spSVD import spSVD_c
from spMoDePy.math_ops import math_ops_c


class POD_c(ModalDecompositionInterface):

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
        self.number_of_decompositions = number_of_decompositions
        
        # Change k to a low value if autoupdate is required
        if params_c.ifautoupdate==True: self.k=self.p 
        self.running_ra=[]

        # Define some buffers
        U_1t = [] 
        D_1t = []
        Vt_1t = []
        for i in range(0,number_of_decompositions):
            U_1t.append(None) 
            D_1t.append(None)
            Vt_1t.append(None)


        self.U_1t = U_1t
        self.D_1t = D_1t
        self.Vt_1t = Vt_1t

        self.svd  = spSVD_c(logger)
        self.math = math_ops_c()

        self.buffer1_index = 0
        self.buffer1_max_index = self.p - 1
        
        self.number_of_updates_buffer1 = 0 
        self.number_of_updates_buffer2 = 0 
    
    def buffer_allocator(self,params,case,rows):
    
        self.buff = [] 
        buff_temp=np.zeros((rows,self.p))
        for i in range(0,(self.number_of_decompositions)):
            self.buff.append(buff_temp.copy())
        
    def firststep(self, params,case,iterator, io_controller,comm):
        if params.ifrestart: 
            self.math.scale_data(self.U_1t[0],io_controller.bm1sqrti,case.mynxyzef,self.U_1t[0].shape[1],'mult')


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

        # Calculate the residual and check if basis needs to be expanded 
        if self.number_of_updates_buffer1 >= 1:
            if params.ifautoupdate==True:
                if params.ifgbl_update == False:
                    ra=self.math.get_perp_ratio(self.U_1t[0],io_controller.Xi.reshape((-1,1)))
                    self.running_ra.append(ra)
                else:
                    ra=self.math.mpi_get_perp_ratio(self.U_1t[0],io_controller.Xi.reshape((-1,1)),comm)
                    self.running_ra.append(ra)
            else:
                ra=0
                self.running_ra.append(ra)
            if params.ifautoupdate==True and ra>=params.minimun_orthogonality_ratio and self.k<params.maximun_number_of_modes: 
                self.k+=1
                print("New k is = " +repr(self.k))


    def update_from_buffer1(self, params, case, iterator, io_controller, comm):
        
        # Get rank info
        rank     = comm.Get_rank()
        size     = comm.Get_size()
                 
        # Perform the update
        if params.ifgbl_update==True:
                self.U_1t[0],self.D_1t[0],self.Vt_1t[0] = self.svd.gblSVD_update_fromBatch(self.U_1t[0],self.D_1t[0],self.Vt_1t[0],self.buff[0][:,:(self.buffer1_index)],self.k, comm)
        else:
                self.U_1t[0],self.D_1t[0],self.Vt_1t[0] = self.svd.lclSVD_update_fromBatch(self.U_1t[0],self.D_1t[0],self.Vt_1t[0],self.buff[0][:,:(self.buffer1_index)],self.k)

        string = 'The shape of the modes after this update is U[%d,%d]' % (self.U_1t[0].shape[0], self.U_1t[0].shape[1])
        self.log.write("info",string)

        self.number_of_updates_buffer1 += 1
        
        self.log.write("info","The total number of updates performed up to now is: "+repr(self.number_of_updates_buffer1))
    
    def load_buffer2(self, *args, **kwargs):
        """Describe later"""
    
    def update_from_buffer2(self, *args, **kwargs):
        """Describe later"""

    def poststream(self, params, case, iterator, io_controller, comm):
        
        # Check if there is information in the buffer that should be taken in case the loop exit without flushing
        if self.buffer1_index > self.buffer1_max_index:
            self.log.write("info","All snapshots where properly included in the updates")
        else: 
            self.log.write("warning","Last loaded snapshot to buffer was: "+repr(self.buffer1_index-1))
            self.log.write("warning","The buffer updates when it is full to position: "+repr(self.buffer1_max_index))
            self.log.write("warning","Data must be updated now to not lose anything,  Performing an update with data in buffer 1")
            self.update_from_buffer1(params, case, iterator, io_controller, comm)

        self.log.write("info","Rescaling the obtained modes")
        # Scale the modes back before gathering them
        self.math.scale_data(self.U_1t[0],io_controller.bm1sqrti,case.mynxyzef,self.k,'div')


        # If local updates where made
        if params.ifgbl_update == False: 
            self.log.write("info","Obtaining global modes from local ones")
            ## Obtain global modes
            self.U_1t[0],self.D_1t[0],self.Vt_1t[0] = self.svd.lclSVD_to_gblSVD(self.U_1t[0],self.D_1t[0],self.Vt_1t[0],params.keep_modes,comm)
            ## Gather the orthogonality record
            ortho = comm.gather(self.running_ra,root=0)
        else:
            ortho = self.running_ra


        # Gather the modes in rank 0
        self.U= [None]
        self.U[0], io_controller.bm1sqrt2 = self.math.gatherModesandMass_atRank0(self.U_1t[0],io_controller.bm1sqrti,case.nxyzef,comm)
    
        # Gather the global element numbers in rank 0 
        comm.Gather([io_controller.lgleli,MPI.INT], [io_controller.lglel, MPI.INT], root=0)
    
