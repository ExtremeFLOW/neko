
import numpy as np
from mpi4py import MPI #equivalent to the use of MPI_init() in C
from spMoDePy.ModalDecompositionInterface import ModalDecompositionInterface
from spMoDePy.mpi_spSVD import spSVD_c
from spMoDePy.math_ops import math_ops_c
import sys


class DMD_c(ModalDecompositionInterface):

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
        U_1t_prm = [] 
        D_1t_prm = []
        Vt_1t_prm = []
        for i in range(0,self.number_of_decompositions):
            U_1t.append(None) 
            D_1t.append(None)
            Vt_1t.append(None)
            U_1t_prm.append(None) 
            D_1t_prm.append(None)
            Vt_1t_prm.append(None)

        self.U_1t = U_1t
        self.D_1t = D_1t
        self.Vt_1t = Vt_1t
        
        self.U_1t_prm = U_1t_prm
        self.D_1t_prm = D_1t_prm
        self.Vt_1t_prm = Vt_1t_prm
         
        self.svd  = spSVD_c(logger)
        self.math = math_ops_c()

        self.buffer1_index = 0
        self.buffer1_max_index = 1 - 1
        
        self.buffer2_index = 0
        self.buffer2_max_index = self.p - 1
        
        self.number_of_updates_buffer1 = 0 
        self.number_of_updates_buffer2 = 0 

    def buffer_allocator(self,params,case,rows):
        
        self.buff = [] 
        buff_temp=np.zeros((rows,1))
        for i in range(0,1):
            self.buff.append(buff_temp.copy())
        
        self.buff2 = [] 
        buff_temp=np.zeros((rows,self.p))
        for i in range(0,(self.number_of_decompositions)):
            self.buff2.append(buff_temp.copy())
    

    def firststep(self, params,case,iterator, io_controller,comm):
        if params.ifrestart:
            self.log.write("critical","Restart not properly implemented for DMD. Finalizing run")
            sys.exit(1)


    def load_buffer1(self, params, case, iterator, io_controller, comm):
        
        if self.buffer1_index > self.buffer1_max_index:
            self.buffer1_index = 0

        ##Scale the data with the mass matrix
        #self.math.scale_data(io_controller.Xi,io_controller.bm1sqrti,case.mynxyzef,1,'mult')
        
        if not iterator.first_iteration: 
            self.log.write("info","Copy data in buffer one position "+repr(self.buffer1_index)+ " to array Xi_previous")
            self.Xi_previous[:,0] = np.copy(self.buff[0][:,self.buffer1_index]) 

        #Fill the buffer
        self.buff[0][:,self.buffer1_index] = np.copy(io_controller.Xi[:,0]) 
        
        if self.buffer1_index == self.buffer1_max_index and iterator.first_iteration:
            self.Xi_previous = np.zeros_like(io_controller.Xi)
            self.log.write("info","Loaded snapshot in buffer 1 in pos: "+repr(self.buffer1_index))
            self.log.write("info", "It is the first step, do not perform any update or buffer load")
            self.buffer1_index += 1
            iterator.update_from_buffer1 = False
            iterator.load_buffer2 = False
        
        elif self.buffer1_index == self.buffer1_max_index:
            self.log.write("info","Loaded snapshot in buffer 1 in pos: "+repr(self.buffer1_index))
            self.log.write("info", "Buffer one is full, proceed to load buffer 2")
            self.buffer1_index += 1
            iterator.update_from_buffer1 = False
            iterator.load_buffer2 = True

    def update_from_buffer1(self, params, case, iterator, io_controller, comm): 
        # Get rank info
        rank     = comm.Get_rank()
        size     = comm.Get_size()
                 
        # Perform the update
        if params.ifgbl_update==True:
            for i in range(0,(self.number_of_decompositions)):
                self.U_1t_prm[i],self.D_1t_prm[i],self.Vt_1t_prm[i] = self.svd.gblSVD_update_fromBatch(self.U_1t[i],self.D_1t[i],self.Vt_1t[i],self.buff[i][:,:(self.buffer1_index)],self.k, comm)
        else:
            for i in range(0,(self.number_of_decompositions)):
                self.U_1t_prm[i],self.D_1t_prm[i],self.Vt_1t_prm[i] = self.svd.lclSVD_update_fromBatch(self.U_1t[i],self.D_1t[i],self.Vt_1t[i],self.buff[i][:,:(self.buffer1_index)],self.k)

        string = 'The shape of the modes after this update is U[%d,%d]' % (self.U_1t_prm[0].shape[0], self.U_1t_prm[0].shape[1])
        self.log.write("info",string)

        # Reset iterator loading sequence
        iterator.update_from_buffer2 = False

        self.number_of_updates_buffer2 += 1
    
    def load_buffer2(self, params, case, iterator, io_controller, comm):
        
        if self.buffer2_index > self.buffer2_max_index:
            self.buffer2_index = 0

        #Fill the buffer
        for i in range(0,(self.number_of_decompositions)):
            self.buff2[i][:,self.buffer2_index] = np.copy(self.Xi_previous[:,0]) 
        
        # Reset the iterator loading sequence
        iterator.load_buffer2 = False

        if self.buffer2_index == self.buffer2_max_index:
            self.log.write("info","Loaded Xi_previous in buffer 2 in pos: "+repr(self.buffer2_index))
            self.log.write("info", "Buffer two is full, proceed to update")
            self.buffer2_index += 1
            iterator.update_from_buffer2 = True
        else: 
            self.log.write("info","Loaded Xi_previous in buffer 2 in pos: "+repr(self.buffer2_index))
            self.buffer2_index += 1
            iterator.update_from_buffer2 = False
        
        # Calculate the residual and check if basis needs to be expanded 
        if self.number_of_updates_buffer1 >= 1:
            if params.ifautoupdate==True:
                self.log.write("critical","Autoupdate not properly implemented for DMD. Finalizing run")
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
            self.log.write("info","The shifted snapshot is properly loaded in the buffer 1")
        else:
            self.log.write("warning","There is no shifted snapshot in the buffer 1")

        
        if self.buffer2_index > self.buffer2_max_index:
            self.log.write("info","All snapshots in buffer 2 were properly included in the updates")
        else: 
            self.log.write("warning","Last loaded snapshot to buffer was: "+repr(self.buffer2_index-1))
            self.log.write("warning","The buffer updates when it is full to position: "+repr(self.buffer2_max_index))
            self.log.write("warning","Data must be updated now to not lose anything,  Performing an update with data in buffer 2")
            self.update_from_buffer2(params, case, iterator, io_controller, comm)

        #self.log.write("info","Rescaling the obtained modes")
        ## Scale the modes back before gathering them
        #for i in range(0,(self.number_of_decompositions)):
        #    self.math.scale_data(self.U_1t[i],io_controller.bm1sqrti,case.mynxyzef,self.k,'div')
 
        self.log.write("info","Last loaded snapshot to buffer 1 was: "+repr(self.buffer1_index-1))
        self.log.write("info","Performing update from buffer 1 to obtain shifted decomposition")
        self.update_from_buffer1(params, case, iterator, io_controller, comm)

        # If local updates were made
        if params.ifgbl_update == False: 
            self.log.write("info","Obtaining global modes from local ones")
            ## Obtain global modes
            for i in range(0,(self.number_of_decompositions)):
                self.U_1t[i],self.D_1t[i],self.Vt_1t[i] = self.svd.lclSVD_to_gblSVD(self.U_1t[i],self.D_1t[i],self.Vt_1t[i],params.keep_modes,comm)
            
            for i in range(0,(self.number_of_decompositions)):
                self.U_1t_prm[i],self.D_1t_prm[i],self.Vt_1t_prm[i] = self.svd.lclSVD_to_gblSVD(self.U_1t_prm[i],self.D_1t_prm[i],self.Vt_1t_prm[i],params.keep_modes,comm)
            
            ## Gather the orthogonality record
            ortho = comm.gather(self.running_ra,root=0)
        else:
            ortho = self.running_ra

        
        self.log.write("info","Obtaining DMD modes from shifted Data sets")
        self.At = [] 
        self.Lambda = []
        self.W = []
        self.U_1t_dmd = []
        for i in range(0,(self.number_of_decompositions)):
            # Calculate the partial matrix in each rank
            partial_At_i = self.U_1t[i].conj().T @ self.U_1t_prm[i]
            # Allocate the total matrix
            partial_At = np.empty_like(partial_At_i)
            # Calculate the total matrix by summing the contents of each rank. Do this in all to all so every rank has it
            comm.Allreduce(partial_At_i, partial_At, op=MPI.SUM)
            # Calculate the actual At matrix
            ## What you do here is: At1=U_str.conj().T @ S @ Vt_str.conj().T @ np.diag(1/D_str)
            ## Further simplified to: U_str.conj().T @ U_prm@np.diag(D_prm)@Vt_prm[:,1:] @ Vt_str.conj().T @ np.diag(1/D_str)
            At = partial_At @np.diag(self.D_1t_prm[i])@self.Vt_1t_prm[i][:,1:] @ self.Vt_1t[i].conj().T @ np.diag(1/self.D_1t[i])

            self.At.append(At.copy())
        
            string = 'The shape of At after this update is At[%d,%d]' % (self.At[0].shape[0], self.At[0].shape[1])
            self.log.write("info",string)

            #Perform eigen decomposition of the advancing matrix
            Lambda,W=np.linalg.eig(At)

            self.Lambda.append(Lambda.copy())
            self.W.append(W.copy())

            #Get the DMD modes from the proyected ones (the eigenvectors of the reduced At)
            self.U_1t_dmd.append(self.U_1t[i]@W)

        # Gather the modes in rank 0
        self.U= []
        for i in range(0,(self.number_of_decompositions)):
            U_, io_controller.bm1sqrt2 = self.math.gatherModesandMass_atRank0(self.U_1t_dmd[i],io_controller.bm1sqrti,case.nxyzef,comm)
            if comm.Get_rank() == 0: self.U.append(U_.copy())

        # Gather the global element numbers in rank 0 
        comm.Gather([io_controller.lgleli,MPI.INT], [io_controller.lglel, MPI.INT], root=0)
