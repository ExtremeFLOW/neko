# Import general modules
import sys
import os
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6



# Import required modules
from mpi4py import MPI #equivalent to the use of MPI_init() in C
import numpy as np
import logging

# Imports from own modules
cwd = os.getcwd()
print(cwd)
sys.path.append(cwd+"/spMoDePy/")
from spMoDePy.params import params_c
from spMoDePy.case import case_c
from spMoDePy.iterator import iterator_c
from spMoDePy.md_logger import logger_c

from spMoDePy.io_controller import io_controller_c
io_dic = {"fromDisk":io_controller_c}

from spMoDePy.md_pod import POD_c
from spMoDePy.md_spod import SPOD_c
from spMoDePy.md_dmd import DMD_c
md_dic = {"POD":POD_c, "SPOD":SPOD_c, "DMD":DMD_c}


# MPI - MPMD
worldcomm = MPI.COMM_WORLD
worldrank = worldcomm.Get_rank()
worldsize = worldcomm.Get_size()
col = 1
comm = worldcomm.Split(col,worldrank)
rank = comm.Get_rank()
size = comm.Get_size()

# Read the input parameters
params = params_c("inputs.json")
if rank == 0: params.write_string()

# Instance the logger
logger = logger_c(level = logging.DEBUG, comm = comm)

# Import adios2 if requested
if params.execution_type == "fromAdios":
    sys.path.append(params.adiosPath)
    import adios2
    logger.write("info","Adios2 has been inported from path: "+params.adiosPath)

    if rank == 0: print("Adios2 has been imported from path: "+params.adiosPath)
    from spMoDePy.io_controller_adios import io_controller_adios_c
    io_dic["fromAdios"] = io_controller_adios_c
        
# Obtain the case info 
case = case_c(params,comm)
if rank == 0: case.write_string()

# Instance the modal decomposition of your choosing
logger.write("info","The chosen decomposition was: "+params.decomposition)
md = md_dic[params.decomposition](params_c=params, logger = logger, comm = comm)

# Instance the controller for io
io_controller = io_dic[params.execution_type](params,case,logger,comm)

# Instance the iterator
iterator = iterator_c(params,case,logger,comm)

while iterator.status == True:

    # Perform IO operations
    if (iterator.first_iteration):
        logger.write("info","Perform actions for the first step")
        io_controller.first_step(params,case,iterator,md,comm)  
    else:
        io_controller.step(params,case,iterator,comm)  
    
    # Perform the streaming modal decomposition
    if iterator.run_update_loop:
        
        if iterator.first_iteration:
            md.firststep(params,case,iterator, io_controller,comm)

        if iterator.load_buffer1:
            md.load_buffer1(params,case,iterator,io_controller,comm)
        
        if iterator.update_from_buffer1:
            md.update_from_buffer1(params,case,iterator,io_controller,comm)
        
        if iterator.load_buffer2:
            md.load_buffer2(params,case,iterator,io_controller,comm)
        
        if iterator.update_from_buffer2:
            md.update_from_buffer2(params,case,iterator,io_controller,comm)


    iterator.mark_iteration()
    iterator.check_status(io_controller)
    logger.write("info","Finished iteration " +repr(iterator.iteration-1))


# Close the io operations
io_controller.close()

# Execute actions post data stream
md.poststream(params,case,iterator,io_controller,comm)

# write data to disk
io_controller.write(params,case,iterator,md,comm)

logger.write("critical","Execution was finalized succesfuly")

