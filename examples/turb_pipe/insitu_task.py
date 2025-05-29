#========================================
# Import and set up general modules
#========================================
import sys
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
# Import MPI
from mpi4py import MPI #equivalent to the use of MPI_init() in C

# Split communicator for MPI - MPMD
worldcomm = MPI.COMM_WORLD
worldrank = worldcomm.Get_rank()
worldsize = worldcomm.Get_size()
col = 1
comm = worldcomm.Split(col,worldrank)
rank = comm.Get_rank()
size = comm.Get_size()

#========================================
# Import modules
#========================================
# general functionality
import numpy as np

# external
import adios2.bindings as adios2

# pysemtools
from pysemtools.io.adios2.stream import DataStreamer
from pysemtools.io.utils import get_fld_from_ndarray
from pysemtools.io.ppymech.neksuite import pynekwrite
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.coef import Coef
from pysemtools.datatypes.field import FieldRegistry
from pysemtools.rom.pod import POD
from pysemtools.rom.io_help import IoHelp

#=========================================
# Define inputs
#=========================================

if comm.Get_rank() == 0:
    print("Python - insitu - reading inputs")

# Read the POD inputs
number_of_pod_fields = 3
pod_batch_size  = 30
pod_keep_modes  = 10
pod_write_modes = 10
dtype_string = "double"
backend = "numpy"
if dtype_string == "single":
    dtype = np.float32
else:
    dtype = np.float64

# Start time
start_time = MPI.Wtime()

#=========================================
# Initialize the streamer and get the mesh
#=========================================

if comm.Get_rank() == 0:
    print("Python - insitu - Initializing streamer")

# Initialize the streamer
ds = DataStreamer(comm)

if comm.Get_rank() == 0:
    print("Python - insitu - Initializing objects")
    
# Recieve the mesh
x = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) # Recieve and reshape x
y = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) # Recieve and reshape y
z = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) # Recieve and reshape z

# Initialize the mesh and mass matrix
msh = Mesh(comm, x = x, y = y, z = z, create_connectivity = False)
coef = Coef(msh, comm)
bm = coef.B

#=========================================
# Initialize the POD
#=========================================

# Instance the POD object
pod = POD(comm, number_of_modes_to_update = pod_keep_modes, global_updates = True, auto_expand = False, bckend = backend)

# Instance io helper that will serve as buffer for the snapshots
ioh = IoHelp(comm, number_of_fields = number_of_pod_fields, batch_size = pod_batch_size, field_size = bm.size, field_data_type=dtype, mass_matrix_data_type=dtype)

# Put the mass matrix in the appropiate format (long 1d array)
mass_list = []
for i in range(0, number_of_pod_fields):
    mass_list.append(np.copy(np.sqrt(bm)))
ioh.copy_fieldlist_to_xi(mass_list)
ioh.bm1sqrt[:,:] = np.copy(ioh.xi[:,:])

#=========================================
# Perform the streaming of data
#=========================================

stream_data = True
while stream_data:

    # Get the data
    u = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) 
    v = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) 
    w = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) 

    # Check if data was recieved or if the stream ended
    if ds.step_status == adios2.StepStatus.OK:
        stream_data = True
    elif ds.step_status == adios2.StepStatus.EndOfStream:
        stream_data = False
        continue

    # Put the snapshot data into a column array
    ioh.copy_fieldlist_to_xi([u, v, w])

    # Load the column array into the buffer
    ioh.load_buffer(scale_snapshot = True)

    # Update POD modes
    if ioh.update_from_buffer:
        pod.update(comm, buff = ioh.buff[:,:(ioh.buffer_index)])

#=========================================
# Perform post-stream operations
#=========================================

# Finilize the streamer
ds.finalize()

# Check if there is information in the buffer that should be taken in case the loop exit without flushing
if ioh.buffer_index > ioh.buffer_max_index:
    ioh.log.write("info","All snapshots where properly included in the updates")
else: 
    ioh.log.write("warning","Last loaded snapshot to buffer was: "+repr(ioh.buffer_index-1))
    ioh.log.write("warning","The buffer updates when it is full to position: "+repr(ioh.buffer_max_index))
    ioh.log.write("warning","Data must be updated now to not lose anything,  Performing an update with data in buffer ")
    pod.update(comm, buff = ioh.buff[:,:(ioh.buffer_index)])

# Scale back the modes
pod.scale_modes(comm, bm1sqrt = ioh.bm1sqrt, op = "div")

# Rotate local modes back to global, This only enters in effect if global_update = false
pod.rotate_local_modes_to_global(comm)

#=========================================
# Write out the modes
#=========================================

# Write the data out
for j in range(0, pod_write_modes):

    if (j+1) < pod.u_1t.shape[1]:

        ## Split the snapshots into the proper fields
        field_list1d = ioh.split_narray_to_1dfields(pod.u_1t[:,j])
        u_mode = get_fld_from_ndarray(field_list1d[0], msh.lx, msh.ly, msh.lz, msh.nelv) 
        v_mode = get_fld_from_ndarray(field_list1d[1], msh.lx, msh.ly, msh.lz, msh.nelv) 
        w_mode = get_fld_from_ndarray(field_list1d[2], msh.lx, msh.ly, msh.lz, msh.nelv) 

        # write the data
        fld = FieldRegistry(comm)        
        fld.add_field(comm, field_name = "u", field = u_mode, dtype = dtype)
        fld.add_field(comm, field_name = "v", field = v_mode, dtype = dtype)
        fld.add_field(comm, field_name = "w", field = w_mode, dtype = dtype)
        pynekwrite(f"./modes0.f{str(j).zfill(5)}", comm=comm, msh=msh, fld=fld, wdsz=4, istep = j) 
        
#=========================================
# Write out singular values and right 
# singular vectors
#=========================================

# Write the singular values and vectors
if comm.Get_rank() == 0:
    np.save("singular_values", pod.d_1t)
    print("Wrote signular values")
    np.save("right_singular_vectors", pod.vt_1t)
    print("Wrote right signular values")

# End time
end_time = MPI.Wtime()
# Print the time
if comm.Get_rank() == 0:
    print("Time to complete: ", end_time - start_time)
