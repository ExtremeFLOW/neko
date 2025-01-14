# Import general modules
import sys
import os
import json
import numpy as np

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

import adios2

# Import data streamer
from pynektools.io.adios2.stream import DataStreamer
from pynektools.io.utils import get_fld_from_ndarray

# Instance the streamer (start the stream)
ds = DataStreamer(comm)
    
# Recieve the data from fortran
x = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv) # Recieve and reshape x
y = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv) # Recieve and reshape y
z = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv) # Recieve and reshape z

# Finalize the stream 
ds.finalize()

# Now that the data is here. Create a mesh object, an empty field and write it to disk
# For this, the pyNek tools are needed
from pynektools.io.ppymech.neksuite import pwritenek
from pynektools.datatypes.msh import Mesh
from pynektools.datatypes.field import Field
from pynektools.datatypes.utils import create_hexadata_from_msh_fld

## Instance the mesh
msh = Mesh(comm, x = x, y = y, z = z)
## Create an empty field and update its metadata
out_fld = Field(comm)
out_fld.fields["scal"].append(np.ones_like(msh.x))
out_fld.update_vars()
## Create the hexadata to write out
out_data = create_hexadata_from_msh_fld(msh = msh, fld = out_fld)
## Write out a file
fname = "python_test0.f"+str(1).zfill(5)
pwritenek("./"+fname,out_data, comm)
