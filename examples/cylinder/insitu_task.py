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
from os.path import join

# external
import adios2.bindings as adios2

# pysemtools
from pysemtools.io.adios2.stream import DataStreamer
from pysemtools.io.utils import get_fld_from_ndarray
from pysemtools.datatypes.msh import Mesh
from pysemtools.interpolation.probes import Probes

# For plotting
import matplotlib.pyplot as plt

import json

#=========================================
# Define some helper functions
#=========================================

def init_mesh_from_adios2(x, y, z, ds: DataStreamer):
    """
    Initialize a Mesh object from coordinates and 
    """
    msh = Mesh(comm, x = x, y = y, z = z)
    msh.lx = ds.lx
    msh.ly = ds.ly
    msh.lz = ds.lz
    msh.gdim = ds.gdim 
    msh.nelv = ds.nelv
    msh.glb_nelv = ds.glb_nelv
    return msh

def log(msg):
    if comm.Get_rank() == 0:
        print("[INFO]", msg)

def get_output_directory(fname):
    with open(fname, 'r') as f:
        data = json.load(f)

    return data["case"]["output_directory"]

#=========================================
# Define some variables/parameters
#=========================================

log("Python - insitu - reading inputs")

dtype_string = "double"
backend = "numpy"
if dtype_string == "single":
    dtype = np.float32
else:
    dtype = np.float64

output_path = get_output_directory("cylinder_insitu.case")
log(f"Outputting insitu snapshots to folder {output_path}")

#=========================================
# Initialize the streamer
#=========================================

log("Python - insitu - Initializing streamer")
ds = DataStreamer(comm)

#=========================================
# Stream the mesh coordinates to build interpolators
#=========================================

log("Receiving mesh...")
x = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) 
y = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) 
z = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) 
log("Data received")

log("Initializing mesh...")
msh = init_mesh_from_adios2(x, y, z, ds)
log("Initializing mesh... done")

#=========================================
# Initialize the structured probes data
#=========================================

if comm.Get_rank() == 0:
    N = 30
    xbounds = [0.6, 4.0]
    ybounds = [-1.0, 1.0]
    z = 2.0

    x = np.linspace(xbounds[0], xbounds[1], N)
    y = np.linspace(ybounds[0], ybounds[1], N)
    X, Y = np.meshgrid(x,y)
    del x,y

    xyz = z*np.ones((N*N,3))
    xyz[:,0] = X.flatten()
    xyz[:,1] = Y.flatten()

else:
    xyz = None

log("Initializing probes...")
pb = Probes(comm, msh=msh, write_coords=False, probes = xyz)

plt.figure(figsize = (6,4))
plt.xlabel("x")
plt.ylabel("y")
plt.title("Waiting to receive data from neko...")
plt.savefig(join(output_path, "probes.png"), dpi = 150)
plt.tight_layout()
log("Initializing probes... done")

#=========================================
# Start streaming data
#=========================================

stream_data = True
j = 0
while stream_data:

    log("Waiting for data from Neko...")
    # Get the data
    u = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) 
    v = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) 
    w = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) 
    log("Data received")

    # Check if data was recieved or if the stream ended
    if ds.step_status == adios2.StepStatus.OK:
        stream_data = True
    elif ds.step_status == adios2.StepStatus.EndOfStream:
        stream_data = False
        continue

    # Interpolate and reshape into a structured grid for plotting
    pb.interpolate_from_field_list(
        j,
        [u,v,w], 
        comm, 
        write_data=False
        )

    if comm.Get_rank() == 0:
        u_probe = pb.interpolated_fields[:,1].reshape((N,N,1)).squeeze()
        v_probe = pb.interpolated_fields[:,2].reshape((N,N,1)).squeeze()
        w_probe = pb.interpolated_fields[:,3].reshape((N,N,1)).squeeze()
        mag = np.sqrt(u_probe**2 + v_probe**2 + w_probe**2)

        plt.title(f"Velocity magnitude (time step {j})")
        log("Plotting probes...")
        plt.contourf(X, Y, mag, levels = 40)
        log("Saving probes...")
        plt.savefig(join(output_path, "probes.png"), dpi = 150)
        log("Saving probes... done")

    j += 1


log("Detected stream ended")

# Finalize the streamer
log("Finalizing streamer...")
ds.finalize()
log("Done...")
