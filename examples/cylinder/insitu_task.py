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
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# external
import adios2.bindings as adios2

# pysemtools
from pysemtools.io.adios2.stream import DataStreamer
from pysemtools.io.utils import get_fld_from_ndarray
from pysemtools.datatypes.msh import Mesh
from pysemtools.interpolation.probes import Probes

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

def init_plot(save_output_path):
    fig, axs = plt.subplots(nrows = 3, figsize = (10,12), sharex = True)

    axs[-1].set_xlabel("x")
    for ax in axs:
        ax.set_aspect('equal')
        ax.set_ylabel("y")
    
    # Set colorbar for x-velocity
    norm = colors.Normalize(vmin=-0.4, vmax=1.4)
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='viridis'),
                         ax = axs[0])
    cbar.set_label("x-velocity")
    
    # Set colorbar for z-vorticity
    norm = colors.Normalize(vmin=-4.0, vmax=4.0)
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='bwr'),
                        ax = axs[1])
    cbar.set_label("z-vorticity")

    # Set colorbar for averaged x-velocity
    norm = colors.Normalize(vmin=-0.4, vmax=1.0)
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='viridis'),
                         ax = axs[2])
    cbar.set_label("avg. x-velocity")

    fig.suptitle("Waiting to receive data from neko...")
    fname = join(save_output_path, "cylinder_insitu_00000.png")
    fig.savefig(fname, dpi = 200)
    fname = "cylinder_insitu_latest.png"
    fig.savefig(fname, dpi = 100)
    log("Initializing probes... done")

    return fig, axs

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

    # Generate a 30x30 grid in the x-y plane at a given z.
    N = 30
    xbounds = [0.6, 4.0]
    ybounds = [-1.0, 1.0]
    z = 2.0

    x = np.linspace(xbounds[0], xbounds[1], N)
    y = np.linspace(ybounds[0], ybounds[1], N)
    X, Y = np.meshgrid(x,y)
    del x,y

    # This is the required format for pysemtools
    xyz = z*np.ones((N*N,3))
    xyz[:,0] = X.flatten()
    xyz[:,1] = Y.flatten()

else:
    xyz = None

log("Initializing probes...")
pb = Probes(comm, msh=msh, write_coords=False, probes = xyz)

#=========================================
# Initialize plots
#=========================================

if comm.Get_rank() == 0:
    fig, axs = init_plot(output_path)

#=========================================
# Start streaming data
#=========================================

stream_data = True
j = 0
while stream_data:

    log("Waiting for data from Neko...")
    # Get the data
    u = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) 
    curl_z = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) 
    avg_u = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype) 
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
        [u,curl_z,avg_u], 
        comm, 
        write_data=False
        )

    if comm.Get_rank() == 0:
        u_probe = pb.interpolated_fields[:,1].reshape((N,N,1)).squeeze()
        curl_probe = pb.interpolated_fields[:,2].reshape((N,N,1)).squeeze()
        avg_probe = pb.interpolated_fields[:,3].reshape((N,N,1)).squeeze()

        fig.suptitle(f"time step {j}")
        log("Plotting probes...")
        axs[0].contourf(X, Y, u_probe, levels = 40, norm=colors.Normalize(vmin=-0.4, vmax=1.4))
        axs[1].contourf(X, Y, curl_probe, levels = 40, cmap = "bwr", norm=colors.Normalize(vmin=-4.0, vmax=4.0))
        axs[2].contourf(X, Y, avg_probe, levels = 40, norm=colors.Normalize(vmin=-0.4, vmax=1.0))
        if j == 0: 
            fig.tight_layout()

        log("Saving probes...")
        fname = join(output_path, f"cylinder_insitu_{j:05d}.png")
        fig.savefig(fname, dpi = 200)
        fname = "cylinder_insitu_latest.png"
        fig.savefig(fname, dpi = 100)
        log("Saving probes... done")

    j += 1


log("Detected stream ended")

# Finalize the streamer
log("Finalizing streamer...")
ds.finalize()
log("Done...")
