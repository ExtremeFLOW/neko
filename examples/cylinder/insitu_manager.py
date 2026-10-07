#========================================
# Import and set up general modules
#========================================
import sys
import os
import argparse
from time import sleep

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

if rank == 0:
    print(f"Python is running with {size} ranks")

#========================================
# Import modules
#========================================
# general functionality
import numpy as np
from os.path import join
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# my classes
from pysemtools_plugins.EnhancedStreamer import EnhancedStreamer

# pysemtools
from pysemtools.monitoring.logger import Logger
log = Logger(comm=comm, module_name="cylinder_insitu_task")

import json

#=========================================
# Define some helper functions
#=========================================

def get_output_directory(fname):
    with open(fname, 'r') as f:
        data = json.load(f)

    return data["case"]["output_directory"]

def get_field_names(fname):
    with open(fname, 'r') as f:
        data = json.load(f)

    for component in data["case"]["simulation_components"]:
        if component.get("type") == "data_streamer":
            return component["fields"]
    raise ValueError("No data_streamer component found")

def init_plot(save_output_path):
    fig, axs = plt.subplots(nrows = 2, figsize = (15,12), sharex = True)

    axs[0].set_ylabel("Pressure [Pa]")
    axs[1].set_ylabel("Velocity magnitude [m/s]")
    axs[1].set_xlabel("x [m]")

    fig.suptitle("Waiting to receive data from neko...")
    fname = join(save_output_path, "cylinder_insitu_00000.png")
    fig.savefig(fname, dpi = 200)

    return fig, axs

#=========================================
# Parse arguments
#=========================================

parser = argparse.ArgumentParser(description="In-situ visualization task for the cylinder case.")
parser.add_argument("--dry-run", action="store_true",
                    help="Import modules and exit without executing any actions.")
parser.add_argument("--timeout", type=int, default=300,
                    help="Timeout in seconds for the ADIOS2 SST stream open (default: 300).")
parser.add_argument("--dtype", choices=["single", "double"], default="double",
                    help="Data type of the streamed fields (default: double).")
args = parser.parse_args()

if args.dry_run:
    log.write("info", "Dry run: all imports succeeded. Exiting.")
    sys.exit(0)

# Remove globalArray_* files
import subprocess
subprocess.run("rm -f globalArray*", shell=True)

#=========================================
# Define some variables/parameters
#=========================================

log.write("info", "Starting insitu task")

dtype_string = args.dtype
backend = "numpy"
if dtype_string == "single":
    dtype = np.float32
else:
    dtype = np.float64

output_path = get_output_directory("cylinder_insitu.case")
log.write("info", f"Outputting insitu snapshots to folder {output_path}")

#=========================================
# Initialize the streamer
#=========================================

streamer_field_names = get_field_names("cylinder_insitu.case")

processor = EnhancedStreamer(
    comm,
    fields = streamer_field_names,
    adios2_timeout = args.timeout,
    create_catalyst_session=True,
    catalyst_pipeline="pipeline_lambda2.py",
    catalyst_channel="field.vtkhdf"
    )

#=========================================
# Stream the mesh coordinates to build interpolators
#=========================================
processor.receive_mesh(dtype_string)

#=========================================
# Define interpolators for the fields
#=========================================

# processor.add_interpolator_from_values(
#     "line",
#     x = np.linspace(-1.0, 3.0, 100),
#     y = 5.05,
#     z = 0.0
#     )

# =========================================
# Initialize the plot
# =========================================

fig, axs = init_plot(output_path)

#=========================================
# Start streaming data
#=========================================

j = 0
while True:

    processor.receive_fields(dtype_string)

    mag = np.sqrt(processor.fields["u"]**2 + \
                  processor.fields["v"]**2 + \
                    processor.fields["w"]**2)

    processor.add_field("mag", mag)

    # Check if data was recieved or if the stream ended
    if not processor.get_adios2_status():
        break

    processor.execute_catalyst_session(j, 0.1*j)
    # processor.update_interpolators(j)

    # p = processor.get_field_from_interpolator("p", "line")
    # m = processor.get_field_from_interpolator("mag", "line")

    # axs[0].plot(processor.interpolators["line"].x, p)
    # axs[1].plot(processor.interpolators["line"].x, m)
    # fig.suptitle(f"Time: {0.1*j:.2f} s")
    # fname = join(output_path, f"cylinder_insitu_00000.png")
    # fig.savefig(fname, dpi=200)
    # axs[0].clear()
    # axs[1].clear()

    j += 1

log.write("info", "Detected stream ended")

processor.finalize()
