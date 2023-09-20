import sys
import os
import copy
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6

from mpi4py import MPI #equivalent to the use of MPI_init() in C
import numpy as np
import scipy.optimize
import time
import pymech as pm
from pymech.neksuite import readnek
from pymech.neksuite import writenek
from tqdm import tqdm
import mpi_spSVD as spSVD
import pandas as pd
    
print("it imports")

# MPI - MPMD
worldcomm = MPI.COMM_WORLD
worldrank = worldcomm.Get_rank()
worldsize = worldcomm.Get_size()
col = 1
comm = worldcomm.Split(col,worldrank)
rank = comm.Get_rank()
size = comm.Get_size()

print("it imports")

# Import adios2
sys.path.append('/home/adalberto/NekDCTB/Nek5000/3rd_party/adios2/lib/python3.8/site-packages/adios2')
import adios2

if __name__ == "__main__":

    # MPI - MPMD
#    worldcomm = MPI.COMM_WORLD
#    worldrank = worldcomm.Get_rank()
#    worldsize = worldcomm.Get_size()
#    col = 1
#    comm = worldcomm.Split(col,worldrank)
#    rank = comm.Get_rank()
#    size = comm.Get_size()
    
    # ADIOS portion
    adios = adios2.ADIOS(comm)

