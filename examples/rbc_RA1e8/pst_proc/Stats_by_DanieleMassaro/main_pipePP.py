"""
Code to post-process statistics of 3D fully-developed turbulent pipe flow.
Adapted from prior Matlab scripts. Main contributors: Daniele Massaro & Saleh Rezaeiravesh. 
KTH Engineering Mechanics, 2020. {dmassaro,salehr}@kth.se
"""

import numpy as np
import os 
import sys
import matplotlib.pyplot as plt
sys.path.append('./')
from MODULES import reader
from MODULES import cartTerms
from MODULES import rotAvg
from MODULES import plotter


#0. Read inputs
Nx,Ny,nu,rho = reader.read_input('input.txt','main')

#1. Read in the int_fld
fname = 'int_fld' 
data = reader.read_int_fld(fname)

#2. Make a dictionary from data 
db = reader.dbMaker(data,Nx,Ny)

#3. Extract the interpolation grid
x,y = reader.makeGrid(data,Nx,Ny)

#4. Compute the terms in Cartersian CS
dbCart = cartTerms.comp(db,Nx,Ny,nu,rho)

#5. Do the rotation and take average over azimuthal direction
dbrotAvg = rotAvg.comp(dbCart,Nx,Ny,nu,rho)

#6. Plot results
plotter.vel(dbrotAvg,nu,rho)
plotter.kBudgets(dbrotAvg,nu,rho)
plotter.Re_stresstensor(dbrotAvg,nu,rho)
plotter.urur(dbrotAvg,nu,rho)
plotter.utut(dbrotAvg,nu,rho)
plotter.uzuz(dbrotAvg,nu,rho)
plotter.uruz(dbrotAvg,nu,rho)

