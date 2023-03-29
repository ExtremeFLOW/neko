#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# helloBPReaderHeatMap2D.py
#
#
#  Created on: Dec 5th, 2017
#      Author: William F Godoy godoywf@ornl.gov
#
import sys
#sys.path.append('/home/adalberto/3rd_party_libs_tar/adios2-installation/lib/python3/dist-packages/adios2/')
sys.path.append('/home/adalberto/NekDCTB/Nek5000/3rd_party/adios2/lib/python3/dist-packages/adios2')

import os
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6

from mpi4py import MPI
import numpy as np
import adios2
import scipy.optimize
import time
import matplotlib.pyplot as plt

#=============================Functions==============================

def distributed_svd(Xi,n,m):
    #Take the partiotioned data Xi in each rank and do the SVD.
    #The outputs are:
        # A partitioned mode matrix U
        # The eigen values D
        # The right orthogonal matrix trasposed Vt


    #Perfrom Svd in all ranks
    tic_in = time.perf_counter()
    Ui,Di,Vti=np.linalg.svd(Xi, full_matrices=False)
    toc_in = time.perf_counter()
    Yi=np.diag(Di)@Vti
    #print(f"Time for SVD of Xi in rank {rank}: {toc_in - tic_in:0.4f} seconds")

    #Gather Yi into Y in rank 0
    #prepare the buffer for recieving
    Y = None
    if rank == 0:
        #Generate the buffer to gather in rank 0
        Y = np.empty((m*size,m))
    comm.Gather(Yi, Y, root=0)

    if rank == 0:
        #If tank is zero, calculate the svd of the combined eigen matrix
        #Perform the svd of the combined eigen matrix
        tic_in = time.perf_counter()
        Uy,Dy,Vty=np.linalg.svd(Y, full_matrices=False)
        toc_in = time.perf_counter()
        #print(f"Time for SVD of Y in rank {rank}: {toc_in - tic_in:0.4f} seconds")
    else:
        #If the rank is not zero, simply create a buffer to recieve the Uy Dy and Vty
        Uy  = np.empty((m*size,m))
        Dy  = np.empty((m))
        Vty = np.empty((m,m))
    comm.Bcast(Uy, root=0)
    comm.Bcast(Dy, root=0)
    comm.Bcast(Vty, root=0)
    #Now matrix multiply each Ui by the corresponding entries in Uy
    U_local=Ui@Uy[rank*m:(rank+1)*m,:]

    return U_local, Dy, Vty

def dist_svd_update(U_1t,D_1t,Vt_1t,Xi,n,j,k): 

    if j==0:
        #Perform the distributed SVD and don't accumulate
        U_1t,D_1t,Vt_1t=distributed_svd(Xi,n,j+1)
    else:
        j1=j
        if j>=k:
            j1=k
        #Find the svd of the new snapshot
        U_tp1,D_tp1,Vt_tp1=distributed_svd(Xi,n,1)
        #2 contruct matrices to Do the updating
        V_tilde=scipy.linalg.block_diag(Vt_1t.T,Vt_tp1.T)
        W=np.append(U_1t@np.diag(D_1t),U_tp1@np.diag(D_tp1),axis=1)
        Uw,Dw,Vtw=distributed_svd(W,n,j1+1)
        #3 Update
        U_1t=Uw
        D_1t=Dw
        Vt_1t=(V_tilde@Vtw.T).T
 
        #Truncate the matrices if needed.
        #Eliminate the lower energy mode, which should be the last one
        if (j+1)>=k:
            U_1t=np.copy(U_1t[:,0:k])
            D_1t=np.copy(D_1t[0:k])
            Vt_1t=np.copy(Vt_1t[0:k,:])

    return U_1t,D_1t,Vt_1t 

def gathermodes(U_1t, n, m):
    U = None #prepare the buffer for recieving
    if rank == 0:
        #Generate the buffer to gather in rank 0
        U = np.empty((n,m))
    comm.Gather(U_1t, U, root=0)
    return U

def gathermodesAndmass(U_1t,BM1, n, m):
    U = None #prepare the buffer for recieving
    bm1sqrt = None #prepare the buffer for recieving
    if rank == 0:
        #Generate the buffer to gather in rank 0
        U = np.empty((n,m))
        bm1sqrt = np.empty((n,1))
    comm.Gather(U_1t, U, root=0)
    comm.Gather(BM1, bm1sqrt, root=0)
    return U, bm1sqrt

#==========================Main code==============================

# MPI - MPMD
worldcomm = MPI.COMM_WORLD
worldrank = worldcomm.Get_rank()
worldsize = worldcomm.Get_size()
col = 1
comm = worldcomm.Split(col,worldrank)
rank = comm.Get_rank()
size = comm.Get_size()

# ADIOS portion
adios = adios2.ADIOS(comm)

# ADIOS IO - Engine
ioRead = adios.DeclareIO("ioReader")
ioRead.SetEngine('InSituMPI')

# Create dummies for the streaming POD
U_1t = None
D_1t = None
Vt_1t = None

#Maximun number of modes to keep
k=100

# Open the insitu array... It is open until stream ends
ibpStream = ioRead.Open("globalArray", adios2.Mode.Read, comm)

# Keep looping waiting for data
while True:

    #=================================
    # Read step in each iteration
    #=================================

    stepStatus = ibpStream.BeginStep()

    if stepStatus == adios2.StepStatus.OK:

        currentStep = ibpStream.CurrentStep() 
        var_inVX = ioRead.InquireVariable("VX")
        var_inBM1 = ioRead.InquireVariable("BM1")
        var_inLGLEL = ioRead.InquireVariable("LGLEL")

        if var_inVX is not None:
            #print('THE SHAPE OF VX IS')
            #print(var_inVX.Shape()) 
            #Set up array sizes in the first step
            if currentStep == 0:
                total_size=var_inVX.Shape()[0]
                my_count= int(total_size/size)
                my_start= int(my_count*rank)
                if total_size % size != 0:
                    if rank < (total_size % size):
                        my_count += 1
                        my_start += rank
                    else:
                        my_start += (total_size%size)


            print(var_inLGLEL.Shape()) 
            #Set up array sizes in the first step
            if currentStep == 0:
                total_size2=var_inLGLEL.Shape()[0]
                my_count2= int(total_size2/size)
                my_start2= int(my_count2*rank)
                if total_size2 % size != 0:
                    if rank < (total_size2 % size):
                        my_count2 += 1
                        my_start2 += rank
                    else:
                        my_start2 += (total_size2%size)
                        
                #Allocate the vector
                inVX = np.zeros((my_count,1), dtype=np.double)
                inBM1 = np.zeros((my_count,1), dtype=np.double)
                inLGLEL = np.zeros((my_count2,1), dtype=np.int32)

            #Select the data in the global array that belongs to me
            var_inVX.SetSelection([[my_start], [my_count]])
            var_inBM1.SetSelection([[my_start], [my_count]])
            var_inLGLEL.SetSelection([[my_start2], [my_count2]])

            #Read the variable into array
            ibpStream.Get(var_inVX, inVX)
            ibpStream.Get(var_inBM1, inBM1)
            ibpStream.Get(var_inLGLEL, inLGLEL)

            ibpStream.EndStep()

    #=================================
    # Post process step in each iteration
    #=================================

            #print the data just to have it
            np.savetxt("elemstep"+repr(currentStep)+"rank"+repr(rank)+".txt", inVX)    
            print("it locked here 5 ") 
            #Scale the variables by their mass matrices
            inBM1sqrt=np.sqrt(inBM1)
            #print('mass matrix shape')
            #print(inBM1sqrt.shape)
            #print('mass matrix')
            #print(inBM1sqrt)
            for i in range(0,my_count):
                inVX[i,0]=inVX[i,0]*inBM1sqrt[i,0]

            #Configure the streaming POD process
            currentsnap=inVX
            n=total_size
            m=currentStep+1

            # Update the svd with each new snapshot
            U_1t,D_1t,Vt_1t = dist_svd_update(U_1t,D_1t,Vt_1t,currentsnap,n,currentStep,k)

            if rank == 0:
                print('Data in rank 0')
                print('Shape of snapshot:')
                print(currentsnap.shape)
                print('Shape of Mode matrix:')
                print(U_1t.shape)
                print('Shape of D:')
                print(D_1t.shape)
                print('Shape of Vt:')
                print(Vt_1t.shape)
                print('current step:')
                print(currentStep)
                print('current m:')
                print(m)
                print('k:')
                print(k)

    elif stepStatus == adios2.StepStatus.EndOfStream:
        print("End of stream")
        break

ibpStream.Close()


#Do this in case you want to keep more data than you actually collected
if k>=m:
    k=m

#Scale the modes back before gathering them
for i in range(0,k):
    for j in range(0,my_count):
        U_1t[j,i]= U_1t[j,i]/inBM1sqrt[j]

#Once the streaming is done. Gather the modes
#U = gathermodes(U_1t,n,m)
U, bm1sqrt = gathermodesAndmass(U_1t,inBM1sqrt,n,k)

#Show the sizes in rank 0
if rank==0:
    print("Shape after SVD")
    print(U.shape)
    print(D_1t.shape)
    print(Vt_1t.shape)

    utot=np.zeros((n,m))

    for i in range(0,m):
        u0=np.loadtxt("elemstep"+repr(i)+"rank"+repr(0)+".txt")    
        u1=np.loadtxt("elemstep"+repr(i)+"rank"+repr(1)+".txt")    
        u01=np.append(u0,u1,axis=0)
        utot[:,i]=u01

    print('shape of snapshot mat')
    print(utot.shape)   
    
    utot_w=np.copy(utot)

    a_file = open("VX.txt", "w")
    for row in utot_w:
        np.savetxt(a_file, row)
    a_file.close()

    for i in range(0,m):
        for j in range(0,n):
            utot_w[j,i]= utot_w[j,i]*bm1sqrt[j]


    a_file = open("scaledVX.txt", "w")
    for row in utot_w:
        np.savetxt(a_file, row)
    a_file.close()


#    print('Calculating the standar svd')
#    # Calculate also with standar svd... this might be slow.
#    Ur2,Dr2,Vtr2=np.linalg.svd(utot_w, full_matrices=False)
#
#    for i in range(0,m):
#        for j in range(0,n):
#            Ur2[j,i]= Ur2[j,i]/bm1sqrt[j]
#
    print('Evaluating results')
    #Evaluate
    urec=U@np.diag(D_1t)@Vt_1t
    #for i in range(0,m):
    #    for j in range(0,n):
    #        urec[j,i]= urec[j,i]/bm1sqrt[j]

#    urec2=Ur2@np.diag(Dr2)@Vtr2
#
#
    print('Writing results')
    a_file = open("U.txt", "w")
    for row in U:
        np.savetxt(a_file, row)
    a_file.close()

    a_file = open("Vt.txt", "w")
    for row in Vt_1t:
        np.savetxt(a_file, row)
    a_file.close()

    np.savetxt("D.txt", D_1t)

    np.savetxt("bm1sqrt.txt", bm1sqrt)

    #np.savetxt("LGLEL.txt", inLGLEL)

    
#    print('Plotting results')
#    #results from streaming
#    A=urec-utot
#    print(urec-utot)
#    plt.figure(figsize=(7,5))
#    ax = plt.subplot(1,1,1)
#    p = ax.pcolormesh(A)
#    plt.colorbar(p)
#    plt.show()
#
#    #results from standar
#    A2=urec2-utot
#    print(urec2-utot)
#    plt.figure(figsize=(7,5))
#    ax = plt.subplot(1,1,1)
#    p = ax.pcolormesh(A2)
#    plt.colorbar(p)
#    plt.show()
    
   
