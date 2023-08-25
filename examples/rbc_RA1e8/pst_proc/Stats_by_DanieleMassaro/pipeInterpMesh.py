#script to generate polar mesh over a circular disk (for pipe post-process)
# Kept similar to pymech
import struct
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import sys 
sys.path.append('./MODULES')
import reader

class point:
    """class defining point variables"""
    def __init__(self,ldim):
        self.pos = np.zeros((ldim))

class pset:
    """class containing data of the point collection"""
    def __init__(self,ldim,npoints):
        self.ldim = ldim
        self.npoints = npoints
        self.pset = [point(ldim) for il in range(npoints)]

def set_pnt_pos(data,il,lpos):
    """set position of the single point"""
    lptn = data.pset[il]
    data_pos = getattr(lptn,'pos')
    for jl in range(data.ldim):
            data_pos[jl] =  lpos[jl]

def write_int_pos(fname,wdsize,emode,data):
    """ write point positions to the file"""
    # open file
    outfile = open(fname, 'wb')
    # word size
    if (wdsize == 4):
        realtype = 'f'
    elif (wdsize == 8):
        realtype = 'd'
    # header
    header = '#iv1 %1i %1i %10i ' %(wdsize,data.ldim,data.npoints)
    header = header.ljust(32)
    outfile.write(header.encode('utf-8'))
    # write tag (to specify endianness)
    #etagb = struct.pack(emode+'f', 6.54321)
    #outfile.write(etagb)
    outfile.write(struct.pack(emode+'f', 6.54321))
    #write point positions
    for il in range(data.npoints):
        lptn = data.pset[il]
        data_pos = getattr(lptn,'pos')
        outfile.write(struct.pack(emode+data.ldim*realtype, *data_pos))

def disk(Rmax,compressedMesh,ReTau,drp1,nR,nTh):
    """ polar mesh over a circular disk """
    if (compressedMesh==int(1)):
       dr1=(drp1/ReTau)*Rmax  #distance from the wall of the first off-wall node
       gam=3.0;  #Grid compression controller >1
       xi=np.linspace(0.0,1.0,nR-1)
       R_=np.zeros(nR) 
       for i in range(0,nR-1):
          R_[i+1] = dr1+(Rmax-dr1)*(1.0+(np.tanh(gam*((xi[i])/(1)-1.0)))/np.tanh(gam))
       R_[0] = 0.0
#      Reverse the order of points: to be from center toward the wall
       R_ = Rmax-R_
       R = np.zeros(nR)
       for i in range(0,nR):
          R[i] = R_[nR-1-i] 
    else:
       R = np.linspace(0.0,Rmax,nR)     
    th = np.linspace(0,2.*mt.pi,nTh) 

    xx = np.zeros((nTh,nR))
    yy = np.zeros((nTh,nR))
    theta = np.zeros((nTh,nR))
    for i in range(0,nR):
       for j in range(0,nTh):
          xx[j,i] = R[i]*mt.cos(th[j])
          yy[j,i] = R[i]*mt.sin(th[j])
          theta[j,i] = th[j]
    

    return xx,yy,theta
    
if __name__ == "__main__":
    # initialise variables
    fname = 'int_pos'
    wdsize = 8
    # little endian
    emode = '<'
    # big endian
    #emode = '<'

#   Read interpolating mesh info
    Nx,Ny,ReTau,Rmax,drp1,compressedMesh = reader.read_input('input.txt','pipeInterpMesh')

    nR = Ny
    nTh = Nx


    #
    # generate polar mesh over a circular disk
    xx,yy,theta = disk(Rmax,compressedMesh,ReTau,drp1,nR,nTh)
    
    # allocate space
    ldim = 2   #2D grid
    npoints = xx.size
    data = pset(ldim,npoints)    
    print('Allocated {0} points'.format(npoints))

    # initialise point positions
    lpos = np.zeros(data.ldim)

    #set the positions in data
    il=0
    
    X = np.transpose(xx)
    Y = np.transpose(yy)

    x = X.flatten()
    y = Y.flatten()
     

    for il in range(0,nR*nTh):
       lpos[0]=x[il]
       lpos[1]=y[il]
       set_pnt_pos(data,il,lpos)
       il = il +1

    # write points to the file
    write_int_pos(fname,wdsize,emode,data)

    # plot the mesh
    plt.plot(x,y,'.b')
    plt.show()
#    

