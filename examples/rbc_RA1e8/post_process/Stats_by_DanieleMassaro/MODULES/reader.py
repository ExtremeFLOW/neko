# Read the binary file int_fld (produced by Nek5000) and input.txt, which contains interpolating mesh parameters.

import struct
import numpy as np

class point:
    """class defining point variables"""

    def __init__(self,ldim,nstat,nderiv):
        self.pos = np.zeros((ldim))
        self.stat = np.zeros((nstat))
        self.deriv = np.zeros((nderiv))

class pset:
    """class containing data of the point collection"""

    def __init__(self,ldim,nstat,nderiv,npoints):
        self.ldim = ldim
        self.nstat = nstat
        self.nderiv = nderiv
        self.npoints = npoints
        self.re = []
        self.bsize = np.zeros((3))
        self.belem = []
        self.pord = []
        self.start_time = []
        self.end_time = []
        self.effav_time = []
        self.dt = []
        self.nrec = []
        self.int_time = []
        self.pset = [point(ldim,nstat,nderiv) for il in range(npoints)]

def read_int(infile,emode,nvar):
    """read integer array"""
    isize = 4
    llist = infile.read(isize*nvar)
    llist = list(struct.unpack(emode+nvar*'i', llist))
    return llist

def read_flt(infile,emode,wdsize,nvar):
    """read real array"""
    if (wdsize == 4):
        realtype = 'f'
    elif (wdsize == 8):
        realtype = 'd'
    llist = infile.read(wdsize*nvar)
    llist = np.frombuffer(llist, dtype=emode+realtype, count=nvar)
    return llist

def read_int_fld(fname):
    """read data from interpolation file"""
    # open file
    infile = open(fname, 'rb')
    # read header
    header = infile.read(460).split()

    # extract word size
    wdsize = int(header[1])

    # identify endian encoding
    etagb = infile.read(4)
    etagL = struct.unpack('<f', etagb)[0]; etagL = int(etagL*1e5)/1e5
    etagB = struct.unpack('>f', etagb)[0]; etagB = int(etagB*1e5)/1e5
    if (etagL == 6.54321):
        emode = '<'
    elif (etagB == 6.54321):
        emode = '>'

    # get simulation parameters
    ren = read_flt(infile,emode,wdsize,1)
    bsize = read_flt(infile,emode,wdsize,3)
    belem = read_int(infile,emode,3)
    pord = read_int(infile,emode,3)
    nstat = read_int(infile,emode,1)
    nderiv = read_int(infile,emode,1)
    stime = read_flt(infile,emode,wdsize,4)
    nrec = read_int(infile,emode,1)
    itime = read_flt(infile,emode,wdsize,1)
    npoints = read_int(infile,emode,1)

    # create main data structure
    ldim = 2
    data = pset(ldim,nstat[0],nderiv[0],npoints[0])

    # fill simulation parameters
    data.re = ren[0]
    data.bsize = bsize
    data.belem = belem
    data.pord = pord
    data.start_time = stime[0]
    data.end_time = stime[1]
    data.effav_time = stime[2]
    data.dt = stime[3]
    data.nrec = nrec[0]
    data.int_time = itime[0]

    # read coordinates
    for il in range(data.npoints):
        lpos = read_flt(infile,emode,wdsize,data.ldim)
        lptn = data.pset[il]
        data_pos = getattr(lptn,'pos')
        for jl in range(data.ldim):
            data_pos[jl] =  lpos[jl]

    # read statistics
    for il in range(data.nstat):
        for jl in range(data.npoints):
            lsts = read_flt(infile,emode,wdsize,1)
            lptn = data.pset[jl]
            data_sts = getattr(lptn,'stat')
            data_sts[il] =  lsts[0]

    # read derivatives
    for il in range(data.nderiv):
        for jl in range(data.npoints):
            ldrv = read_flt(infile,emode,wdsize,1)
            lptn = data.pset[jl]
            data_drv = getattr(lptn,'deriv')
            data_drv[il] =  ldrv[0]

    return data

def print_sim_data(data):
    """print simulation data"""
    print('Simulation data:')
    print('    Re = {0}'.format(data.re))
    print('    Lx = {0}'.format(data.bsize[0]))
    print('    Ly = {0}'.format(data.bsize[1]))
    print('    Lz = {0}'.format(data.bsize[2]))
    print('    nelx = {0}'.format(data.belem[0]))
    print('    nely = {0}'.format(data.belem[1]))
    print('    nelz = {0}'.format(data.belem[2]))
    print('    lx1 = {0}'.format(data.pord[0]))
    print('    ly1 = {0}'.format(data.pord[1]))
    print('    lz1 = {0}'.format(data.pord[2]))
    print('    nstat = {0}'.format(data.nstat))
    print('    nderiv = {0}'.format(data.nderiv))
    print('    start time = {0}'.format(data.start_time))
    print('    end time = {0}'.format(data.end_time))
    print('    effective average time = {0}'.format(data.effav_time))
    print('    time step = {0}'.format(data.dt))
    print('    nrec = {0}'.format(data.nrec))
    print('    time interval = {0}'.format(data.int_time))
    print('    npoints = {0}'.format(data.npoints))

def print_point_data(data,il):
    """print data related to single point"""
    print('Point data, npt = {0}'.format(il+1))
    lptn = data.pset[il]
    data_pos = getattr(lptn,'pos')
    print('x,y',data_pos)
    data_sts = getattr(lptn,'stat')
    print(data_sts)
    data_drv = getattr(lptn,'deriv')
    print(data_drv.size)

def makeGrid(data,Nx,Ny):
    """
    Make the interpolation grid from the position of the points
    """
    npts=data.npoints
    x=np.zeros(npts)
    y=np.zeros(npts)
    for i in range(npts):
        lptn = data.pset[i]
        data_pos = getattr(lptn,'pos')
        x[i]=data_pos[0]
        y[i]=data_pos[1]
    x=np.reshape(x,[Nx,Ny],'F')    
    y=np.reshape(y,[Nx,Ny],'F')    
    return x,y

def dbMaker(data,Nx,Ny):
    """
    Making database for fields 'F1', 'F2, ...,'D64'
    Each value is an array of shape (Nx,Ny) 
    """
    nStats=data.nstat
    nDerivs=data.nderiv
    nPts=data.npoints
    f=np.zeros((nPts,nStats))
    g=np.zeros((nPts,nDerivs))
    db={}

    for j in range(nPts):
        lptn = data.pset[j]
        f[j,:] = getattr(lptn,'stat')
        g[j,:] = getattr(lptn,'deriv')

    for i in range(nStats):
        fName='F'+str(i+1)
        f2=f[:,i]
        f2=np.reshape(f2,[Nx,Ny],'F')    
        db.update({fName:f2})

    for i in range(nDerivs):
        fName='D'+str(i+1)
        f2=g[:,i]
        f2=np.reshape(f2,[Nx,Ny],'F')    
        db.update({fName:f2})
    return db    

def read_input(filename,s):
   """read parameters of the interpolating mesh from input.txt""" 
   # 2D interpolating mesh parameters
   # Read param from .md file into a dict
   params = {}
   with open(filename) as fpar:
       for line in fpar:
          line = line.strip()
          key_val = line.split("= ")
          if len(key_val) == 2:
             params[key_val[0].strip()] = key_val[1]
   # Convert dict
   Nx  = int(params['Nx'])
   Ny  = int(params['Ny'])
   nu  = 1/float(params['1/nu'])
   rho = float(params['rho'])

   ReTau          = float(params['target_Retau'])
   Rmax           = float(params['Rmax'])
   drp1           = float(params['drp1'])
   compressedMesh = int(params['compressedMesh'])

   if (s=='main'):
      return Nx,Ny,nu,rho 
   elif (s=='pipeInterpMesh'):
      return Nx,Ny,ReTau,Rmax,drp1,compressedMesh
   elif (s=='rotAvg'):
      return ReTau,Rmax,drp1,compressedMesh

   return

