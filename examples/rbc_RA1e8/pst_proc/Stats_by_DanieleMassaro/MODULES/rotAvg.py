# Rotation to cyl. coordinate and averaging over azimuthal direction \theta. 
# The main statistical quantities are also saved in .txt files.

import numpy as np
import math
import sys
sys.path.append('../')
from pipeInterpMesh import disk
sys.path.append('MODULES/')
import reader

def comp(dbCart,Nx,Ny,nu,rho):

    nR = Ny
    nTh = Nx

#   Read interpolating mesh info
    ReTau,Rmax,drp1,compressedMesh = reader.read_input('input.txt','rotAvg')

    xx,yy,theta = disk(Rmax,compressedMesh,ReTau,drp1,nR,nTh)
  

#   re-order x,y
    X = np.transpose(xx)
    Y = np.transpose(yy)

    k = 0
    x1 = np.zeros(nTh*nR)
    y1 = np.zeros(nTh*nR)
    angle = np.array([])
    for i in range(0,nR):
       a = theta[:,i]
       angle = np.concatenate((angle,a), axis=None)
 
    for j in range(0,nTh):
       for i in range(0,nR):
          x1[k] = X[i,j]
          y1[k] = Y[i,j]
          k = k +1


#   Transpose the imported data tensors
    U = np.transpose(dbCart['U']); V = np.transpose(dbCart['V']); W = np.transpose(dbCart['W'])
    uu = np.transpose(dbCart['uu']); uv = np.transpose(dbCart['uv']); uw = np.transpose(dbCart['uw']); vv = np.transpose(dbCart['vv']); vw = np.transpose(dbCart['vw'])
    ww = np.transpose(dbCart['ww']) 
    dUdx = np.transpose(dbCart['dUdx']); dUdy = np.transpose(dbCart['dUdy']); dUdz = np.transpose(dbCart['dUdz'])
    dVdx = np.transpose(dbCart['dVdx']); dVdy = np.transpose(dbCart['dVdy']); dVdz = np.transpose(dbCart['dVdz'])
    dWdx = np.transpose(dbCart['dWdx']); dWdy = np.transpose(dbCart['dWdy']); dWdz = np.transpose(dbCart['dWdz']) 
    Pxx = np.transpose(dbCart['Pxx']); Pxy = np.transpose(dbCart['Pxy']); Pxz = np.transpose(dbCart['Pxz']); Pyy = np.transpose(dbCart['Pyy']); Pyz = np.transpose(dbCart['Pyz'])
    Pzz = np.transpose(dbCart['Pzz'])
    Dxx = np.transpose(dbCart['Dxx']); Dxy = np.transpose(dbCart['Dxy']); Dxz = np.transpose(dbCart['Dxz']); Dyy = np.transpose(dbCart['Dyy']); Dyz = np.transpose(dbCart['Dyz'])
    Dzz = np.transpose(dbCart['Dzz'])
    Cxx = np.transpose(dbCart['Cxx']); Cxy = np.transpose(dbCart['Cxy']); Cxz = np.transpose(dbCart['Cxz']); Cyy = np.transpose(dbCart['Cyy']); Cyz = np.transpose(dbCart['Cyz'])
    Czz = np.transpose(dbCart['Czz'])
    Txx = np.transpose(dbCart['Txx']); Txy = np.transpose(dbCart['Txy']); Txz = np.transpose(dbCart['Txz']); Tyy = np.transpose(dbCart['Tyy']); Tyz = np.transpose(dbCart['Tyz'])
    Tzz = np.transpose(dbCart['Tzz'])
    VDxx = np.transpose(dbCart['VDxx']); VDxy = np.transpose(dbCart['VDxy']); VDxz = np.transpose(dbCart['VDxz']); VDyy = np.transpose(dbCart['VDyy']); VDyz = np.transpose(dbCart['VDyz'])
    VDzz = np.transpose(dbCart['VDzz'])
    PTxx = np.transpose(dbCart['PTxx']); PTxy = np.transpose(dbCart['PTxy']); PTxz = np.transpose(dbCart['PTxz']); PTyy = np.transpose(dbCart['PTyy']); PTyz = np.transpose(dbCart['PTyz'])
    PTzz = np.transpose(dbCart['PTzz'])
    PSxx = np.transpose(dbCart['PSxx']); PSxy = np.transpose(dbCart['PSxy']); PSxz = np.transpose(dbCart['PSxz']); PSyy = np.transpose(dbCart['PSyy']); PSyz = np.transpose(dbCart['PSyz'])
    PSzz = np.transpose(dbCart['PSzz'])
    uuu = np.transpose(dbCart['uuu']); uvv = np.transpose(dbCart['uvv']); uuw = np.transpose(dbCart['uuw']); uuv = np.transpose(dbCart['uuv']); uvw = np.transpose(dbCart['uvw'])
    uww = np.transpose(dbCart['uww']); vvv = np.transpose(dbCart['vvv']); vvw = np.transpose(dbCart['vvw']); vww = np.transpose(dbCart['vww']); www = np.transpose(dbCart['www'])
    skew_tensor = np.zeros((3,3,3))



#   Rotation to cylindrical coord.
    n2 = -1
    n3 = 0


#   Initialisation rotated quantities. The variables zz_ are zz after rotation
    r1 = np.zeros((nTh,nR))
    Ur = np.zeros((nTh,nR)); Ut = np.zeros((nTh,nR)); Uz = np.zeros((nTh,nR))
    urur = np.zeros((nTh,nR)); urut = np.zeros((nTh,nR)); uruz = np.zeros((nTh,nR)); utut = np.zeros((nTh,nR)); uzuz = np.zeros((nTh,nR)); utuz = np.zeros((nTh,nR))
    dUrdr = np.zeros((nTh,nR)); dUrdt = np.zeros((nTh,nR)); dUrdz = np.zeros((nTh,nR)); dUtdr = np.zeros((nTh,nR)); dUtdt = np.zeros((nTh,nR)); dUtdz = np.zeros((nTh,nR))
    dUzdr = np.zeros((nTh,nR)); dUzdt = np.zeros((nTh,nR)); dUzdz = np.zeros((nTh,nR))
    Prr = np.zeros((nTh,nR)); Ptt = np.zeros((nTh,nR)); Pzz_ = np.zeros((nTh,nR)); Prt = np.zeros((nTh,nR)); Prz = np.zeros((nTh,nR)); Ptz = np.zeros((nTh,nR))
    Drr = np.zeros((nTh,nR)); Dtt = np.zeros((nTh,nR)); Dzz_ = np.zeros((nTh,nR)); Drt = np.zeros((nTh,nR)); Drz = np.zeros((nTh,nR)); Dtz = np.zeros((nTh,nR))
    Crr = np.zeros((nTh,nR)); Ctt = np.zeros((nTh,nR)); Czz_ = np.zeros((nTh,nR)); Crt = np.zeros((nTh,nR)); Crz = np.zeros((nTh,nR)); Ctz = np.zeros((nTh,nR))
    Trr = np.zeros((nTh,nR)); Ttt = np.zeros((nTh,nR)); Tzz_ = np.zeros((nTh,nR)); Trt = np.zeros((nTh,nR)); Trz = np.zeros((nTh,nR)); Ttz = np.zeros((nTh,nR))
    VDrr = np.zeros((nTh,nR)); VDtt = np.zeros((nTh,nR)); VDzz_ = np.zeros((nTh,nR)); VDrt = np.zeros((nTh,nR)); VDrz = np.zeros((nTh,nR)); VDtz = np.zeros((nTh,nR))
    PTrr = np.zeros((nTh,nR)); PTtt = np.zeros((nTh,nR)); PTzz_ = np.zeros((nTh,nR)); PTrt = np.zeros((nTh,nR)); PTrz = np.zeros((nTh,nR)); PTtz = np.zeros((nTh,nR))
    PSrr = np.zeros((nTh,nR)); PStt = np.zeros((nTh,nR)); PSzz_ = np.zeros((nTh,nR)); PSrt = np.zeros((nTh,nR)); PSrz = np.zeros((nTh,nR)); PStz = np.zeros((nTh,nR))
    ururur = np.zeros((nTh,nR)); ututut = np.zeros((nTh,nR)); uzuzuz = np.zeros((nTh,nR)); ururut = np.zeros((nTh,nR)); ururuz = np.zeros((nTh,nR)); urutut = np.zeros((nTh,nR))
    ututuz = np.zeros((nTh,nR)); uruzuz = np.zeros((nTh,nR)); utuzuz = np.zeros((nTh,nR)); urutuz = np.zeros((nTh,nR))


    for i in range(0,nTh):        

       n1 = n2 +1 
       if (i==0 ):
          n2 = n1+nR-1
       else:         
          n2 = n1+nR-1
       
       r1[i,:] = np.sqrt(np.power(x1[n1:n2+1],2) + np.power(y1[n1:n2+1],2)) 

       # Rotation matrix
       R = np.matrix([[math.cos(angle[i]), math.sin(angle[i]), 0], [-math.sin(angle[i]), math.cos(angle[i]), 0], [0,             0,            1]])


       for jj in range(0,nR):

          # Mean velocities. Tensors of Rank 1.
          U_tens = np.matrix([[U[jj,i]], [V[jj,i]], [W[jj,i]]])
          prod = R*U_tens 
          Ur[i,jj] = prod[0]
          Ut[i,jj] = prod[1]
          Uz[i,jj] = prod[2]  

          # Reynolds stress tensor. Tensors of Rank 2.
          S = np.matrix([[uu[jj,i], uv[jj,i], uw[jj,i]], [uv[jj,i], vv[jj,i], vw[jj,i]], [uw[jj,i], vw[jj,i], ww[jj,i]]])          
          prod = R*S*np.transpose(R)
          urur[i,jj] = prod[0,0]
          urut[i,jj] = prod[0,1]
          uruz[i,jj] = prod[0,2]
          utut[i,jj] = prod[1,1]
          uzuz[i,jj] = prod[2,2]
          utuz[i,jj] = prod[1,2]

          # Velocity gradient tensor. Tensor of Rank 2.     
          S = np.matrix([[dUdx[jj,i], dUdy[jj,i], dUdz[jj,i]], [dVdx[jj,i], dVdy[jj,i], dVdz[jj,i]], [dWdx[jj,i], dWdy[jj,i], dWdz[jj,i]]])  
          prod = R*S*np.transpose(R)
          dUrdr[i,jj] = prod[0,0]
          dUrdt[i,jj] = prod[0,1]
          dUrdz[i,jj] = prod[0,2]
          dUtdr[i,jj] = prod[1,0]
          dUtdt[i,jj] = prod[1,1]
          dUtdz[i,jj] = prod[1,2]
          dUzdr[i,jj] = prod[2,0]
          dUzdt[i,jj] = prod[2,1]  
          dUzdz[i,jj] = prod[2,2]  
       
          # Production tensor. Tensor of Rank 2.
          S = np.matrix([[Pxx[jj,i], Pxy[jj,i], Pxz[jj,i]], [Pxy[jj,i], Pyy[jj,i], Pyz[jj,i]], [Pxz[jj,i], Pyz[jj,i], Pzz[jj,i]]])  
          prod = R*S*np.transpose(R)
          Prr[i,jj] = prod[0,0]
          Ptt[i,jj] = prod[1,1]
          Pzz_[i,jj] = prod[2,2]
          Prt[i,jj] = prod[0,1]
          Prz[i,jj] = prod[0,2]
          Ptz[i,jj] = prod[1,2]

          # Dissipation tensor. Tensor of Rank 2.
          S = np.matrix([[Dxx[jj,i], Dxy[jj,i], Dxz[jj,i]], [Dxy[jj,i], Dyy[jj,i], Dyz[jj,i]], [Dxz[jj,i], Dyz[jj,i], Dzz[jj,i]]]) 
          prod = R*S*np.transpose(R) 
          Drr[i,jj] = prod[0,0]
          Dtt[i,jj] = prod[1,1]
          Dzz_[i,jj] = prod[2,2]
          Drt[i,jj] = prod[0,1]
          Drz[i,jj] = prod[0,2]
          Dtz[i,jj] = prod[1,2]

          # Mean convection tensor. Tensor of Rank 2.
          S = np.matrix([[Cxx[jj,i], Cxy[jj,i], Cxz[jj,i]], [Cxy[jj,i], Cyy[jj,i], Cyz[jj,i]], [Cxz[jj,i], Cyz[jj,i], Czz[jj,i]]]) 
          prod = R*S*np.transpose(R) 
          Crr[i,jj] = prod[0,0]
          Ctt[i,jj] = prod[1,1]
          Czz_[i,jj] = prod[2,2]
          Crt[i,jj] = prod[0,1]
          Crz[i,jj] = prod[0,2]
          Ctz[i,jj] = prod[1,2]


          # Turbulent transport tensor. Tensor of Rank 2.
          S = np.matrix([[Txx[jj,i], Txy[jj,i], Txz[jj,i]], [Txy[jj,i], Tyy[jj,i], Tyz[jj,i]], [Txz[jj,i], Tyz[jj,i], Tzz[jj,i]]]) 
          prod = R*S*np.transpose(R) 
          Trr[i,jj] = prod[0,0]
          Ttt[i,jj] = prod[1,1]
          Tzz_[i,jj] = prod[2,2]
          Trt[i,jj] = prod[0,1]
          Trz[i,jj] = prod[0,2]
          Ttz[i,jj] = prod[1,2]

          # Viscous diffusion tensor. Tensor of Rank 2.
          S = np.matrix([[VDxx[jj,i], VDxy[jj,i], VDxz[jj,i]], [VDxy[jj,i], VDyy[jj,i], VDyz[jj,i]], [VDxz[jj,i], VDyz[jj,i], VDzz[jj,i]]]) 
          prod = R*S*np.transpose(R) 
          VDrr[i,jj] = prod[0,0]
          VDtt[i,jj] = prod[1,1]
          VDzz_[i,jj] = prod[2,2]
          VDrt[i,jj] = prod[0,1]
          VDrz[i,jj] = prod[0,2]
          VDtz[i,jj] = prod[1,2]

          # Pressure transport tensor. Tensor of Rank 2.   
          S = np.matrix([[PTxx[jj,i], PTxy[jj,i], PTxz[jj,i]], [PTxy[jj,i], PTyy[jj,i], PTyz[jj,i]], [PTxz[jj,i], PTyz[jj,i], PTzz[jj,i]]]) 
          prod = R*S*np.transpose(R) 
          PTrr[i,jj] = prod[0,0]
          PTtt[i,jj] = prod[1,1]
          PTzz_[i,jj] = prod[2,2]
          PTrt[i,jj] = prod[0,1]
          PTrz[i,jj] = prod[0,2]
          PTtz[i,jj] = prod[1,2]

          # Pressure strain tensor. Tensor of Rank 2.
          S = np.matrix([[PSxx[jj,i], PSxy[jj,i], PSxz[jj,i]], [PSxy[jj,i], PSyy[jj,i], PSyz[jj,i]], [PSxz[jj,i], PSyz[jj,i], PSzz[jj,i]]]) 
          prod = R*S*np.transpose(R) 
          PSrr[i,jj] = prod[0,0]
          PStt[i,jj] = prod[1,1]
          PSzz_[i,jj] = prod[2,2]
          PSrt[i,jj] = prod[0,1]
          PSrz[i,jj] = prod[0,2]
          PStz[i,jj] = prod[1,2]

          # Skewness tensor. Tensor of Rank 3.  
          skew_tensor[:,:,0] = np.matrix([[uuu[jj,i], uvv[jj,i], uuw[jj,i]], [uuv[jj,i], uvv[jj,i], uvw[jj,i]], [uuw[jj,i], uvw[jj,i], uww[jj,i]]]) 
          skew_tensor[:,:,1] = np.matrix([[uuv[jj,i], uvv[jj,i], uvw[jj,i]], [uvv[jj,i], vvv[jj,i], vvw[jj,i]], [uvw[jj,i], vvw[jj,i], vww[jj,i]]]) 
          skew_tensor[:,:,2] = np.matrix([[uuw[jj,i], uvw[jj,i], uww[jj,i]], [uvw[jj,i], vvw[jj,i], vww[jj,i]], [uww[jj,i], vww[jj,i], www[jj,i]]]) 

          aabc = np.zeros((3,3,3))
          adef = skew_tensor[:,:,:]
          for aa in range(0,3):
             for bb in range(0,3):
                for cc in range(0,3):
                   for dd in range(0,3):
                      for ee in range(0,3):
                         for ff in range(0,3):
                            aabc[aa,bb,cc] = aabc[aa,bb,cc]+R[aa,dd]*R[bb,ee]*R[cc,ff]*adef[dd,ee,ff]

          ururur[i,jj] = aabc[0,0,0]
          ututut[i,jj] = aabc[1,1,1]
          uzuzuz[i,jj] = aabc[2,2,2]
          ururut[i,jj] = aabc[0,1,0]
          ururuz[i,jj] = aabc[0,2,0]
          urutut[i,jj] = aabc[1,1,0]
          ututuz[i,jj] = aabc[1,2,1]
          uruzuz[i,jj] = aabc[2,2,0]
          utuzuz[i,jj] = aabc[2,2,1]
          urutuz[i,jj] = aabc[1,2,0] 


#   Average over the azimuthal direction
#   Initialisation averaged quantities.
    Ur_m = np.zeros(nR); Ut_m = np.zeros(nR); Uz_m = np.zeros(nR)
    urur_m = np.zeros(nR); urut_m = np.zeros(nR); uruz_m = np.zeros(nR); utut_m = np.zeros(nR); uzuz_m = np.zeros(nR); utuz_m = np.zeros(nR)
    dUrdr_m = np.zeros(nR); dUrdt_m = np.zeros(nR); dUrdz_m = np.zeros(nR); dUtdr_m = np.zeros(nR); dUtdt_m = np.zeros(nR); dUtdz_m = np.zeros(nR)
    dUzdr_m = np.zeros(nR); dUzdt_m = np.zeros(nR); dUzdz_m = np.zeros(nR)
    Prr_m = np.zeros(nR); Ptt_m = np.zeros(nR); Pzz_m = np.zeros(nR); Prt_m = np.zeros(nR); Prz_m = np.zeros(nR); Ptz_m = np.zeros(nR)
    Drr_m = np.zeros(nR); Dtt_m = np.zeros(nR); Dzz_m = np.zeros(nR); Drt_m = np.zeros(nR); Drz_m = np.zeros(nR); Dtz_m = np.zeros(nR)
    Crr_m = np.zeros(nR); Ctt_m = np.zeros(nR); Czz_m = np.zeros(nR); Crt_m = np.zeros(nR); Crz_m = np.zeros(nR); Ctz_m = np.zeros(nR)
    Trr_m = np.zeros(nR); Ttt_m = np.zeros(nR); Tzz_m = np.zeros(nR); Trt_m = np.zeros(nR); Trz_m = np.zeros(nR); Ttz_m = np.zeros(nR)
    VDrr_m = np.zeros(nR); VDtt_m = np.zeros(nR); VDzz_m = np.zeros(nR); VDrt_m = np.zeros(nR); VDrz_m = np.zeros(nR); VDtz_m = np.zeros(nR)
    PTrr_m = np.zeros(nR); PTtt_m = np.zeros(nR); PTzz_m = np.zeros(nR); PTrt_m = np.zeros(nR); PTrz_m = np.zeros(nR); PTtz_m = np.zeros(nR)
    PSrr_m = np.zeros(nR); PStt_m = np.zeros(nR); PSzz_m = np.zeros(nR); PSrt_m = np.zeros(nR); PSrz_m = np.zeros(nR); PStz_m = np.zeros(nR)
    ururur_m = np.zeros(nR); ututut_m = np.zeros(nR); uzuzuz_m = np.zeros(nR); ururut_m = np.zeros(nR); ururuz_m = np.zeros(nR); urutut_m = np.zeros(nR)
    ututuz_m = np.zeros(nR); uruzuz_m = np.zeros(nR); utuzuz_m = np.zeros(nR); urutuz_m = np.zeros(nR)


    for j in range(0,nR):
          
       # Mean velocities. Tensors of Rank 1.        
       Ur_m[j] = np.mean(Ur[:,j])
       Ut_m[j] = np.mean(Ut[:,j])       
       Uz_m[j] = np.mean(Uz[:,j])

       # Reynolds stress tensor. Tensors of Rank 2.
       urur_m[j] = np.mean(urur[:,j])
       urut_m[j] = np.mean(urut[:,j])
       uruz_m[j] = np.mean(uruz[:,j])
       utut_m[j] = np.mean(utut[:,j])
       uzuz_m[j] = np.mean(uzuz[:,j])
       utuz_m[j] = np.mean(utuz[:,j])

       # Velocity gradient tensor. Tensor of Rank 2.
       dUrdr_m[j] = np.mean(dUrdr[:,j])
       dUrdt_m[j] = np.mean(dUrdt[:,j])
       dUrdz_m[j] = np.mean(dUrdz[:,j])
       dUtdr_m[j] = np.mean(dUtdr[:,j])
       dUtdt_m[j] = np.mean(dUtdt[:,j])
       dUtdz_m[j] = np.mean(dUtdz[:,j])
       dUzdr_m[j] = np.mean(dUzdr[:,j])
       dUzdt_m[j] = np.mean(dUzdt[:,j]) 
       dUzdz_m[j] = np.mean(dUzdz[:,j])  

         
       # Production tensor. Tensor of Rank 2.     
       Prr_m[j] = np.mean(Prr[:,j])
       Ptt_m[j] = np.mean(Ptt[:,j])
       Pzz_m[j] = np.mean(Pzz_[:,j])
       Prt_m[j] = np.mean(Prt[:,j])
       Prz_m[j] = np.mean(Prz[:,j])
       Ptz_m[j] = np.mean(Ptz[:,j])

       # Dissipation tensor. Tensor of Rank 2.
       Drr_m[j] = np.mean(Drr[:,j])
       Dtt_m[j] = np.mean(Dtt[:,j])
       Dzz_m[j] = np.mean(Dzz_[:,j])
       Drt_m[j] = np.mean(Drt[:,j])
       Drz_m[j] = np.mean(Drz[:,j])
       Dtz_m[j] = np.mean(Dtz[:,j])
  
       # Mean convection tensor. Tensor of Rank 2.
       Crr_m[j] = np.mean(Crr[:,j])
       Ctt_m[j] = np.mean(Ctt[:,j])
       Czz_m[j] = np.mean(Czz_[:,j])
       Crt_m[j] = np.mean(Crt[:,j])
       Crz_m[j] = np.mean(Crz[:,j])
       Ctz_m[j] = np.mean(Ctz[:,j])
  
       # Turbulent transport tensor. Tensor of Rank 2.
       Trr_m[j] = np.mean(Trr[:,j])
       Ttt_m[j] = np.mean(Ttt[:,j])
       Tzz_m[j] = np.mean(Tzz_[:,j])
       Trt_m[j] = np.mean(Trt[:,j])
       Trz_m[j] = np.mean(Trz[:,j])
       Ttz_m[j] = np.mean(Ttz[:,j])

       # Viscous diffusion tensor. Tensor of Rank 2.
       VDrr_m[j] = np.mean(VDrr[:,j])
       VDtt_m[j] = np.mean(VDtt[:,j])
       VDzz_m[j] = np.mean(VDzz_[:,j])
       VDrt_m[j] = np.mean(VDrt[:,j])
       VDrz_m[j] = np.mean(VDrz[:,j])
       VDtz_m[j] = np.mean(VDtz[:,j])
 
       # Pressure transport tensor. Tensor of Rank 2.   
       PTrr_m[j] = np.mean(PTrr[:,j])
       PTtt_m[j] = np.mean(PTtt[:,j])
       PTzz_m[j] = np.mean(PTzz_[:,j])
       PTrt_m[j] = np.mean(PTrt[:,j])
       PTrz_m[j] = np.mean(PTrz[:,j])
       PTtz_m[j] = np.mean(PTtz[:,j])

       # Pressure strain tensor. Tensor of Rank 2.
       PSrr_m[j] = np.mean(PSrr[:,j])
       PStt_m[j] = np.mean(PStt[:,j])
       PSzz_m[j] = np.mean(PSzz_[:,j])
       PSrt_m[j] = np.mean(PSrt[:,j])
       PSrz_m[j] = np.mean(PSrz[:,j])
       PStz_m[j] = np.mean(PStz[:,j])


       # Skewness tensor. Tensor of Rank 3. 
       ururur_m[j] = np.mean(ururur[:,j])
       ututut_m[j] = np.mean(ututut[:,j])
       uzuzuz_m[j] = np.mean(uzuzuz[:,j])
       ururut_m[j] = np.mean(ururut[:,j])
       ururuz_m[j] = np.mean(ururuz[:,j])
       urutut_m[j] = np.mean(urutut[:,j])
       ututuz_m[j] = np.mean(ututuz[:,j])
       uruzuz_m[j] = np.mean(uruzuz[:,j])
       utuzuz_m[j] = np.mean(utuzuz[:,j])

    # Compute TKE budget terms
    P_k  = 0.5*(Prr_m  + Ptt_m  + Pzz_m)   # production
    T_k  = 0.5*(Trr_m  + Ttt_m  + Tzz_m)   # turbulent transport
    PS_k = 0.5*(PSrr_m + PStt_m + PSzz_m)  # pressure-strain
    PT_k = 0.5*(PTrr_m + PTtt_m + PTzz_m)  # pressure-transport
    Pi_k =     (PT_k - PS_k)
    VD_k = 0.5*(VDrr_m + VDtt_m + VDzz_m)  # viscous diffusion
    D_k  = 0.5*(Drr_m  + Dtt_m  + Dzz_m)   # dissipation
    C_k  = 0.5*(Crr_m  + Ctt_m  + Czz_m)   # mean convection

    # Write postprocessed results in a file
    # make an array of radii
    r_ = np.zeros(nR)
    for i in range(0,nR):
       r_[i] = np.sqrt(np.power(X[i,0],2) + np.power(Y[i,0],2))

    rMax = np.max(r_)
    r_ = rMax - r_    #Changing the discrete radius from wall toward center

    # Compute uTau
    uTau = np.sqrt(-dUzdr_m[-1]*nu)
    # method 2
    tauW = (Uz_m[-2] - Uz_m[-1])/(r_[-2]-r_[-1])
    uTau2 = np.sqrt(nu*tauW)

    ReTau = uTau*rMax/nu
    ReTau2 = uTau2*rMax/nu
    print('\n')
    print('To double check results u_tau and Re_tau are computed with two different methods. ')
    print('Method 1: as square root of the product between the averaged field dUzdr at the wall and nu. Retau is obtained from the uTau*rMax/nu.')
    print('Method 2: as square root of nu*tauW, where tauW is computed as finite difference of the averaged field Uz at the wall. Retau is obtained from the uTau*rMax/nu. \n')
    print('uTau (method1)= ',uTau)
    print('uTau (method2)= ',uTau2)
    print('Re_tau (method1)= ',ReTau)
    print('Re_tau (method2)= ',ReTau2)



#   Saving Dict. of rotated and averaged quantities(not scaled) 
    dbRotAvg = {'Ur':Ur_m, 'Ut':Ut_m, 'Uz':Uz_m, \
                'urur':urur_m, 'urut':urut_m, 'uruz':uruz_m, 'utut':utut_m, 'uzuz':uzuz_m, 'utuz':utuz_m, \
                'dUrdr':dUrdr_m, 'dUrdt':dUrdt_m, 'dUrdz':dUrdz_m, 'dUtdr':dUtdr_m, 'dUtdt':dUtdt_m, 'dUtdz':dUtdz_m, 'dUzdr':dUzdr_m, 'dUzdt':dUzdt_m, 'dUzdz_m':dUzdz, \
                'Prr':Prr_m, 'Ptt':Ptt_m, 'Pzz':Pzz_m, 'Prt':Prt_m, 'Prz':Prz_m, 'Ptz':Ptz_m, \
                'Drr':Drr_m, 'Dtt':Dtt_m, 'Dzz':Dzz_m, 'Drt':Drt_m, 'Drz':Drz_m, 'Dtz':Dtz_m, \
                'Crr':Crr_m, 'Ctt':Ctt_m, 'Czz':Czz_m, 'Crt':Crt_m, 'Crz':Crz_m, 'Ctz':Ctz_m, \
                'Trr':Trr_m, 'Ttt':Ttt_m, 'Tzz':Tzz_m, 'Trt':Trt_m, 'Trz':Trz_m, 'Ttz':Ttz_m, \
                'VDrr':VDrr_m, 'VDtt':VDtt_m, 'VDzz':VDzz_m, 'VDrt':VDrt_m, 'VDrz':VDrz_m, 'VDtz':VDtz_m, \
                'PTrr':PTrr_m, 'PTtt':PTtt_m, 'PTzz':PTzz_m, 'PTrt':PTrt_m, 'PTrz':PTrz_m, 'PTtz':PTtz_m, \
                'PSrr':PSrr_m, 'PStt':PStt_m, 'PSzz':PSzz_m, 'PSrt':PSrt_m, 'PSrz':PSrz_m, 'PStz':PStz_m, \
                'ururur':ururur_m, 'ututut':ututut_m, 'uzuzuz':uzuzuz_m, 'ururut':ururut_m, 'ururuz':ururuz_m, 'urutut':urutut_m, 'ututuz':ututuz_m, 'uruzuz':uruzuz_m, 'utuzuz':utuzuz_m, \
                'P_k':P_k, 'T_k':T_k, 'PS_k':PS_k, 'PT_k':PT_k, 'Pi_k':Pi_k, 'VD_k':VD_k, 'D_k':D_k, 'C_k':C_k, \
                'r':r_, 'nu':nu, 'u_tau': uTau, 'Re_tau':ReTau \
                }

#   Saving txt file of rotated and averaged quantities

    #velocity profiles + derivatives
    fileName = 'OUTPUT/turbPipe_meanVelDer.txt'
    F = open(fileName,'w') 
    F.write("# Postprocessed results of turbulent pipe flow \n") 
    F.write("# Mean velocity and their derivatives \n") 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# uTau = %g \n" % uTau) 
    F.write("# nu = %g \n"   % nu) 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# r\t Ur\t Ut\t Uz\t dUrdr\t dUrdt\t dUrdz\t dUtdr\t dUtdt\t dUtdz\t dUzdr\t dUzdt\t dUzdz \n")            
    F.write("# ------------------------------------------------------------------\n") 
    for i in range(0,nR):
       F.write("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t   \n" % \
               (r_[i], Ur_m[i], Ut_m[i], Uz_m[i], dUrdr_m[i], dUrdt_m[i], dUrdz_m[i], dUtdr_m[i], dUtdt_m[i], dUtdz_m[i], dUzdr_m[i], dUzdt_m[i], dUzdz_m[i]))
    F.close()


    #Reynolds stress components
    fileName = 'OUTPUT/turbPipe_ReStressTens.txt'
    F = open(fileName,'w') 
    F.write("# Postprocessed results of turbulent pipe flow \n") 
    F.write("# Reynolds stress components \n") 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# uTau = %g \n" % uTau) 
    F.write("# nu = %g \n"   % nu) 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# r\t urp\t urut\t uruz\t utp\t uzp\t utuz \n") 
    F.write("# ------------------------------------------------------------------\n") 
    for i in range(0,nR):
       F.write("%g\t%g\t%g\t%g\t %g\t%g\t%g\t \n" % \
               (r_[i], np.sqrt(urur_m[i]-Ur_m[i]*Ur_m[i]), urut_m[i]-Ur_m[i]*Ut_m[i], uruz_m[i]-Ur_m[i]*Uz_m[i], np.sqrt(utut_m[i]-Ut_m[i]*Ut_m[i]), \
               np.sqrt(abs(uzuz_m[i]-Uz_m[i]*Uz_m[i])), utuz_m[i]-Ut_m[i]*Uz_m[i]) )
    F.close()

    #TKE budget terms
    fileName = 'OUTPUT/turbPipe_kBudgets.txt'
    F = open(fileName,'w') 
    F.write("# Postprocessed results of turbulent pipe flow \n") 
    F.write("# TKE budget terms: nondim by multiplying by nu/uTau^4 \n") 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# uTau = %g \n" % uTau) 
    F.write("# nu = %g \n"   % nu) 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# r+\t P_k+ \t T_k+ \t PS_k+ \t PT_k+ \t VD_k+ \t D_k+ \t C_k+ \n") 
    F.write("# ------------------------------------------------------------------\n") 
    fac = nu/np.power(uTau,4)
    for i in range(0,nR):
       F.write("%g\t%g\t%g\t%g\t %g\t%g\t%g\t%g\t \n" % \
               (r_[i]*ReTau, P_k[i]*fac, T_k[i]*fac, PS_k[i]*fac, PT_k[i]*fac, VD_k[i]*fac, D_k[i]*fac, C_k[i]*fac) )
    F.close()

    #budget terms of <urur>
    fileName = 'OUTPUT/turbPipe_rrBudgets.txt'
    F = open(fileName,'w') 
    F.write("# Postprocessed results of turbulent pipe flow \n") 
    F.write("# budget terms of <urur>: nondim by multiplying by nu/uTau^4 \n") 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# uTau = %g \n" % uTau) 
    F.write("# nu = %g \n"   % nu) 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# r+\t P_rr+ \t T_rr+ \t PS_rr+ \t PT_rr+ \t VD_rr+ \t D_rr+ \t C_rr+ \n") 
    F.write("# ------------------------------------------------------------------\n") 
    for i in range(0,nR):
       F.write("%g\t%g\t%g\t%g\t %g\t%g\t%g\t%g \n" % \
               (r_[i]*ReTau, Prr_m[i]*fac, Trr_m[i]*fac, PSrr_m[i]*fac, PTrr_m[i]*fac, VDrr_m[i]*fac, Drr_m[i]*fac, Crr_m[i]*fac) )
    F.close()

    #budget terms of <utut>
    fileName = 'OUTPUT/turbPipe_ttBudgets.txt'
    F = open(fileName,'w') 
    F.write("# Postprocessed results of turbulent pipe flow \n") 
    F.write("# budget terms of <utut>: nondim by multiplying by nu/uTau^4 \n") 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# uTau = %g \n" % uTau) 
    F.write("# nu = %g \n"   % nu) 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# r+\t P_tt+ \t T_tt+ \t PS_tt+ \t PT_tt+ \t VD_tt+ \t D_tt+ \t C_tt+ \n") 
    F.write("# ------------------------------------------------------------------\n") 
    for i in range(0,nR):
       F.write("%g\t%g\t%g\t%g\t %g\t%g\t%g\t%g \n" % \
               (r_[i]*ReTau, Ptt_m[i]*fac, Ttt_m[i]*fac, PStt_m[i]*fac, PTtt_m[i]*fac, VDtt_m[i]*fac, Dtt_m[i]*fac, Ctt_m[i]*fac) )
    F.close()

    #budget terms of <uzuz>
    fileName = 'OUTPUT/turbPipe_zzBudgets.txt'
    F = open(fileName,'w') 
    F.write("# Postprocessed results of turbulent pipe flow \n") 
    F.write("# budget terms of <uzuz>: nondim by multiplying by nu/uTau^4 \n") 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# uTau = %g \n" % uTau) 
    F.write("# nu = %g \n"   % nu) 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# r+\t P_zz+ \t T_zz+ \t PS_zz+ \t PT_zz+ \t VD_zz+ \t D_zz+ \t C_zz+ \n") 
    F.write("# ------------------------------------------------------------------\n") 
    for i in range(0,nR):
       F.write("%g\t%g\t%g\t%g\t %g\t%g\t%g\t%g \n" % \
               (r_[i]*ReTau, Pzz_m[i]*fac, Tzz_m[i]*fac, PSzz_m[i]*fac, PTzz_m[i]*fac, VDzz_m[i]*fac, Dzz_m[i]*fac, Czz_m[i]*fac) )
    F.close()

    #budget terms of <uruz>
    fileName = 'OUTPUT/turbPipe_rzBudgets.txt'
    F = open(fileName,'w') 
    F.write("# Postprocessed results of turbulent pipe flow \n") 
    F.write("# budget terms of <uruz>: nondim by multiplying by nu/uTau^4 \n") 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# uTau = %g \n" % uTau) 
    F.write("# nu = %g \n"   % nu) 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# r+\t P_rz+ \t T_rz+ \t PS_rz+ \t PT_rz+ \t VD_rz+ \t D_rz+ \t C_rz+ \n") 
    F.write("# ------------------------------------------------------------------\n") 
    for i in range(0,nR):
       F.write("%g\t%g\t%g\t%g\t %g\t%g\t%g\t%g \n" % \
               (r_[i]*ReTau, Prz_m[i]*fac, Trz_m[i]*fac, PSrz_m[i]*fac, PTrz_m[i]*fac, VDrz_m[i]*fac, Drz_m[i]*fac, Crz_m[i]*fac) )
    F.close()


    return dbRotAvg
