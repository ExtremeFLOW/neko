# Extract statistical fields over the Cartesian interpolating mesh. Look at nom_fields.txt to see fields nomenclature.
import numpy as np


def comp(db,Nx,Ny,nu,rho):

#   Mean velocities. Tensors of Rank 1.
    U = db['F1']
    V = db['F2']
    W = db['F3']


#   Reynolds stress tensor. Tensors of Rank 2.
    uu = db['F5']
    vv = db['F6']
    ww = db['F7']
    uv = db['F9']
    uw = db['F11']
    vw = db['F10']


#   Mean, RMS, skewness and flatness of pressure
    P = db['F4']
    pp = db['F8'] - P*P
    ppp = db['F27'] - 3*P*pp - P*P*P
    pppp = db['F38'] - 4*P*ppp - 6*P*P*pp - P*P*P*P
#   Normalize pressure
    prms = np.power(pp,0.5)
    pskew = np.divide(ppp,np.power(pp,1.5))
    pflat = np.divide(pppp,np.power(pp,2))


#   Skewness tensor. Tensor of Rank 3.
#   Form of the tensor.
#   [ uuu uuv uuw ] [ uuv uvv uvw ] [ uuw uvw uww ] 
#   [ uuv uvv uvw ] [ uvv vvv vvw ] [ uvw vvw vww ]
#   [ uuw uvw uww ] [ uvw vvw vww ] [ uww vww www ]
    uuu = db['F24'] - 3*U*uu - U*U*U 
    vvv = db['F25'] - 3*V*vv - V*V*V
    www = db['F26'] - 3*W*ww - W*W*W
    uuv = db['F28'] - 2*U*uv - V*uu - U*U*V
    uuw = db['F29'] - 2*U*uw - W*uu - U*U*W
    uvv = db['F30'] - 2*V*uv - U*vv - V*V*U
    vvw = db['F31'] - 2*V*vw - W*vv - V*V*W
    uww = db['F32'] - 2*W*uw - U*ww - W*W*U
    vww = db['F33'] - 2*W*vw - V*ww - W*W*V
    uvw = db['F34'] - U*vw - V*uw - W*uv - U*V*W


#   Normalized skewness
#    uskew = np.divide(uuu,np.real(np.power(uu,3/2)))
#    vskew = np.divide(vvv,np.real(np.power(vv,3/2)))
#    wskew = np.divide(www,np.real(np.power(ww,3/2)))

#   Velocity gradient tensor. Tensor of Rank 2.
    dUdx = db['D1']
    dVdx = db['D3']
    dWdx = db['D5']
    dUdy = db['D2']
    dVdy = db['D4']
    dWdy = db['D6']
    dUdz = np.zeros((Nx,Ny))
    dVdz = np.zeros((Nx,Ny))
    dWdz = np.zeros((Nx,Ny))


#   Production tensor. Tensor of Rank 2.
#   Definition of the production tensor assuming fully-developed flow, i.e.,
#   d()/dz=0, as a function of x and y.
#   Pij=-(<uiuk>*dUj/dxk+<ujuk>*dUi/dxk)
#   P11=-2*(<uu>*dU/dx+<uv>*dU/dy) 
#   P12=-(<uu>*dV/dx +<uv>*dV/dy+<uv>*dU/dx+<vv>*dU/dy)
#   P13=-(<uu>*dW/dx +<uv>*dW/dy+<uw>*dU/dx+<vw>*dU/dy)
#   P22=-2*(<uv>*dV/dx+<vv>*dV/dy) 
#   P23=-(<uv>*dW/dx+<vv>*dW/dy+<uw>*dV/dx+<vw>*dV/dy)
#   P33=-2*(<uw>*dW/dx+<vw>*dW/dy)
    Pxx = -2.0*(uu*dUdx + uv*dUdy)
    Pyy = -2.0*(uv*dVdx + vv*dVdy)
    Pzz = -2.0*(uw*dWdx + vw*dWdy)
    Pxy = -(uu*dVdx + uv*dVdy + uv*dUdx + vv*dUdy)
    Pyz = -(uv*dWdx + vv*dWdy + uw*dVdx + vw*dVdy)
    Pxz = -(uu*dWdx + uv*dWdy + uw*dUdx + vw*dUdy)


#   Dissipation tensor. Tensor of Rank 2.
#   Definition of the dissipation tensor assuming fully-developed flow, i.e.,
#   d()/dz=0, as a function of x and y.
#   Dij=-2*nu*<dui/dxk*duj/dxk>

#   e11_tot=<((du/dx)^2+(du/dy)^2+(du/dz)^2)>        % F39 
#   e22_tot=<((dv/dx)^2+(dv/dy)^2+(dv/dz)^2)>        % F40 
#   e33_tot=<((dw/dx)^2+(dw/dy)^2+(dw/dz)^2)>        % F41 
#   e12_tot=<(du/dx*dv/dx+du/dy*dv/dy+du/dz*dv/dz)>  % F42  
#   e13_tot=<(du/dx*dw/dx+du/dy*dw/dy+du/dz*dw/dz)>  % F43 
#   e23_tot=<(dv/dx*dw/dx+dv/dy*dw/dy+dv/dz*dw/dz)>  % F44 

#   e11=e11_tot-(dU/dx)^2-(dU/dy)^2 
#   e22=e22_tot-(dV/dx)^2-(dV/dy)^2 
#   e33=e33_tot-(dW/dx)^2-(dW/dy)^2 
#   e12=e12_tot-(dU/dx*dV/x)-(dU/dy*dV/dy)
#   e13=e13_tot-(dU/dx*dW/x)-(dU/dy*dW/dy)
#   e23=e23_tot-(dV/dx*dW/x)-(dV/dy*dW/dy)

    e11_tot = db['F39']
    e22_tot = db['F40']
    e33_tot = db['F41']
    e12_tot = db['F42']
    e13_tot = db['F43']  
    e23_tot = db['F44'] 

    e11 = e11_tot - dUdx*dUdx - dUdy*dUdy
    e22 = e22_tot - dVdx*dVdx - dVdy*dVdy
    e33 = e33_tot - dWdx*dWdx - dWdy*dWdy
    e12 = e12_tot - dUdx*dVdx - dUdy*dVdy
    e13 = e13_tot - dUdx*dWdx - dUdy*dWdy
    e23 = e23_tot - dVdx*dWdx - dVdy*dWdy

    Dxx = -2*nu*e11
    Dyy = -2*nu*e22
    Dzz = -2*nu*e33
    Dxy = -2*nu*e12 
    Dyz = -2*nu*e23
    Dxz = -2*nu*e13


#   Mean convection tensor. Tensor of Rank 2.
#   Definition of the mean convection tensor assuming 
#   fully-developed flow, i.e., d()/dz=0, as a function of x and y.
#   Cij=Uk*d<uiuj>/dxk
#   Note that under this definition: Production + Dissipation - Convection ...
#   C11=U*d(uu)/dx+V*d(uu)/dy
#   C22=U*d(vv)/dx+V*d(vv)/dy
#   C33=U*d(ww)/dx+V*d(ww)/dy
#   C12=U*d(uv)/dx+V*d(uv)/dy
#   C13=U*d(uw)/dx+V*d(uw)/dy
#   C23=U*d(vw)/dx+V*d(vw)/dy

    duudx = db['D9'] - 2*U*dUdx
    dvvdx = db['D11'] - 2*V*dVdx
    dwwdx = db['D13'] - 2*W*dWdx
    duvdx = db['D17'] - U*dVdx - V*dUdx
    duwdx = db['D21'] - U*dWdx - W*dUdx
    dvwdx = db['D19'] - V*dWdx - W*dVdx
    duudy = db['D10'] - 2*U*dUdy
    dvvdy = db['D12'] - 2*V*dVdy
    dwwdy = db['D14'] - 2*W*dWdy
    duvdy = db['D18'] - U*dVdy - V*dUdy
    duwdy = db['D22'] - U*dWdy - W*dUdy
    dvwdy = db['D20'] - V*dWdy - W*dVdy

    Cxx = U*duudx + V*duudy
    Cyy = U*dvvdx + V*dvvdy
    Czz = U*dwwdx + V*dwwdy
    Cxy = U*duvdx + V*duvdy
    Cxz = U*duwdx + V*duwdy
    Cyz = U*dvwdx + V*dvwdy


#   Turbulent transport tensor. Tensor of Rank 2.
#   Definition of the turbulent transport tensor assuming fully-developed 
#   flow, i.e., d()/dz=0, as a function of x and y.
#   Tij=-d<uiujuk>/dxk
#   T11=-(d<uuuk>)/dxk, K=1,3=-[d<uuu>/dx+d<uuv>/dy+d<uuw>/dz]
#   T22=-(d<vvuk>)/dxk, K=1,3=-[d<vvu>/dx+d<vvv>/dy+d<vvw>/dz]
#   T33=-(d<wwuk>)/dxk, K=1,3=-[d<wwu>/dx+d<wwv>/dy+d<www>/dz]
#   T12=-(d<uvuk>)/dxk, K=1,3=-[d<uvu>/dx+d<uvv>/dy+d<uvw>/dz]
#   T13=-(d<uwuk>)/dxk, K=1,3=-[d<uwu>/dx+d<uwv>/dy+d<uww>/dz]
#   T23=-(d<vwuk>)/dxk, K=1,3=-[d<vwu>/dx+d<vwv>/dy+d<vww>/dz]

    duuudx = db['D23'] - 3*U*U*dUdx - 3*(U*duudx+uu*dUdx)
    dvvudx = db['D35'] - 2*(V*duvdx+uv*dVdx) - (U*dvvdx+vv*dUdx) - (V*V*dUdx+2*U*V*dVdx)
    dwwudx = db['D39'] - 2*(W*duwdx+uw*dWdx) - (U*dwwdx+ww*dUdx) - (W*W*dUdx+2*U*W*dWdx)
    duvudx = db['D31'] - 2*(U*duvdx+uv*dUdx) - (V*duudx+uu*dVdx) - (U*U*dVdx+2*U*V*dUdx)
    duwudx = db['D33'] - 2*(U*duwdx+uw*dUdx) - (W*duudx+uu*dWdx) - (U*U*dWdx+2*U*W*dUdx)
    dvwudx = db['D43'] - (U*dvwdx+vw*dUdx) - (V*duwdx+uw*dVdx) - (W*duvdx+uv*dWdx) - (U*V*dWdx+U*W*dVdx+V*W*dUdx)
      
    duuvdy = db['D32'] - 2*(U*duvdy+uv*dUdy) - (V*duudy+uu*dVdy) - (U*U*dVdy+2*U*V*dUdy)
    dvvvdy = db['D26'] - 3*(V*dvvdy+vv*dVdy) - 3*V*V*dVdy
    dwwvdy = db['D42'] - 2*(W*dvwdy+vw*dWdy) - (V*dwwdy+ww*dVdy) - (W*W*dVdy+2*V*W*dWdy)
    duvvdy = db['D36'] - 2*(V*duvdy+uv*dVdy) - (U*dvvdy+vv*dUdy) - (V*V*dUdy+2*U*V*dVdy)
    duwvdy = db['D44'] - (U*dvwdy+vw*dUdy) - (V*duwdy+uw*dVdy) - (W*duvdy+uv*dWdy) - (U*V*dWdy+U*W*dVdy+V*W*dUdy)
    dvwvdy = db['D38'] - 2*(V*dvwdy+vw*dVdy) - (W*dvvdy+vv*dWdy) - (V*V*dWdy+2*V*W*dVdy)

    duuwdz = np.zeros((Nx,Ny))
    dvvwdz = np.zeros((Nx,Ny))
    dwwwdz = np.zeros((Nx,Ny))
    duvwdz = np.zeros((Nx,Ny))
    duwwdz = np.zeros((Nx,Ny))
    dvwwdz = np.zeros((Nx,Ny))

    Txx = - (duuudx + duuvdy + duuwdz)
    Tyy = - (dvvudx + dvvvdy + dvvwdz)
    Tzz = - (dwwudx + dwwvdy + dwwwdz)
    Txy = - (duvudx + duvvdy + duvwdz)
    Txz = - (duwudx + duwvdy + duwwdz)
    Tyz = - (dvwudx + dvwvdy + dvwwdz)


#   Viscous diffusion tensor. Tensor of Rank 2.
#   Definition of the viscous diffusion tensor assuming fully-developed 
#   flow, i.e., d()/dz=0, as a function of x and y.
#   VDij=nu*d2(uiuj)/dxk2
#   VD11=nu*(d2(uu)/dx2+d2(uu)/dy2+d2(uu)/dz2)
#   VD22=nu*(d2(vv)/dx2+d2(vv)/dy2+d2(vv)/dz2)
#   VD33=nu*(d2(ww)/dx2+d2(ww)/dy2+d2(ww)/dz2)
#   VD12=nu*(d2(uv)/dx2+d2(uv)/dy2+d2(uv)/dz2)
#   VD13=nu*(d2(uw)/dx2+d2(uw)/dy2+d2(uw)/dz2)
#   VD23=nu*(d2(vw)/dx2+d2(vw)/dy2+d2(vw)/dz2)

    d2uudx2 = db['D51'] - 2*(U*db['D45']+dUdx*dUdx)
    d2vvdx2 = db['D53'] - 2*(V*db['D47']+dVdx*dVdx)
    d2wwdx2 = db['D55'] - 2*(W*db['D49']+dWdx*dWdx)
    d2uvdx2 = db['D57'] - (V*db['D45']+U*db['D47']+2*dUdx*dVdx)
    d2uwdx2 = db['D59'] - (U*db['D49']+W*db['D45']+2*dUdx*dWdx)
    d2vwdx2 = db['D61'] - (V*db['D49']+W*db['D47']+2*dVdx*dWdx)

    d2uudy2 = db['D52'] - 2*(U*db['D46']+dUdy*dUdy)
    d2vvdy2 = db['D54'] - 2*(V*db['D48']+dVdy*dVdy)
    d2wwdy2 = db['D56'] - 2*(W*db['D50']+dWdy*dWdy)
    d2uvdy2 = db['D58'] - (V*db['D46']+U*db['D48']+2*dUdy*dVdy)
    d2uwdy2 = db['D60'] - (U*db['D50']+W*db['D46']+2*dUdy*dWdy)
    d2vwdy2 = db['D62'] - (V*db['D50']+W*db['D48']+2*dVdy*dWdy)

    d2uudz2 = np.zeros((Nx,Ny))
    d2vvdz2 = np.zeros((Nx,Ny))
    d2wwdz2 = np.zeros((Nx,Ny))
    d2uvdz2 = np.zeros((Nx,Ny))
    d2uwdz2 = np.zeros((Nx,Ny))
    d2vwdz2 = np.zeros((Nx,Ny))

    VDxx = nu*(d2uudx2 + d2uudy2 + d2uudz2)
    VDyy = nu*(d2vvdx2 + d2vvdy2 + d2vvdz2)
    VDzz = nu*(d2wwdx2 + d2wwdy2 + d2wwdz2)
    VDxy = nu*(d2uvdx2 + d2uvdy2 + d2uvdz2)
    VDxz = nu*(d2uwdx2 + d2uwdy2 + d2uwdz2)
    VDyz = nu*(d2vwdx2 + d2vwdy2 + d2vwdz2)


#   Velocity-pressure-gradient tensor. Tensor of Rank 2.
#   Definition of the velocity-pressure-gradient tensor assuming 
#   fully-developed flow, i.e., d()/dz=0, as a function of x and y.
#   Piij=-1/rho*(<ui*dp/dxj>+<uj*dp/dxi>)
#   Pi11=-2/rho*<u*dp/dx>
#   Pi22=-2/rho*<v*dp/dy>
#   Pi33=-2/rho*<w*dp/dz>
#   Pi12=-1/rho*(<u*dp/dy>+<v*dp/dx>)
#   Pi13=-1/rho*(<u*dp/dz>+<w*dp/dx>)
#   Pi23=-1/rho*(<v*dp/dz>+<w*dp/dy>)

#   Now, since we don't compute <ui*dp/dxj>, we use the following to express
#   these terms as a function of velocity gradients.
#   <ui*dp/dxj>=d(<p*ui>)/dxj-<p*dui/dxj>
# =>We define the pressure transport and pressure strain tensors.


#   Pressure transport tensor. Tensor of Rank 2.
#   Definition of the pressure transport tensor assuming 
#   fully-developed flow, i.e., d()/dz=0, as a function of x and y.
#   PTij=-1/rho*(<d(p*ui)/dxj>+<d(p*uj)/dxi>)
#   PT11=-2/rho*<d(p*u)/dx>
#   PT22=-2/rho*<d(p*v)/dy>
#   PT33=-2/rho*<d(p*w)/dz>
#   PT12=-1/rho*(<d(p*u)/dy>+<d(p*v)/dx>)
#   PT13=-1/rho*(<d(p*u)/dz>+<d(p*w)/dx>)
#   PT23=-1/rho*(<d(p*v)/dz>+<d(p*w)/dy>)

    dpudx = db['D63'] - P*dUdx - U*db['D7']
    dpvdx = db['D65'] - P*dVdx - V*db['D7']
    dpwdx = db['D67'] - P*dWdx - W*db['D7']

    dpudy = db['D64'] - P*dUdy - U*db['D8']
    dpvdy = db['D66'] - P*dVdy - V*db['D8']
    dpwdy = db['D68'] - P*dWdy - W*db['D8']

    dpudz = np.zeros((Nx,Ny))
    dpvdz = np.zeros((Nx,Ny))
    dpwdz = np.zeros((Nx,Ny))

    PTxx = 2.0/rho*dpudx
    PTyy = 2.0/rho*dpvdy
    PTzz = 2.0/rho*dpwdz
    PTxy = 1.0/rho*(dpudy+dpvdx)
    PTxz = 1.0/rho*(dpudz+dpwdx)
    PTyz = 1.0/rho*(dpvdz+dpwdy)


#   Pressure strain tensor. Tensor of Rank 2.
#   Definition of the pressure strain tensor assuming 
#   fully-developed flow, i.e., d()/dz=0, as a function of x and y.
#   PSij=-1/rho*(<p*dui/dxj>+<p*duj/dxi>)
#   PS11=-2/rho*<p*du/dx>
#   PS22=-2/rho*<p*dv/dy>
#   PS33=-2/rho*<p*dw/dz>
#   PS12=-1/rho*(<p*du/dy>+<p*dv/dx>)
#   PS13=-1/rho*(<p*du/dz>+<p*dw/dx>)
#   PS23=-1/rho*(<p*dv/dz>+<p*dw/dy>)

#   <pdudx>       % F15
#   <pdudy>       % F16
#   <pdudz>       % F17
#   <pdvdx>       % F18
#   <pdvdy>       % F19
#   <pdvdz>       % F20
#   <pdwdx>       % F21
#   <pdwdy>       % F22
#   <pdwdz>       % F23

    pdudx = db['F15'] - P*dUdx
    pdudy = db['F16'] - P*dUdy
    pdudz = db['F17'] - P*dUdz

    pdvdx = db['F18'] - P*dVdx
    pdvdy = db['F19'] - P*dVdy
    pdvdz = db['F20'] - P*dVdz

    pdwdx = db['F21'] - P*dWdx
    pdwdy = db['F22'] - P*dWdy
    pdwdz = db['F23'] - P*dWdz

    PSxx = 2.0/rho*pdudx
    PSyy = 2.0/rho*pdvdy
    PSzz = 2.0/rho*pdwdz
    PSxy = 1.0/rho*(pdudy+pdvdx)
    PSxz = 1.0/rho*(pdudz+pdwdx)
    PSyz = 1.0/rho*(pdvdz+pdwdy)

#   Construct velocity-pressure-gradient-tensor
    Pixx = -1.0/rho*(PTxx - PSxx)
    Piyy = -1.0/rho*(PTyy - PSyy)
    Pizz = -1.0/rho*(PTzz - PSzz)
    Pixy = -1.0/rho*(PTxy - PSxy)
    Pixz = -1.0/rho*(PTxz - PSxz)
    Piyz = -1.0/rho*(PTyz - PSyz)

#   Budget for each component of the Reynolds stress tensor 
#   Without mean convection
    Sxx = Pxx+Dxx+Txx+VDxx+Pixx
    Syy = Pyy+Dyy+Tyy+VDyy+Piyy
    Szz = Pzz+Dzz+Tzz+VDzz+Pizz
    Sxy = Pxy+Dxy+Txy+VDxy+Pixy
    Sxz = Pxz+Dxz+Txz+VDxz+Pixz
    Syz = Pyz+Dyz+Tyz+VDyz+Piyz

#   With mean convection
    Scxx = Pxx+Dxx+Txx+VDxx+Pixx-Cxx
    Scyy = Pyy+Dyy+Tyy+VDyy+Piyy-Cyy
    Sczz = Pzz+Dzz+Tzz+VDzz+Pizz-Czz
    Scxy = Pxy+Dxy+Txy+VDxy+Pixy-Cxy
    Scxz = Pxz+Dxz+Txz+VDxz+Pixz-Cxz
    Scyz = Pyz+Dyz+Tyz+VDyz+Piyz-Cyz




#   Save averaged-interpolated fields in cartesian coordinate(database dbCS)
    dbCS = {'U': U, 'V': V,'W': W, \
            'uu':uu, 'vv':vv, 'ww':ww, 'uv':uv, 'uw':uw, 'vw':vw, \
            'P':P, 'pp':pp, 'ppp':ppp, 'pppp':pppp, \
            'uuu':uuu, 'vvv':vvv, 'www':www, 'uuv':uuv, 'uuw':uuw, 'uvv':uvv, 'vvw':vvw, 'uww':uww, 'vww':vww, 'uvw':uvw, \
            'Pxx':Pxx, 'Pyy':Pyy, 'Pzz':Pzz, 'Pxy':Pxy, 'Pxz':Pxz, 'Pyz':Pyz, \
            'Dxx':Dxx, 'Dyy':Dyy, 'Dzz':Dzz, 'Dxy':Dxy, 'Dxz':Dxz, 'Dyz':Dyz, \
            'Txx':Txx, 'Tyy':Tyy, 'Tzz':Tzz, 'Txy':Txy, 'Txz':Txz, 'Tyz':Tyz, \
            'VDxx':VDxx, 'VDyy':VDyy, 'VDzz':VDzz, 'VDxy':VDxy, 'VDxz':VDxz, 'VDyz':VDyz, \
            'Pixx':Pixx, 'Piyy':Piyy, 'Pizz':Pizz, 'Pixy':Pixy, 'Pixz':Pixz, 'Piyz':Piyz, \
            'Cxx':Cxx, 'Cyy':Cyy, 'Czz':Czz, 'Cxy':Cxy, 'Cxz':Cxz, 'Cyz':Cyz, \
            'PTxx':PTxx, 'PTyy':PTyy, 'PTzz':PTzz, 'PTxy':PTxy, 'PTxz':PTxz, 'PTyz':PTyz, \
            'PSxx':PSxx, 'PSyy':PSyy, 'PSzz':PSzz, 'PSxy':PSxy, 'PSxz':PSxz, 'PSyz':PSyz, \
            'dUdx':dUdx, 'dUdy':dUdy, 'dUdz':dUdz, 'dVdx':dVdx, 'dVdy':dVdy, 'dVdz':dVdz, 'dWdx':dWdx, 'dWdy':dWdy, 'dWdz':dWdz, 'nu':nu \
           }


    return dbCS
