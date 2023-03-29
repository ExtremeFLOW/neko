#!/usr/bin/python
import numpy as np
from math import pi

#-----------------------------------------------------------------------
## Define functions for the calculation of the quadrature points (Taken from the lecture notes)
def GLC_pwts(n):
    """ 
    Gauss-Lobatto-Chebyshev (GLC) points and weights over [-1,1]    
    Args: 
      `n`: int, number of nodes
    Returns 
       `x`: 1D numpy array of size `n`, nodes         
       `w`: 1D numpy array of size `n`, weights
    """
    def delt(i,n):
        del_=1.
        if i==0 or i==n-1:
           del_=0.5 
        return del_
    x=np.cos(np.arange(n)*pi/(n-1))
    w=np.zeros(n)    
    for i in range(n):
        tmp_=0.0
        for k in range(int((n-1)/2)):
            tmp_+=delt(2*k,n)/(1-4.*k**2)*np.cos(2*i*pi*k/(n-1))
        w[i]=tmp_*delt(i,n)*4/float(n-1)
    return x,w 

def GLL_pwts(n,eps=10**-8,maxIter=1000):
    """
    Generating `n `Gauss-Lobatto-Legendre (GLL) nodes and weights using the 
    Newton-Raphson iteration.
    Args:    
      `n`: int
         Number of GLL nodes
      `eps`: float (optional) 
         Min error to keep the iteration running
      `maxIter`: float (optional)
         Max number of iterations
    Outputs:
      `xi`: 1D numpy array of size `n`
         GLL nodes
      `w`: 1D numpy array of size `n`
         GLL weights
    Reference:
       Canuto C., Hussaini M. Y., Quarteroni A., Tang T. A., 
       "Spectral Methods in Fluid Dynamics," Section 2.3. Springer-Verlag 1987.
       https://link.springer.com/book/10.1007/978-3-642-84108-8
    """
    V=np.zeros((n,n))  #Legendre Vandermonde Matrix
    #Initial guess for the nodes: GLC points
    xi,w_=GLC_pwts(n)
    iter_=0
    err=1000
    xi_old=xi
    while iter_<maxIter and err>eps:
        iter_+=1
        #Update the Legendre-Vandermonde matrix
        V[:,0]=1.
        V[:,1]=xi
        for j in range(2,n):
            V[:,j]=((2.*j-1)*xi*V[:,j-1] - (j-1)*V[:,j-2])/float(j)
        #Newton-Raphson iteration 
        xi=xi_old-(xi*V[:,n-1]-V[:,n-2])/(n*V[:,n-1])
        err=max(abs(xi-xi_old).flatten())
        xi_old=xi
    if (iter_>maxIter and err>eps):
       print('gllPts(): max iterations reached without convergence!')
    #Weights
    w=2./(n*(n-1)*V[:,n-1]**2.)
    return xi,w

def get_transform_matrix(n,dim):
    # Get the quadrature nodes
    x,w_=GLL_pwts(n) # The outputs of this functions are not exactly in the order we want (start from 1 not -1)

    # Reorder the quadrature nodes
    x=np.flip(x)
    w=np.flip(w_)

    # Create a diagonal matrix
    WW=np.eye(n)
    for i in range(0,n):
        WW[i,i]=w[i]

    ## First we need the legendre polynomials
    # order of the polynomials
    p=n
    # Create a counter for the loops
    p_v=np.arange(p)

    # get the legendre polynomial matrix

    # The polynomials are stored in a matrix with the following structure:
    #  |  pi_0(x0)  ... pi_0(x_n)
    #  |  pi_1(x0)  ... pi_1(x_n)
    #  |  ...
    #  |  pi_p(x0)  ... pi_p(x_n)
    #  The acending rows represent accending polynomial order, the different columns represent different x_i

    # Allocate space
    Leg=np.zeros((p,p))
    #First row is filled with 1 according to recursive formula
    Leg[0,:]=np.ones((1,p))
    #Second row is filled with x according to recursive formula
    Leg[1,:]=np.multiply(np.ones((1,p)),x)    

    # Apply the recursive formula for all x_i
    for j in range(1,len(p_v)-1):
        for k_ in p_v:
            Leg[j+1,k_]=((2*j+1)*x[k_]*Leg[j,k_]-j*Leg[j-1,k_])/(j+1)

    Leg=Leg.T # nek and I transpose it for transform

    #Scaling factor as in books
    delta=np.ones(n)
    for i in range(0,n):
        #delta[i]=2/(2*i+1)       #it is the same both ways
        delta[i]=2/(2*(i+1)-1)
    delta[n-1]=2/(n-1)
    #print(delta)
    #Scaling factor to normalize
    for i in range(0,n):
        delta[i]=np.sqrt(1/delta[i])

    #apply the scaling factor
    for i in range(0,n):
        for j in range(0,n):
            Leg[i,j]=Leg[i,j]*delta[j]

    AA=np.matmul(Leg.T,np.matmul(WW,Leg))

    #2d transformation matrix
    v=Leg
    V2d=np.kron(v,v)
    vinv=Leg.T@WW
    Vinv2d=np.kron(vinv,vinv)

    #3d transformation matrix
    v=Leg
    V=np.kron(v,np.kron(v,v))
    vinv=Leg.T@WW
    Vinv=np.kron(vinv,np.kron(vinv,vinv))
    
    if dim==2:
        V=V2d
        Vinv=Vinv2d
    if dim==3:
        V=V
        Vinv=Vinv
        
    return V,Vinv
        
            
def get_derivative_matrix(n,dim):
    # Get the quadrature nodes
    x,w_=GLL_pwts(n) # The outputs of this functions are not exactly in the order we want (start from 1 not -1)

    # Reorder the quadrature nodes
    x=np.flip(x)
    w=np.flip(w_)

    # Create a diagonal matrix
    WW=np.eye(n)
    for i in range(0,n):
        WW[i,i]=w[i]

    ## First we need the legendre polynomials
    # order of the polynomials
    p=n
    # Create a counter for the loops
    p_v=np.arange(p)

    # get the legendre polynomial matrix

    # The polynomials are stored in a matrix with the following structure:
    #  |  pi_0(x0)  ... pi_0(x_n)
    #  |  pi_1(x0)  ... pi_1(x_n)
    #  |  ...
    #  |  pi_p(x0)  ... pi_p(x_n)
    #  The acending rows represent accending polynomial order, the different columns represent different x_i

    # Allocate space
    Leg=np.zeros((p,p))
    #First row is filled with 1 according to recursive formula
    Leg[0,:]=np.ones((1,p))
    #Second row is filled with x according to recursive formula
    Leg[1,:]=np.multiply(np.ones((1,p)),x)    

    # Apply the recursive formula for all x_i
    for j in range(1,len(p_v)-1):
        for k_ in p_v:
            Leg[j+1,k_]=((2*j+1)*x[k_]*Leg[j,k_]-j*Leg[j-1,k_])/(j+1)


    D_N=np.zeros((p,p))

    # Simply apply the values as given in the book
    for i in range(0,len(p_v)):
        for j in range(0,len(p_v)):
            if (i != j):
                D_N[i,j]=(Leg[p-1,i]/Leg[p-1,j])*(1/(x[i]-x[j]))
            if (i==0 and j==0):
                D_N[i,j]=-(((p-1)+1)*(p-1))/4    
            if (i==(p-1) and j==(p-1)):
                D_N[i,j]=(((p-1)+1)*(p-1))/4
            if (i==j and i != 0 and i!=(p-1)):
                D_N[i,j]=0

    DX2D=np.kron(D_N,np.eye(p))
    DY2D=np.kron(np.eye(p),D_N)
 
    if dim==2:
        DX=DX2D
        DY=DY2D
        
    if dim==3:
        a=1
        
    return DX,DY,D_N

def local_grad2(y,DX,DY,nel):

    dur = np.copy(y)
    dus = np.copy(y)

    for e in range(0,nel):
        dur[:,e]= DX@y[:,e]
        dus[:,e]= DY@y[:,e]

    return dur, dus


def grad2(y,rx,sx,ry,sy,jac,n,dim,nxyz,nel):

    DR,DS,D_N = get_derivative_matrix(n,dim)
    
    dr,ds = local_grad2(y,DR,DS,nel)
    
    dx = np.zeros(y.shape)
    dy = np.zeros(y.shape)

    for e in range(0,nel):
        for i in range(0,nxyz):
            dx[i,e]=jac[i,e]*(dr[i,e]*rx[i,e]+ds[i,e]*sx[i,e]) 
            dy[i,e]=jac[i,e]*(dr[i,e]*ry[i,e]+ds[i,e]*sy[i,e])    
    
    return dx,dy
