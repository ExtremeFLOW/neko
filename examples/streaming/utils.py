import sys
import os
import copy

import numpy as np
import scipy.optimize
import time
import matplotlib.pyplot as plt
import pymech as pm
from pymech.neksuite import readnek
from pymech.neksuite import writenek
from tqdm import tqdm
import collections

def reshape_data(y_in, nxyz,nel):
    #Reshape the arrays for simplicity
    y_in=y_in.reshape(-1,1)
    y_out=np.zeros((nxyz,nel))
    e=0
    for e in range(0, nel):
        y_out[:,e]=y_in[e*nxyz:(e+1)*nxyz,0]

    return y_out

def get_mod_info(path,casename,prefix):
        
    path=path
    casename=casename

    filename=path+'PODmod'+'_'+prefix+'_'+casename+'0.f'+str(0+1).zfill(5)
    data=readnek(filename)
    # simulation parameters
    buff  = np.zeros((10))
    buff[0]=data.nel #nel
    buff[1]=data.lr1[0] #lx1
    buff[2]=data.lr1[1] #ly1
    buff[3]=data.lr1[2] #lz1
    buff[4]=np.prod(data.lr1) #nxyz
    buff[5]= data.nel*data.lr1[0]*data.lr1[1]*data.lr1[2] #nxyze
    if data.lr1[2]==1:
        buff[6]=2
    else:
        buff[6]=3


    nel=int(buff[0])
    lx1=int(buff[1])
    ly1=int(buff[2])
    lz1=int(buff[3])
    nxyz=int(buff[4])
    nxyze=int(buff[5])
    dim=int(buff[6])

    return [nel, lx1, ly1, lz1, nxyz, nxyze, dim]


def get_data(path, casename, prefix,qoi):

    filename_sig=path+'PODsig'+'_'+prefix+'_'+casename+'.npy'
    filename_vtr=path+'PODvtr'+'_'+prefix+'_'+casename+'.npy'

    ## Load the time coefficient information
    D = np.load(filename_sig)
    Vt = np.load(filename_vtr)
    m = Vt.shape[0] 

    ## Get the information of the modes
    [nel, lx1, ly1, lz1, nxyz, nxyze, dim] = get_mod_info(path,casename,prefix)
    ## Allocate space for the modes
    U = np.zeros((nxyze,m))
    ## Read the data from the modes one by one
    pbar= tqdm(total=m)
    for j in range(0,m):
        filename=path+'PODmod'+'_'+prefix+'_'+casename+'0.f'+str(j+1).zfill(5)
        data=readnek(filename)
        for e in range(0,nel):
            #Rearange 
            ve=data.elem[e].vel[qoi,:,:,:].reshape((nxyz,1),order='F')
            #Copy into a bigger vector to scatter
            U[e*nxyz:e*nxyz+nxyz,j]=np.copy((ve[:,0]))
        pbar.update(1)
    pbar.close()
    
    return U,D,Vt

def get_error_metrics(UU,DD,VTVT,pref,mesh,ifplot,ref,exp,path,casename):

    print('Error metrics')
    pbar= tqdm(total=len(UU))

    # Mesh data
    xx_mesh=mesh[0]
    yy_mesh=mesh[1]
    zz_mesh=mesh[2]
    elemplot=mesh[3]
    # Reference data
    U_ref=ref[0]
    D_ref=ref[1]
    Vt_ref=ref[2]
    p_ref=500
    k_ref=500
    S_ref = U_ref@np.diag(D_ref)@Vt_ref

    # Error lists
    DDerr=[]
    UUerr=[]
    VTerr=[]
    SSerr=[]

    for i in range(0,len(UU)):

        U=UU[i]
        D=DD[i]
        Vt=VTVT[i]
        p=exp[i][0]
        k=exp[i][1]
        mr=exp[i][2]
        prefix=pref[i]

        gr=1
        first_mode=0
        num_plots=5
        lev=100
        xx=xx_mesh[:,elemplot].flatten()
        yy=yy_mesh[:,elemplot].flatten()
        zz=zz_mesh[:,elemplot].flatten()

        
        # Reconstruction and reconstruction error
        S = U@np.diag(D)@Vt
        SErr = np.zeros((S.shape[1]))
        for j in range(0,S.shape[1]):
            SErr[j]=np.sqrt(np.mean((abs(S[:,j])-abs(S_ref[:,j]))**2)) 
        
        #Check difference in singular values
        DErr = np.sqrt((D_ref[0:k]-D)**2)

        #Check the root mean squared error of each of the modes
        UErr = np.zeros((U.shape[1]))
        for j in range(0,U.shape[1]):
            UErr[j]=np.sqrt(np.mean((abs(U[:,j])-abs(U_ref[:,j]))**2)) 
         
        #Check the root mean squared error of Vt
        VtErr = np.zeros((Vt.shape[1]))
        for j in range(0,Vt.shape[1]):
            VtErr[j]=np.sqrt(np.mean((abs(Vt[:,j])-abs(Vt_ref[:k,j]))**2)) 

        DDerr.append(DErr)
        UUerr.append(UErr)
        VTerr.append(VtErr)
        SSerr.append(SErr)

        if ifplot == True:
            
            # Plot info for singular values
            fig, ax = plt.subplots(1,1)
            ax.plot(D_ref[0:k],'-k')
            ax.plot(D,'--b')
            plt.yscale('log')        
            filename=path+'sig_err'+'_'+prefix+'_'+casename+'.png'
            plt.savefig(filename,format="png") 
        
            # Plot info for modes
            ## Plot root mean square error per mode
            fig, ax = plt.subplots(1,1)
            ax.plot(UErr)
            plt.yscale('log')        
            filename=path+'mod_err'+'_'+prefix+'_'+casename+'.png'
            plt.savefig(filename,format="png")
        
            ## Plot slice of the mode to understand
            ### See which part of the mode to plot
            U_plot_ref=[]
            U_plot=[]
            for j in range(first_mode,first_mode+num_plots):
                temp=reshape_data(U_ref[:,j],xx_mesh.shape[0],xx_mesh.shape[1])
                temp1=temp[:,elemplot].flatten()
                U_plot_ref.append(temp1)
                temp=reshape_data(U[:,j],xx_mesh.shape[0],xx_mesh.shape[1])
                temp1=temp[:,elemplot].flatten()
                U_plot.append(temp1)
            
            ### Perform the plot
            fig, ax = plt.subplots(num_plots, 2)
            cc=[]
            cc2=[]
            for j in range(0,num_plots):
                c = ax[j,0].tricontourf(xx,yy,U_plot_ref[j],levels=lev,cmap='RdBu_r')
                cc.append(c)
                c = ax[j,1].tricontourf(xx,yy,U_plot[j],levels=lev,cmap='RdBu_r')
            cc2.append(c)

            for j in range(0,2):
                for f in range(0,num_plots):
                    ax[f,j].tick_params(left = False, right = False , labelleft = False , labelbottom = False, bottom = False)
       
            #for j in range(0,num_plots):
            #    cbar = fig.colorbar(cc[j],ax=ax[j,0])
            #    cbar = fig.colorbar(cc2[j],ax=ax[j,1])

            filename=path+'mod_dif'+'_'+prefix+'_'+casename+'.png'
            plt.savefig(filename,format="png",bbox_inches="tight") 

            # Plot the info on the right subspace
            fig, ax = plt.subplots(1,1)
            for f in range(0,k,int(k/gr)):
                ax.plot(Vt_ref[f,:],'-k')
                ax.plot(Vt[f,:],'--b')
            filename=path+'vtr_err'+'_'+prefix+'_'+casename+'.png'
            plt.savefig(filename,format="png")

            # Reconstrution error
            fig, ax = plt.subplots(1,1)
            ax.plot(SErr)
            plt.yscale('log')        
            filename=path+'rec_err'+'_'+prefix+'_'+casename+'.png'
            plt.savefig(filename,format="png")

        pbar.update(1)
    pbar.close()

    return DDerr, UUerr, VTerr, SSerr


def read_mesh(path,casename):
    
    filename=path+casename+'0.f'+str(0+1).zfill(5)
    data=readnek(filename)

    # simulation parameters
    nel=data.nel
    lx1=data.lr1[0]
    ly1=data.lr1[1]
    lz1=data.lr1[2]
    nxyz=np.prod(data.lr1)
    nxyze= nel*lx1*ly1*lz1

    n=lx1
    if lz1==1:
        dim=2
    else:
        dim=3

    # Get the transformation matrix
    #[V,Vinv]=fun.get_transform_matrix(n,dim)

    x_glb=np.zeros((nxyze,1))
    y_glb=np.zeros((nxyze,1))
    z_glb=np.zeros((nxyze,1))

    pbar= tqdm(total=nel)
    for e in range(0,nel): 

        # Read and process the element data
        xe=data.elem[e].pos[0,:,:,:].reshape((nxyz,1),order='F')
        ye=data.elem[e].pos[1,:,:,:].reshape((nxyz,1),order='F')

        if dim==3:
            ze=data.elem[e].pos[2,:,:,:].reshape((nxyz,1),order='F')
                
        x_glb[e*nxyz:(e+1)*nxyz,:]=xe
        y_glb[e*nxyz:(e+1)*nxyz,:]=ye

        if dim==3:
            z_glb[e*nxyz:(e+1)*nxyz,:]=ze

        pbar.update(1)

    pbar.close()
    return x_glb,y_glb,z_glb,nel,nxyz, lx1, ly1, lz1,n,dim


def prepare_mesh(path,casename,xl:float=0,xr:float=20,yd:float=0,yu:float=14,zd:float=0):

    # Read the mesh
    x_mesh,y_mesh,z_mesh, nel,nxyz, _,_,_,_,_ = read_mesh(path,casename)
    xx_mesh=reshape_data(x_mesh,nxyz,nel)
    yy_mesh=reshape_data(y_mesh,nxyz,nel)
    zz_mesh=reshape_data(z_mesh,nxyz,nel)
    
    # Select where to perform the plots
    ifplot = True
    #xl=0
    #xr=20
    #yd=0
    #yu=14
    #zd=0
    zu=zd+0.01

    # Calculate elements to be plotted
    elemplot=[]
    for e in range(0,nel):
        xm=np.mean(xx_mesh[:,e])
        ym=np.mean(yy_mesh[:,e])
        zm=np.mean(zz_mesh[:,e])
        if xm>=xl and xm<=xr and ym>=yd and ym <=yu and zm>=zd and zm <zu:
            elemplot.append(e)
    
    return xx_mesh,yy_mesh,zz_mesh,elemplot,nel,nxyz

def create_sobol_matrix(exp,params,SSerr,which_snaps,which_q):

    '''
    which_comp alludes to which were the comparisons made.
    which_comp = 0 --> pp vs kk
    which_comp = 1 --> pp vs mrmr
    '''

    if which_q == "PvsK":
        exp_ind1=0
        exp_ind2=1
        q1=params[exp_ind1]
        q2=params[exp_ind2]
    elif which_q == "PvsMR":
        exp_ind1=0
        exp_ind2=2
        q1=params[exp_ind1]
        q2=params[exp_ind2]

    #allocate the vector
    n=[]
    n.append(len(q1))
    n.append(len(q2))
    fEx = np.zeros(n)
    
    # Construct the quantity list
    q=[]
    q.append(np.array(q1))
    q.append(np.array(q2))
    print(q)

    pdf=[]
    for i in range(len(q)):
        pdf.append(np.ones(n[i])/(q[i][-1]-q[i][0]))

    # Construct the function output
    for i in range(0,n[0]): 
        for j in range(0,n[1]):
            # Loop to search in the list
            for ii in range(0,len(exp)):
                if q[0][i] == exp[ii][exp_ind1] and q[1][j] == exp[ii][exp_ind2]:
                    fEx[i,j]=1

    print(fEx)
    # Construct the function output
    for i in range(0,n[0]): 
        for j in range(0,n[1]):
            # Loop to search in the list
            for ii in range(0,len(exp)):
                if q[0][i] == exp[ii][exp_ind1] and q[1][j] == exp[ii][exp_ind2]:
                    fEx[i,j]=SSerr[ii][which_snaps]

    print(fEx)


    fval_pce = np.zeros((n[0]*n[1]))
    iii=0
    for j in range(0,n[1]):
        for i in range(0,n[0]):
            # Loop to search in the list
            for ii in range(0,len(exp)):
                if q[0][i] == exp[ii][exp_ind1] and q[1][j] == exp[ii][exp_ind2]:
                    fval_pce[iii]=SSerr[ii][which_snaps]
                    iii+=1


    return q, pdf,  fEx, fval_pce

def get_prefix(exp, gb, au,i):

    p=int(exp[i][0])
    k=int(exp[i][1])
    min_ratio=float(exp[i][2])
    str_mr = str(min_ratio).replace('.','')
    ifgbl_update=bool(int(gb))
    ifautoupdate=bool(int(au))

    if ifgbl_update==True: 
        if ifautoupdate==True:
            prefix='auto_gbl_p'+repr(p)+'_k'+repr(k)+'_mr'+str_mr
        else:
            prefix='fixd_gbl_p'+repr(p)+'_k'+repr(k)+'_mr'+str_mr
    else:
        if ifautoupdate==True:
            prefix='auto_lcl_p'+repr(p)+'_k'+repr(k)+'_mr'+str_mr
        else:
            prefix='fixd_lcl_p'+repr(p)+'_k'+repr(k)+'_mr'+str_mr
    
    return prefix


















