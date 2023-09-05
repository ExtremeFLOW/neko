cc-----------------------------------------------------------------------
C
C  USER SPECIFIED ROUTINES:
C
C     - boundary conditions
C     - initial conditions
C     - variable properties
C     - forcing function for fluid (f)
C     - forcing function for passive scalar (q)
C     - general purpose routine for checking errors etc.
C
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'


      if (ifield.eq.1) then 
        UTRANS = param(1)      ! rho=1.000
        UDIFF  = param(2)      ! nue=sqrt(Pr/Ra)
      else          
        UTRANS = param(7)      ! rho*cp=1.000
        UDIFF  = param(8)      ! kappa=1/sqrt(Ra*Pr)
      endif

  
      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      FFX = 0.0
      FFY = 0.0
      FFZ = 0.0

      if (if3d) then
         ffx = 0.0   
         ffy = 0.0
         ffz = temp   ! Buoyancy force in 3d 
      elseif (ifaxis) then
         ffx = -temp
         ffy = 0.0
         ffz = 0.0
      else   ! 2D
         ffx = 0.0
         ffy = temp
         ffz = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      QVOL   = 0.0
      SOURCE = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'

c-----JS 11/02/2018 

      parameter (lt=lx1*ly1*lz1*lelt)

      real dtdx(lt),dtdy(lt),dtdz(lt)
      real  uzt(lt),epst(lt),epsv(lt)
      real dudx(lt),dudy(lt),dudz(lt)
      real dvdx(lt),dvdy(lt),dvdz(lt)
      real dwdx(lt),dwdy(lt),dwdz(lt)
      real grad_valv,grad_valt

c-----Comment out the line below if not restarting
C      if (istep.gt.0.and.istep.lt.nbdinp) call my_full_restart_load


      call my_full_restart_save  ! save add'l files for full-restart      


      n = nx1*ny1*nz1*nelv
      if (istep.le.100.or.mod(istep,20).eq.0) then
         wmin = glmin(vz,n)
         wmax = glmax(vz,n)
         tmin = glmin(t ,n)
         tmax = glmax(t ,n)
         if (nid.eq.0) write(6,1) istep,time,tmin,tmax,wmin,wmax
    1    format(i8,1p5e12.4,' uz_t_mx')
      endif


c-----The first data file in each job dumps the spectral element grid
c     together with the fields (u,T,p, ...) 
c---------------------------------------------------------------------
      ifxyo = .true.
      if (istep.gt.iostep) ifxyo = .false.

      if ((mod(istep,iostep).eq.0).and.(istep.ge.iostep)) then
!-----Temperature gradient  
        call gradm1(dtdx,dtdy,dtdz,t)

!-----Velocity gradients
        call gradm1(dudx,dudy,dudz,vx)
        call gradm1(dvdx,dvdy,dvdz,vy)
        call gradm1(dwdx,dwdy,dwdz,vz)
 
        do i=1,n

!-----Convective current
          uzt(i) = t(i,1,1,1,1)*vz(i,1,1,1)

!-----Thermal dissipation epsT=(grad T)**2/sqrt(RaPr)
          grad_valt=0.

          grad_valt=dtdx(i)**2+dtdy(i)**2+dtdz(i)**2
          epst(i)=vdiff(i,1,1,1,2)*grad_valt

!-----Energy dissipation epsv=0.5*(du_i/dx_j+du_j/dx_i)**2*sqrt(Pr/Ra)
          grad_valv=0.
 
          grad_valv=grad_valv+(dudx(i)+dudx(i))**2
          grad_valv=grad_valv+(dudy(i)+dvdx(i))**2
          grad_valv=grad_valv+(dudz(i)+dwdx(i))**2
 
          grad_valv=grad_valv+(dvdx(i)+dudy(i))**2
          grad_valv=grad_valv+(dvdy(i)+dvdy(i))**2
          grad_valv=grad_valv+(dvdz(i)+dwdy(i))**2

          grad_valv=grad_valv+(dwdx(i)+dudz(i))**2
          grad_valv=grad_valv+(dwdy(i)+dvdz(i))**2
          grad_valv=grad_valv+(dwdz(i)+dwdz(i))**2
          epsv(i)=0.5*vdiff(i,1,1,1,1)*grad_valv
        enddo

        call vertical_mean(uzt,1)
        call vertical_mean(epst,2)
        call vertical_mean(epsv,3)
        call vertical_mean(t,4)
        call vertical_mean(dtdz,5)

	call volume_mean(epst,1)
	call volume_mean(epsv,2)
        call volume_mean(uzt,3)

      endif
 
      call userchk_flux


      return
      end
c-----------------------------------------------------------------------
      subroutine vertical_mean(ff,i_name)
      include 'SIZE'
      include 'TOTAL'

      parameter(nelz=300)

      integer i_name

      real ff(lx1,ly1,lz1,lelt)
      real fbar(lz1,nelz)
      real zbar(lz1,nelz)
      real wght(lz1,nelz)
      real work(lz1,nelz)

      integer e,eg,ex,ey,ez,f

      nelxy = nelgv/nelz 

!-----Set both arrays to zero
      call rzero(fbar,lz1*nelz)
      call rzero(zbar,lz1*nelz)
      call rzero(wght,lz1*nelz)


!-----Pick face 5 to evaluate surface Jacobian
      f = 5 

      do e=1,nelv

         eg = lglel(e)
         call get_exyz(ex,ey,ez,eg,nelxy,1,nelz)

         do k=1,nz1
           do i=1,nx1*ny1
             fbar(k,ez) = fbar(k,ez)+area(i,1,f,e)*ff(i,1,k,e)
             zbar(k,ez) = zbar(k,ez)+area(i,1,f,e)*zm1(i,1,k,e)
             wght(k,ez) = wght(k,ez)+area(i,1,f,e)
           enddo
         enddo

      enddo

!-----Gather over all processes (-> mpi_allreduce)
      call gop(fbar,work,'+  ',lz1*nelz)
      call gop(zbar,work,'+  ',lz1*nelz)
      call gop(wght,work,'+  ',lz1*nelz)

!-----Area average
      do i=1,lz1*nelz
        fbar(i,1)=fbar(i,1)/wght(i,1) 
        zbar(i,1)=zbar(i,1)/wght(i,1)
      enddo

!-----Output of the vertical profile of plane averages 
!
!     Example with 32 (=nelz) vertical elements and lx1=ly1=lz1=4:
!      -> 2 internal vertical GLL points per element
!      -> 97 vertical planes, here output of 96
!
!------------------------------------------------------------ 
      if(nid.eq.0)then

        if(i_name.eq.1)OPEN(10,file="ver_uzte.dat",position="append")
        if(i_name.eq.2)OPEN(10,file="ver_epst.dat",position="append")
        if(i_name.eq.3)OPEN(10,file="ver_epsv.dat",position="append")
        if(i_name.eq.4)OPEN(10,file="ver_temp.dat",position="append")
        if(i_name.eq.5)OPEN(10,file="ver_dtdz.dat",position="append")

        do k=1,nelz
        do i=1,lz1-1
          write(10,*)zbar(i,k),fbar(i,k)
        enddo   
        enddo
        CLOSE(10)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine pdf_calc(ff,step,i_offset,i_name)
      include 'SIZE'
      include 'TOTAL'

      parameter(npdf=3601)

      real ff(lx1,ly1,lz1,lelt)
      real pdf(npdf)
      real work(npdf)
      real val, offset, wght, norm
      integer e 

!-----Set arrays to zero
      call rzero(pdf,npdf)
      call rzero(work,npdf)

!-----Offset
      if(i_offset==0) offset=0.0 
      if(i_offset==1) offset=int(npdf/2)*step 

      do e=1,nelv
        do k=1,nz1
          do j=1,ny1
            do i=1,nx1

               wght= bm1(i,j,k,e) 
               val= ff(i,j,k,e)

               do ipdf=1,npdf
                dm1=(ipdf -1 )*step-offset
                dm2= ipdf     *step-offset

                if((val.ge.dm1).and.(val.lt.dm2))then
                  pdf(ipdf)=pdf(ipdf)+wght
                endif 

               enddo    
            enddo
          enddo
        enddo
      enddo

!-----Normalization     
      do i=1,npdf 
        pdf(i)=pdf(i)/volvm1
      enddo

!-----Gather over all processes (-> mpi_allreduce)
      call gop(pdf,work,'+  ',npdf)

!------------------------------------------------------------ 
      if(nid.eq.0)then

        if(i_name.eq.1)OPEN(10,file="pdf_uzte.dat",position="append")
        if(i_name.eq.2)OPEN(10,file="pdf_epst.dat",position="append")
        if(i_name.eq.3)OPEN(10,file="pdf_epsv.dat",position="append")
        if(i_name.eq.4)OPEN(10,file="pdf_temp.dat",position="append")
        if(i_name.eq.5)OPEN(10,file="pdf_etak.dat",position="append")

        do ipdf=1,npdf
          write(10,*) (ipdf-1)*step-offset, pdf(ipdf)
        enddo
        CLOSE(10)

!-------Normalization check
        norm=0.0
        do ipdf=1,npdf
          norm=norm+pdf(ipdf)
        enddo
        if(i_name.eq.1)write(77,*)'uz*T',norm
        if(i_name.eq.2)write(77,*)'epsT',norm
        if(i_name.eq.3)write(77,*)'epsv',norm
        if(i_name.eq.4)write(77,*)'Temp',norm
        if(i_name.eq.5)write(77,*)'etaK',norm

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine volume_mean(ff,i_name)
      include 'SIZE'
      include 'TOTAL'

     
      real ff(lx1,ly1,lz1,lelt)
      real work1
      real vol_check, mean_val1, mean_val2, mean_val3, mean_val4
      integer e 

!-----Set arrays to zero
      call rzero(mean_val1,1)
      call rzero(mean_val2,1)
      call rzero(mean_val3,1)
      call rzero(mean_val4,1)

      call rzero(vol_check,1)
      call rzero(work1,1)

      do e=1,nelv
        do k=1,nz1
          do j=1,ny1
            do i=1,nx1

               vol_check=vol_check + bm1(i,j,k,e)
	       mean_val1=mean_val1 + ff(i,j,k,e)*bm1(i,j,k,e)
	       mean_val2=mean_val2 + (ff(i,j,k,e)**2)*bm1(i,j,k,e)
	       mean_val3=mean_val3 + (ff(i,j,k,e)**3)*bm1(i,j,k,e)	       
	       mean_val4=mean_val4 + (ff(i,j,k,e)**4)*bm1(i,j,k,e)

            enddo
          enddo
        enddo
      enddo

!-----Normalization      
      mean_val1=mean_val1/volvm1
      mean_val2=mean_val2/volvm1
      mean_val3=mean_val3/volvm1
      mean_val4=mean_val4/volvm1

!-----Gather over all processes (-> mpi_allreduce)
      call gop(mean_val1,work1,'+  ',1)
      call gop(mean_val2,work1,'+  ',1)
      call gop(mean_val3,work1,'+  ',1)
      call gop(mean_val4,work1,'+  ',1)
      call gop(vol_check,work1,'+  ',1)

      if(nid.eq.0)  write(*,*)'  0',volvm1
      if(nid.eq.100)write(*,*)'100',volvm1

!------------------------------------------------------------ 
      if(nid.eq.0)then

100   format(5(1x,e17.9))
 
      if(i_name.eq.1)OPEN(10,file="mean_epst.dat",position="append")
      if(i_name.eq.2)OPEN(10,file="mean_epsv.dat",position="append")
      if(i_name.eq.3)OPEN(10,file="mean_uzte.dat",position="append")

      write(10,100)mean_val1,mean_val2,mean_val3,mean_val4,vol_check

      CLOSE(10)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'


      ux=0.0
      uy=0.0
      uz=0.0

!-----Job 0: add a little bias to get convection started
!      arg  = -istep/.005
!      bias = x*0.2*exp(arg)   ! use this to break axisymmetry
!      temp = 1-z + bias       ! early in the simulation, only

!-----Job 1, 2, ...
      temp=1-z

      return
      end
c-----------------------------------------------------------------------
c
c this makes an initial condition of random thermal perturbations
c

      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      real harvest
 
      ux = 0.0
      uy = 0.0
      uz = 0.0
      call random_number(harvest)
      temp = 1-z + 0.001*harvest

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat
      include 'SIZE'
      include 'TOTAL'

      integer e

      n = 8*nelv
      do i=1,n            ! Project edge vertices onto circle
         x=xc(i,1)
         y=yc(i,1)
         rr = x*x+y*y
         if (rr.gt.0) rr = sqrt(rr)
         if (rr.gt.0.049999) then
            xc(i,1) = xc(i,1)*0.05/rr
            yc(i,1) = yc(i,1)*0.05/rr
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      include 'SIZE'
      include 'TOTAL'

      return
      end
c-----------------------------------------------------------------------
c
c     Automatically added by makenek
c
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk_flux
      include 'SIZE'  
      include 'TOTAL' 

c     JDS DEBUG 6/6/2008
c     This routine is for post-processing data, such as computing
c     integrated heat fluxes, etc.
c     JS 9/20/2012
c     reducing output and change flux to Nusselt at both plates 

      parameter (lt=lx1*ly1*lz1*lelt)
      real dtdx(lt),dtdy(lt),dtdz(lt),w(lt)

      common /chkcmnr/ atime,timel,flux(2)
      common /chkcmni/ iesurf(0:2*lelt,2),ifsurf(0:2*lelt,2)

      integer icalld
      save    icalld
      data    icalld  /0/
      CHARACTER*3 CB1B

      if (icalld.eq.0) then
         icalld = icalld + 1
         atime  = 0.
         timel  = time

         call find_srf(iesurf(0,1),ifsurf(0,1),5)    ! bottom face
         call find_srf(iesurf(0,2),ifsurf(0,2),6)    ! top    face

      endif
      if (istep.le.0) return


!-----Compute temperature gradients for surface fluxes
      call gradm1(dtdx,dtdy,dtdz,t)

      do k = 1,2
         flux(k)  = 0.0
         do isurf = 1,iesurf(0,k)
            ie    = iesurf(isurf,k)
            iface = ifsurf(isurf,k) 
            CB1B = CBC(IFACE,IE,1)
            if (CB1B.EQ.'E  ') then
!-------------do nothing
            else
                 call surface_flux(dq,dtdx,dtdy,dtdz,ie,iface,w)
            flux(k)  = flux(k)  + dq
            end if
         enddo
      enddo

!-----Sum over all processors
      call gop(flux,w,'+  ',2)

!-----Output of Nusselt number at top/bottom plate
      dtime = time  - timel
      atime = atime + dtime
      pi4=4.0*atan(1.0)*(0.05**2)

      if (nid.eq.0) then
        write(6,1) istep,time,atime,flux(1)/pi4,-flux(2)/pi4
      endif
    1 format(i6,' Nusselt',1p4e14.6)

      timel = time

      return
      end
c-----------------------------------------------------------------------
      subroutine find_srf(iesurf,ifsurf,inface)
c
c     Find list of surfaces over which the flux should be computed
c
c     Number of such surfaces is returned in iesurf(0).
c
c
      include 'SIZE'
      include 'TOTAL'

      integer iesurf(0:lelt),ifsurf(0:lelt)

      nsurf = 0
      nfaces = 2*ndim
      do ie=1,nelv
         nsurf         = nsurf+1
         iesurf(nsurf) = ie
         ifsurf(nsurf) = inface
      enddo
      iesurf(0) = nsurf
      ifsurf(0) = nsurf
      return
      end
c-----------------------------------------------------------------------
      function ran2(delta)
c
      real*8 F7,T
      PARAMETER ( F7=78125., T30=1073741824.)
c
      save t
c
      integer icalld
      save    icalld
      data    icalld/0/
C
C    INITIALIZATION:
C    Generate PSEUDO-RANDOM  (0., 1.)  DATA
C    USING A RANDOM NUMBER GENERATOR BASED ON THE RECURSION
C           X(N+1) = 5**7 * X(N)  (MOD 2**30)
C    THIS RECURSION WILL GENERATE 2**28 (APPROX. 268 MILLION) NUMBERS
C    BEFORE REPEATING.  FOR THIS SCHEME TO WORK PROPERLY, THE HARDWARE
C    MULTIPLY OPERATION MUST BE CORRECT TO 47 BITS OF PRECISION.
C
      if (icalld.eq.0) T = F7 / T30
      icalld = 1
      T   = MOD (F7 * T, 1.)
      ran2 =  delta*t
      return
      end
c-----------------------------------------------------------------------
c
c     .usr file with full-restart routines, called from userchk
c
c-----------------------------------------------------------------------
      subroutine my_full_restart_save    ! Call this from userchk

c     Saves files for next full restart


      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'  ! max_rst


c     This is the full-restart save part:

      max_rst = 2*(nbdinp-1)    ! max # of rst files saved

 	      
      nps1 = 0
      if (ifheat) nps1 = 1 + npscal

      mostep = mod1(istep,iostep)
      if (istep.gt.iostep.and.mostep.lt.nbdinp) 
     $    call outpost2(vx,vy,vz,pr,t,nps1,'rst')


      return
      end
c-----------------------------------------------------------------------
      subroutine my_full_restart_load    ! Call this from userchk
      include 'SIZE'
      include 'TOTAL'

      character*80 s80

      call blank(s80,80)

 
      if (istep.eq.1) s80 ='rstG0p1_3000.f00001'  



      call bcast(s80,80)
      call chcopy(initc,s80,80)

c     time_curr = time

      nfiles = 1
      call restart(nfiles)  ! Note -- time is reset.

      if (nid.ne.0) time=0
      time = glmax(time,1)  ! Synchronize time across all processors

c     time = time_curr      ! Preserve current simulation time

      return
      end

c-----------------------------------------------------------------------
      subroutine set_obj  ! define objects for surface integrals
c
      include 'SIZE'
      include 'TOTAL'
c
      integer e,f
c
c     Define new objects

      nobj = 1
      iobj = 0
      do ii=nhis+1,nhis+nobj
         iobj = iobj+1
         hcode(10,ii) = 'I'
         hcode( 1,ii) = 'F' ! 'F'
         hcode( 2,ii) = 'F' ! 'F'
         hcode( 3,ii) = 'F' ! 'F'
         lochis(1,ii) = iobj
      enddo
      nhis = nhis + nobj

      if (maxobj.lt.nobj) write(6,*) 'increase maxobj in SIZE'
      if (maxobj.lt.nobj) call exitt

      do e=1,nelv
      do f=1,2*ndim
         if (cbc(f,e,1).eq.'W  ') then
            iobj = 0
c           if (f.eq.1) iobj=1  ! lower wall
c           if (f.eq.3) iobj=2  ! upper wall
            iobj=1              ! cylinder wall
            if (iobj.gt.0) then
               nmember(iobj) = nmember(iobj) + 1
               mem = nmember(iobj)
               ieg = lglel(e)
               object(iobj,mem,1) = ieg
               object(iobj,mem,2) = f
c              write(6,1) iobj,mem,f,ieg,e,nid,' OBJ'
    1          format(6i9,a4)
            endif
c
         endif
      enddo
      enddo
c     write(6,*) 'number',(nmember(k),k=1,4)
c
      return
      end
c-----------------------------------------------------------------------
