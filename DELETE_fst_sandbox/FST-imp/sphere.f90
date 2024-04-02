module sphere
  use neko
  use global_params
  ! use global_params
  implicit none
  
contains

  !     this subroutine computes a set of Np point which are 
  !     (more or less) uniformly distributed on a sphere with
  !     radius R.
  !
  !      input:    Np      Number of points
  !                        (will be changed in case another number of
  !                        points needs to be used)
  !                rad     Radius of sphere
  !                x,y,z   Arrays holding the coordinates (1..Np)
  !                rotx,roty,rotz  Rotation angles
  !                file    if file=true
  !                        --> wire.dat will be written
  subroutine compute_sphere(Np,x,y,z,rad,rotx,roty,rotz,file)
    
    integer :: Np
    real(kind=rp) :: x(Np), y(Np), z(Np), rotx, roty, rotz
    logical :: file

    integer :: i, j ,k , Nn ,kk
    real(kind=rp) :: dphi, dtheta  ,theta ,phi
    integer :: Nphi,num ,Npn
    
    real(kind=rp) :: theta1, Nphir
    logical new
    
    real(kind=rp) :: w,d,rad
    integer :: N,Nl,Npp,Npn1
    integer :: l(1000,2)
    character(len=LOG_SIZE) :: log_buf
    
    N = 1000


    !      real pi
    !      parameter (pi=3.141592653589793)

    new = .true.
    !     this flag sets the algorithm for uniformly distribution


    if (Np.gt.N) then
      write(log_buf,*) "number of points too big. Please change in sphere.f90."
      call neko_log%message(log_buf)
      stop
    end if

    if (Np.le.2) then
      call asp(x,y,z,1,0._rp,0._rp,0._rp)
      Nl=0.
      Np=1.
    else if (Np.eq.4) then
      !let's have a tetraeder
      Np=4
      call asp(x,y,z,1,-1.d0/6.d0*sqrt(3.d0),-0.5d0,0._rp)
      call asp(x,y,z,2,-1.d0/6.d0*sqrt(3.d0),0.5d0,0.d0)
      call asp(x,y,z,3,1.d0/3.d0*sqrt(3.d0),0.d0,0.d0)
      call asp(x,y,z,4,0.d0,0.d0,1.d0/3.d0*sqrt(6.d0))

      call trans(x,y,z,Np,0._rp,0._rp,-sqrt(6.d0)/12.d0)
      call scale1(x,y,z,Np,sqrt(6.d0)*2d0/3.d0)

      Nl=6
      call asl(l,1,1,2)
      call asl(l,2,2,3)
      call asl(l,3,3,1)
      call asl(l,4,4,1)
      call asl(l,5,4,2)
      call asl(l,6,4,3)

    else if (Np.eq.6) then
      !let's have a octaeder

      Np=6
      call asp(x,y,z,1,0.d0,0.d0,sqrt(2.d0)/2.d0)
      call asp(x,y,z,2,0.d0,1.d0,sqrt(2.d0)/2.d0)
      call asp(x,y,z,3,1.d0,1.d0,sqrt(2.d0)/2.d0)
      call asp(x,y,z,4,1.d0,0.d0,sqrt(2.d0)/2.d0)
      call asp(x,y,z,5,0.5d0,0.5d0,0.d0)
      call asp(x,y,z,6,0.5d0,0.5d0,sqrt(2.d0))

      call trans(x,y,z,Np,-0.5d0,-0.5d0,-sqrt(2._rp)/2._rp)
      call scale1(x,y,z,Np,sqrt(2.d0))


      Nl=12
      call asl(l,1,1,2)
      call asl(l,2,2,3)
      call asl(l,3,3,4)
      call asl(l,4,4,1)
      call asl(l,5,1,5)
      call asl(l,6,2,5)
      call asl(l,7,3,5)
      call asl(l,8,4,5)
      call asl(l,9,1,6)
      call asl(l,10,2,6)
      call asl(l,11,3,6)
      call asl(l,12,4,6)

    else if (Np.eq.8) then
      !let's have a cube
      Np=8
      call asp(x,y,z,1,-sqrt(3.d0)/3.d0,sqrt(3.d0)/3.d0,-sqrt(3.d0)/3.d0)
      call asp(x,y,z,2,-sqrt(3.d0)/3.d0,-sqrt(3.d0)/3.d0,-sqrt(3.d0)/3.d0)
      call asp(x,y,z,3,sqrt(3.d0)/3.d0,-sqrt(3.d0)/3.d0,-sqrt(3.d0)/3.d0)
      call asp(x,y,z,4,sqrt(3.d0)/3.d0,sqrt(3.d0)/3.d0,-sqrt(3.d0)/3.d0)
      call asp(x,y,z,5,-sqrt(3.d0)/3.d0,sqrt(3.d0)/3.d0,sqrt(3.d0)/3.d0)
      call asp(x,y,z,6,-sqrt(3.d0)/3.d0,-sqrt(3.d0)/3.d0,sqrt(3.d0)/3.d0)
      call asp(x,y,z,7,sqrt(3.d0)/3.d0,-sqrt(3.d0)/3.d0,sqrt(3.d0)/3.d0)
      call asp(x,y,z,8,sqrt(3.d0)/3.d0,sqrt(3.d0)/3.d0,sqrt(3.d0)/3.d0)

      Nl=12
      call asl(l,1,1,2)
      call asl(l,2,2,3)
      call asl(l,3,3,4)
      call asl(l,4,4,1)
      call asl(l,5,5,6)
      call asl(l,6,6,7)
      call asl(l,7,7,8)
      call asl(l,8,8,5)
      call asl(l,9,1,5)
      call asl(l,10,2,6)
      call asl(l,11,3,7)
      call asl(l,12,4,8)

    else if (Np.eq.12) then
      !     let's have an ikosaeder
      Np=12
      w=0.5*(sqrt(5.)+1.d0)
      call asp(x,y,z,1,w/2.d0, 0.d0, 0.5d0*(w-1))
      call asp(x,y,z,2,w/2.d0, 0.d0, 0.5d0*(w+1))
      call asp(x,y,z,3,0.d0,0.5d0*(w-1.d0),0.5d0*w)
      call asp(x,y,z,4,0.d0,0.5d0*(w+1.d0),0.5d0*w)
      call asp(x,y,z,5,0.5d0*w,w,0.5d0*(w-1.d0))
      call asp(x,y,z,6,0.5d0*w,w,0.5d0*(w+1.d0))
      call asp(x,y,z,7,w,(w+1.d0)/2.,w/2.)
      call asp(x,y,z,8,w,(w-1.d0)/2.,w/2.)
      call asp(x,y,z,9,0.5*(w+1.d0), 0.5*w,w)
      call asp(x,y,z,10,0.5*(w-1),0.5*w,w)
      call asp(x,y,z,11,0.5*(w+1),0.5*w,0.d0)
      call asp(x,y,z,12,0.5*(w-1),0.5*w,0.d0)

      call trans(x,y,z,Np,-0.5*w,-0.5*w,-0.5*w)
      call scale1(x,y,z,Np,2./(sqrt(w**2+1.d0)))



      Nl=30
      call asl(l,1,1,2)
      call asl(l,2,3,4)
      call asl(l,3,5,6)
      call asl(l,4,7,8)
      call asl(l,5,9,10)
      call asl(l,6,11,12)
      call asl(l,7,11,5)
      call asl(l,8,11,7)
      call asl(l,9,11,8)
      call asl(l,10,11,1)
      call asl(l,11,12,1)
      call asl(l,12,12,3)
      call asl(l,13,12,4)
      call asl(l,14,12,5)
      call asl(l,15,9,6)
      call asl(l,16,9,7)
      call asl(l,17,9,8)
      call asl(l,18,9,2)
      call asl(l,19,10,2)
      call asl(l,20,10,3)
      call asl(l,21,10,4)
      call asl(l,22,10,6)
      call asl(l,23,8,1)
      call asl(l,24,8,2)
      call asl(l,25,3,1)
      call asl(l,26,3,2)
      call asl(l,27,4,5)
      call asl(l,28,4,6)
      call asl(l,29,7,5)
      call asl(l,30,7,6)

    else if (Np.eq.20) then
      !        let's have a dodekaeder

      w=0.5*(sqrt(5.)+3.)
      Np=20
      call asp(x,y,z,1,0.5*w,0.5*(w-1.d0),0.d0)
      call asp(x,y,z,2,.5*w,.5*(w+1.d0),0.d0)
      call asp(x,y,z,3,w-0.5,w-.5,.5d0)
      call asp(x,y,z,4,w,.5*w,.5*(w-1.d0))
      call asp(x,y,z,5,w-.5,0.5d0,.5d0)
      call asp(x,y,z,6,.5*(w+1.d0),0.d0,.5*w)
      call asp(x,y,z,7,.5*(w-1.d0),0.d0,.5*w)
      call asp(x,y,z,8,0.5d0,0.5d0,0.5d0)
      call asp(x,y,z,9,0.d0,.5*w,(w-1.d0)*.5d0)
      call asp(x,y,z,10,.5d0,w-.5d0,.5d0)
      call asp(x,y,z,11,.5*(w-1.d0),w,.5*w)
      call asp(x,y,z,12,.5*(w+1.d0),w,.5*w)
      call asp(x,y,z,13,w-.5,w-.5,w-.5)
      call asp(x,y,z,14,w,.5*w,.5*(w+1.d0))
      call asp(x,y,z,15,w-.5,.5d0,w-.5)
      call asp(x,y,z,16,.5*w,.5*(w-1.d0),w)
      call asp(x,y,z,17,.5d0,.5d0,w-.5d0)
      call asp(x,y,z,18,0.d0,.5*w,.5*(w+1.d0))
      call asp(x,y,z,19,.5d0,w-.5,w-.5)
      call asp(x,y,z,20,.5d0*w,.5*(w+1.d0),w)

      call trans(x,y,z,Np,-0.5*w,-0.5*w,-0.5*w)
      call scale1(x,y,z,Np,2./(sqrt(w**2+1.d0)))

      Nl=30
      call asl(l,1,1,2)
      call asl(l,2,2,3)
      call asl(l,3,3,4)
      call asl(l,4,4,5)
      call asl(l,5,5,6)
      call asl(l,6,6,7)
      call asl(l,7,7,8)
      call asl(l,8,8,9)
      call asl(l,9,9,10)
      call asl(l,10,10,11)
      call asl(l,11,11,12)
      call asl(l,12,12,13)
      call asl(l,13,13,14)
      call asl(l,14,14,15)
      call asl(l,15,15,16)
      call asl(l,16,16,17)
      call asl(l,17,17,18)
      call asl(l,18,18,19)
      call asl(l,19,19,20)
      call asl(l,20,3,12)
      call asl(l,21,19,11)
      call asl(l,22,13,20)
      call asl(l,23,16,20)
      call asl(l,24,7,17)
      call asl(l,25,6,15)
      call asl(l,26,1,5)
      call asl(l,27,1,8)
      call asl(l,28,2,10)
      call asl(l,29,4,14)
      call asl(l,30,9,18)


    else if (.not.(new)) then

      !       having uniformly distributed sphere

      Npp=Np

      Nn=(1.d0+sqrt(-3.+2.*Npp))/2.
      Nn = Nn*2
      Npn =  real(Nn, kind=rp)**2./2.-real(Nn, kind=rp)+2.
      if (file) then
        if (Np.ne.Npn) then
          call neko_log%message('number of points is not exactly the same...')
          call print_param('Np ', real(Np, kind=rp))
          call print_param('Npn', real(Npn, kind=rp))
        else
          call print_param('Np', real(Np, kind=rp))
        end if
      end if

      Np=Npn

      dphi= 2*pi / real(Nn, kind=rp)
      dtheta = 2*pi / real(Nn, kind=rp)

      k=1
      x(k)=0.
      y(k)=0.
      z(k)=1.d0
      do j=1, Nn/2-1
        do i=1, (Nn)
          k=k+1
          x(k) = cos(i*dphi)*sin(j*dtheta)
          y(k) = sin(i*dphi)*sin(j*dtheta)
          z(k) = cos ( j*dtheta )
        end do
      end do
      k=k+1
      x(k)=0.
      y(k)=0.
      z(k)=-1.d0

      k=0
      do i=2, Nn+1
        k=k+1
        call asl(l,k,1,i)
        k=k+1
        call asl(l,k,Npn,Npn-i+1)
      end do

      do j=1,Nn/2-1
        do i=1, Nn  -1
          k=k+1
          call asl(l,k,(j-1)*Nn+i+1,(j-1)*Nn+i+2)
          if (j.ne.1) then
            k=k+1
            call asl(l,k,(j-1)*Nn+i+1,(j-2)*Nn+i+1)
          end if
        end do
        k=k+1
        call asl(l,k,(j-1)*Nn+i+1,(j-1)*Nn+2)
        if (j.ne.1) then
          k=k+1
          call asl(l,k,(j-1)*Nn+i+1,(j-2)*Nn+i+1)
        end if
      end do
      Nl=k

    else if (new)  then
      !     having it differently
      k=0
      kk=0

      Nn=1
      500   Npn=0
      dtheta = 2*pi / real(Nn, kind=rp)
      do j=0, (Nn)/2
        theta=j*dtheta
        if (sin(theta).eq.0.) then
          dphi = 99999999.
        else
          dphi = dtheta / sin(theta)
        end if
        Nphir = max(2*pi / dphi ,1.d0)
        Nphi = Nphir+0.5
        dphi = 2*pi / real(Nphi, kind=rp)
        Npn=Npn+Nphi
      end do
      if (Npn.LE.Np) then
      !         write(*,*) Npn
        Npn1=Npn
        Nn=Nn+1
        goto 500
      end if
      Npn=Npn1
      Nn=Nn-1
      if (file) then
          if (Np.ne.Npn) then
            write(log_buf, *) 'number of points is not exactly the same...'
            call neko_log%message(log_buf)
            call print_param('Np ', real(Np, kind=rp))
            call print_param('Npn', real(Npn, kind=rp))
          else
            call print_param('Np', real(Np, kind=rp))
          end if
      end if


      dtheta = 2*pi / real(Nn, kind=rp)
      do j=0, (Nn)/2
        theta=j*dtheta
        if (sin(theta).eq.0.) then
          dphi = 99999999.
        else
          dphi = dtheta / sin(theta)
        end if
        Nphir = max(2*pi / dphi ,1.d0)
        Nphi = Nphir+0.5
        dphi = 2*pi / real(Nphi, kind=rp)
        do i=1, (Nphi)
          phi=i*dphi
          k=k+1
          x(k) = cos(phi)*sin(theta)
          y(k) = sin(phi)*sin(theta)
          z(k) = cos ( theta )
          kk=kk+1
          call asl(l,kk,k,k+1)
        end do
      end do
      Np=k
      Nl=kk-1

      end if


        !     do the scaling to rad
      call scale1(x,y,z,Np,rad)

        !     do the rotation
      call rot3d(Np,x,y,z,rotx,roty,rotz)


        !      check...
        !      do i=1, Np
        !        d=x(i)**2.+y(i)**2.+z(i)**2.
        !        write(*,*) d,d**2.
        !      end do




      if (file) then
      write(log_buf, *) 'writing to file wire.dat'
        call neko_log%message(log_buf)
        open(file='wire.dat', status='unknown', form='formatted', unit=10)
        write(10,'(2i5)') Np,Nl
        do i=1, Np
          write(10,'(3E18.9)') x(i),y(i),z(i)
        end do
        do i=1, Nl
          write(10,'(2i5)') l(i,1),l(i,2)
        end do
        close(10)
      end if


  end subroutine compute_sphere

  subroutine asp(x,y,z,i,xx,yy,zz)
    implicit none
    ! integer N
    ! N=1000
    integer i
    real(kind=rp) :: xx,yy,zz
    real(kind=rp) :: x(1000),y(1000),z(1000)
    x(i)=xx
    y(i)=yy
    z(i)=zz
  end subroutine asp
  
  subroutine asl(l,i,i1,i2)
    implicit none
    integer N
    parameter (N=1000)
    integer i,i1,i2
    integer l(N,2)
    l(i,1)=i1
    l(i,2)=i2
  end subroutine asl
  
  subroutine trans(x,y,z,Np,xx,yy,zz)
    implicit none
    integer Np,i
    real(kind=rp) :: xx,yy,zz
    real(kind=rp) :: x(Np),y(Np),z(Np)
    do i=1,Np
      x(i)=x(i)+xx
      y(i)=y(i)+yy
      z(i)=z(i)+zz
    end do
  end subroutine trans
  
  subroutine scale1(x,y,z,Np,r)
    implicit none
    integer Np,i
    real(kind=rp) :: r
    real(kind=rp) :: x(Np),y(Np),z(Np)
    do i=1,Np
      x(i)=x(i)*r
      y(i)=y(i)*r
      z(i)=z(i)*r
    end do
  end subroutine scale1

  subroutine rot3d(Np,x,y,z,rotx,roty,rotz)
    implicit none
    integer Np,i
    real(kind=rp) :: rotx,roty,rotz
    real(kind=rp) :: x(Np),y(Np),z(Np)
    real(kind=rp) :: cx,cy,cz,sx,sy,sz
    real(kind=rp) :: xx,yy,zz ,d

    cx=cos(rotx)
    sx=sin(rotx)
    cy=cos(roty)
    sy=sin(roty)
    cz=cos(rotz)
    sz=sin(rotz)

    do i=1,Np
      xx=x(i)
      yy=y(i)
      zz=z(i)
  
      x(i)=xx*cy*cz+yy*cy*sz-zz*sy
      y(i)=xx*(cz*sx*sy-cx*sz)+yy*(cx*cz+sx*sy*sz)+zz*cy*sx
      z(i)=xx*(cx*cz*sy+sx*sz)-yy*(cz*sx-cx*sy*sz)+zz*cx*cy
  
    end do
  end subroutine rot3d
  
end module sphere
