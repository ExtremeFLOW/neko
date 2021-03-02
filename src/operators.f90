module operators
  use num_types
  use space  
  use math
  use mesh
  use mxm_wrapper
  use coefs
  use field
  use gather_scatter
  use mathops
  implicit none
contains
  subroutine dudxyz (du,u,dr,ds,dt,coef)
!--------------------------------------------------------------
!
!     du  - dU/dx or dU/dy or dU/dz
!     u   - a field variable defined on mesh 1
!     dr  - dr/dx or dr/dy or dr/dz  
!     ds  - ds/dx or ds/dy or ds/dz
!     dt  - dt/dx or dt/dy or dt/dz
    type(coef_t), intent(in), target :: coef
    real(kind=dp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), intent(inout) ::  du
    real(kind=dp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), intent(inout) ::  u, dr, ds, dt
    real(kind=dp) :: drst(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz)
    type(space_t), pointer :: Xh 
    type(mesh_t), pointer :: msh
     integer :: e, k, lxy, lyz, lxyz
     Xh => coef%Xh
     msh => coef%msh 
     lxy  = Xh%lx*Xh%ly
     lyz  = Xh%ly*Xh%lz
     lxyz = Xh%lx*Xh%ly*Xh%lz
      
     do e=1,msh%nelv
      if (msh%nelv .eq. 2) then
         call mxm     (Xh%dx,Xh%lx,u(1,1,1,e),Xh%lx,du(1,1,1,e),lyz)
         call col2    (du(1,1,1,e),dr(1,1,1,e),lxyz)
         call mxm     (U(1,1,1,e),Xh%lx,Xh%dyt,Xh%ly,drst,Xh%ly)
         call addcol3 (du(1,1,1,e),drst,ds(1,1,1,e),lxyz)
      else
         call mxm   (Xh%dx,Xh%lx,U(1,1,1,e),Xh%lx,du(1,1,1,e),lyz)
         call col2  (du(1,1,1,e),dr(1,1,1,e),lxyz)
         do k=1,Xh%lz
            call mxm  (u(1,1,k,e),Xh%lx,Xh%dyt,Xh%ly,drst(1,1,k),Xh%ly)
         end do
         call addcol3 (du(1,1,1,e),drst,ds(1,1,1,e),lxyz)
         call mxm     (U(1,1,1,e),lxy,Xh%dzt,Xh%lz,drst,Xh%lz)
         call addcol3 (du(1,1,1,e),drst,dt(1,1,1,e),lxyz)
      end if
   end do
   call col2 (du,coef%jacinv,coef%dof%n_dofs)
  end subroutine dudxyz

  !> Equals wgradm1 in nek5000. Gradient of velocity vectors.
  subroutine opgrad(ux,uy,uz,u,coef) ! weak form of grad 

  !Compute gradient of T -- mesh 1 to mesh 1 (vel. to vel.)

    type(coef_t) :: coef  
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: ux,uy,uz,u
    real(kind=dp) :: ur(coef%Xh%lxyz)
    real(kind=dp) :: us(coef%Xh%lxyz)
    real(kind=dp) :: ut(coef%Xh%lxyz)
    integer :: e, i, N

    N = coef%Xh%lx - 1
    do e=1,coef%msh%nelv
       if(coef%msh%gdim .eq. 3) then
         call local_grad3(ur,us,ut,u(1,e),N,coef%Xh%dx,coef%Xh%dxt)
            do i=1,coef%Xh%lxyz
               ux(i,e) = coef%Xh%w3(i,1,1)*(ur(i)*coef%drdx(i,1,1,e) &
                       + us(i)*coef%dsdx(i,1,1,e) &
                       + ut(i)*coef%dtdx(i,1,1,e) )
               uy(i,e) = coef%Xh%w3(i,1,1)*(ur(i)*coef%drdy(i,1,1,e) &
                       + us(i)*coef%dsdy(i,1,1,e) &
                       + ut(i)*coef%dtdy(i,1,1,e) )
               uz(i,e) = coef%Xh%w3(i,1,1)*(ur(i)*coef%drdz(i,1,1,e) &
                       + us(i)*coef%dsdz(i,1,1,e) &
                       + ut(i)*coef%dtdz(i,1,1,e) )
          enddo
       else

            call local_grad2(ur,us,u(1,e),N,coef%Xh%dx,coef%Xh%dyt)

            do i=1,coef%Xh%lxyz
               ux(i,e) = coef%Xh%w3(i,1,1)*(ur(i)*coef%drdx(i,1,1,e) &
                       + us(i)*coef%dsdx(i,1,1,e) )
               uy(i,e) = coef%Xh%w3(i,1,1)*(ur(i)*coef%drdy(i,1,1,e) &
                       + us(i)*coef%dsdy(i,1,1,e) )
          enddo
       endif
    enddo
  end subroutine opgrad
  
  subroutine local_grad3(ur, us, ut, u, n, D, Dt)
    integer, intent(inout) :: n
    real(kind=dp), intent(inout) :: ur(0:n, 0:n, 0:n)
    real(kind=dp), intent(inout) :: us(0:n, 0:n, 0:n)
    real(kind=dp), intent(inout) :: ut(0:n, 0:n, 0:n)
    real(kind=dp), intent(inout) :: u(0:n, 0:n, 0:n)
    real(kind=dp), intent(inout) :: D(0:n, 0:n)
    real(kind=dp), intent(inout) :: Dt(0:n, 0:n)
    integer :: m1, m2, k
  
    m1 = n + 1
    m2 = m1*m1
  
    call mxm(D ,m1,u,m1,ur,m2)
    do k=0,n
       call mxm(u(0,0,k),m1,Dt,m1,us(0,0,k),m1)
    enddo
    call mxm(u,m2,Dt,m1,ut,m1)
    
  end subroutine local_grad3

  subroutine local_grad2(ur, us, u, n, D, Dt)
    integer, intent(inout) :: n    
    real(kind=dp), intent(inout) :: ur(0:n, 0:n)
    real(kind=dp), intent(inout) :: us(0:n, 0:n)
    real(kind=dp), intent(inout) :: u(0:n, 0:n)
    real(kind=dp), intent(inout) :: D(0:n, 0:n)
    real(kind=dp), intent(inout) :: Dt(0:n, 0:n)
    integer :: m1, m2, k
  
    m1 = n + 1
  
    call mxm(D ,m1,u,m1,ur,m1)
    call mxm(u,m1,Dt,m1,us,m1)
    
  end subroutine local_grad2

  !> Othogonalize with regard to vector (1,1,1,1,1,1...,1)^T.
  subroutine ortho(x,n ,glb_n)
    integer, intent(in) :: n
    integer, intent(in) :: glb_n
    real(kind=dp), dimension(n), intent(inout) :: x
    real(kind=dp) :: rlam

    rlam = glsum(x,n)/glb_n
    call cadd(x,-rlam,n)

  end subroutine ortho
  !> Compute DT*X (entire field)
  !> This needs to be revised... the loop over n1,n2 is probably unesccssary
  subroutine cdtp (dtx,x,dr,ds,dt, coef)
    type(coef_t) :: coef
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dtx, x, dr, ds, dt
    real(kind=dp) :: wx(coef%Xh%lxyz)
    real(kind=dp) :: ta1(coef%Xh%lxyz)
    real(kind=dp) :: ta2(coef%Xh%lxyz)
    real(kind=dp) :: ta3(coef%Xh%lxyz)
    integer :: e, i1, i2, n1, n2, iz
    type(space_t), pointer :: Xh 
 
    Xh => coef%Xh
    n1 = Xh%lx*Xh%ly
    n2 = Xh%lx*Xh%ly

    do e=1,coef%msh%nelv
       call col3 (wx,coef%B(1,1,1,e),x(1,e),Xh%lxyz)
       call invcol2(wx,coef%jac(1,1,1,e),Xh%lxyz)
       call col3 (ta1,wx,dr(1,e),Xh%lxyz)
       call mxm  (Xh%dxt,Xh%lx,ta1,Xh%lx,dtx(1,e),Xh%lyz)
       call col3 (ta1,wx,ds(1,e),Xh%lxyz)
       i1 = 1
       i2 = 1
       do iz=1,Xh%lz
          call mxm  (ta1(i2),Xh%lx,Xh%dy,Xh%ly,ta2(i1),Xh%ly)
          i1 = i1 + n1
          i2 = i2 + n2
       enddo
       call add2 (dtx(1,e),ta2,Xh%lxyz)
       call col3 (ta1,wx,dt(1,e),Xh%lxyz)
       call mxm  (ta1,Xh%lxy,Xh%dz,Xh%lz,ta2,Xh%lz)
       call add2 (dtx(1,e),ta2,Xh%lxyz)
    enddo
  end subroutine cdtp
   
  subroutine conv1(du,u, vx, vy, vz, Xh, coef, nelv, gdim)  ! used to be conv1n
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: nelv, gdim
    real(kind=dp), intent(inout) ::  du(Xh%lxyz,nelv)
    real(kind=dp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  u, vx, vy, vz
!   Store the inverse jacobian to speed this operation up
    real(kind=dp), dimension(Xh%lx,Xh%ly,Xh%lz) :: dudr, duds, dudt
    integer :: ie, iz, i
!   Compute vel.grad(u)
    do ie=1,nelv
      if (gdim .eq. 3) then
         call mxm   (Xh%dx,Xh%lx,u(1,1,1,ie),Xh%lx,dudr,Xh%lxy)
         do iz=1,Xh%lz
           call mxm (u(1,1,iz,ie),Xh%lx,Xh%dyt,Xh%ly,duds(1,1,iz),Xh%ly)
         enddo
         call mxm   (u(1,1,1,ie),Xh%lxy,Xh%dzt,Xh%lz,dudt,Xh%lz)
         do i=1,Xh%lxyz
            du(i,ie) = coef%jacinv(i,1,1,ie)*( &
                       vx(i,1,1,ie)*( &
                       coef%drdx(i,1,1,ie)*dudr(i,1,1) &
                     + coef%dsdx(i,1,1,ie)*duds(i,1,1) &
                     + coef%dtdx(i,1,1,ie)*dudt(i,1,1)) &
                     + vy(i,1,1,ie)*( &
                       coef%drdy(i,1,1,ie)*dudr(i,1,1) &
                     + coef%dsdy(i,1,1,ie)*duds(i,1,1) &
                     + coef%dtdy(i,1,1,ie)*dudt(i,1,1)) &
                     + vz(i,1,1,ie)*( &
                       coef%drdz(i,1,1,ie)*dudr(i,1,1) &
                     + coef%dsdz(i,1,1,ie)*duds(i,1,1) &
                     + coef%dtdz(i,1,1,ie)*dudt(i,1,1)))
         enddo
       else
!        2D
         call mxm (Xh%dx,Xh%lx,u(1,1,1,ie),Xh%lx,dudr,Xh%lyz)
         call mxm (u(1,1,1,ie),Xh%lx,Xh%dyt,Xh%ly,duds,Xh%ly)
         do i=1,Xh%lxyz
            du(i,ie) = coef%jacinv(i,1,1,ie)*( &
                       vx(i,1,1,ie)*( &
                       coef%drdx(i,1,1,ie)*dudr(i,1,1) &
                     + coef%dsdx(i,1,1,ie)*duds(i,1,1)) &
                     + vy(i,1,1,ie)*( &
                       coef%drdy(i,1,1,ie)*dudr(i,1,1) &
                     + coef%dsdy(i,1,1,ie)*duds(i,1,1)))
          enddo
        endif
     enddo
  end subroutine conv1

  subroutine curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    type(field_t), intent(inout) :: w1
    type(field_t), intent(inout) :: w2
    type(field_t), intent(inout) :: w3
    type(field_t), intent(inout) :: u1
    type(field_t), intent(inout) :: u2
    type(field_t), intent(inout) :: u3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(coef_t), intent(inout)  :: c_Xh 
    integer :: gdim, n

    n = w1%dof%size()
    gdim = c_Xh%msh%gdim

    !     this%work1=dw/dy ; this%work2=dv/dz
    call dudxyz(work1%x, u3%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    if (gdim .eq. 3) then
       call dudxyz(work2%x, u2%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       call sub3(w1%x, work1%x, work2%x, n)
    else
       call copy(w1%x, work1%x, n)
    endif
    !     this%work1=du/dz ; this%work2=dw/dx
    if (gdim .eq. 3) then
       call dudxyz(work1%x, u1%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       call dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    else
       call rzero (work1%x, n)
       call dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    endif
    !     this%work1=dv/dx ; this%work2=du/dy
    call dudxyz(work1%x, u2%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call dudxyz(work2%x, u1%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call sub3(w3%x, work1%x, work2%x, n)
    !!    BC dependent, Needs to change if cyclic

    call opcolv(w1%x,w2%x,w3%x,c_Xh%B, gdim, n)
    call gs_op(c_Xh%gs_h, w1, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w2, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w3, GS_OP_ADD) 
    call opcolv  (w1%x,w2%x,w3%x,c_Xh%Binv, gdim, n)

  end subroutine curl

  function cfl(dt, u, v, w, Xh, coef, nelv, gdim)
    type(space_t) :: Xh
    type(coef_t) :: coef
    integer :: nelv, gdim
    real(kind=dp) :: dt
    real(kind=dp) ::  du(Xh%lxyz,nelv)
    real(kind=dp), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  u, v, w
    real(kind=dp) :: cflr, cfls, cflt, cflm, cfl_temp(1)
    real(kind=dp) :: ur, us, ut
    real(kind=dp) :: cfl
    integer :: i,j,k,e
    cfl_temp(1) = 0d0
    if (gdim .eq. 3) then
       do e=1,nelv
          do k=1,Xh%lz
          do j=1,Xh%ly
          do i=1,Xh%lx
             ur = ( u(i,j,k,e)*coef%drdx(i,j,k,e) &
                +   v(i,j,k,e)*coef%drdy(i,j,k,e) &
                +   w(i,j,k,e)*coef%drdz(i,j,k,e) ) * coef%jacinv(i,j,k,e)
             us = ( u(i,j,k,e)*coef%dsdx(i,j,k,e) &
                +   v(i,j,k,e)*coef%dsdy(i,j,k,e) &
                +   w(i,j,k,e)*coef%dsdz(i,j,k,e) ) * coef%jacinv(i,j,k,e)
             ut = ( u(i,j,k,e)*coef%dtdx(i,j,k,e) &
                +   v(i,j,k,e)*coef%dtdy(i,j,k,e) &
                +   w(i,j,k,e)*coef%dtdz(i,j,k,e) ) * coef%jacinv(i,j,k,e)
 
             cflr = abs(dt*ur*Xh%dr_inv(i))
             cfls = abs(dt*us*Xh%ds_inv(j))
             cflt = abs(dt*ut*Xh%dt_inv(k))
 
             cflm = cflr + cfls + cflt
             cfl_temp(1)  = max(cfl_temp(1),cflm)
          enddo
          enddo
          enddo
       enddo
    else
       do e=1,nelv
          do j=1,Xh%ly
          do i=1,Xh%lx
             ur = ( u(i,j,1,e)*coef%drdx(i,j,1,e) &
                +   v(i,j,1,e)*coef%drdy(i,j,1,e) ) * coef%jacinv(i,j,k,e)
             us = ( u(i,j,1,e)*coef%dsdx(i,j,1,e) &
                +   v(i,j,1,e)*coef%dsdy(i,j,1,e) ) * coef%jacinv(i,j,k,e)

             cflr = abs(dt*ur*Xh%dr_inv(i))
             cfls = abs(dt*us*Xh%ds_inv(j))

             cflm = cflr + cfls
             cfl_temp(1)  = max(cfl_temp(1),cflm)

          enddo
          enddo
       enddo
    endif
    cfl = glmax(cfl_temp,1)
  end function cfl
  
end module operators
