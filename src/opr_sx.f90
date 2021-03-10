!> Operators SX-Aurora backend
module opr_sx
  use sx_dudxyz
  use gather_scatter
  use num_types
  use space
  use coefs
  use math
  use mesh
  use field

  use mathops
  implicit none
contains

    subroutine opr_sx_dudxyz(du, u, dr, ds, dt, coef)
    type(coef_t), intent(in), target :: coef
    real(kind=dp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), intent(inout) ::  du
    real(kind=dp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), intent(inout) ::  u, dr, ds, dt


    select case(coef%Xh%lx)
    case(12)
       call sx_dudxyz_lx12(du, u, dr, ds, dt, coef)
    case(11)
       call sx_dudxyz_lx11(du, u, dr, ds, dt, coef)
    case(10)
       call sx_dudxyz_lx10(du, u, dr, ds, dt, coef)
    case(9)
       call sx_dudxyz_lx9(du, u, dr, ds, dt, coef)
    case(8)
       call sx_dudxyz_lx8(du, u, dr, ds, dt, coef)
    case(6)
       call sx_dudxyz_lx6(du, u, dr, ds, dt, coef)
    end select
    
   end subroutine opr_sx_dudxyz

   subroutine opr_sx_opgrad(ux,uy,uz,u,coef) 
     type(coef_t), intent(in) :: coef  
     real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: ux
     real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uy
     real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uz
     real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: u

  end subroutine opr_sx_opgrad

  subroutine opr_sx_cdtp(dtx,x,dr,ds,dt, coef)
    type(coef_t), intent(in) :: coef
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dtx
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: x
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dr
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: ds
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dt

  end subroutine opr_sx_cdtp

  subroutine opr_sx_conv1(du,u, vx, vy, vz, Xh, coef, nelv, gdim)  
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: nelv, gdim
    real(kind=dp), intent(inout) ::  du(Xh%lxyz,nelv)
    real(kind=dp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  u
    real(kind=dp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vx
    real(kind=dp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vy
    real(kind=dp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vz

   end subroutine opr_sx_conv1

   subroutine opr_sx_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
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
     call opr_sx_dudxyz(work1%x, u3%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
     if (gdim .eq. 3) then
        call opr_sx_dudxyz(work2%x, u2%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
        call sub3(w1%x, work1%x, work2%x, n)
     else
        call copy(w1%x, work1%x, n)
     endif
     !     this%work1=du/dz ; this%work2=dw/dx
     if (gdim .eq. 3) then
        call opr_sx_dudxyz(work1%x, u1%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
        call opr_sx_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
        call sub3(w2%x, work1%x, work2%x, n)
     else
        call rzero (work1%x, n)
        call opr_sx_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
        call sub3(w2%x, work1%x, work2%x, n)
     endif
     !     this%work1=dv/dx ; this%work2=du/dy
     call opr_sx_dudxyz(work1%x, u2%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
     call opr_sx_dudxyz(work2%x, u1%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
     call sub3(w3%x, work1%x, work2%x, n)
     !!    BC dependent, Needs to change if cyclic

     call opcolv(w1%x,w2%x,w3%x,c_Xh%B, gdim, n)
     call gs_op(c_Xh%gs_h, w1, GS_OP_ADD) 
     call gs_op(c_Xh%gs_h, w2, GS_OP_ADD) 
     call gs_op(c_Xh%gs_h, w3, GS_OP_ADD) 
     call opcolv  (w1%x,w2%x,w3%x,c_Xh%Binv, gdim, n)
     
   end subroutine opr_sx_curl
  
end module opr_sx
