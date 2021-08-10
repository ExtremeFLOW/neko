!> Operators accelerator backends
module opr_device
  use num_types
  use device    
  use space
  use coefs
  use math
  use mesh
  use field
  use gather_scatter
  use mathops
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  interface
     subroutine hip_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
          dx_d, dy_d, dz_d, jacinv_d, nel, lx) &
          bind(c, name='hip_dudxyz')
          use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, jacinv_d
       integer(c_int) :: nel, lx
     end subroutine hip_dudxyz
  end interface
#endif
  
contains

  subroutine opr_device_dudxyz(du, u, dr, ds, dt, coef)
    type(coef_t), intent(in), target :: coef
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly, &
         coef%Xh%lz,coef%msh%nelv), intent(inout) ::  du
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly, &
         coef%Xh%lz,coef%msh%nelv), intent(in) ::  u, dr, ds, dt
    type(c_ptr) :: du_d, u_d, dr_d, ds_d, dt_d

    du_d = device_get_ptr(du, size(du))
    u_d = device_get_ptr(u, size(u))

    dr_d = device_get_ptr(dr, size(dr))
    ds_d = device_get_ptr(ds, size(ds))
    dt_d = device_get_ptr(dt, size(dt))

    associate(Xh => coef%Xh, msh => coef%msh, dof => coef%dof)    
#ifdef HAVE_HIP
      call hip_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
           Xh%dx_d, Xh%dy_d, Xh%dz_d, coef%jacinv_d, &
           msh%nelv, Xh%lx)
#else
      call neko_error('No device backend configured')
#endif
    end associate
  
  end subroutine opr_device_dudxyz

  subroutine opr_device_opgrad(ux,uy,uz,u,coef) 
    type(coef_t), intent(in) :: coef  
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: ux
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uy
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uz
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: u

#ifdef HAVE_HIP
#else
    call neko_error('No device backend configured')
#endif
    
  end subroutine opr_device_opgrad


  subroutine opr_device_cdtp(dtx,x,dr,ds,dt, coef)
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dtx
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: x
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: dr
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: ds
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: dt

#ifdef HAVE_HIP
#else
    call neko_error('No device backend configured')
#endif

  end subroutine opr_device_cdtp

  subroutine opr_device_conv1(du,u, vx, vy, vz, Xh, coef, nelv, gdim)  
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: nelv, gdim
    real(kind=rp), intent(inout) ::  du(Xh%lxyz,nelv)
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  u
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vx
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vy
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vz

#ifdef HAVE_HIP
#else
    call neko_error('No device backend configured')
#endif
    
  end subroutine opr_device_conv1

  subroutine opr_device_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    type(field_t), intent(inout) :: w1
    type(field_t), intent(inout) :: w2
    type(field_t), intent(inout) :: w3
    type(field_t), intent(inout) :: u1
    type(field_t), intent(inout) :: u2
    type(field_t), intent(inout) :: u3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(coef_t), intent(in)  :: c_Xh
    integer :: gdim, n

    n = w1%dof%size()
    gdim = c_Xh%msh%gdim

    !     this%work1=dw/dy ; this%work2=dv/dz
    call opr_device_dudxyz(work1%x, u3%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    if (gdim .eq. 3) then
       call opr_device_dudxyz(work2%x, u2%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       call sub3(w1%x, work1%x, work2%x, n)
    else
       call copy(w1%x, work1%x, n)
    endif
    !     this%work1=du/dz ; this%work2=dw/dx
    if (gdim .eq. 3) then
       call opr_device_dudxyz(work1%x, u1%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       call opr_device_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    else
       call rzero (work1%x, n)
       call opr_device_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    endif
    !     this%work1=dv/dx ; this%work2=du/dy
    call opr_device_dudxyz(work1%x, u2%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call opr_device_dudxyz(work2%x, u1%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call sub3(w3%x, work1%x, work2%x, n)
    !!    BC dependent, Needs to change if cyclic

    call opcolv(w1%x,w2%x,w3%x,c_Xh%B, gdim, n)
    call gs_op(c_Xh%gs_h, w1, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w2, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w3, GS_OP_ADD) 
    call opcolv  (w1%x,w2%x,w3%x,c_Xh%Binv, gdim, n)

  end subroutine opr_device_curl



end module opr_device
