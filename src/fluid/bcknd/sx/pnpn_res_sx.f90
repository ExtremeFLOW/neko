!> Residuals in the Pn-Pn formulation (CPU version)
module pnpn_res_sx
  use gather_scatter
  use pnpn_residual
  use operators
  implicit none
  private
  
  type, public, extends(pnpn_prs_res_t) :: pnpn_prs_res_sx_t
   contains
     procedure, nopass :: compute => pnpn_prs_res_sx_compute
  end type pnpn_prs_res_sx_t

  type, public, extends(pnpn_vel_res_t) :: pnpn_vel_res_sx_t
   contains
     procedure, nopass :: compute => pnpn_vel_res_sx_compute
  end type pnpn_vel_res_sx_t

contains

  subroutine pnpn_prs_res_sx_compute(p, p_res, u, v, w, u_e, v_e, w_e, &
       ta1, ta2, ta3, wa1, wa2, wa3, work1, work2, f_Xh, c_Xh, gs_Xh, &
       bc_prs_surface, bc_sym_surface, Ax, bd, dt, Re, rho)
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(inout) :: u_e, v_e, w_e
    type(field_t), intent(inout) :: ta1, ta2, ta3
    type(field_t), intent(inout) :: wa1, wa2, wa3
    type(field_t), intent(inout) :: work1, work2
    type(field_t), intent(inout) :: p_res
    type(source_t), intent(inout) :: f_Xh
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(facet_normal_t), intent(inout) :: bc_prs_surface
    type(facet_normal_t), intent(inout) :: bc_sym_surface
    class(Ax_t), intent(inout) :: Ax
    real(kind=rp), intent(inout) :: bd
    real(kind=rp), intent(in) :: dt
    real(kind=rp), intent(in) :: Re
    real(kind=rp), intent(in) :: rho
    real(kind=rp) :: dtbd
    integer :: n
    integer :: i

    n = c_Xh%dof%size()

    do i = 1, n
       c_Xh%h1(i,1,1,1) = 1.0_rp / rho
       c_Xh%h2(i,1,1,1) = 0.0_rp
    end do
    c_Xh%ifh2 = .false.
    
    call curl(ta1, ta2, ta3, u_e, v_e, w_e, work1, work2, c_Xh)
    call curl(wa1, wa2, wa3, ta1, ta2, ta3, work1, work2, c_Xh)

    do i = 1, n
       wa1%x(i,1,1,1) = (wa1%x(i,1,1,1) * (1.0_rp / Re) /rho) * c_Xh%B(i,1,1,1)
       wa2%x(i,1,1,1) = (wa2%x(i,1,1,1) * (1.0_rp / Re) /rho) * c_Xh%B(i,1,1,1)
       wa3%x(i,1,1,1) = (wa3%x(i,1,1,1) * (1.0_rp / Re) /rho) * c_Xh%B(i,1,1,1)
    end do

    do i = 1, n
       ta1%x(i,1,1,1) = f_Xh%u(i,1,1,1) / rho - wa1%x(i,1,1,1)
       ta2%x(i,1,1,1) = f_Xh%v(i,1,1,1) / rho - wa2%x(i,1,1,1)
       ta3%x(i,1,1,1) = f_Xh%w(i,1,1,1) / rho - wa3%x(i,1,1,1)
    end do
     
    call gs_op(gs_Xh, ta1, GS_OP_ADD) 
    call gs_op(gs_Xh, ta2, GS_OP_ADD) 
    call gs_op(gs_Xh, ta3, GS_OP_ADD) 

    do i = 1, n
       ta1%x(i,1,1,1) = ta1%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
       ta2%x(i,1,1,1) = ta2%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
       ta3%x(i,1,1,1) = ta3%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
    end do
    
    call Ax%compute(p_res%x,p%x,c_Xh,p%msh,p%Xh)

    call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call cdtp(wa3%x, ta3%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
    do i = 1, n
       p_res%x(i,1,1,1) = (-p_res%x(i,1,1,1)) &
                        + wa1%x(i,1,1,1) + wa2%x(i,1,1,1) + wa3%x(i,1,1,1)
    end do

    !
    ! Surface velocity terms
    !
    do i = 1, n
       wa1%x(i,1,1,1) = 0.0_rp
       wa2%x(i,1,1,1) = 0.0_rp
       wa3%x(i,1,1,1) = 0.0_rp
    end do
    
    call bc_sym_surface%apply_surfvec(wa1%x,wa2%x,wa3%x,ta1%x, ta2%x, ta3%x, n)
    
    dtbd = bd / dt
    do i = 1, n
       ta1%x(i,1,1,1) = 0.0_rp
       ta2%x(i,1,1,1) = 0.0_rp
       ta3%x(i,1,1,1) = 0.0_rp
    end do
    
    call bc_prs_surface%apply_surfvec(ta1%x, ta2%x, ta3%x, u%x, v%x, w%x, n)

    do i = 1, n
       p_res%x(i,1,1,1) = p_res%x(i,1,1,1) &
            - (dtbd * (ta1%x(i,1,1,1) + ta2%x(i,1,1,1) + ta3%x(i,1,1,1)))&
            - (wa1%x(i,1,1,1) + wa2%x(i,1,1,1) + wa3%x(i,1,1,1))
    end do

  end subroutine pnpn_prs_res_sx_compute

  subroutine pnpn_vel_res_sx_compute(Ax, u, v, w, u_res, v_res, w_res, &
       p, ta1, ta2, ta3, f_Xh, c_Xh, msh, Xh, Re, rho, bd, dt, n)
    class(ax_t), intent(in) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh    
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(inout) :: u_res, v_res, w_res
    type(field_t), intent(inout) :: ta1, ta2, ta3
    type(source_t), intent(inout) :: f_Xh
    type(coef_t), intent(inout) :: c_Xh
    real(kind=rp), intent(in) :: Re
    real(kind=rp), intent(in) :: rho
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n
    integer :: i

    do i = 1, n
       c_Xh%h1(i,1,1,1) = (1.0_rp / Re)
       c_Xh%h2(i,1,1,1) = rho * (bd / dt)
    end do
    c_Xh%ifh2 = .true.

    call Ax%compute(u_res%x, u%x, c_Xh, msh, Xh)
    call Ax%compute(v_res%x, v%x, c_Xh, msh, Xh)
    call Ax%compute(w_res%x, w%x, c_Xh, msh, Xh)

    call opgrad(ta1%x, ta2%x, ta3%x, p%x, c_Xh)

    do i = 1, n
       u_res%x(i,1,1,1) = (-u_res%x(i,1,1,1)) - ta1%x(i,1,1,1) + f_Xh%u(i,1,1,1)
       v_res%x(i,1,1,1) = (-v_res%x(i,1,1,1)) - ta2%x(i,1,1,1) + f_Xh%v(i,1,1,1)
       w_res%x(i,1,1,1) = (-w_res%x(i,1,1,1)) - ta3%x(i,1,1,1) + f_Xh%w(i,1,1,1)
    end do
    
  end subroutine pnpn_vel_res_sx_compute  
     
end module pnpn_res_sx
