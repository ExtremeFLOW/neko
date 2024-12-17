module euler_res_cpu
  use euler_residual, only : euler_rhs_t
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use num_types, only : rp
  use operators, only: div
  use math, only: subcol3, copy, sub2, add2, add3, col2, col3, addcol3, cmult, cfill, invcol3
  use gs_ops, only : GS_OP_ADD, GS_OP_MIN_ABS

  type, public, extends(euler_rhs_t) :: euler_res_cpu_t
   contains
     procedure, nopass :: compute => euler_res_cpu_compute
  end type euler_res_cpu_t

contains
  subroutine euler_res_cpu_compute(rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, &
                rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
                c_Xh, gs_Xh)
    type(field_t), intent(inout) :: rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    real(kind=rp), allocatable :: temp(:), f_x(:), f_y(:), f_z(:)
    integer :: n
    real(kind=rp) :: h, c_avisc

    h = 0.01_rp / 1.0_rp ! grid size / polynomial degreedm
    c_avisc = 0.5_rp*h
    n = c_Xh%dof%size()
    allocate(temp(n))
    allocate(f_x(n))
    allocate(f_y(n))
    allocate(f_z(n))

    !> rho = rho - dt * div(m)
    call div(rhs_rho_field%x, m_x%x, m_y%x, m_z%x, c_Xh)
    ! artificial diffusion for rho
    call Ax%compute(temp, rho_field%x, c_Xh, p%msh, p%Xh)
    call gs_Xh%op(temp, n, GS_OP_ADD)
    call col2(temp, c_Xh%Binv, n)
    call cmult(temp, c_avisc, n) ! first-order viscosity
    call add2(rhs_rho_field%x, temp, n)

    ! m = m - dt * div(rho * u * u^T + p*I)
    !> m_x
    call copy(f_x, p%x, n)
    call addcol3(f_x, m_x%x, u%x, n)
    call col3(f_y, m_x%x, v%x, n)
    call col3(f_z, m_x%x, w%x, n)
    call div(rhs_m_x%x, f_x, f_y, f_z, c_Xh)
    ! artificial diffusion for m_x
    call Ax%compute(temp, m_x%x, c_Xh, p%msh, p%Xh)
    call gs_Xh%op(temp, n, GS_OP_ADD)
    call col2(temp, c_Xh%Binv, n)
    call cmult(temp, c_avisc, n) ! first-order viscosity
    call add2(rhs_m_x%x, temp, n)

    !> m_y
    call col3(f_x, m_y%x, u%x, n)
    call copy(f_y, p%x, n)
    call addcol3(f_y, m_y%x, v%x, n)
    call col3(f_z, m_y%x, w%x, n)
    call div(rhs_m_y%x, f_x, f_y, f_z, c_Xh)
    ! artificial diffusion for m_y
    call Ax%compute(temp, m_y%x, c_Xh, p%msh, p%Xh)
    call gs_Xh%op(temp, n, GS_OP_ADD)
    call col2(temp, c_Xh%Binv, n)
    call cmult(temp, c_avisc, n) ! first-order viscosity
    call add2(rhs_m_y%x, temp, n)

    !> m_z
    call col3(f_x, m_z%x, u%x, n)
    call col3(f_y, m_z%x, v%x, n)
    call copy(f_z, p%x, n)
    call addcol3(f_z, m_z%x, w%x, n)
    call div(rhs_m_z%x, f_x, f_y, f_z, c_Xh)
    ! artificial diffusion for m_z
    call Ax%compute(temp, m_z%x, c_Xh, p%msh, p%Xh)
    call gs_Xh%op(temp, n, GS_OP_ADD)
    call col2(temp, c_Xh%Binv, n)
    call cmult(temp, c_avisc, n) ! first-order viscosity
    call add2(rhs_m_z%x, temp, n)

    ! E = E - dt * div(u * (E + p))
    call add3(temp, E%x, p%x, n)
    call col3(f_x, u%x, temp, n)
    call col3(f_y, v%x, temp, n)
    call col3(f_z, w%x, temp, n)
    call div(rhs_E%x, f_x, f_y, f_z, c_Xh)
    ! artificial diffusion for E
    call Ax%compute(temp, E%x, c_Xh, p%msh, p%Xh)
    call gs_Xh%op(temp, n, GS_OP_ADD)
    call col2(temp, c_Xh%Binv, n)
    call cmult(temp, c_avisc, n) ! first-order viscosity
    call add2(rhs_E%x, temp, n)

    ! call gs_Xh%op(rhs_rho_field%x, n, GS_OP_MIN_ABS)
    ! call gs_Xh%op(rhs_m_x%x, n, GS_OP_MIN_ABS)
    ! call gs_Xh%op(rhs_m_y%x, n, GS_OP_MIN_ABS)
    ! call gs_Xh%op(rhs_m_z%x, n, GS_OP_MIN_ABS)
    ! call gs_Xh%op(rhs_E%x, n, GS_OP_MIN_ABS)

    deallocate(temp)
    deallocate(f_x)
    deallocate(f_y)
    deallocate(f_z)

  end subroutine euler_res_cpu_compute

end module euler_res_cpu