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
  use scratch_registry, only: neko_scratch_registry

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
    integer :: n
    real(kind=rp) :: h, c_avisc
    type(field_t), pointer :: temp, f_x, f_y, f_z
    integer :: temp_indices(4)

    h = 0.01_rp / 1.0_rp ! grid size / polynomial degreedm
    c_avisc = 0.5_rp*h
    n = c_Xh%dof%size()
    call neko_scratch_registry%request_field(temp, temp_indices(1))
    call neko_scratch_registry%request_field(f_x, temp_indices(2))
    call neko_scratch_registry%request_field(f_y, temp_indices(3))
    call neko_scratch_registry%request_field(f_z, temp_indices(4))

    !> rho = rho - dt * div(m)
    call div(rhs_rho_field%x, m_x%x, m_y%x, m_z%x, c_Xh)
    ! artificial diffusion for rho
    call Ax%compute(temp%x, rho_field%x, c_Xh, p%msh, p%Xh)
    call gs_Xh%op(temp, GS_OP_ADD)

    do concurrent (i = 1:n)
      rhs_rho_field%x(i,1,1,1) = rhs_rho_field%x(i,1,1,1) &
        + c_avisc * c_Xh%Binv(i,1,1,1) * temp%x(i,1,1,1)
    end do

    ! m = m - dt * div(rho * u * u^T + p*I)
    !> m_x
    do concurrent (i = 1:n)
      f_x%x(i,1,1,1) = m_x%x(i,1,1,1) * u%x(i,1,1,1) &
                        + p%x(i,1,1,1)
      f_y%x(i,1,1,1) = m_x%x(i,1,1,1) * v%x(i,1,1,1)
      f_z%x(i,1,1,1) = m_x%x(i,1,1,1) * w%x(i,1,1,1)
    end do
    call div(rhs_m_x%x, f_x%x, f_y%x, f_z%x, c_Xh)
    ! artificial diffusion for m_x
    call Ax%compute(temp%x, m_x%x, c_Xh, p%msh, p%Xh)
    call gs_Xh%op(temp, GS_OP_ADD)
    do concurrent (i = 1:n)
      rhs_m_x%x(i,1,1,1) = rhs_m_x%x(i,1,1,1) &
        + c_avisc * c_Xh%Binv(i,1,1,1) * temp%x(i,1,1,1)
    end do

    !> m_y
    do concurrent (i = 1:n)
      f_x%x(i,1,1,1) = m_y%x(i,1,1,1) * u%x(i,1,1,1)
      f_y%x(i,1,1,1) = m_y%x(i,1,1,1) * v%x(i,1,1,1) &
                        + p%x(i,1,1,1)
      f_z%x(i,1,1,1) = m_y%x(i,1,1,1) * w%x(i,1,1,1)
    end do
    call div(rhs_m_y%x, f_x%x, f_y%x, f_z%x, c_Xh)
    ! artificial diffusion for m_y
    call Ax%compute(temp%x, m_y%x, c_Xh, p%msh, p%Xh)
    call gs_Xh%op(temp, GS_OP_ADD)
    do concurrent (i = 1:n)
      rhs_m_y%x(i,1,1,1) = rhs_m_y%x(i,1,1,1) &
        + c_avisc * c_Xh%Binv(i,1,1,1) * temp%x(i,1,1,1)
    end do

    !> m_z
    do concurrent (i = 1:n)
      f_x%x(i,1,1,1) = m_z%x(i,1,1,1) * u%x(i,1,1,1)
      f_y%x(i,1,1,1) = m_z%x(i,1,1,1) * v%x(i,1,1,1)
      f_z%x(i,1,1,1) = m_z%x(i,1,1,1) * w%x(i,1,1,1) &
                        + p%x(i,1,1,1)
    end do
    call div(rhs_m_z%x, f_x%x, f_y%x, f_z%x, c_Xh)
    ! artificial diffusion for m_z
    call Ax%compute(temp%x, m_z%x, c_Xh, p%msh, p%Xh)
    call gs_Xh%op(temp, GS_OP_ADD)
    do concurrent (i = 1:n)
      rhs_m_z%x(i,1,1,1) = rhs_m_z%x(i,1,1,1) &
        + c_avisc * c_Xh%Binv(i,1,1,1) * temp%x(i,1,1,1)
    end do

    ! E = E - dt * div(u * (E + p))
    do concurrent (i = 1:n)
      f_x%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) &
                        * u%x(i,1,1,1)
      f_y%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) &
                        * v%x(i,1,1,1)
      f_z%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) &
                        * w%x(i,1,1,1)
    end do
    call div(rhs_E%x, f_x%x, f_y%x, f_z%x, c_Xh)
    ! artificial diffusion for E
    call Ax%compute(temp%x, E%x, c_Xh, p%msh, p%Xh)
    call gs_Xh%op(temp, GS_OP_ADD)
    do concurrent (i = 1:n)
      rhs_E%x(i,1,1,1) = rhs_E%x(i,1,1,1) &
        + c_avisc * c_Xh%Binv(i,1,1,1) * temp%x(i,1,1,1)
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine euler_res_cpu_compute

end module euler_res_cpu