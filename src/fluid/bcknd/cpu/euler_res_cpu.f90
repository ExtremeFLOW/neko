module euler_res_cpu
  use euler_residual, only : euler_rhs_t
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use num_types, only : rp
  use operators, only: div
  use math, only: subcol3, copy, sub2, add2, add3, col2, col3, addcol3, cmult, cfill, invcol3
  use gs_ops, only : GS_OP_ADD
  use scratch_registry, only: neko_scratch_registry

  type, public, extends(euler_rhs_t) :: euler_res_cpu_t
   contains
     procedure, nopass :: compute => euler_res_cpu_compute
  end type euler_res_cpu_t

contains
  subroutine euler_res_cpu_compute(rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, &
                rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
                coef, gs, h, c_avisc_low)
    type(field_t), intent(inout) :: rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    real(kind=rp) :: c_avisc_low
    integer :: n
    type(field_t), pointer :: temp, f_x, f_y, f_z, &
                              visc_rho, visc_m_x, visc_m_y, visc_m_z, visc_E
    integer :: temp_indices(9)

    n = coef%dof%size()
    call neko_scratch_registry%request_field(temp, temp_indices(1))
    call neko_scratch_registry%request_field(f_x, temp_indices(2))
    call neko_scratch_registry%request_field(f_y, temp_indices(3))
    call neko_scratch_registry%request_field(f_z, temp_indices(4))

    !> rho = rho - dt * div(m)
    call div(rhs_rho_field%x, m_x%x, m_y%x, m_z%x, coef)

    !> m = m - dt * div(rho * u * u^T + p*I)
    ! m_x
    do concurrent (i = 1:n)
      f_x%x(i,1,1,1) = m_x%x(i,1,1,1) * m_x%x(i,1,1,1) / rho_field%x(i, 1, 1, 1) &
                        + p%x(i,1,1,1)
      f_y%x(i,1,1,1) = m_x%x(i,1,1,1) * m_y%x(i,1,1,1) / rho_field%x(i, 1, 1, 1)
      f_z%x(i,1,1,1) = m_x%x(i,1,1,1) * m_z%x(i,1,1,1) / rho_field%x(i, 1, 1, 1)
    end do
    call div(rhs_m_x%x, f_x%x, f_y%x, f_z%x, coef)
    ! m_y
    do concurrent (i = 1:n)
      f_x%x(i,1,1,1) = m_y%x(i,1,1,1) * m_x%x(i,1,1,1) / rho_field%x(i, 1, 1, 1)
      f_y%x(i,1,1,1) = m_y%x(i,1,1,1) * m_y%x(i,1,1,1) / rho_field%x(i, 1, 1, 1) &
                        + p%x(i,1,1,1)
      f_z%x(i,1,1,1) = m_y%x(i,1,1,1) * m_z%x(i,1,1,1) / rho_field%x(i, 1, 1, 1)
    end do
    call div(rhs_m_y%x, f_x%x, f_y%x, f_z%x, coef)
    ! m_z
    do concurrent (i = 1:n)
      f_x%x(i,1,1,1) = m_z%x(i,1,1,1) * m_x%x(i,1,1,1) / rho_field%x(i, 1, 1, 1)
      f_y%x(i,1,1,1) = m_z%x(i,1,1,1) * m_y%x(i,1,1,1) / rho_field%x(i, 1, 1, 1)
      f_z%x(i,1,1,1) = m_z%x(i,1,1,1) * m_z%x(i,1,1,1) / rho_field%x(i, 1, 1, 1) &
                        + p%x(i,1,1,1)
    end do
    call div(rhs_m_z%x, f_x%x, f_y%x, f_z%x, coef)

    !> E = E - dt * div(u * (E + p))
    do concurrent (i = 1:n)
      f_x%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) &
                        * u%x(i,1,1,1)
      f_y%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) &
                        * v%x(i,1,1,1)
      f_z%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) &
                        * w%x(i,1,1,1)
    end do
    call div(rhs_E%x, f_x%x, f_y%x, f_z%x, coef)

    call gs%op(rhs_rho_field, GS_OP_ADD)
    call gs%op(rhs_m_x, GS_OP_ADD)
    call gs%op(rhs_m_y, GS_OP_ADD)
    call gs%op(rhs_m_z, GS_OP_ADD)
    call gs%op(rhs_E, GS_OP_ADD)

    do concurrent (i = 1:rhs_E%dof%size())
      rhs_rho_field%x(i,1,1,1) = rhs_rho_field%x(i,1,1,1) * coef%mult(i,1,1,1)
      rhs_m_x%x(i,1,1,1) = rhs_m_x%x(i,1,1,1) * coef%mult(i,1,1,1)
      rhs_m_y%x(i,1,1,1) = rhs_m_y%x(i,1,1,1) * coef%mult(i,1,1,1)
      rhs_m_z%x(i,1,1,1) = rhs_m_z%x(i,1,1,1) * coef%mult(i,1,1,1)
      rhs_E%x(i,1,1,1) = rhs_E%x(i,1,1,1) * coef%mult(i,1,1,1)
    end do

    call neko_scratch_registry%request_field(visc_rho, temp_indices(5))
    call neko_scratch_registry%request_field(visc_m_x, temp_indices(6))
    call neko_scratch_registry%request_field(visc_m_y, temp_indices(7))
    call neko_scratch_registry%request_field(visc_m_z, temp_indices(8))
    call neko_scratch_registry%request_field(visc_E, temp_indices(9))

    ! Calculate artificial diffusion
    call Ax%compute(visc_rho%x, rho_field%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_m_x%x, m_x%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_m_y%x, m_y%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_m_z%x, m_z%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_E%x, E%x, coef, p%msh, p%Xh)

    call gs%op(visc_rho, GS_OP_ADD)
    call gs%op(visc_m_x, GS_OP_ADD)
    call gs%op(visc_m_y, GS_OP_ADD)
    call gs%op(visc_m_z, GS_OP_ADD)
    call gs%op(visc_E, GS_OP_ADD)

    do concurrent (i = 1:n)
      rhs_rho_field%x(i,1,1,1) = rhs_rho_field%x(i,1,1,1) &
        + c_avisc_low * h%x(i,1,1,1) * coef%Binv(i,1,1,1) * visc_rho%x(i,1,1,1)
      rhs_m_x%x(i,1,1,1) = rhs_m_x%x(i,1,1,1) &
        + c_avisc_low * h%x(i,1,1,1) * coef%Binv(i,1,1,1) * visc_m_x%x(i,1,1,1)
      rhs_m_y%x(i,1,1,1) = rhs_m_y%x(i,1,1,1) &
        + c_avisc_low * h%x(i,1,1,1) * coef%Binv(i,1,1,1) * visc_m_y%x(i,1,1,1)
      rhs_m_z%x(i,1,1,1) = rhs_m_z%x(i,1,1,1) &
        + c_avisc_low * h%x(i,1,1,1) * coef%Binv(i,1,1,1) * visc_m_z%x(i,1,1,1)
      rhs_E%x(i,1,1,1) = rhs_E%x(i,1,1,1) &
        + c_avisc_low * h%x(i,1,1,1) * coef%Binv(i,1,1,1) * visc_E%x(i,1,1,1)
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine euler_res_cpu_compute

end module euler_res_cpu