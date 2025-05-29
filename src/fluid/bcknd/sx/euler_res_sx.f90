module euler_res_sx
  use euler_residual, only : euler_rhs_t
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use num_types, only : rp
  use operators, only: div
  use math, only: subcol3, copy, sub2, add2, add3, &
       col2, col3, addcol3, cmult, cfill, invcol3
  use gs_ops, only : GS_OP_ADD
  use scratch_registry, only: neko_scratch_registry
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  use field_list, only : field_list_t
  implicit none
  private

  type, public, extends(euler_rhs_t) :: euler_res_sx_t
   contains
     procedure, nopass :: step => advance_primitive_variables_sx
     procedure, nopass :: evaluate_rhs_sx
  end type euler_res_sx_t

contains
  subroutine advance_primitive_variables_sx(rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
       coef, gs, h, c_avisc_low, rk_scheme, dt)
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    real(kind=rp) :: c_avisc_low
    class(runge_kutta_time_scheme_t), intent(in) :: rk_scheme
    real(kind=rp), intent(in) :: dt
    integer :: n, s, i, j, k
    real(kind=rp) :: t, c
    type(field_t), pointer :: k_rho_1, k_rho_2, k_rho_3, k_rho_4, &
         k_m_x_1, k_m_x_2, k_m_x_3, k_m_x_4, &
         k_m_y_1, k_m_y_2, k_m_y_3, k_m_y_4, &
         k_m_z_1, k_m_z_2, k_m_z_3, k_m_z_4, &
         k_E_1, k_E_2, k_E_3, k_E_4, &
         temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E
    integer :: temp_indices(25)
    type(field_list_t) :: k_rho, k_m_x, k_m_y, k_m_z, k_E

    n = p%dof%size()
    s = rk_scheme%order
    call neko_scratch_registry%request_field(k_rho_1, temp_indices(1))
    call neko_scratch_registry%request_field(k_rho_2, temp_indices(2))
    call neko_scratch_registry%request_field(k_rho_3, temp_indices(3))
    call neko_scratch_registry%request_field(k_rho_4, temp_indices(4))
    call neko_scratch_registry%request_field(k_m_x_1, temp_indices(5))
    call neko_scratch_registry%request_field(k_m_x_2, temp_indices(6))
    call neko_scratch_registry%request_field(k_m_x_3, temp_indices(7))
    call neko_scratch_registry%request_field(k_m_x_4, temp_indices(8))
    call neko_scratch_registry%request_field(k_m_y_1, temp_indices(9))
    call neko_scratch_registry%request_field(k_m_y_2, temp_indices(10))
    call neko_scratch_registry%request_field(k_m_y_3, temp_indices(11))
    call neko_scratch_registry%request_field(k_m_y_4, temp_indices(12))
    call neko_scratch_registry%request_field(k_m_z_1, temp_indices(13))
    call neko_scratch_registry%request_field(k_m_z_2, temp_indices(14))
    call neko_scratch_registry%request_field(k_m_z_3, temp_indices(15))
    call neko_scratch_registry%request_field(k_m_z_4, temp_indices(16))
    call neko_scratch_registry%request_field(k_E_1, temp_indices(17))
    call neko_scratch_registry%request_field(k_E_2, temp_indices(18))
    call neko_scratch_registry%request_field(k_E_3, temp_indices(19))
    call neko_scratch_registry%request_field(k_E_4, temp_indices(20))
    call neko_scratch_registry%request_field(temp_rho, temp_indices(21))
    call neko_scratch_registry%request_field(temp_m_x, temp_indices(22))
    call neko_scratch_registry%request_field(temp_m_y, temp_indices(23))
    call neko_scratch_registry%request_field(temp_m_z, temp_indices(24))
    call neko_scratch_registry%request_field(temp_E, temp_indices(25))

    call k_rho%init(4)
    call k_rho%assign(1, k_rho_1)
    call k_rho%assign(2, k_rho_2)
    call k_rho%assign(3, k_rho_3)
    call k_rho%assign(4, k_rho_4)
    call k_m_x%init(4)
    call k_m_x%assign(1, k_m_x_1)
    call k_m_x%assign(2, k_m_x_2)
    call k_m_x%assign(3, k_m_x_3)
    call k_m_x%assign(4, k_m_x_4)
    call k_m_y%init(4)
    call k_m_y%assign(1, k_m_y_1)
    call k_m_y%assign(2, k_m_y_2)
    call k_m_y%assign(3, k_m_y_3)
    call k_m_y%assign(4, k_m_y_4)
    call k_m_z%init(4)
    call k_m_z%assign(1, k_m_z_1)
    call k_m_z%assign(2, k_m_z_2)
    call k_m_z%assign(3, k_m_z_3)
    call k_m_z%assign(4, k_m_z_4)
    call k_E%init(4)
    call k_E%assign(1, k_E_1)
    call k_E%assign(2, k_E_2)
    call k_E%assign(3, k_E_3)
    call k_E%assign(4, k_E_4)

    ! Runge-Kutta stages
    do i = 1, s
       call copy(temp_rho%x, rho_field%x, n)
       call copy(temp_m_x%x, m_x%x, n)
       call copy(temp_m_y%x, m_y%x, n)
       call copy(temp_m_z%x, m_z%x, n)
       call copy(temp_E%x, E%x, n)

       do j = 1, i-1
          do concurrent (k = 1:n)
             temp_rho%x(k,1,1,1) = temp_rho%x(k,1,1,1) &
                  + dt * rk_scheme%coeffs_A(i, j) * k_rho%items(j)%ptr%x(k,1,1,1)
             temp_m_x%x(k,1,1,1) = temp_m_x%x(k,1,1,1) &
                  + dt * rk_scheme%coeffs_A(i, j) * k_m_x%items(j)%ptr%x(k,1,1,1)
             temp_m_y%x(k,1,1,1) = temp_m_y%x(k,1,1,1) &
                  + dt * rk_scheme%coeffs_A(i, j) * k_m_y%items(j)%ptr%x(k,1,1,1)
             temp_m_z%x(k,1,1,1) = temp_m_z%x(k,1,1,1) &
                  + dt * rk_scheme%coeffs_A(i, j) * k_m_z%items(j)%ptr%x(k,1,1,1)
             temp_E%x(k,1,1,1) = temp_E%x(k,1,1,1) &
                  + dt * rk_scheme%coeffs_A(i, j) * k_E%items(j)%ptr%x(k,1,1,1)
          end do
       end do

       ! Compute f(U) = rhs(U) for the intermediate values
       call evaluate_rhs_sx(k_rho%items(i)%ptr, k_m_x%items(i)%ptr, &
            k_m_y%items(i)%ptr, k_m_z%items(i)%ptr, &
            k_E%items(i)%ptr, &
            temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E, &
            p, u, v, w, Ax, &
            coef, gs, h, c_avisc_low)
    end do

    ! Update the solution
    do i = 1, s
       do concurrent (k = 1:n)
          rho_field%x(k,1,1,1) = rho_field%x(k,1,1,1) &
               + dt * rk_scheme%coeffs_b(i) * k_rho%items(i)%ptr%x(k,1,1,1)
          m_x%x(k,1,1,1) = m_x%x(k,1,1,1) &
               + dt * rk_scheme%coeffs_b(i) * k_m_x%items(i)%ptr%x(k,1,1,1)
          m_y%x(k,1,1,1) = m_y%x(k,1,1,1) &
               + dt * rk_scheme%coeffs_b(i) * k_m_y%items(i)%ptr%x(k,1,1,1)
          m_z%x(k,1,1,1) = m_z%x(k,1,1,1) &
               + dt * rk_scheme%coeffs_b(i) * k_m_z%items(i)%ptr%x(k,1,1,1)
          E%x(k,1,1,1) = E%x(k,1,1,1) &
               + dt * rk_scheme%coeffs_b(i) * k_E%items(i)%ptr%x(k,1,1,1)
       end do
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine advance_primitive_variables_sx

  subroutine evaluate_rhs_sx(rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, &
       rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
       coef, gs, h, c_avisc_low)
    type(field_t), intent(inout) :: rhs_rho_field, rhs_m_x, &
         rhs_m_y, rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    real(kind=rp) :: c_avisc_low
    integer :: i, n
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

    ! Move div to the rhs and apply the artificial viscosity
    do concurrent (i = 1:n)
       rhs_rho_field%x(i,1,1,1) = -rhs_rho_field%x(i,1,1,1) &
            - c_avisc_low * h%x(i,1,1,1) * coef%Binv(i,1,1,1) * visc_rho%x(i,1,1,1)
       rhs_m_x%x(i,1,1,1) = -rhs_m_x%x(i,1,1,1) &
            - c_avisc_low * h%x(i,1,1,1) * coef%Binv(i,1,1,1) * visc_m_x%x(i,1,1,1)
       rhs_m_y%x(i,1,1,1) = -rhs_m_y%x(i,1,1,1) &
            - c_avisc_low * h%x(i,1,1,1) * coef%Binv(i,1,1,1) * visc_m_y%x(i,1,1,1)
       rhs_m_z%x(i,1,1,1) = -rhs_m_z%x(i,1,1,1) &
            - c_avisc_low * h%x(i,1,1,1) * coef%Binv(i,1,1,1) * visc_m_z%x(i,1,1,1)
       rhs_E%x(i,1,1,1) = -rhs_E%x(i,1,1,1) &
            - c_avisc_low * h%x(i,1,1,1) * coef%Binv(i,1,1,1) * visc_E%x(i,1,1,1)
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine evaluate_rhs_sx

end module euler_res_sx
