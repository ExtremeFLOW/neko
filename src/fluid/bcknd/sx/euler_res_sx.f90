module euler_res_sx
  use euler_residual, only : euler_rhs_t
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use num_types, only : rp
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  use utils, only : neko_error

  type, public, extends(euler_rhs_t) :: euler_res_sx_t
   contains
     procedure, nopass :: step => advance_primitive_variables_sx
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
    call neko_error("SX not supported")
  end subroutine advance_primitive_variables_sx

end module euler_res_sx