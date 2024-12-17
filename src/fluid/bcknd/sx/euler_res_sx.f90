module euler_res_sx
  use euler_residual, only : euler_rhs_t
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t

  type, public, extends(euler_rhs_t) :: euler_res_sx_t
   contains
     procedure, nopass :: compute => euler_res_sx_compute
  end type euler_res_sx_t

contains
  subroutine euler_res_sx_compute(rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, &
                rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
                c_Xh, gs_Xh)
    type(field_t), intent(inout) :: rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh

  end subroutine euler_res_sx_compute

end module euler_res_sx