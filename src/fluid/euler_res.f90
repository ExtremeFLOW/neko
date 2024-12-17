module euler_residual
  use gather_scatter, only : gs_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use facet_normal, only : facet_normal_t
  use space, only : space_t
  use mesh, only : mesh_t
  use num_types, only : rp
  implicit none
  private

  !> Abstract type to compute rhs
  type, public, abstract :: euler_rhs_t
   contains
     procedure(euler_rhs), nopass, deferred :: compute
  end type euler_rhs_t

  abstract interface
     subroutine euler_rhs(rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, &
                          rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
                          c_Xh, gs_Xh)

       import field_t
       import Ax_t
       import gs_t
       import coef_t
       import rp
       type(field_t), intent(inout) :: rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
       type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
       type(field_t), intent(in) :: p, u, v, w
       class(Ax_t), intent(inout) :: Ax
       type(coef_t), intent(inout) :: c_Xh
       type(gs_t), intent(inout) :: gs_Xh
     end subroutine euler_rhs
  end interface

  interface
     module subroutine euler_rhs_factory(object)
       class(euler_rhs_t), allocatable, intent(inout) :: object
     end subroutine euler_rhs_factory
  end interface

  public :: euler_rhs_factory

end module euler_residual