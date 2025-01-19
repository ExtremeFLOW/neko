module euler_residual
  use gather_scatter, only : gs_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use facet_normal, only : facet_normal_t
  use space, only : space_t
  use mesh, only : mesh_t
  use num_types, only : rp
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  implicit none
  private

  !> Abstract type to compute rhs
  type, public, abstract :: euler_rhs_t
   contains
     procedure(euler_rhs), nopass, deferred :: step
  end type euler_rhs_t

  abstract interface
     subroutine euler_rhs(rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
                          coef, gs, h, c_avisc_low, rk_scheme, dt)
       import field_t
       import Ax_t
       import gs_t
       import coef_t
       import rp
       import runge_kutta_time_scheme_t
       type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
       type(field_t), intent(in) :: p, u, v, w, h
       class(Ax_t), intent(inout) :: Ax
       type(coef_t), intent(inout) :: coef
       type(gs_t), intent(inout) :: gs
       real(kind=rp) :: c_avisc_low
       class(runge_kutta_time_scheme_t), intent(in) :: rk_scheme
       real(kind=rp), intent(in) :: dt
     end subroutine euler_rhs
  end interface

  interface
     module subroutine euler_rhs_factory(object)
       class(euler_rhs_t), allocatable, intent(inout) :: object
     end subroutine euler_rhs_factory
  end interface

  public :: euler_rhs_factory

end module euler_residual