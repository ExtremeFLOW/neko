!> Residuals in the Pn-Pn formulation (SX version)
module scalar_residual_sx
  use scalar_residual, only : scalar_residual_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use num_types, only : rp    
  use math, only : copy, cfill
  use field, only : field_t
  use mesh, only : mesh_t
  use ax_product, only : ax_t
  use space, only : space_t
  use coefs, only : coef_t
  implicit none
  private

  type, public, extends(scalar_residual_t) :: scalar_residual_sx_t
   contains
     procedure, nopass :: compute => scalar_residual_sx_compute
  end type scalar_residual_sx_t

contains

  subroutine scalar_residual_sx_compute(Ax, s, s_res, f_Xh, c_Xh, msh, Xh, &
      lambda, rhocp, bd, dt, n)
    class(ax_t), intent(in) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(field_t), intent(inout) :: s
    type(field_t), intent(inout) :: s_res
    type(field_t), intent(inout) :: f_Xh
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: lambda
    real(kind=rp), intent(in) :: rhocp
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n
    integer :: i

    call copy(c_Xh%h1, lambda%x, n)
    call cfill(c_Xh%h2, rhocp * bd / dt, n)
    c_Xh%ifh2 = .true.

    call Ax%compute(s_res%x, s%x, c_Xh, msh, Xh)

    do i = 1, n
       s_res%x(i,1,1,1) = (-s_res%x(i,1,1,1)) + f_Xh%x(i,1,1,1)
    end do

  end subroutine scalar_residual_sx_compute

end module scalar_residual_sx
