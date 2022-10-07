!> Residuals in the Pn-Pn formulation (SX version)
module scalar_residual_sx
  use gather_scatter
  use scalar_residual
  use operators
  implicit none
  private
  
  type, public, extends(scalar_residual_t) :: scalar_residual_sx_t
   contains
     procedure, nopass :: compute => scalar_residual_sx_compute
  end type scalar_residual_sx_t

contains

  subroutine scalar_residual_sx_compute(Ax, s, s_res, f_Xh, c_Xh, msh, Xh, &
      Pr, Re, rho, bd, dt, n)
    class(ax_t), intent(in) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh    
    type(field_t), intent(inout) :: s
    type(field_t), intent(inout) :: s_res
    type(source_scalar_t), intent(inout) :: f_Xh
    type(coef_t), intent(inout) :: c_Xh
    real(kind=rp), intent(in) :: Pr
    real(kind=rp), intent(in) :: Re
    real(kind=rp), intent(in) :: rho
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n
    integer :: i

    do i = 1, n
       c_Xh%h1(i,1,1,1) = (1.0_rp / (Pr * Re))
       ! todo :should not be just rho here.
       c_Xh%h2(i,1,1,1) = rho * (bd / dt)
    end do
    c_Xh%ifh2 = .true.

    call Ax%compute(s_res%x, s%x, c_Xh, msh, Xh)

    do i = 1, n
       s_res%x(i,1,1,1) = (-s_res%x(i,1,1,1)) + f_Xh%s(i,1,1,1)
    end do
    
  end subroutine scalar_residual_sx_compute  
     
end module scalar_residual_sx
