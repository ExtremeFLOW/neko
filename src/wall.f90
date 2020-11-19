!> Defines wall boundary conditions
module wall
  use num_types
  use dirichlet
  implicit none
  private

  !> No-slip Wall boundary condition
  type, public, extends(dirichlet_t) :: no_slip_wall_t
   contains
     procedure, pass(this) :: apply => no_slip_wall_apply
  end type no_slip_wall_t

contains

  !> Boundary condition apply for a no-slip wall condition
  subroutine no_slip_wall_apply(this, x, n)
    class(no_slip_wall_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout),  dimension(n) :: x
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = 0d0
    end do
  end subroutine no_slip_wall_apply
  
end module wall
