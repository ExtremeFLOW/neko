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
     procedure, pass(this) :: apply_mult => no_slip_wall_apply_mult
  end type no_slip_wall_t

contains

  !> Boundary condition apply for a no-slip wall condition
  !! to a vector @a x
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

  !> Boundary condition apply for a no-slip wall condition
  !! to vectors @a x, @a y and @a z
  subroutine no_slip_wall_apply_mult(this, x, y, z, n)
    class(no_slip_wall_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout),  dimension(n) :: x
    real(kind=dp), intent(inout),  dimension(n) :: y
    real(kind=dp), intent(inout),  dimension(n) :: z
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = 0d0
       y(k) = 0d0
       z(k) = 0d0
    end do
    
  end subroutine no_slip_wall_apply_mult
  
end module wall
