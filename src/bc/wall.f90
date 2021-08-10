!> Defines wall boundary conditions
module wall
  use num_types
  use dirichlet
  implicit none
  private

  !> No-slip Wall boundary condition
  type, public, extends(dirichlet_t) :: no_slip_wall_t
   contains
     procedure, pass(this) :: apply_scalar => no_slip_wall_apply_scalar
     procedure, pass(this) :: apply_vector => no_slip_wall_apply_vector
  end type no_slip_wall_t

contains

  !> Boundary condition apply for a no-slip wall condition
  !! to a vector @a x
  subroutine no_slip_wall_apply_scalar(this, x, n)
    class(no_slip_wall_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = 0d0
    end do
    
  end subroutine no_slip_wall_apply_scalar
  
  !> Boundary condition apply for a no-slip wall condition
  !! to vectors @a x, @a y and @a z
  subroutine no_slip_wall_apply_vector(this, x, y, z, n)
    class(no_slip_wall_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = 0d0
       y(k) = 0d0
       z(k) = 0d0
    end do
    
  end subroutine no_slip_wall_apply_vector
  
end module wall
