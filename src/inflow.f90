  !> Defines inflow dirichlet conditions
module inflow
  use num_types
  use dirichlet
  implicit none
  private
  
  !> Dirichlet condition for inlet (vector valued)
  type, public, extends(dirichlet_t) :: inflow_t
     real(kind=dp), dimension(3) :: x = (/0d0, 0d0, 0d0 /)
   contains
     procedure, pass(this) :: apply_scalar => inflow_apply_scalar
     procedure, pass(this) :: apply_vector => inflow_apply_vector
     procedure, pass(this) :: set_inflow => inflow_set_vector
  end type inflow_t

contains

  !> No-op scalar apply
  subroutine inflow_apply_scalar(this, x, n)
    class(inflow_t) , intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout),  dimension(n) :: x
  end subroutine inflow_apply_scalar

  !> Apply inflow conditions (vector valued)
  subroutine inflow_apply_vector(this, x, y, z, n)
    class(inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout),  dimension(n) :: x
    real(kind=dp), intent(inout),  dimension(n) :: y
    real(kind=dp), intent(inout),  dimension(n) :: z
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = this%x(1)
       y(k) = this%x(2)
       z(k) = this%x(3)
    end do
  end subroutine inflow_apply_vector

  !> Set inflow vector
  subroutine inflow_set_vector(this, x)
    class(inflow_t), intent(inout) :: this
    real(kind=dp), dimension(3), intent(inout) :: x
    this%x = x
  end subroutine inflow_set_vector
    
  
end module inflow
