!> Defines inflow dirichlet conditions
module usr_inflow
  use num_types
  use inflow
  implicit none
  private
  
  !> User defined dirichlet condition for inlet (vector valued)
  type, public, extends(inflow_t) :: usr_inflow_t
   contains
     procedure, pass(this) :: apply_scalar => usr_inflow_apply_scalar
     procedure, pass(this) :: apply_vector => usr_inflow_apply_vector
  end type usr_inflow_t

contains

  !> No-op scalar apply
  subroutine usr_inflow_apply_scalar(this, x, n)
    class(usr_inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout),  dimension(n) :: x
  end subroutine usr_inflow_apply_scalar

  !> Apply inflow conditions (vector valued)
  subroutine usr_inflow_apply_vector(this, x, y, z, n)
    class(usr_inflow_t), intent(inout) :: this
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
  end subroutine usr_inflow_apply_vector
  
end module usr_inflow
