!> Defines inflow dirichlet conditions
module inflow
  use device_inflow
  use num_types
  use dirichlet
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  !> Dirichlet condition for inlet (vector valued)
  type, public, extends(dirichlet_t) :: inflow_t
     real(kind=rp), dimension(3) :: x = (/0d0, 0d0, 0d0 /)
   contains
     procedure, pass(this) :: apply_scalar => inflow_apply_scalar
     procedure, pass(this) :: apply_vector => inflow_apply_vector
     procedure, pass(this) :: apply_scalar_dev => inflow_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => inflow_apply_vector_dev
     procedure, pass(this) :: set_inflow => inflow_set_vector
  end type inflow_t
  
contains

  !> No-op scalar apply
  subroutine inflow_apply_scalar(this, x, n)
    class(inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
  end subroutine inflow_apply_scalar

  !> No-op scalar apply (device version)
  subroutine inflow_apply_scalar_dev(this, x_d)
    class(inflow_t), intent(inout), target :: this
    type(c_ptr) :: x_d
  end subroutine inflow_apply_scalar_dev
  
  !> Apply inflow conditions (vector valued)
  subroutine inflow_apply_vector(this, x, y, z, n)
    class(inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = this%x(1)
       y(k) = this%x(2)
       z(k) = this%x(3)
    end do
  end subroutine inflow_apply_vector

  !> Apply inflow conditions (vector valued) (device version)
  subroutine inflow_apply_vector_dev(this, x_d, y_d, z_d)
    class(inflow_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    call device_inflow_apply_vector(this%msk_d, x_d, y_d, z_d, &
                                    c_loc(this%x), size(this%msk))
  end subroutine inflow_apply_vector_dev

  !> Set inflow vector
  subroutine inflow_set_vector(this, x)
    class(inflow_t), intent(inout) :: this
    real(kind=rp), dimension(3), intent(inout) :: x
    this%x = x
  end subroutine inflow_set_vector
    
  
end module inflow
