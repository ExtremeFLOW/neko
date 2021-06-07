!> Defines wall boundary conditions
module wall
  use num_types
  use dirichlet
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> No-slip Wall boundary condition
  type, public, extends(dirichlet_t) :: no_slip_wall_t
   contains
     procedure, pass(this) :: apply_scalar => no_slip_wall_apply_scalar
     procedure, pass(this) :: apply_vector => no_slip_wall_apply_vector
     procedure, pass(this) :: apply_scalar_dev => no_slip_wall_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => no_slip_wall_apply_vector_dev
  end type no_slip_wall_t

  interface
     subroutine hip_no_slip_wall_apply_scalar(msk, x, m) &
          bind(c, name='hip_no_slip_wall_apply_scalar')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: msk, x
     end subroutine hip_no_slip_wall_apply_scalar
  end interface
  
  interface
     subroutine hip_no_slip_wall_apply_vector(msk, x, y, z, m) &
          bind(c, name='hip_no_slip_wall_apply_vector')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: msk, x, y, z
     end subroutine hip_no_slip_wall_apply_vector
  end interface

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

  !> Boundary condition apply for a no-slip wall condition
  !! to a vector @a x (device version)
  subroutine no_slip_wall_apply_scalar_dev(this, x_d)
    class(no_slip_wall_t), intent(inout) :: this
    type(c_ptr) :: x_d

    call hip_no_slip_wall_apply_scalar(this%msk_d, x_d, size(this%msk))
    
  end subroutine no_slip_wall_apply_scalar_dev
  
  !> Boundary condition apply for a no-slip wall condition
  !! to vectors @a x, @a y and @a z (device version)
  subroutine no_slip_wall_apply_vector_dev(this, x_d, y_d, z_d)
    class(no_slip_wall_t), intent(inout) :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d

    call hip_no_slip_wall_apply_vector(this%msk_d, x_d, y_d, z_d, size(this%msk))
    
  end subroutine no_slip_wall_apply_vector_dev
  
end module wall
