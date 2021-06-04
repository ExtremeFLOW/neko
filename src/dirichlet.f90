!> Defines a dirichlet boundary condition
module dirichlet
  use num_types
  use utils
  use bc
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Generic Dirichlet boundary condition
  !! \f$ x = g \f$ on \f$\partial \Omega\f$
  type, public, extends(bc_t) :: dirichlet_t
     real(kind=rp), private :: g
   contains
     procedure, pass(this) :: apply_scalar => dirichlet_apply_scalar
     procedure, pass(this) :: apply_vector => dirichlet_apply_vector
     procedure, pass(this) :: apply_scalar_dev => dirichlet_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => dirichlet_apply_vector_dev
     procedure, pass(this) :: set_g => dirichlet_set_g
  end type dirichlet_t

  interface
     subroutine hip_dirichlet_apply_scalar(msk, x, g, m) &
          bind(c, name='hip_dirichlet_apply_scalar')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double) :: g
       integer(c_int) :: m
       type(c_ptr), value :: msk, x
     end subroutine hip_dirichlet_apply_scalar
  end interface
  
  interface
     subroutine hip_dirichlet_apply_vector(msk, x, y, z, g, m) &
          bind(c, name='hip_dirichlet_apply_vector')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double) :: g
       integer(c_int) :: m
       type(c_ptr), value :: msk, x, y, z
     end subroutine hip_dirichlet_apply_vector
  end interface
     
contains

  !> Boundary condition apply for a generic Dirichlet condition
  !! to a vector @a x
  subroutine dirichlet_apply_scalar(this, x, n)
    class(dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = this%g
    end do
  end subroutine dirichlet_apply_scalar

  !> Boundary condition apply for a generic Dirichlet condition
  !! to vectors @a x, @a y and @a z
  subroutine dirichlet_apply_vector(this, x, y, z, n)
    class(dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = this%g
       y(k) = this%g
       z(k) = this%g
    end do
    
  end subroutine dirichlet_apply_vector

  !> Boundary condition apply for a generic Dirichlet condition
  !! to a vector @a x (device version)
  subroutine dirichlet_apply_scalar_dev(this, x_d)
    class(dirichlet_t), intent(inout) :: this
    type(c_ptr) :: x_d

    call hip_dirichlet_apply_scalar(this%msk_d, x_d, &
                                    this%g, this%msk(0))
    
  end subroutine dirichlet_apply_scalar_dev

  !> Boundary condition apply for a generic Dirichlet condition 
  !! to vectors @a x, @a y and @a z (device version)
  subroutine dirichlet_apply_vector_dev(this, x_d, y_d, z_d)
    class(dirichlet_t), intent(inout) :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d

    call hip_dirichlet_apply_vector(this%msk_d, x_d, y_d, z_d, &
                                    this%g, this%msk(0))
    
  end subroutine dirichlet_apply_vector_dev

  !> Set value of \f$ g \f$
  subroutine dirichlet_set_g(this, g)
    class(dirichlet_t), intent(inout) :: this
    real(kind=rp), intent(in) :: g

    this%g = g
    
  end subroutine dirichlet_set_g
  
end module dirichlet
