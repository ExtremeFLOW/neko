!> Defines a dirichlet boundary condition
module dirichlet
  use num_types
  use bc  
  implicit none
  private

  !> Generic Dirichlet boundary condition
  !! \f$ x = g \f$ on \f$\partial \Omega\f$
  type, public, extends(bc_t) :: dirichlet_t
     real(kind=rp), private :: g
   contains
     procedure, pass(this) :: apply_scalar => dirichlet_apply_scalar
     procedure, pass(this) :: apply_vector => dirichlet_apply_vector
     procedure, pass(this) :: set_g => dirichlet_set_g
  end type dirichlet_t
     
contains

  !> Boundary condition apply for a generic Dirichlet condition
  !! to a vector @a x
  subroutine dirichlet_apply_scalar(this, x, n)
    class(dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    integer :: i, m, k

    m = this%msk(0)
    print *, 'yes', this%g
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

  !> Set value of \f$ g \f$
  subroutine dirichlet_set_g(this, g)
    class(dirichlet_t), intent(inout) :: this
    real(kind=rp), intent(in) :: g

    this%g = g
    
  end subroutine dirichlet_set_g
  
end module dirichlet
