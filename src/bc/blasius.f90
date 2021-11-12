 !> Defines a Blasius profile dirichlet condition
module blasius
  use num_types
  use coefs
  use utils
  use inflow
  use flow_profile
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Blasius profile for inlet (vector valued)
  type, public, extends(inflow_t) :: blasius_t
     type(coef_t), pointer :: c => null()
     real(kind=rp) :: delta
     procedure(blasius_profile), nopass, pointer :: bla => null()
   contains
     procedure, pass(this) :: apply_scalar => blasius_apply_scalar
     procedure, pass(this) :: apply_vector => blasius_apply_vector
     procedure, pass(this) :: apply_scalar_dev => blasius_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => blasius_apply_vector_dev
     procedure, pass(this) :: set_params => blasius_set_params
     procedure, pass(this) :: set_coef => blasius_set_coef
  end type blasius_t

contains
  
  !> No-op scalar apply
  subroutine blasius_apply_scalar(this, x, n)
    class(blasius_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
  end subroutine blasius_apply_scalar

  !> No-op scalar apply (device version)
  subroutine blasius_apply_scalar_dev(this, x_d)
    class(blasius_t), intent(inout), target :: this
    type(c_ptr) :: x_d
  end subroutine blasius_apply_scalar_dev
  
  !> Apply blasius conditions (vector valued)
  subroutine blasius_apply_vector(this, x, y, z, n)
    class(blasius_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    integer :: i, m, k, idx(4), facet

    associate(xc => this%c%dof%x, yc => this%c%dof%y, zc => this%c%dof%z, &
         nx => this%c%nx, ny => this%c%ny, nz => this%c%nz, &
         lx => this%c%Xh%lx)
      m = this%msk(0)
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, lx, lx, lx)
         select case(facet)
         case(1,2)
            x(k) = this%bla(zc(idx(1), idx(2), idx(3), idx(4)), &
                 this%delta, this%x(1))
            y(k) = 0.0_rp 
            z(k) = 0.0_rp
         case(3,4)
            x(k) = 0.0_rp 
            y(k) = this%bla(xc(idx(1), idx(2), idx(3), idx(4)), &
                 this%delta, this%x(2))
            z(k) = 0.0_rp
         case(5,6)
            x(k) = 0.0_rp
            y(k) = 0.0_rp
            z(k) = this%bla(yc(idx(1), idx(2), idx(3), idx(4)), &
                 this%delta, this%x(3))
         end select            
      end do
    end associate
  end subroutine blasius_apply_vector

  !> Apply blasius conditions (vector valued) (device version)
  subroutine blasius_apply_vector_dev(this, x_d, y_d, z_d)
    class(blasius_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
  end subroutine blasius_apply_vector_dev

  !> Set Blasius parameters
  subroutine blasius_set_params(this, delta)
    class(blasius_t), intent(inout) :: this
    real(kind=rp) :: delta
    this%delta = delta
    this%bla => blasius_sin
  end subroutine blasius_set_params

  !> Assign coefficients (facet normals etc)
  subroutine blasius_set_coef(this, c)
    class(blasius_t), intent(inout) :: this
    type(coef_t), target, intent(inout) :: c
    this%c => c
  end subroutine blasius_set_coef
     
end module blasius
