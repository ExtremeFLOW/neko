!> Dirichlet condition applied in the facet normal direction
module facet_normal
  use num_types
  use dirichlet
  use coefs
  use math
  use utils
  implicit none
  private

  !> Dirichlet condition in facet normal direction
  type, public, extends(dirichlet_t) :: facet_normal_t
     type(coef_t), pointer :: c => null()
   contains
     procedure, pass(this) :: apply_scalar => facet_normal_apply_scalar
     procedure, pass(this) :: apply_vector => facet_normal_apply_vector
     procedure, pass(this) :: apply_surfvec => facet_normal_apply_surfvec
     procedure, pass(this) :: set_coef => facet_normal_set_coef
  end type facet_normal_t

contains

  !> No-op scalar apply
  subroutine facet_normal_apply_scalar(this, x, n)
    class(facet_normal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout), dimension(n) :: x
  end subroutine facet_normal_apply_scalar

  !> No-op vector apply
  subroutine facet_normal_apply_vector(this, x, y, z, n)
    class(facet_normal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout), dimension(n) :: x
    real(kind=dp), intent(inout), dimension(n) :: y
    real(kind=dp), intent(inout), dimension(n) :: z
  end subroutine facet_normal_apply_vector

  !> Apply in facet normal direction (vector valued)
  subroutine facet_normal_apply_surfvec(this, x, y, z, u, v, w, n)
    class(facet_normal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout), dimension(n) :: x
    real(kind=dp), intent(inout), dimension(n) :: y
    real(kind=dp), intent(inout), dimension(n) :: z
    real(kind=dp), intent(inout), dimension(n) :: u
    real(kind=dp), intent(inout), dimension(n) :: v
    real(kind=dp), intent(inout), dimension(n) :: w
    integer :: i, m, k, idx(4), facet
    
    if (.not. associated(this%c)) then
       call neko_error('No coefficients assigned')
    end if
    associate(c => this%c)
      m = this%msk(0)
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, c%Xh%lx, c%Xh%lx, c%Xh%lx)
         select case(facet)
         case(1,2)          
            x(k) = u(k) * c%nx(idx(2), idx(3), facet, idx(4)) &
                 * c%area(idx(2), idx(3), facet, idx(4))
            y(k) = v(k) * c%ny(idx(2), idx(3), facet, idx(4)) &
                 * c%area(idx(2), idx(3), facet, idx(4))
            z(k) = w(k) * c%nz(idx(2), idx(3), facet, idx(4)) &
                 * c%area(idx(2), idx(3), facet, idx(4))
         case(3,4)
            x(k) = u(k) * c%nx(idx(1), idx(3), facet, idx(4)) &
                 * c%area(idx(1), idx(3), facet, idx(4))
            y(k) = v(k) * c%ny(idx(1), idx(3), facet, idx(4)) &
                 * c%area(idx(1), idx(3), facet, idx(4))
            z(k) = w(k) * c%nz(idx(1), idx(3), facet, idx(4)) &
                 * c%area(idx(1), idx(3), facet, idx(4))          
         case(5,6)
            x(k) = u(k) * c%nx(idx(1), idx(2), facet, idx(4)) &
                 * c%area(idx(1), idx(2), facet, idx(4))
            y(k) = v(k) * c%ny(idx(1), idx(2), facet, idx(4)) &
                 * c%area(idx(1), idx(2), facet, idx(4))
            z(k) = w(k) * c%nz(idx(1), idx(2), facet, idx(4)) &
                 * c%area(idx(1), idx(2), facet, idx(4))
         end select
      end do
    end associate
    
  end subroutine facet_normal_apply_surfvec
  
  !> Assign coefficients (facet normals etc)
  subroutine facet_normal_set_coef(this, c)
    class(facet_normal_t), intent(inout) :: this
    type(coef_t), target, intent(inout) :: c
    this%c => c
  end subroutine facet_normal_set_coef
  
end module facet_normal
