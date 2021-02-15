!> Defines inflow dirichlet conditions
module usr_inflow
  use num_types
  use coefs
  use inflow
  use utils
  implicit none
  private
  
  !> User defined dirichlet condition for inlet (vector valued)
  type, public, extends(inflow_t) :: usr_inflow_t
     type(coef_t), pointer :: c => null()
     procedure(usr_inflow_eval), nopass, pointer :: eval => null()
   contains
     procedure, pass(this) :: apply_scalar => usr_inflow_apply_scalar
     procedure, pass(this) :: apply_vector => usr_inflow_apply_vector
     procedure, pass(this) :: validate => usr_inflow_validate
     procedure, pass(this) :: set_coef => usr_inflow_set_coef
     procedure, pass(this) :: set_eval => usr_inflow_set_eval
  end type usr_inflow_t

  !> Abstract interface defining a user defined inflow condition (pointwise)
  abstract interface
     subroutine usr_inflow_eval(u, v, w, x, y, z, nx, ny, nz)
       import dp
       real(kind=dp), intent(inout) :: u
       real(kind=dp), intent(inout) :: v
       real(kind=dp), intent(inout) :: w
       real(kind=dp), intent(in) :: x
       real(kind=dp), intent(in) :: y
       real(kind=dp), intent(in) :: z
       real(kind=dp), intent(in) :: nx
       real(kind=dp), intent(in) :: ny
       real(kind=dp), intent(in) :: nz
     end subroutine usr_inflow_eval
  end interface

  public :: usr_inflow_eval

contains
     
  !> No-op scalar apply
  subroutine usr_inflow_apply_scalar(this, x, n)
    class(usr_inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout),  dimension(n) :: x
  end subroutine usr_inflow_apply_scalar

  !> Apply user defined inflow conditions (vector valued)
  subroutine usr_inflow_apply_vector(this, x, y, z, n)
    class(usr_inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout),  dimension(n) :: x
    real(kind=dp), intent(inout),  dimension(n) :: y
    real(kind=dp), intent(inout),  dimension(n) :: z
    integer :: i, m, k, idx(4), facet

    associate(xc => this%c%dof%x, yc => this%c%dof%y, zc => this%c%dof%z, &
         nx => this%c%nx, ny => this%c%ny, nz => this%c%nz, &
         lx => this%c%Xh%lx)
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, lx, lx, lx)
         select case(facet)
         case(1,2)          
            call this%eval(x(k), y(k), z(k), &
                 xc(idx(1), idx(2), idx(3), idx(4)), &
                 yc(idx(1), idx(2), idx(3), idx(4)), &
                 zc(idx(1), idx(2), idx(3), idx(4)), &
                 nx(idx(2), idx(3), facet, idx(4)), &
                 ny(idx(2), idx(3), facet, idx(4)), &
                 nz(idx(2), idx(3), facet, idx(4)))
         case(3,4)
            call this%eval(x(k), y(k), z(k), &
                 xc(idx(1), idx(2), idx(3), idx(4)), &
                 yc(idx(1), idx(2), idx(3), idx(4)), &
                 zc(idx(1), idx(2), idx(3), idx(4)), &       
                 nx(idx(1), idx(3), facet, idx(4)), &
                 ny(idx(1), idx(3), facet, idx(4)), &
                 nz(idx(1), idx(3), facet, idx(4)))
         case(5,6)
            call this%eval(x(k), y(k), z(k), &
                 xc(idx(1), idx(2), idx(3), idx(4)), &
                 yc(idx(1), idx(2), idx(3), idx(4)), &
                 zc(idx(1), idx(2), idx(3), idx(4)), &                     
                 nx(idx(1), idx(2), facet, idx(4)), &
                 ny(idx(1), idx(2), facet, idx(4)), &
                 nz(idx(1), idx(2), facet, idx(4)))
         end select
      end do
    end associate
    
  end subroutine usr_inflow_apply_vector

  !> Assign coefficients (facet normals etc)
  subroutine usr_inflow_set_coef(this, c)
    class(usr_inflow_t), intent(inout) :: this
    type(coef_t), target, intent(inout) :: c
    this%c => c
  end subroutine usr_inflow_set_coef

  !> Assign user provided eval function
  subroutine usr_inflow_set_eval(this, usr_eval)
    class(usr_inflow_t), intent(inout) :: this
    procedure(usr_inflow_eval) :: usr_eval
    this%eval => usr_eval
  end subroutine usr_inflow_set_eval

  !> Validate user inflow condition
  subroutine usr_inflow_validate(this)
    class(usr_inflow_t), intent(inout) :: this
    logical :: valid

    valid = .true. ! Assert it's going to be ok...    
    if (.not. associated(this%c)) then
       call neko_warning('Missing coefficients')
       valid = .false.       
    end if

    if (.not. associated(this%eval)) then
       call neko_warning('Missing eval function')
       valid = .false.
    end if

    if (.not. valid) then
       call neko_error('Invalid user defined inflow condition')
    end if
    
  end subroutine usr_inflow_validate
  
end module usr_inflow
