!> Defines a field
!
module field
  use num_types
  use math
  use mesh
  use space
  implicit none
  
  type field_t
     !> @todo Is x really a good name for a field?
     real(kind=dp), allocatable :: x(:,:,:,:,:)
     
     type(space_t), pointer :: Xh !< Function space \f$ X_h \f$
     type(mesh_t), pointer :: msh !< Mesh
     character(len=80) :: name
  end type field_t

  interface assignment(=)
     module procedure field_assign_field
  end interface assignment(=)

  interface field_add
     module procedure field_add_field, field_add_scalar
  end interface field_add

contains

  !> Initialize a field @a f on the mesh @a msh
  subroutine field_init(f, msh, space, fld_name)
    type(field_t), intent(inout) :: f       !< Field to be initialized
    type(mesh_t), target, intent(in) :: msh !< underlying mesh of the field
    type(space_t), target, intent(in) :: space !< Function space for the field
    character(len=*), optional :: fld_name     !< Name of the field
    integer :: ierr
    integer :: lx, ly, lz, nelv, ndim

    call field_free(f)

    f%Xh => space
    f%msh => msh

    lx = f%Xh%lx
    ly = f%Xh%ly
    lz = f%Xh%lz
    nelv = f%msh%nelv
    ndim = f%Xh%ndim
    
    if (.not. allocated(f%x)) then
       allocate(f%x(lx, ly, lz, nelv, ndim), stat = ierr)        
       f%x = 0d0
    end if

    if (present(fld_name)) then
       f%name = fld_name
    else
       f%name = "Field"
    end if
    
  end subroutine field_init

  !> Deallocate a field @a f
  subroutine field_free(f)
    type(field_t), intent(inout) :: f
    
    if (allocated(f%x)) then
       deallocate(f%x)
    end if

    nullify(f%msh)
    nullify(f%Xh)

  end subroutine field_free

  subroutine field_assign_field(this_f, f)
    type(field_t), intent(inout) :: this_f
    type(field_t), intent(in) :: f
    integer :: n

    n = f%msh%nelv * f%Xh%lx * f%Xh%ly * f%Xh%lz * f%Xh%ndim
    call copy(this_f%x, f%x, n)
    
  end subroutine field_assign_field

  !> Add \f$ F(u_1, u_2, ... , u_n) =
  !! F(u_1, u_2, ... , u_n) + G(u_1, u_2, ... , u_n) \f$
  !! @note Component wise
  subroutine field_add_field(f, g)
    type(field_t), intent(inout) :: f
    type(field_t), intent(inout) :: g
    integer :: n

    n = f%msh%nelv * f%Xh%lx * f%Xh%ly * f%Xh%lz * f%Xh%ndim    
    call add2(f%x, g%x, n)

  end subroutine field_add_field


  !> Add \f$ F(u_1, u_2, ... , u_n) =
  !! F(u_1, u_2, ... , u_n) + a \f$
  subroutine field_add_scalar(f, a)
    type(field_t), intent(inout) :: f
    real(kind=dp), intent(inout) :: a
    integer :: n

    n = f%msh%nelv * f%Xh%lx * f%Xh%ly * f%Xh%lz * f%Xh%ndim    
    call cadd(f%x, a, n)

  end subroutine field_add_scalar

end module field

