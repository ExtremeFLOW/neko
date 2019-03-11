!> Defines a field
!
module field
  use num_types
  use math
  use mesh
  use space
  implicit none
  
  type field_t
     real(kind=dp), allocatable :: x(:,:,:,:)
     real(kind=dp), allocatable :: y(:,:,:,:)
     real(kind=dp), allocatable :: z(:,:,:,:)     

     type(space_t), pointer :: Vh !< Function space \f$ V_h \f$
     type(mesh_t), pointer :: msh
  end type field_t

  interface assignment(=)
     module procedure field_assign_field
  end interface assignment(=)

  interface field_add
     module procedure field_add_field, field_add_scalar
  end interface field_add

contains

  !> Initialize a field @a f on the mesh @a msh
  subroutine field_init(f, msh, space)
    type(field_t), intent(inout) :: f       !< Field to be initialized
    type(mesh_t), target, intent(in) :: msh !< Underlying mesh of the field
    type(space_t), target, intent(in) :: space !< Function space for the field
    integer :: ierr
    integer :: lx, ly, lz, nelv

    call field_free(f)

    f%Vh => space
    f%msh => msh

    lx = f%Vh%lx
    ly = f%Vh%ly
    lz = f%Vh%lz
    nelv = f%msh%nelv
        
     if (.not. allocated(f%x)) then
        allocate(f%x(lx, ly, lz, nelv), stat = ierr)        
        f%x = 0d0
     end if

     if (.not. allocated(f%y)) then
        allocate(f%y(lx, ly, lz, nelv), stat = ierr)
        f%y = 0d0
     end if

     if (.not. allocated(f%z)) then
        allocate(f%z(lx, ly, lz, nelv), stat = ierr)
        f%z = 0d0
     end if

  end subroutine field_init

  !> Deallocate a field @a f
  subroutine field_free(f)
    type(field_t), intent(inout) :: f
    
    if (allocated(f%x)) then
       deallocate(f%x)
    end if

    if (allocated(f%y)) then
       deallocate(f%y)
    end if

    if (allocated(f%z)) then
       deallocate(f%z)
    end if

    nullify(f%msh)
    nullify(f%Vh)

  end subroutine field_free

  subroutine field_assign_field(this_f, f)
    type(field_t), intent(inout) :: this_f
    type(field_t), intent(in) :: f
    integer :: n

    n = f%msh%nelv * f%Vh%lx * f%Vh%ly * f%Vh%lz
    call copy(this_f%x, f%x, n)
    call copy(this_f%y, f%y, n)
    call copy(this_f%z, f%z, n)

  end subroutine field_assign_field

  !> Add \f$ F(u_1, u_2, ... , u_n) =
  !! F(u_1, u_2, ... , u_n) + G(u_1, u_2, ... , u_n) \f$
  !! @note Component wise
  subroutine field_add_field(f, g)
    type(field_t), intent(inout) :: f
    type(field_t), intent(inout) :: g
    integer :: n

    n = f%msh%nelv * f%Vh%lx * f%Vh%ly * f%Vh%lz
    call add2(f%x, g%x, n)
    call add2(f%y, g%y, n)
    call add2(f%z, g%z, n)

  end subroutine field_add_field


  !> Add \f$ F(u_1, u_2, ... , u_n) =
  !! F(u_1, u_2, ... , u_n) + a \f$
  subroutine field_add_scalar(f, a)
    type(field_t), intent(inout) :: f
    real(kind=dp), intent(inout) :: a
    integer :: n

    n = f%msh%nelv * f%Vh%lx * f%Vh%ly * f%Vh%lz
    call cadd(f%x, a, n)
    call cadd(f%y, a, n)
    call cadd(f%z, a, n)

  end subroutine field_add_scalar

end module field

