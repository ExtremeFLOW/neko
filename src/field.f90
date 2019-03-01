!> Defines a field
!
module field
  use num_types
  use math
  use mesh
  implicit none
  
  type field_t
     real(kind=dp), allocatable :: x(:,:,:,:)
     real(kind=dp), allocatable :: y(:,:,:,:)
     real(kind=dp), allocatable :: z(:,:,:,:)     

     integer :: lx1
     integer :: ly1
     integer :: lz1

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
  subroutine field_init(f, msh, lx1, ly1, lz1)
    type(field_t), intent(inout) :: f !< Field to be initialized
    type(mesh_t), target, intent(in) :: msh !< Underlying mesh of the field
    integer, intent(in) :: lx1  !< Polynomial dimension in x-direction
    integer, intent(in) :: ly1  !< Polynomial dimension in y-direction
    integer, intent(in) :: lz1  !< Polynomial dimension in z-direction
    integer :: ierr

    f%lx1 = lx1
    f%ly1 = ly1
    f%lz1 = lz1
    f%msh => msh
    
     if (.not. allocated(f%x)) then
        allocate(f%x(f%lx1, f%ly1, f%lz1, f%msh%nelv), stat = ierr)        
        f%x = 0d0
     end if

     if (.not. allocated(f%y)) then
        allocate(f%y(f%lx1, f%ly1, f%lz1, f%msh%nelv), stat = ierr)
        f%y = 0d0
     end if

     if (.not. allocated(f%z)) then
        allocate(f%z(f%lx1, f%ly1, f%lz1, f%msh%nelv), stat = ierr)
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

  end subroutine field_free

  subroutine field_assign_field(this_f, f)
    type(field_t), intent(inout) :: this_f
    type(field_t), intent(inout) :: f
    integer :: n

    n = f%msh%nelv * f%lx1 * f%ly1 * f%lz1
    call copy(this_f%x, f%x, n)
    call copy(this_f%y, f%y, n)
    call copy(this_F%z, f%z, n)

  end subroutine field_assign_field

  !> Add \f$ F(u_1, u_2, ... , u_n) =
  !! F(u_1, u_2, ... , u_n) + G(u_1, u_2, ... , u_n) \f$
  !! @note Component wise
  subroutine field_add_field(f, g)
    type(field_t), intent(inout) :: f
    type(field_t), intent(inout) :: g
    integer :: n

    n = f%msh%nelv * f%lx1 * f%ly1 * f%lz1
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


    n = f%msh%nelv * f%lx1 * f%ly1 * f%lz1
    call cadd(f%x, a, n)
    call cadd(f%y, a, n)
    call cadd(f%z, a, n)

  end subroutine field_add_scalar

end module field

