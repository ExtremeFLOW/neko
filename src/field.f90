!> Defines a field
!
module field
  use num_types
  use mesh
  implicit none
  
  type field_t
     real(kind=dp), allocatable :: x(:,:,:,:)
     real(kind=dp), allocatable :: y(:,:,:,:)
     real(kind=dp), allocatable :: z(:,:,:,:)     
     type(mesh_t), pointer :: msh
  end type field_t

  interface assignment(=)
     module procedure field_assign_field
  end interface assignment(=)

  interface field_add
     module procedure field_add_field, field_add_scalar
  end interface field_add

contains

  !> Initialize a field @a f
  subroutine field_init(f, mesh)
    type(field_t), intent(inout) :: f !< Field
    type(mesh_t), target, intent(in) :: mesh
    integer :: ierr

    f%msh => mesh
    
    if (.not. allocated(f%x)) then
       allocate(f%x(f%msh%lx1, f%msh%ly1, f%msh%lz1, f%msh%lelv), stat = ierr)
       f%x = 0d0
    end if

    if (.not. allocated(f%y)) then
       allocate(f%y(f%msh%lx1, f%msh%ly1, f%msh%lz1, f%msh%lelv), stat = ierr)
       f%y = 0d0
    end if

    if (.not. allocated(f%z)) then
       allocate(f%x(f%msh%lx1, f%msh%ly1, f%msh%lz1, f%msh%lelv), stat = ierr)
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
    type(field_t), intent(in) :: f
    integer :: i, j, k, l

    do i = 1, f%msh%lelv
       do l = 1, f%msh%lz1
          do k = 1, f%msh%ly1
             do j = 1, f%msh%lz1
                this_f%x(j, k, l, i) = f%x(j, k, l, i)
                this_f%y(j, k, l, i) = f%y(j, k, l, i)
                this_f%z(j, k, l, i) = f%z(j, k, l, i)
             end do
          end do
       end do
    end do
    
  end subroutine field_assign_field

  !> Add \f$ F(u_1, u_2, ... , u_n) =
  !! F(u_1, u_2, ... , u_n) + G(u_1, u_2, ... , u_n) \f$
  !! @note Component wise
  subroutine field_add_field(f, g)
    type(field_t), intent(inout) :: f
    type(field_t), intent(in) :: g
    integer i, j, k, l

    do i = 1, f%msh%lelv
       do l = 1, f%msh%lz1
          do k = 1, f%msh%ly1
             do j = 1, f%msh%lz1
                f%x(j, k, l, i) = f%x(j, k, l, i) + g%x(j, k, l, i)
                f%y(j, k, l, i) = f%y(j, k, l, i) + g%y(j, k, l, i)
                f%z(j, k, l, i) = f%z(j, k, l, i) + g%z(j, k, l, i)
             end do
          end do
       end do
    end do
    
  end subroutine field_add_field


  !> Add \f$ F(u_1, u_2, ... , u_n) =
  !! F(u_1, u_2, ... , u_n) + a \f$
  subroutine field_add_scalar(f, a)
    type(field_t), intent(inout) :: f
    real(kind=dp), intent(in) :: a
    integer i, j, k, l

    do i = 1, f%msh%lelv
       do l = 1, f%msh%lz1
          do k = 1, f%msh%ly1
             do j = 1, f%msh%lz1
                f%x(j, k, l, i) = f%x(j, k, l, i) + a
                f%y(j, k, l, i) = f%y(j, k, l, i) + a
                f%z(j, k, l, i) = f%z(j, k, l, i) + a
             end do
          end do
       end do
    end do
    
  end subroutine field_add_scalar

end module field

