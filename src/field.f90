!> Defines a field
!
module field
  use num_types

  
  type field_t
     real(kind=dp), allocatable :: x(:,:,:,:)
     real(kind=dp), allocatable :: y(:,:,:,:)
     real(kind=dp), allocatable :: z(:,:,:,:)     
  end type field_t

  interface assignment(=)
     module procedure field_assign_field
  end interface assignment(=)

contains

  !> Initialize a field @a f
  subroutine field_init(f, nx, ny, nz, ne)
    type(field_t), intent(inout) :: f !< Field
    integer, intent(in) :: nx         !< Points in x-dir
    integer, intent(in) :: ny         !< Points in y-dir
    integer, intent(in) :: nz         !< Points in z-dir
    integer, intent(in) :: ne         !< Number of elements
    integer :: ierr
    
    if (.not. allocated(f%x)) then
       allocate(f%x(nx, ny, nz, ne), stat = ierr)
       f%x = 0d0
    end if

    if (.not. allocated(f%y)) then
       allocate(f%y(nx, ny, nz, ne), stat = ierr)
       f%y = 0d0
    end if

    if (.not. allocated(f%z)) then
       allocate(f%z(nx, ny, nz, ne), stat = ierr)
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
    integer :: nj, nk, nl, ne

    nj = size(f%x, 1)
    nk = size(f%x, 2)
    nl = size(f%x, 3)
    ne = size(f%x, 4)

    do i = 1, n
       do l = 1, nl
          do k = 1, nk
             do j = 1, nj
                this_f%x(j, k, l, i) = f%x(j, k, l, i)
                this_f%y(j, k, l, i) = f%y(j, k, l, i)
                this_f%z(j, k, l, i) = f%z(j, k, l, i)
             end do
          end do
       end do
    end do
    

  end subroutine field_assign_field

end module field

