!> Defines a function space
module space
  use num_types
  use speclib
  implicit none

  type space_t
     integer :: lx1             !< Polynomial dimension in x-direction
     integer :: ly1             !< Polynomial dimension in y-direction
     integer :: lz1             !< Polynomial dimension in z-direction
     
     real(kind=sp) , allocatable :: zgm1(:, :)
     real(kind=sp) , allocatable :: zgm2(:, :)
     real(kind=sp) , allocatable :: zgm3(:, :)
     
     real(kind=sp), allocatable :: wxm1(:)
     real(kind=sp), allocatable :: wym1(:)
     real(kind=sp), allocatable :: wzm1(:)

     !> @todo Store gll points etc in the space
  end type space_t
  
contains

  !> Initialize a function space @a s with given polynomial dimensions
  subroutine space_init(s, lx1, ly1, lz1)
    type(space_t), intent(inout) :: s
    integer, intent(in) :: lx1
    integer, intent(in) :: ly1
    integer, intent(in) :: lz1

    call space_free(s)

    s%lx1 = lx1
    s%ly1 = ly1
    s%lz1 = lz1

    allocate(s%zgm1(lx1, 3))
    allocate(s%zgm2(ly1, 3))
    allocate(s%zgm3(lz1, 3))

    allocate(s%wxm1(lx1))
    allocate(s%wym1(ly1))
    allocate(s%wzm1(lz1))
    

    call zwgll(s%zgm1(1,1), s%wxm1, lx1)
    call zwgll(s%zgm1(1,2), s%wym1, ly1)
    call zwgll(s%zgm1(1,3), s%wzm1, lz1)

  end subroutine space_init

  subroutine space_free(s)
    type(space_t), intent(in) :: s
    
    if (allocated(s%zgm1)) then
       deallocate(s%zgm1)
    end if

    if (allocated(s%zgm2)) then
       deallocate(s%zgm2)
    end if

    if (allocated(s%zgm3)) then
       deallocate(s%zgm3)
    end if

    if (allocated(s%wxm1)) then
       deallocate(s%wxm1)
    end if

    if (allocated(s%wym1)) then
       deallocate(s%wym1)
    end if

    if (allocated(s%wzm1)) then
       deallocate(s%wzm1)
    end if

  end subroutine space_free


end module space
