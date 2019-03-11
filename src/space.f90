!> Defines a function space
module space
  use num_types
  implicit none

  type space_t
     integer :: lx1             !< Polynomial dimension in x-direction
     integer :: ly1             !< Polynomial dimension in y-direction
     integer :: lz1             !< Polynomial dimension in z-direction

     !> @todo Store gll points etc in the space
  end type space_t
  
  interface space_t
     module procedure space_init
  end interface space_t
contains

  subroutine space_init(this, lx1, ly1, lz1)
    type(space_t), intent(inout) :: this
    integer, intent(in) :: lx1
    integer, intent(in) :: ly1
    integer, intent(in) :: lz1

    this%lx1 = lx1
    this%ly1 = ly1
    this%lz1 = lz1

  end subroutine space_init


end module space
