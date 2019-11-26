!> Implements an unordered set ADT
!! @details A unordered set storing a fixed data-type @a data 
module uset
  use num_types
  use htable
  implicit none
  private

  !> Integer based unordered set
  type, public :: uset_i4_t
     type(htable_i4_t) :: t
  end type uset_i4_t

contains

  subroutine uset_i4_init(s, n)
    type(uset_i4_t), intent(inout) :: s
    integer, value, optional :: n
    integer :: key

    if (present(n)) then
       call s%t%init(n)
    else
       call s%t%init(64)
    end if    
  end subroutine uset_i4_init
  
  subroutine uset_i4_free(s)
    type(uset_i4_t), intent(inout) :: s

    call s%t%free()
    
  end subroutine uset_i4_free

  subroutine uset_i4_add(s, key)
    type(uset_i4_t), intent(inout) :: s
    integer :: key

    call s%t%set(key, 1)
  end subroutine uset_i4_add


end module uset
