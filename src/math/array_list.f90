module vector_list
  use vector, only : vector_t
  use num_types, only : rp
  implicit none
  private

  !> An array of vectors
  type, public :: vector_list_t
     real(kind=rp), allocatable :: items(:)
   contains
     !> Append a field to the list.
     procedure, pass(this) :: append => vector_list_append
     !> Destructor
     procedure, pass(this) :: free => vector_list_free
  end type vector_list_t

contains
  !> Append an vector to the list.
  !! @param item The vector to append.
  subroutine vector_list_append(this, item)
    class(vector_list_t), intent(inout) :: this
    class(vector_t), intent(in), target :: item
    real(kind=rp), allocatable :: tmp(:)
    integer :: len

    len = size(this%fields)

    allocate(tmp(len+1))
    tmp(1:len) = this%fields
    call move_alloc(tmp, this%fields)
    this%fields(len+1)%f => f

  end subroutine vector_list_append

  !> Destructor.
  subroutine vector_list_free(this)
    class(field_list_t), intent(inout) :: this
    integer :: i, n_fields

    if (allocated(this%fields)) then
       n_fields = size(this%fields)
       do i=1, n_fields
          call this%fields(i)%f%free()
          nullify(this%fields(i)%f)
       end do
       deallocate(this%fields)
    end if
  end subroutine vector_list_free



end module vector_list
