! Copyright (c) 2018-2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
! Defines a vector list
module vector_list
  use, intrinsic :: iso_fortran_env, only : error_unit
  use vector, only : vector_ptr_t, vector_t
  use iso_c_binding, only : c_ptr
  use num_types, only : rp
  use utils, only : neko_error, NEKO_VARNAME_LEN
  use comm, only : pe_rank
  implicit none
  private

  !> vector_list_t, To be able to group vectors together
  type, public :: vector_list_t
     type(vector_ptr_t), allocatable :: items(:)
   contains
     !> Constructor. Allocates array and pointers.
     procedure, pass(this) :: init => vector_list_init
     !> Destructor.
     procedure, pass(this) :: free => vector_list_free
     !> Append a vector to the list.
     procedure, pass(this) :: append => vector_list_append
     generic :: get => get_by_index, get_by_name
     !> Get an item pointer by array index
     procedure, pass(this) :: get_by_index => vector_list_get_by_index
     !> Get an item pointer by vector name
     procedure, pass(this) :: get_by_name => vector_list_get_by_name
     !> Point item at given index.
     generic :: assign => assign_to_ptr, assign_to_vector_ptr
     procedure, pass(this) :: assign_to_ptr => vector_list_assign_to_ptr
     procedure, pass(this) :: assign_to_vector_ptr => &
          vector_list_assign_to_vector_ptr
     procedure, pass(this) :: assign_to_vector => vector_list_assign_to_vector

     !> Get device pointer for a given index.
     procedure, pass(this) :: x_d => vector_list_x_d
     !> Get pointer to the raw data array for a given index.
     procedure, pass(this) :: x => vector_list_x
     !> Get number of items in the list.
     procedure, pass(this) :: size => vector_list_size
     !> Get the size of an item in the list.
     procedure, pass(this) :: item_size => vector_list_item_size
     !> Get the name for an item in the list.
     procedure, pass(this) :: name => vector_list_name
  end type vector_list_t

contains
  !> Constructor. Just allocates the array.
  !! @param size The size of the list to preallocate
  subroutine vector_list_init(this, size)
    class(vector_list_t), intent(inout) :: this
    integer, intent(in) :: size

    call this%free()

    allocate(this%items(size))
  end subroutine vector_list_init

  !> Get number of items in the list.
  pure function vector_list_size(this) result(n)
    class(vector_list_t), intent(in) :: this
    integer :: n
    n = size(this%items)
  end function vector_list_size

  !> Get an item pointer by array index
  !! @param i The index of the item.
  function vector_list_get_by_index(this, i) result(f)
    class(vector_list_t), target, intent(inout) :: this
    type(vector_t), pointer :: f
    integer, intent(in) :: i
    f => this%items(i)%ptr
  end function vector_list_get_by_index

  !> Get an item pointer by array index
  !! @param name The name of the item.
  function vector_list_get_by_name(this, name) result(f)
    class(vector_list_t), target, intent(inout) :: this
    type(vector_t), pointer :: f
    character(len=*), intent(in) :: name
    integer :: i

    nullify(f)

    do i = 1, this%size()
       if (this%name(i) .eq. trim(name)) then
          f => this%items(i)%ptr
          return
       end if
    end do

    if (pe_rank .eq. 0) then
       write(error_unit,*) "Current vector list contents:"

       do i = 1, this%size()
          write(error_unit,*) "- ", this%name(i)
       end do
    end if

    call neko_error("No vector with name " // trim(name) // " found in list")
  end function vector_list_get_by_name

  !> Append a vector to the list.
  !! @param f The vector to append.
  subroutine vector_list_append(this, f)
    class(vector_list_t), intent(inout) :: this
    class(vector_t), intent(in), target :: f
    type(vector_ptr_t), allocatable :: tmp(:)
    integer :: len

    if (.not. allocated(this%items)) then
       allocate(this%items(1))
       call this%items(1)%init(f)
       return
    end if

    len = size(this%items)

    allocate(tmp(len+1))
    tmp(1:len) = this%items
    call move_alloc(tmp, this%items)
    call this%items(len+1)%init(f)

  end subroutine vector_list_append

  !> Destructor.
  subroutine vector_list_free(this)
    class(vector_list_t), intent(inout) :: this
    integer :: i, n_fields

    if (allocated(this%items)) then
       n_fields = this%size()
       do i = 1, n_fields
          call this%items(i)%free()
       end do
       deallocate(this%items)
    end if

  end subroutine vector_list_free

  !> Get device pointer for a given index
  !! @param i The index of the item.
  function vector_list_x_d(this, i) result(ptr)
    class(vector_list_t), intent(in) :: this
    integer, intent(in) :: i
    type(c_ptr) :: ptr

    ptr = this%items(i)%ptr%x_d
  end function vector_list_x_d

  function vector_list_x(this, i) result(x)
    class(vector_list_t), target, intent(in) :: this
    real(kind=rp), pointer, contiguous :: x(:)
    integer, intent(in) :: i
    x => this%items(i)%ptr%x
  end function vector_list_x

  !> Get the size of item `i`.
  !! @param i The index of the item.
  function vector_list_item_size(this, i) result(size)
    class(vector_list_t), target, intent(in) :: this
    integer, intent(in) :: i
    integer :: size

    size = this%items(i)%ptr%size()

  end function vector_list_item_size

  !> Point item at a given index.
  !! @param i The index of the item.
  !! @param ptr A vector pointer to point the item to.
  subroutine vector_list_assign_to_ptr(this, i, ptr)
    class(vector_list_t), intent(inout) :: this
    integer, intent(in) :: i
    type(vector_t), pointer, intent(in) :: ptr

    call this%items(i)%init(ptr)

  end subroutine vector_list_assign_to_ptr

  !> Point item at a given index.
  !! @param i The index of the item.
  !! @param ptr An encapsulated vector pointer to point the item to.
  subroutine vector_list_assign_to_vector_ptr(this, i, ptr)
    class(vector_list_t), intent(inout) :: this
    integer, intent(in) :: i
    type(vector_ptr_t), target, intent(in) :: ptr

    call this%items(i)%init(ptr%ptr)
  end subroutine vector_list_assign_to_vector_ptr

  !> Point item at a given index.
  !! @param i The index of the item.
  !! @param vec A vector to point the item to.
  subroutine vector_list_assign_to_vector(this, i, vec)
    class(vector_list_t), intent(inout) :: this
    integer, intent(in) :: i
    type(vector_t), target, intent(in) :: vec

    call this%items(i)%init(vec)
  end subroutine vector_list_assign_to_vector

  !> Get the name for an item in the list.
  !! @param i The index of the item.
  function vector_list_name(this, i) result(result)
    class(vector_list_t), target, intent(in) :: this
    integer, intent(in) :: i
    character(len=NEKO_VARNAME_LEN) :: result

    result = this%items(i)%ptr%name
  end function vector_list_name


end module vector_list
