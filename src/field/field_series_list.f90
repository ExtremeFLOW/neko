! Copyright (c) 2024, The Neko Authors
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
!     disclaimer in the documentation and/or materials provided
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
module field_series_list
  use field_series, only : field_series_t, field_series_ptr_t
  use field, only : field_t
  use utils, only : neko_error
  implicit none
  private

  !> field_series_list_t, To be able to group field series together
  type, public :: field_series_list_t
     type(field_series_ptr_t), allocatable :: items(:)
   contains
     !> Constructor. Allocates array and pointers.
     procedure, pass(this) :: init => field_series_list_init
     !> Destructor.
     procedure, pass(this) :: free => field_series_list_free
     !> Append a field series to the list.
     procedure, pass(this) :: append => field_series_list_append
     !> Get an item pointer by array index
     procedure, pass(this) :: get => field_series_list_get_by_index
     !> Point item at given index.
     procedure, pass(this) :: assign => field_series_list_assign_to_ptr
     !> Get number of items in the list.
     procedure, pass(this) :: size => field_series_list_size
     !> Get the size of the field series for item `i`.
     procedure, pass(this) :: series_size => field_series_list_series_size
  end type field_series_list_t

contains
  !> Constructor. Just allocates the array.
  !! @param size The size of the list to preallocate
  subroutine field_series_list_init(this, size)
    class(field_series_list_t), intent(inout) :: this
    integer, intent(in) :: size

    call this%free()

    allocate(this%items(size))
  end subroutine field_series_list_init

  !> Get number of items in the list.
  pure function field_series_list_size(this) result(n)
    class(field_series_list_t), intent(in) :: this
    integer :: n
    if (allocated(this%items)) then
       n = size(this%items)
    else
       n = 0
    end if
  end function field_series_list_size

  !> Get an item pointer by array index
  !! @param i The index of the item.
  function field_series_list_get_by_index(this, i) result(fs)
    class(field_series_list_t), intent(in) :: this
    type(field_series_t), pointer :: fs
    integer, intent(in) :: i
    fs => this%items(i)%ptr
  end function field_series_list_get_by_index

  !> Append a field series to the list.
  !! @param fs The field series to append.
  subroutine field_series_list_append(this, fs)
    class(field_series_list_t), intent(inout) :: this
    class(field_series_t), intent(in), target :: fs
    type(field_series_ptr_t), allocatable :: tmp(:)
    integer :: len

    if (.not. allocated(this%items)) then
       allocate(this%items(1))
       this%items(1)%ptr => fs
       return
    end if

    len = size(this%items)

    allocate(tmp(len+1))
    tmp(1:len) = this%items
    call move_alloc(tmp, this%items)
    this%items(len+1)%ptr => fs

  end subroutine field_series_list_append

  !> Destructor.
  subroutine field_series_list_free(this)
    class(field_series_list_t), intent(inout) :: this
    integer :: i, n_series

    if (allocated(this%items)) then
       n_series = this%size()
       do i=1, n_series
          if (associated(this%items(i)%ptr)) then
             call this%items(i)%ptr%free()
          end if
          nullify(this%items(i)%ptr)
       end do
       deallocate(this%items)
    end if

  end subroutine field_series_list_free

  !> Get the size of the field series for item `i`.
  !! @param i The index of the item.
  function field_series_list_series_size(this, i) result(size)
    class(field_series_list_t), intent(in) :: this
    integer, intent(in) :: i
    integer :: size

    size = this%items(i)%ptr%size()

  end function field_series_list_series_size

  !> Point item at a given index.
  !! @param i The index of the item.
  !! @param ptr A field series pointer to point the item to.
  subroutine field_series_list_assign_to_ptr(this, i, ptr)
    class(field_series_list_t), intent(inout) :: this
    integer, intent(in) :: i
    type(field_series_t), pointer, intent(in) :: ptr

    this%items(i)%ptr => ptr
  end subroutine field_series_list_assign_to_ptr

end module field_series_list