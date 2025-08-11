! Copyright (c) 2025, The Neko Authors
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
!> Contains the `field_series_list_t` type for managing multiple field series.
module field_series_list
  use field_series, only : field_series_t, field_series_ptr_t
  use utils, only : neko_error
  implicit none
  private

  !> A list of field series pointers, used for managing multiple scalar lag fields
  type, public :: field_series_list_t
     type(field_series_ptr_t), allocatable :: items(:)
     integer, private :: n_items = 0
   contains
     !> Constructor.
     procedure, pass(this) :: init => field_series_list_init
     !> Destructor.
     procedure, pass(this) :: free => field_series_list_free
     !> Add a field series to the list.
     procedure, pass(this) :: append => field_series_list_append
     !> Get a field series by index.
     procedure, pass(this) :: get => field_series_list_get
     !> Get the number of items in the list.
     procedure, pass(this) :: size => field_series_list_size
  end type field_series_list_t

contains

  !> Initialize a field series list with a given capacity
  subroutine field_series_list_init(this, capacity)
    class(field_series_list_t), intent(inout) :: this
    integer, intent(in) :: capacity

    if (capacity <= 0) then
       call neko_error('Field series list capacity must be positive')
    end if

    allocate(this%items(capacity))
    this%n_items = 0

  end subroutine field_series_list_init

  !> Free the field series list
  subroutine field_series_list_free(this)
    class(field_series_list_t), intent(inout) :: this

    if (allocated(this%items)) then
       deallocate(this%items)
    end if
    this%n_items = 0

  end subroutine field_series_list_free

  !> Add a field series to the list
  subroutine field_series_list_append(this, fld_series)
    class(field_series_list_t), intent(inout) :: this
    type(field_series_t), target, intent(in) :: fld_series
    type(field_series_ptr_t), allocatable :: tmp(:)
    integer :: len

    if (.not. allocated(this%items)) then
       allocate(this%items(1))
       this%n_items = 0
    end if

    if (this%n_items >= size(this%items)) then
       len = size(this%items)
       allocate(tmp(len + 1))
       tmp(1:len) = this%items
       call move_alloc(tmp, this%items)
    end if

    this%n_items = this%n_items + 1
    this%items(this%n_items)%ptr => fld_series

  end subroutine field_series_list_append

  !> Get a field series by index
  function field_series_list_get(this, index) result(field_series_ptr)
    class(field_series_list_t), intent(in) :: this
    integer, intent(in) :: index
    type(field_series_t), pointer :: field_series_ptr

    if (index < 1 .or. index > this%n_items) then
       call neko_error('Field series list index out of bounds')
    end if

    field_series_ptr => this%items(index)%ptr

  end function field_series_list_get

  !> Get the number of items in the list
  pure function field_series_list_size(this) result(n)
    class(field_series_list_t), intent(in) :: this
    integer :: n

    n = this%n_items

  end function field_series_list_size

end module field_series_list
