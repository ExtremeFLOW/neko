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
!> Defines a list of `bc_t`.
module bc_list
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use device, only : device_get_ptr
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr
  use bc, only : bc_t, bc_ptr_t
  implicit none
  private

  !> A list of allocatable @ref `bc_t`.
  !! Follows the standard interface of lists.
  type, public :: bc_list_t
     ! The items of the list.
     class(bc_ptr_t), allocatable, private :: items(:)
     !> Number of items in the list that are themselves allocated.
     integer, private :: size_
     !> Capacity, i.e. the size_ of the items list. Some items may themselves be
     !! unallocated.
     integer, private :: capacity
   contains
     !> Constructor.
     procedure, pass(this) :: init => bc_list_init
     !> Destructor.
     procedure, pass(this) :: free => bc_list_free
     !> Append an item to the end of the list.
     procedure, pass(this) :: append => bc_list_append
     !> Get the item at the given index.
     procedure, pass(this) :: get => bc_list_get
     !> Apply all boundary conditions in the list.
     generic :: apply => apply_scalar, apply_vector
     !> Appply the boundary conditions to a scalar field.
     procedure, pass(this) :: apply_scalar => bc_list_apply_scalar
     !> Appply the boundary conditions to a vector field.
     procedure, pass(this) :: apply_vector => bc_list_apply_vector
     !> Check wether the list is empty
     procedure, pass(this) :: is_empty => bc_list_is_empty
     !> Return the number of items in the list.
     procedure :: size => bc_list_size
  end type bc_list_t

contains

  !> Constructor.
  !! @param size The size of the list to allocate.
  subroutine bc_list_init(this, size)
    class(bc_list_t), intent(inout), target :: this
    integer, optional :: size
    integer :: n

    call this%free()

    n = 1
    if (present(size)) n = size

    allocate(this%items(n))

    this%size_ = 0
    this%capacity = n

  end subroutine bc_list_init

  !> Destructor.
  !! @note This will only nullify all pointers, not deallocate any
  !! conditions pointed to by the list
  subroutine bc_list_free(this)
    class(bc_list_t), intent(inout) :: this
    integer :: i

    if (allocated(this%items)) then
       do i = 1, this%size_
          this%items(i)%ptr => null()
       end do

       deallocate(this%items)
    end if

    this%size_ = 0
    this%capacity = 0
  end subroutine bc_list_free

  !> Append a condition to the end of the list.
  !! @param bc The boundary condition to add.
  subroutine bc_list_append(this, bc)
    class(bc_list_t), intent(inout) :: this
    class(bc_t), intent(inout), target :: bc
    type(bc_ptr_t), allocatable :: tmp(:)

    !> Do not add if bc is empty
    if(bc%marked_facet%size() .eq. 0) return

    if (this%size_ .ge. this%capacity) then
       this%capacity = this%capacity * 2
       allocate(tmp(this%capacity))
       tmp(1:this%size_) = this%items
       call move_alloc(tmp, this%items)
    end if

    this%size_ = this%size_ + 1
    this%items(this%size_)%ptr => bc

  end subroutine bc_list_append

  !> Get the item at the given index.
  !! @param i The index of the item to get.
  !! @return The item at the given index.
  function bc_list_get(this, i) result(bc)
    class(bc_list_t), intent(in) :: this
    class(bc_t), pointer :: bc
    integer, intent(in) :: i

    if (i .lt. 1 .or. i .gt. this%size_) then
       call neko_error("Index out of bounds in bc_list_get")
    end if

    bc => this%items(i)%ptr

  end function bc_list_get

  !> Apply a list of boundary conditions to a scalar field
  !! @param x The field to apply the boundary conditions to.
  !! @param n The size of x.
  !! @param t Current time.
  !! @param tstep Current time-step.
  subroutine bc_list_apply_scalar(this, x, n, t, tstep)
    class(bc_list_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    type(c_ptr) :: x_d
    integer :: i

    if (NEKO_BCKND_DEVICE .eq. 1) then
       x_d = device_get_ptr(x)
       do i = 1, this%size_
          call this%items(i)%ptr%apply_scalar_dev(x_d, t, tstep)
       end do
    else
       do i = 1, this%size_
          call this%items(i)%ptr%apply_scalar(x, n, t, tstep)
       end do
    end if
  end subroutine bc_list_apply_scalar

  !> Apply a list of boundary conditions to a vector field.
  !! @param x The x comp of the field for which to apply the bcs.
  !! @param y The y comp of the field for which to apply the bcs.
  !! @param z The z comp of the field for which to apply the bcs.
  !! @param n The size of x, y, z.
  !! @param t Current time.
  !! @param tstep Current time-step.
  subroutine bc_list_apply_vector(this, x, y, z, n, t, tstep)
    class(bc_list_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    integer :: i

    if (NEKO_BCKND_DEVICE .eq. 1) then
       x_d = device_get_ptr(x)
       y_d = device_get_ptr(y)
       z_d = device_get_ptr(z)
       do i = 1, this%size_
          call this%items(i)%ptr%apply_vector_dev(x_d, y_d, z_d, t, tstep)
       end do
    else
       do i = 1, this%size_
          call this%items(i)%ptr%apply_vector(x, y, z, n, t, tstep)
       end do
    end if

  end subroutine bc_list_apply_vector

  !> Return whether the list is empty.
  function bc_list_is_empty(this) result(is_empty)
    class(bc_list_t), intent(in), target :: this
    logical :: is_empty
    integer :: i

    is_empty = .true.
    do i = 1, this%size_

       if (.not. allocated(this%items(i)%ptr%msk)) then
          call neko_error("bc not finalized, error in bc_list%is_empty")
       end if

       if (this%items(i)%ptr%msk(0) > 0) is_empty = .false.

    end do
  end function bc_list_is_empty

  !> Return the number of items in the list.
  pure function bc_list_size(this) result(size)
    class(bc_list_t), intent(in), target :: this
    integer :: size

    size = this%size_
  end function bc_list_size

end module bc_list
