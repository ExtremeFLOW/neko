! Copyright (c) 2020-2024, The Neko Authors
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
  use neko_config
  use num_types
  use device
  use, intrinsic :: iso_c_binding, only : c_ptr
  use bc, only : bc_t, bc_ptr_t
  implicit none
  private

  !> A list of allocatable @ref `bc_t`.
  !! Follows the standard interface of lists.
  type, public :: bc_list_t
     type(bc_ptr_t), allocatable :: items(:)
     !> Number of items.
     integer, private :: size_
     !> Capacity.
     integer :: capacity
   contains
     !> Constructor.
     procedure, pass(this) :: init => bc_list_init
     !> Destructor.
     procedure, pass(this) :: free => bc_list_free
     !> Append an item to the end of the list.
     procedure, pass(this) :: append => bc_list_append
     !> Apply all boundary conditions in the list.
     generic :: apply => apply_scalar, apply_vector
     !> Appply the boundary conditions to a scalar field.
     procedure, pass(this) :: apply_scalar => bc_list_apply_scalar
     !> Appply the boundary conditions to a vector field.
     procedure, pass(this) :: apply_vector => bc_list_apply_vector
     !> Return the number of items in the list.
     procedure, pass(this) :: size => bc_list_size
  end type bc_list_t

  interface bc_list_apply
     module procedure bc_list_apply_scalar, bc_list_apply_vector
  end interface bc_list_apply

contains

  !> Constructor.
  !! @param size The size of the list to allocate.
  subroutine bc_list_init(this, size)
    class(bc_list_t), intent(inout), target :: this
    integer, optional :: size
    integer :: n, i

    call this%free()

    if (present(size)) then
       n = size
    else
       n = 1
    end if

    allocate(this%items(n))

    do i = 1, n
       this%items(i)%ptr => null()
    end do

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
       do i =1, this%size()
!         call this%items(i)%ptr%free()
         nullify(this%items(i)%ptr)
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
       if (present(t) .and. present(tstep)) then
          do i = 1, this%size()
             call this%items(i)%ptr%apply_scalar_dev(x_d, t=t, tstep=tstep)
          end do
       else if (present(t)) then
          do i = 1, this%size()
             call this%items(i)%ptr%apply_scalar_dev(x_d, t=t)
          end do
       else if (present(tstep)) then
          do i = 1, this%size()
             call this%items(i)%ptr%apply_scalar_dev(x_d, tstep=tstep)
          end do
       else
          do i = 1, this%size()
             call this%items(i)%ptr%apply_scalar_dev(x_d)
          end do
       end if
    else
       if (present(t) .and. present(tstep)) then
          do i = 1, this%size()
             call this%items(i)%ptr%apply_scalar(x, n, t, tstep)
          end do
       else if (present(t)) then
          do i = 1, this%size()
             call this%items(i)%ptr%apply_scalar(x, n, t=t)
          end do
       else if (present(tstep)) then
          do i = 1, this%size()
             call this%items(i)%ptr%apply_scalar(x, n, tstep=tstep)
          end do
       else
          do i = 1, this%size()
             call this%items(i)%ptr%apply_scalar(x, n)
          end do
       end if
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
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
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
       if (present(t) .and. present(tstep)) then
          do i = 1, this%size()
             call this%items(i)%ptr%apply_vector_dev(x_d, y_d, z_d, t, tstep)
          end do
       else if (present(t)) then
          do i = 1, this%size()
             call this%items(i)%ptr%apply_vector_dev(x_d, y_d, z_d, t=t)
          end do
       else if (present(tstep)) then
          do i = 1, this%size()
             call this%items(i)%ptr%apply_vector_dev(x_d, y_d, z_d, tstep=tstep)
          end do
       else
          do i = 1, this%size()
             call this%items(i)%ptr%apply_vector_dev(x_d, y_d, z_d)
          end do
       end if
    else
       if (present(t) .and. present(tstep)) then
          do i = 1, this%size()
             call this%items(i)%ptr%apply_vector(x, y, z, n, t, tstep)
          end do
       else if (present(t)) then
          do i = 1, this%size()
             call this%items(i)%ptr%apply_vector(x, y, z, n, t=t)
          end do
       else if (present(tstep)) then
          do i = 1, this%size()
             call this%items(i)%ptr%apply_vector(x, y, z, n, tstep=tstep)
          end do
       else
          do i = 1, this%size()
             call this%items(i)%ptr%apply_vector(x, y, z, n)
          end do
       end if
    end if

  end subroutine bc_list_apply_vector

  !> Return the number of items in the list.
  pure function bc_list_size(this) result(size)
    class(bc_list_t), intent(in), target :: this
    integer :: size

    size = this%size_
  end function bc_list_size


end module bc_list
