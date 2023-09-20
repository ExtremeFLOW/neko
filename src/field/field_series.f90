! Copyright (c) 2021-2023, The Neko Authors
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
!> Stores a series fields
module field_series
  use num_types
  use field
  implicit none
  private
  
  type, public :: field_series_t
     type(field_t), pointer :: f => null()
     type(field_t), allocatable :: lf(:) 
     integer, private :: len = 0
   contains
     procedure, pass(this) :: init => field_series_init
     procedure, pass(this) :: free => field_series_free
     procedure, pass(this) :: update => field_series_update
     procedure, pass(this) :: set => field_series_set
     procedure, pass(this) :: size => field_series_size
  end type field_series_t

contains

  !> Initialize a field series of length @a len for a field @a f
  subroutine field_series_init(this, f, len)
    class(field_series_t), intent(inout) :: this
    type(field_t), intent(inout), target :: f
    integer :: len
    character(len=80) :: name
    integer :: i

    call this%free()

    this%f => f
    this%len = len

    allocate(this%lf(len))

    do i = 1, this%len
       name = trim(f%name)//'_lag'//char(i)
       call this%lf(i)%init(this%f%dof, name)
    end do

  end subroutine field_series_init

  !> Deallocates a field series
  subroutine field_series_free(this)
    class(field_series_t), intent(inout) :: this
    integer :: i

    if (associated(this%f)) then
       nullify(this%f)
    end if

    do i = 1, this%len
       call this%lf(i)%free()
    end do
    
  end subroutine field_series_free

  !> Return the size of the field series
  function field_series_size(this) result(len)
    class(field_series_t), intent(in) :: this
    integer :: len
    len = this%len
  end function field_series_size

  !> Update a field series (evict oldest entry)
  subroutine field_series_update(this)
    class(field_series_t), intent(inout) :: this
    integer :: i
    !$omp parallel private(i)
    do i = this%len, 2, -1
       this%lf(i) = this%lf(i-1)
    end do

    this%lf(1) = this%f
    !$omp end parallel
    
  end subroutine field_series_update

  !> Set all fields in a series to @a g
  subroutine field_series_set(this, g)
    class(field_series_t), intent(inout) :: this
    type(field_t), intent(in) :: g
    integer :: i
    !$omp parallel
    do i = 1, this%len
       this%lf(i) = g
    end do
    !$omp end parallel
  end subroutine field_series_set
  
end module field_series
