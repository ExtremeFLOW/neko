! Copyright (c) 2021-2026, The Neko Authors
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
!> Contains the `vector_series_t` type.
module vector_series
  use vector, only : vector_t
  implicit none
  private

  !> Stores a series (sequence) of vectors, logically connected to a base vector,
  !! and arranged according to some ordering.
  !! Currently used to store time-lagged values of solution vectors.
  type, public :: vector_series_t
     type(vector_t), pointer :: v => null()
     type(vector_t), allocatable :: lv(:)
     integer, private :: len = 0
   contains
     !> Constructor.
     procedure, pass(this) :: init => vector_series_init
     !> Destructor.
     procedure, pass(this) :: free => vector_series_free
     procedure, pass(this) :: update => vector_series_update
     procedure, pass(this) :: set => vector_series_set
     !> Return the size of the vector series.
     procedure, pass(this) :: size => vector_series_size
  end type vector_series_t

  !> A wrapper for a pointer to a `vector_series_t`.
  type, public :: vector_series_ptr_t
     type(vector_series_t), pointer :: ptr => null()
  end type vector_series_ptr_t

contains

  !> Initialize a vector series of length @a len for a vector @a v
  subroutine vector_series_init(this, v, len)
    class(vector_series_t), intent(inout) :: this
    type(vector_t), intent(inout), target :: v
    integer :: len
    character(len=80) :: name
    character(len=5) :: id_str
    integer :: i

    call this%free()

    this%v => v
    this%len = len

    allocate(this%lv(len))

    do i = 1, this%len
       write(id_str, '(I0)') i
       name = trim(this%v%name)//'_lag'//id_str
       call this%lv(i)%init(this%v%size(), name)
    end do

  end subroutine vector_series_init

  !> Deallocates a vector series
  subroutine vector_series_free(this)
    class(vector_series_t), intent(inout) :: this
    integer :: i

    if (associated(this%v)) then
       nullify(this%v)
    end if


    if (allocated(this%lv)) then
       do i = 1, this%len
          call this%lv(i)%free()
       end do

       deallocate(this%lv)
    end if

  end subroutine vector_series_free

  !> Return the size of the vector series
  function vector_series_size(this) result(len)
    class(vector_series_t), intent(in) :: this
    integer :: len
    len = this%len
  end function vector_series_size

  !> Update a vector series (evict oldest entry)
  subroutine vector_series_update(this)
    class(vector_series_t), intent(inout) :: this
    integer :: i

    do i = this%len, 2, -1
       this%lv(i) = this%lv(i-1)
    end do

    this%lv(1) = this%v

  end subroutine vector_series_update

  !> Set all vectors in a series to @a g
  subroutine vector_series_set(this, g)
    class(vector_series_t), intent(inout) :: this
    type(vector_t), intent(in) :: g
    integer :: i

    do i = 1, this%len
       this%lv(i) = g
    end do

  end subroutine vector_series_set

end module vector_series
