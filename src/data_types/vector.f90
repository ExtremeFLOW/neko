! Copyright (c) 2022-2025, The Neko Authors
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
!> Defines a vector
module vector
  use array, only : array_t
  use num_types, only : rp
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  implicit none
  private

  type, public, extends(array_t) :: vector_t
     !> Vector shaped pointer.
     real(kind=rp), pointer, dimension(:) :: x => null()
   contains
     !> Initialise a vector of size `n`.
     procedure, pass(this) :: init => vector_init
     !> Free the vector
     procedure, pass(this) :: free => vector_free

  end type vector_t

  type, public :: vector_ptr_t
     type(vector_t), pointer :: ptr => null()
   contains
     !> Constructor. Just assigns the pointer
     procedure, pass(this) :: init => vector_ptr_init
     !> Destructor. Just nullifies the pointer.
     procedure, pass(this) :: free => vector_ptr_free
  end type vector_ptr_t

contains

  !> Initialise a vector of size @a n.
  subroutine vector_init(this, n, name)
    class(vector_t), intent(inout), target :: this
    integer, intent(in) :: n
    character(len=*), intent(in), optional :: name

    call this%init_base(n, name)

    this%x => this%data(:)

  end subroutine vector_init

  !> Free the vector
  subroutine vector_free(this)
    class(vector_t), intent(inout) :: this

    call this%free_base()
    nullify(this%x)

  end subroutine vector_free

  ! ========================================================================== !
  ! vector pointer type subroutines

  subroutine vector_ptr_init(this, ptr)
    class(vector_ptr_t), intent(inout) :: this
    type(vector_t), target, intent(in) :: ptr

    call this%free()
    this%ptr => ptr
  end subroutine vector_ptr_init

  subroutine vector_ptr_free(this)
    class(vector_ptr_t), intent(inout) :: this

    if (associated(this%ptr)) then
       nullify(this%ptr)
    end if

  end subroutine vector_ptr_free

end module vector
