! Copyright (c) 2022, The Neko Authors
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
  use neko_config
  use num_types
  use device
  use device_math
  use utils
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  type, public ::  vector_t
     real(kind=rp), allocatable :: x(:) !< Vector entries
     type(c_ptr) :: x_d = C_NULL_PTR    !< Device pointer
     integer :: n  = 0                  !< Size of vector
   contains
     procedure, pass(v) :: init => vector_init
     procedure, pass(v) :: free => vector_free
     procedure, pass(v) :: size => vector_size
     procedure, pass(v) :: vector_assign_vector
     procedure, pass(v) :: vector_assign_scalar
     generic :: assignment(=) => vector_assign_vector, &
          vector_assign_scalar
  end type vector_t

  type, public :: vector_ptr_t
     type(vector_t), pointer :: v
  end type

contains

  !> Initialise a vector of size @a n
  subroutine vector_init(v, n)
    class(vector_t), intent(inout) :: v
    integer, intent(in) :: n

    call v%free()

    allocate(v%x(n))
    v%x = 0.0_rp
   
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(v%x, v%x_d, n)
       call device_cfill(v%x_d, 0.0_rp, n)
    end if

    v%n = n
    
  end subroutine vector_init

  !> Deallocate a vector
  subroutine vector_free(v)
    class(vector_t), intent(inout) :: v

    if (allocated(v%x)) then
       deallocate(v%x)
    end if

    if (c_associated(v%x_d)) then
       call device_free(v%x_d)
    end if

    v%n = 0
        
  end subroutine vector_free

  !> Return the number of entries in the vector
  function vector_size(v) result(s)
    class(vector_t), intent(inout) :: v
    integer :: s
    s = v%n    
  end function vector_size

  !> Assignment \f$ v = w \f$
  subroutine vector_assign_vector(v, w)
    class(vector_t), intent(inout) :: v
    type(vector_t), intent(in) :: w

    if (allocated(v%x)) then
       call v%free()
    end if

    if (.not. allocated(v%x)) then

       v%n = w%n
       allocate(v%x(v%n))
       
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(v%x, v%x_d, v%n)
       end if
       
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(v%x_d, w%x_d, v%n)
    else
       v%x = w%x
    end if

  end subroutine vector_assign_vector

  !> Assignment \f$ v = s \f$
  subroutine vector_assign_scalar(v, s)
    class(vector_t), intent(inout) :: v
    real(kind=rp), intent(in) :: s

    if (.not. allocated(v%x)) then
       call neko_error('Vector not allocated')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(v%x_d, s, v%n)
    else
       v%x = s
    end if

  end subroutine vector_assign_scalar

  
end module vector
