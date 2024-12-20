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
  use neko_config, only: NEKO_BCKND_DEVICE
  use math, only: sub3, chsign, add3, cmult2, cadd2, cfill, copy
  use num_types, only: rp
  use device, only: device_map, device_free, c_ptr, C_NULL_PTR
  use device_math, only: device_copy, device_cfill, device_cmult, &
       device_sub3, device_cmult2, device_add3, device_cadd2
  use utils, only: neko_error
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public ::  vector_t
     !> Vector entries.
     real(kind=rp), allocatable :: x(:)
     !> Device pointer.
     type(c_ptr) :: x_d = C_NULL_PTR
     !> Size of vector.
     integer :: n  = 0
   contains
     !> Initialise a vector of size `n`.
     procedure, pass(v) :: init => vector_init
     !> Deallocate a vector.
     procedure, pass(v) :: free => vector_free
     !> Returns the number of entries in the vector.
     procedure, pass(v) :: size => vector_size
     !> Assignment \f$ v = w \f$
     procedure, pass(v) :: vector_assign_vector
     !> Assignment \f$ v = s \f$.
     procedure, pass(v) :: vector_assign_scalar
     !> Vector-vector addition \f$ v = a + b \f$.
     procedure, pass(a) :: vector_add_vector
     !> Vector-scalar addition \f$ v = a + c \f$.
     procedure, pass(a) :: vector_add_scalar_left
     !> Scalar-vector addition \f$ v = c + a \f$.
     procedure, pass(a) :: vector_add_scalar_right
     !> Vector-vector subtraction \f$ v = a - b \f$.
     procedure, pass(a) :: vector_sub_vector
     !> Vector-scalar subtraction \f$ v = a - c \f$.
     procedure, pass(a) :: vector_sub_scalar_left
     !> Scalar-vector subtraction \f$ v = c - a \f$.
     procedure, pass(a) :: vector_sub_scalar_right
     !> Vector-scalar multiplication \f$ v = a*c \f$.
     procedure, pass(a) :: vector_cmult_left
     !> Scalar-vector multiplication \f$ v = c*a \f$.
     procedure, pass(a) :: vector_cmult_right

     generic :: assignment(=) => vector_assign_vector, &
          vector_assign_scalar
     generic :: operator(+) => vector_add_vector, &
          vector_add_scalar_left, vector_add_scalar_right
     generic :: operator(-) => vector_sub_vector, &
          vector_sub_scalar_left, vector_sub_scalar_right
     generic :: operator(*) => vector_cmult_left, vector_cmult_right
  end type vector_t

  type, public :: vector_ptr_t
     type(vector_t), pointer :: ptr
  end type vector_ptr_t

contains

  !> Initialise a vector of size @a n.
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

  !> Deallocate a vector.
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

  !> Return the number of entries in the vector.
  function vector_size(v) result(s)
    class(vector_t), intent(inout) :: v
    integer :: s
    s = v%n
  end function vector_size

  !> Assignment \f$ v = w \f$.
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

  !> Assignment \f$ v = s \f$.
  subroutine vector_assign_scalar(v, s)
    class(vector_t), intent(inout) :: v
    real(kind=rp), intent(in) :: s

    if (.not. allocated(v%x)) then
       call neko_error('Vector not allocated')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(v%x_d, s, v%n)
    else
       call cfill(v%x, s, v%n)
    end if

  end subroutine vector_assign_scalar

  !> Vector-vector addition \f$ v = a + b \f$.
  function vector_add_vector(a, b) result(v)
    class(vector_t), intent(in) :: a, b
    type(vector_t) :: v

    if (a%n .ne. b%n) call neko_error("Vectors must be the same length!")

    v%n = a%n
    allocate(v%x(v%n))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(v%x, v%x_d, v%n)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add3(v%x_d, a%x_d, b%x_d, v%n)
    else
       call add3(v%x, a%x, b%x, v%n)
    end if

  end function vector_add_vector

  !> Vector-scalar addition \f$ v = a + c \f$.
  function vector_add_scalar_left(a, c) result(v)
    class(vector_t), intent(in) :: a
    real(kind=rp), intent(in) :: c
    type(vector_t) :: v

    v%n = a%n
    allocate(v%x(v%n))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(v%x, v%x_d, v%n)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cadd2(v%x_d, a%x_d, c, v%n)
    else
       call cadd2(v%x, a%x, c, v%n)
    end if

  end function vector_add_scalar_left

  !> Scalar-vector addition \f$ v = c + a \f$.
  function vector_add_scalar_right(c, a) result(v)
    real(kind=rp), intent(in) :: c
    class(vector_t), intent(in) :: a
    type(vector_t) :: v

    v = vector_add_scalar_left(a, c)

  end function vector_add_scalar_right

  !> Vector-vector subtraction \f$ v = a - b \f$.
  function vector_sub_vector(a, b) result(v)
    class(vector_t), intent(in) :: a, b
    type(vector_t) :: v

    if (a%n .ne. b%n) call neko_error("Vectors must be the same length!")

    v%n = a%n
    allocate(v%x(v%n))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(v%x, v%x_d, v%n)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_sub3(v%x_d, a%x_d, b%x_d, v%n)
    else
       call sub3(v%x, a%x, b%x, v%n)
    end if

  end function vector_sub_vector

  !> Vector-scalar subtraction \f$ v = a - c \f$.
  function vector_sub_scalar_left(a, c) result(v)
    class(vector_t), intent(in) :: a
    real(kind=rp), intent(in) :: c
    type(vector_t) :: v

    v%n = a%n
    allocate(v%x(v%n))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(v%x, v%x_d, v%n)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cadd2(v%x_d, a%x_d, -1.0_rp*c, v%n)
    else
       call cadd2(v%x, a%x, -1.0_rp*c, a%n)
    end if

  end function vector_sub_scalar_left

  !> Scalar-vector subtraction \f$ v = c - a \f$.
  function vector_sub_scalar_right(c, a) result(v)
    real(kind=rp), intent(in) :: c
    class(vector_t), intent(in) :: a
    type(vector_t) :: v

    v = vector_sub_scalar_left(a, c)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(v%x_d, -1.0_rp, v%n)
    else
       v%x = -v%x
       !call chsign(v%x, v%n)
    end if

  end function vector_sub_scalar_right

  !> Vector-scalar multiplication \f$ v = a*c \f$.
  function vector_cmult_left(a, c) result(v)
    class(vector_t), intent(in) :: a
    real(kind=rp), intent(in) :: c
    type(vector_t) :: v

    v%n = a%n
    allocate(v%x(v%n))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(v%x, v%x_d, v%n)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult2(v%x_d, a%x_d, c, v%n)
    else
       call cmult2(v%x, a%x, c, v%n)
    end if

  end function vector_cmult_left

  !> Scalar-vector multiplication \f$ v = c*a \f$.
  function vector_cmult_right(c, a) result(v)
    real(kind=rp), intent(in) :: c
    class(vector_t), intent(in) :: a
    type(vector_t) :: v

    v = vector_cmult_left(a, c)

  end function vector_cmult_right

end module vector
