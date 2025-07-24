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
  use math, only: sub3, chsign, add3, cmult2, cadd2, cfill, copy, col3, cdiv2, &
       col2, invcol3
  use num_types, only: rp
  use device, only: device_map, device_free
  use device_math, only: device_copy, device_cfill, device_cmult, &
       device_sub3, device_cmult2, device_add3, device_cadd2, device_col3, &
       device_col2, device_invcol3, device_cdiv2
  use utils, only: neko_error
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public :: vector_t
     !> Vector entries.
     real(kind=rp), allocatable :: x(:)
     !> Device pointer.
     type(c_ptr) :: x_d = C_NULL_PTR
     !> Size of vector.
     integer, private :: n = 0
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
     !> Pointwise vector multiplication \f$ v = a*b \f$.
     procedure, pass(a) :: vector_pointwise_mult
     !> Pointwise vector power \f$ v = a*b \f$.
     procedure, pass(a) :: vector_pointwise_power
     !> Scalar-vector division \f$ v = c / a \f$.
     procedure, pass(a) :: vector_cdiv_left
     !> Vector-scalar division \f$ v = a / c \f$.
     procedure, pass(a) :: vector_cdiv_right
     !> Pointwise vector division \f$ v = a / b \f$.
     procedure, pass(a) :: vector_pointwise_div
     !> Change the sign of the vector.
     procedure, pass(a) :: vector_chsign

     generic :: assignment(=) => vector_assign_vector, &
          vector_assign_scalar
     generic :: operator(+) => vector_add_vector, &
          vector_add_scalar_left, vector_add_scalar_right
     generic :: operator(-) => vector_sub_vector, &
          vector_sub_scalar_left, vector_sub_scalar_right, &
          vector_chsign
     generic :: operator(*) => vector_cmult_left, vector_cmult_right, &
          vector_pointwise_mult
     generic :: operator(/) => vector_cdiv_left, vector_cdiv_right, &
          vector_pointwise_div
     ! Seems to crash the cray compiler
     !generic :: operator(**) => vector_pointwise_power

     ! Private interfaces
     procedure, pass(a), private :: alloc => vector_allocate

  end type vector_t

  type, public :: vector_ptr_t
     type(vector_t), pointer :: ptr
  end type vector_ptr_t

contains

  !> Initialise a vector of size @a n.
  subroutine vector_init(v, n)
    class(vector_t), intent(inout) :: v
    integer, intent(in) :: n

    call v%alloc(n)
    v%x = 0.0_rp
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(v%x_d, 0.0_rp, n)
    end if

  end subroutine vector_init

  !> Vector allocation without initialisation.
  subroutine vector_allocate(a, n)
    class(vector_t), intent(inout) :: a
    integer, intent(in) :: n

    if (n .eq. 0) call neko_error('Vector cannot have size 0')
    call a%free()

    a%n = n
    allocate(a%x(n))
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(a%x, a%x_d, n)
    end if

  end subroutine vector_allocate

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
  pure function vector_size(v) result(s)
    class(vector_t), intent(in) :: v
    integer :: s
    s = v%n
  end function vector_size

  !> Assignment \f$ v = w \f$.
  subroutine vector_assign_vector(v, w)
    class(vector_t), intent(inout) :: v
    type(vector_t), intent(in) :: w

    if (v%n .ne. w%n) call v%alloc(w%n)

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

    if (v%n .eq. 0) call neko_error('Vector not allocated')

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

    if (a%n .ne. b%n) call neko_error("Vectors must be the same length")

    call v%alloc(a%n)

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

    call v%alloc(a%n)

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

    if (a%n .ne. b%n) call neko_error("Vectors must be the same length")

    call v%alloc(a%n)

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

    call v%alloc(a%n)

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
    end if

  end function vector_sub_scalar_right

  !> Change the sign of the vector.
  function vector_chsign(a) result(v)
    class(vector_t), intent(in) :: a
    type(vector_t) :: v

    call v%alloc(a%n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       v = a
       call device_cmult(v%x_d, -1.0_rp, v%n)
    else
       v%x = -a%x
    end if

  end function vector_chsign

  !> Vector-scalar multiplication \f$ v = a*c \f$.
  function vector_cmult_left(a, c) result(v)
    class(vector_t), intent(in) :: a
    real(kind=rp), intent(in) :: c
    type(vector_t) :: v

    call v%alloc(a%n)

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

  !> Pointwise vector multiplication \f$ v = a*b \f$.
  function vector_pointwise_mult(a, b) result(v)
    class(vector_t), intent(in) :: a, b
    type(vector_t) :: v

    if (a%n .ne. b%n) call neko_error("Vectors must be the same length")

    call v%alloc(a%n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col3(v%x_d, a%x_d, b%x_d, v%n)
    else
       call col3(v%x, a%x, b%x, v%n)
    end if

  end function vector_pointwise_mult

  !> Pointwise power \f$ v = a^b \f$. OBS integer b
  !! Todo: Incredibly poor performance, needs to be optimized.
  function vector_pointwise_power(a, b) result(v)
    class(vector_t), intent(in) :: a
    integer, intent(in) :: b
    type(vector_t) :: v
    integer :: i

    call v%alloc(a%n)
    v = 1.0_rp
    if (b .eq. 0) then
       return
    end if

    do i = 1, b
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(v%x_d, a%x_d, v%n)
       else
          call col2(v%x, a%x, v%n)
       end if
    end do

  end function vector_pointwise_power

  !> Scalar-vector division \f$ v = c / a \f$.
  function vector_cdiv_left(c, a) result(v)
    real(kind=rp), intent(in) :: c
    class(vector_t), intent(in) :: a
    type(vector_t) :: v

    call v%alloc(a%n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cdiv2(v%x_d, a%x_d, c, v%n)
    else
       call cdiv2(v%x, a%x, c, v%n)
    end if

  end function vector_cdiv_left

  !> Vector-scalar division \f$ v = a / c \f$.
  function vector_cdiv_right(a, c) result(v)
    class(vector_t), intent(in) :: a
    real(kind=rp), intent(in) :: c
    type(vector_t) :: v

    call v%alloc(a%n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult2(v%x_d, a%x_d, 1.0_rp / c, v%n)
    else
       call cmult2(v%x, a%x, 1.0_rp / c, v%n)
    end if

  end function vector_cdiv_right

  !> Pointwise vector division \f$ v = a / b \f$.
  function vector_pointwise_div(a, b) result(v)
    class(vector_t), intent(in) :: a, b
    type(vector_t) :: v

    if (a%n .ne. b%n) call neko_error("Vectors must be the same length")

    call v%alloc(a%n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_invcol3(v%x_d, a%x_d, b%x_d, v%n)
    else
       call invcol3(v%x, a%x, b%x, v%n)
    end if

  end function vector_pointwise_div

end module vector
