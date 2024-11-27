! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC.
!
! The UChicago Argonne, LLC as Operator of Argonne National
! Laboratory holds copyright in the Software. The copyright holder
! reserves all rights except those expressly granted to licensees,
! and U.S. Government license rights.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the disclaimer below.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the disclaimer (as noted below)
! in the documentation and/or other materials provided with the
! distribution.
!
! 3. Neither the name of ANL nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
! UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF
! ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Additional BSD Notice
! ---------------------
! 1. This notice is required to be provided under our contract with
! the U.S. Department of Energy (DOE). This work was produced at
! Argonne National Laboratory under Contract
! No. DE-AC02-06CH11357 with the DOE.
!
! 2. Neither the United States Government nor UCHICAGO ARGONNE,
! LLC nor any of their employees, makes any warranty,
! express or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.
!
! 3. Also, reference herein to any specific commercial products, process,
! or services by trade name, trademark, manufacturer or otherwise does
! not necessarily constitute or imply its endorsement, recommendation,
! or favoring by the United States Government or UCHICAGO ARGONNE LLC.
! The views and opinions of authors expressed
! herein do not necessarily state or reflect those of the United States
! Government or UCHICAGO ARGONNE, LLC, and shall
! not be used for advertising or product endorsement purposes.
!
module field_math
  use neko_config, only: NEKO_BCKND_DEVICE
  use num_types, only: rp
  use field, only: field_t
  use math, only: rzero, rone, copy, cmult, cadd, cfill, invcol1, vdot3, add2, &
       add3, add4, sub2, sub3, add2s1, add2s2, addsqr2s2, cmult2, invcol2, &
       col2, col3, subcol3, add3s2, addcol3, addcol4, glsum, glsc2, glsc3
  use device_math, only: device_rzero, device_rone, device_copy, device_cmult, &
       device_cadd, device_cfill, device_invcol1, device_vdot3, device_add2, &
       device_add3, device_add4, device_sub2, device_sub3, device_add2s1, &
       device_add2s2, device_addsqr2s2, device_cmult2, device_invcol2, &
       device_col2, device_col3, device_subcol3, device_add3s2, &
       device_addcol3, device_addcol4, device_glsum, device_glsc2, device_glsc3
  implicit none
  private

  public :: field_rzero, field_rone, field_copy, field_cmult, &
       field_cadd, field_cfill, field_invcol1, field_vdot3, &
       field_add2, field_sub2, field_sub3, field_add2s1, &
       field_add2s2, field_addsqr2s2, field_cmult2, &
       field_invcol2, field_col2, field_col3, field_subcol3, &
       field_add3s2, field_addcol3, field_addcol4, field_glsum, &
       field_glsc2, field_glsc3
  
contains

  !> Zero a real vector
  subroutine field_rzero(a, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(a%x_d, size)
    else
       call rzero(a%x, size)
    end if
  end subroutine field_rzero

  !> Set all elements to one
  subroutine field_rone(a, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rone(a%x_d, size)
    else
       call rone(a%x, size)
    end if
  end subroutine field_rone

  !> Copy a vector \f$ a = b \f$
  subroutine field_copy(a, b, n)
    integer, intent(in), optional :: n
    type(field_t), intent(in) :: b
    type(field_t), intent(inout) :: a
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(a%x_d, b%x_d, size)
    else
       call copy(a%x, b%x, size)
    end if
  end subroutine field_copy

  !> Multiplication by constant c \f$ a = c \cdot a \f$
  subroutine field_cmult(a, c, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(a%x_d, c, size)
    else
       call cmult(a%x, c, size)
    end if
  end subroutine field_cmult

  !> Add a scalar to vector \f$ a = \sum a_i + s \f$
  subroutine field_cadd(a, s, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    real(kind=rp), intent(in) :: s
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cadd(a%x_d, s, size)
    else
       call cadd(a%x, s, size)
    end if
  end subroutine field_cadd

  !> Set all elements to a constant c \f$ a = c \f$
  subroutine field_cfill(a, c, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(a%x_d, c, size)
    else
       call cfill(a%x, c, size)
    end if
  end subroutine field_cfill

  !> Invert a vector \f$ a = 1 / a \f$
  subroutine field_invcol1(a, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_invcol1(a%x_d, size)
    else
       call invcol1(a%x, size)
    end if

  end subroutine field_invcol1

  !> Compute a dot product \f$ dot = u \cdot v \f$ (3-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  subroutine field_vdot3(dot, u1, u2, u3, v1, v2, v3, n)
    integer, intent(in), optional :: n
    type(field_t), intent(in) :: u1, u2, u3
    type(field_t), intent(in) :: v1, v2, v3
    type(field_t), intent(out) :: dot
    integer :: size

    if (present(n)) then
       size = n
    else
       size = dot%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_vdot3(dot%x_d, &
            u1%x_d, u2%x_d, u3%x_d, &
            v1%x_d, v2%x_d, v3%x_d, &
            size)
    else
       call vdot3(dot%x, &
            u1%x, u2%x, u3%x, &
            v1%x, v2%x, v3%x, &
            size)
    end if

  end subroutine field_vdot3

  !> Vector addition \f$ a = a + b \f$
  subroutine field_add2(a, b, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(in) :: b
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(a%x_d, b%x_d, size)
    else
       call add2(a%x, b%x, size)
    end if

  end subroutine field_add2

  !> Vector addition \f$ a = b + c \f$
  subroutine field_add3(a, b, c, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(in) :: b, c
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add3(a%x_d, b%x_d, c%x_d, size)
    else
       call add3(a%x, b%x, c%x, size)
    end if

  end subroutine field_add3

  !> Vector addition \f$ a = b + c + d \f$
  subroutine field_add4(a, b, c, d, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(in) :: b, c, d
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add4(a%x_d, b%x_d, c%x_d, d%x_d, size)
    else
       call add4(a%x, b%x, c%x, d%x, size)
    end if

  end subroutine field_add4

  !> Vector substraction \f$ a = a - b \f$
  subroutine field_sub2(a, b, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(inout) :: b
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_sub2(a%x_d, b%x_d, size)
    else
       call sub2(a%x, b%x, size)
    end if

  end subroutine field_sub2

  !> Vector subtraction \f$ a = b - c \f$
  subroutine field_sub3(a, b, c, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: c
    type(field_t), intent(inout) :: b
    type(field_t), intent(out) :: a
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_sub3(a%x_d, b%x_d, c%x_d, size)
    else
       call sub3(a%x, b%x, c%x, size)
    end if

  end subroutine field_sub3


  !> Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
  !! (multiplication on first argument)
  subroutine field_add2s1(a, b, c1, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(inout) :: b
    real(kind=rp), intent(in) :: c1
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2s1(a%x_d, b%x_d, c1, size)
    else
       call add2s1(a%x, b%x, c1, size)
    end if

  end subroutine field_add2s1

  !> Vector addition with scalar multiplication  \f$ a = a + c_1 b \f$
  !! (multiplication on second argument)
  subroutine field_add2s2(a, b, c1, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(inout) :: b
    real(kind=rp), intent(in) :: c1
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2s2(a%x_d, b%x_d, c1, size)
    else
       call add2s2(a%x, b%x, c1, size)
    end if

  end subroutine field_add2s2

  !> Returns \f$ a = a + c1 * (b * b )\f$
  subroutine field_addsqr2s2(a, b, c1, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(in) :: b
    real(kind=rp), intent(in) :: c1
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_addsqr2s2(a%x_d, b%x_d, c1, size)
    else
       call addsqr2s2(a%x, b%x, c1, size)
    end if

  end subroutine field_addsqr2s2

  !> Multiplication by constant c \f$ a = c \cdot b \f$
  subroutine field_cmult2(a, b, c, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(in) :: b
    real(kind=rp), intent(in) :: c
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult2(a%x_d, b%x_d, c, size)
    else
       call cmult2(a%x, b%x, c, size)
    end if

  end subroutine field_cmult2

  !> Vector division \f$ a = a / b \f$
  subroutine field_invcol2(a, b, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(in) :: b
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_invcol2(a%x_d, b%x_d, size)
    else
       call invcol2(a%x, b%x, size)
    end if

  end subroutine field_invcol2


  !> Vector multiplication \f$ a = a \cdot b \f$
  subroutine field_col2(a, b, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(in) :: b
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(a%x_d, b%x_d, size)
    else
       call col2(a%x, b%x, size)
    end if

  end subroutine field_col2

  !> Vector multiplication with 3 vectors \f$ a =  b \cdot c \f$
  subroutine field_col3(a, b, c, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(in) :: b
    type(field_t), intent(in) :: c
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col3(a%x_d, b%x_d, c%x_d, size)
    else
       call col3(a%x, b%x, c%x, size)
    end if

  end subroutine field_col3

  !> Returns \f$ a = a - b*c \f$
  subroutine field_subcol3(a, b, c, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(in) :: b
    type(field_t), intent(in) :: c
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_subcol3(a%x_d, b%x_d, c%x_d, size)
    else
       call subcol3(a%x, b%x, c%x, size)
    end if

  end subroutine field_subcol3

  !> Returns \f$ a = c1 * b + c2 * c \f$
  subroutine field_add3s2(a, b, c, c1, c2, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(in) :: b
    type(field_t), intent(in) :: c
    real(kind=rp), intent(in) :: c1, c2
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add3s2(a%x_d, b%x_d, c%x_d, c1, c2, size)
    else
       call add3s2(a%x, b%x, c%x, c1, c2, size)
    end if

  end subroutine field_add3s2

  !> Returns \f$ a = a + b*c \f$
  subroutine field_addcol3(a, b, c, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(in) :: b
    type(field_t), intent(in) :: c
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_addcol3(a%x_d, b%x_d, c%x_d, size)
    else
       call addcol3(a%x, b%x, c%x, size)
    end if

  end subroutine field_addcol3

  !> Returns \f$ a = a + b*c*d \f$
  subroutine field_addcol4(a, b, c, d, n)
    integer, intent(in), optional :: n
    type(field_t), intent(inout) :: a
    type(field_t), intent(in) :: b
    type(field_t), intent(in) :: c
    type(field_t), intent(in) :: d
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_addcol4(a%x_d, b%x_d, c%x_d, d%x_d, size)
    else
       call addcol4(a%x, b%x, c%x, d%x, size)
    end if

  end subroutine field_addcol4

  function field_glsum(a, n) result(sum)
    integer, intent(in), optional :: n
    type(field_t), intent(in) :: a
    real(kind=rp) :: sum
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       sum = device_glsum(a%x_d, size)
    else
       sum = glsum(a%x, size)
    end if

  end function field_glsum

  function field_glsc2(a, b, n) result(norm)
    integer, intent(in), optional :: n
    type(field_t), intent(in) :: a, b
    real(kind=rp) :: norm
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       norm = device_glsc2(a%x_d, b%x_d, size)
    else
       norm = glsc2(a%x, b%x, size)
    end if

  end function field_glsc2

  function field_glsc3(a, b, c, n) result(norm)
    integer, intent(in), optional :: n
    type(field_t), intent(in) :: a, b, c
    real(kind=rp) :: norm
    integer :: size

    if (present(n)) then
       size = n
    else
       size = a%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       norm = device_glsc3(a%x_d, b%x_d, c%x_d, size)
    else
       norm = glsc3(a%x, b%x, c%x, size)
    end if

  end function field_glsc3

end module field_math
