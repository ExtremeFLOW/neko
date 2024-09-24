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
module host_math
  use num_types, only : rp, dp, sp, qp, i4
  use comm
  implicit none

  !> Machine epsilon \f$ \epsilon \f$
  real(kind=rp), private, parameter :: NEKO_EPS = epsilon(1.0_rp)

  !> \f$ ln(2) \f$
  real(kind=rp), private, parameter :: NEKO_M_LN2 = log(2.0_rp)

  !> \f$ \pi \f$
  real(kind=rp), private, parameter :: pi = 4._rp*atan(1._rp)

  interface host_abscmp
     module procedure host_sabscmp, host_dabscmp, host_qabscmp
  end interface host_abscmp

  interface host_host_sort
     module procedure host_sortrp, host_sorti4
  end interface host_host_sort

  interface host_swap
     module procedure host_swapdp, host_swapi4
  end interface host_swap

  interface host_reord
     module procedure host_reorddp, host_reordi4
  end interface host_reord

  interface host_flipv
     module procedure host_flipvdp, host_flipvi4
  end interface host_flipv

  interface host_relcmp
     module procedure host_srelcmp, host_drelcmp, host_qrelcmp
  end interface host_relcmp

contains

  !> Return single precision absolute comparison \f$ | x - y | < \epsilon \f$
  pure function host_sabscmp(x, y)
    real(kind=sp), intent(in) :: x
    real(kind=sp), intent(in) :: y
    logical :: host_sabscmp

    host_sabscmp = abs(x - y) .lt. NEKO_EPS

  end function host_sabscmp

  !> Return double precision absolute comparison \f$ | x - y | < \epsilon \f$
  pure function host_dabscmp(x, y)
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(in) :: y
    logical :: host_dabscmp

    host_dabscmp = abs(x - y) .lt. NEKO_EPS

  end function host_dabscmp

  !> Return double precision absolute comparison \f$ | x - y | < \epsilon \f$
  pure function host_qabscmp(x, y)
    real(kind=qp), intent(in) :: x
    real(kind=qp), intent(in) :: y
    logical :: host_qabscmp

    host_qabscmp = abs(x - y) .lt. NEKO_EPS

  end function host_qabscmp

  !> Return single precision relative comparison \f$ | x - y |<= \epsilon*|y| \f$
  pure function host_srelcmp(x, y, eps)
    real(kind=sp), intent(in) :: x
    real(kind=sp), intent(in) :: y
    real(kind=sp), intent(in), optional :: eps
    logical :: host_srelcmp

    if (present(eps)) then
       host_srelcmp = abs(x - y) .le. abs(y) * eps
    else
       host_srelcmp = abs(x - y) .le. abs(y) * NEKO_EPS
    end if

  end function host_srelcmp

  !> Return double precision relative comparison \f$ | x - y |/|y| < \epsilon \f$
  pure function host_drelcmp(x, y, eps)
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(in) :: y
    real(kind=dp), intent(in), optional :: eps
    logical :: host_drelcmp

    if (present(eps)) then
       host_drelcmp = abs(x - y) .le. abs(y) * eps
    else
       host_drelcmp = abs(x - y) .le. abs(y) * NEKO_EPS
    end if

  end function host_drelcmp


  !> Return quad precision relative comparison \f$ | x - y |/|y| < \epsilon \f$
  pure function host_qrelcmp(x, y, eps)
    real(kind=qp), intent(in) :: x
    real(kind=qp), intent(in) :: y
    real(kind=qp), intent(in), optional :: eps
    logical :: host_qrelcmp
    if (present(eps)) then
       host_qrelcmp = abs(x - y) .le. abs(y) * eps
    else
       host_qrelcmp = abs(x - y) .le. abs(y) * NEKO_EPS
    end if

  end function host_qrelcmp

  !> Zero a real vector
  subroutine host_rzero(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a

    a = 0.0_rp

  end subroutine host_rzero

  !> Zero an integer vector
  subroutine host_izero(a, n)
    integer, intent(in) :: n
    integer, dimension(n), intent(inout) :: a

    a = 0

  end subroutine host_izero

  !> Sets row e to 0 in matrix a
  subroutine host_row_zero(a, m, n, e)
    integer, intent(in) :: m, n, e
    real(kind=rp), intent(inout) :: a(m,n)
    integer :: j

    do j = 1, n
       a(e,j) = 0.0_rp
    end do
  end subroutine host_row_zero

  !> Set all elements to one
  subroutine host_rone(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a

    a = 1.0_rp

  end subroutine host_rone

  !> Copy a vector \f$ a = b \f$
  subroutine host_copy(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(out) :: a

    a = b

  end subroutine host_copy

  !> Copy a masked vector \f$ a(mask) = b(mask) \f$.
  !! @param a Destination array of size `n`.
  !! @param b Source array of size `n`.
  !! @param mask Mask array of length m+1, where `mask(0) =m`
  !! the length of the mask array.
  !! @param n Size of the arrays `a` and `b`.
  !! @param m Size of the mask array `mask`.
  subroutine host_masked_copy(a, b, mask, n, m)
    integer, intent(in) :: n, m
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(out) :: a
    integer, dimension(0:m) :: mask
    integer :: i, j

    do i = 1, m
       j = mask(i)
       a(j) = b(j)
    end do

  end subroutine host_masked_copy

  !> Copy a masked vector to reduced contigous vector
  !! \f$ a = b(mask) \f$.
  !! @param a Destination array of size `m`.
  !! @param b Source array of size `n`.
  !! @param mask Mask array of length m+1, where `mask(0) =m`
  !! the length of the mask array.
  !! @param n Size of the array `b`.
  !! @param m Size of the mask array `mask` and `a`.
  subroutine host_masked_red_copy(a, b, mask, n, m)
    integer, intent(in) :: n, m
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(m), intent(inout) :: a
    integer, dimension(0:m) :: mask
    integer :: i, j

    do i = 1, m
       j = mask(i)
       a(i) = b(j)
    end do

  end subroutine host_masked_red_copy

  !> @brief Fill a constant to a masked vector.
  !! \f$ a_i = c, for i in mask \f$
  subroutine host_cfill_mask(a, c, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = c
    end do

  end subroutine host_cfill_mask

  !> Multiplication by constant c \f$ a = c \cdot a \f$
  subroutine host_cmult(a, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), intent(in) :: c

    a = c * a

  end subroutine host_cmult

  !> Add a scalar to vector \f$ a_i = a_i + s \f$
  subroutine host_cadd(a, s, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), intent(in) :: s

    a = a + s

  end subroutine host_cadd

  !> Add a scalar to vector \f$ a_i = b_i + s \f$
  subroutine host_cadd2(a, b, s, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: s
    integer :: i

    do i = 1, n
       a(i) = b(i) + s
    end do
  end subroutine host_cadd2

  !> Set all elements to a constant c \f$ a = c \f$
  subroutine host_cfill(a, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: a
    real(kind=rp), intent(in) :: c

    a = c

  end subroutine host_cfill

  !> Sum a vector of length n
  function host_glsum(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: a
    real(kind=rp) :: tmp, host_glsum
    integer :: ierr

    tmp = sum(a)
    call MPI_Allreduce(tmp, host_glsum, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

  end function host_glsum

  !> Max of a vector of length n
  function host_glmax(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: a
    real(kind=rp) :: tmp, host_glmax
    integer :: ierr

    tmp = maxval(a)
    call MPI_Allreduce(tmp, host_glmax, 1, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
  end function host_glmax

  !> Max of an integer vector of length n
  function host_glimax(a, n)
    integer, intent(in) :: n
    integer, dimension(n), intent(in) :: a
    integer :: tmp, host_glimax
    integer :: ierr

    tmp = maxval(a)
    call MPI_Allreduce(tmp, host_glimax, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)
  end function host_glimax

  !> Min of a vector of length n
  function host_glmin(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: a
    real(kind=rp) :: tmp, host_glmin
    integer :: ierr

    tmp = minval(a)
    call MPI_Allreduce(tmp, host_glmin, 1, &
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
  end function host_glmin

  !> Min of an integer vector of length n
  function host_glimin(a, n)
    integer, intent(in) :: n
    integer, dimension(n), intent(in) :: a
    integer :: tmp, host_glimin
    integer :: ierr

    tmp = minval(a)
    call MPI_Allreduce(tmp, host_glimin, 1, &
         MPI_INTEGER, MPI_MIN, NEKO_COMM, ierr)

  end function host_glimin




  !> Change sign of vector \f$ a = -a \f$
  subroutine host_chsign(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    integer :: i

    a = -a

  end subroutine host_chsign

  !> maximum value of a vector of length @a n
  function host_vlmax(vec, n) result(tmax)
    integer :: n, i
    real(kind=rp), intent(in) :: vec(n)
    real(kind=rp) :: tmax

    tmax = maxval(vec)

  end function host_vlmax

  !> minimun value of a vector of length @a n
  function host_vlmin(vec, n) result(tmin)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: vec(n)
    real(kind=rp) :: tmin

    tmin = minval(vec)

  end function host_vlmin

  !> Invert a vector \f$ a = 1 / a \f$
  subroutine host_invcol1(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    integer :: i

    a = 1.0_rp / a

  end subroutine host_invcol1

  !> Invert a vector \f$ a = b / c \f$
  subroutine host_invcol3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c

    a = b / c

  end subroutine host_invcol3

  !> Compute inverted vector \f$ a = 1 / b \f$
  subroutine host_invers2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: a
    real(kind=rp), dimension(n), intent(in) :: b

    a = 1.0_rp / b

  end subroutine host_invers2

  !> Compute a cross product \f$ u = v \times w \f$
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  subroutine host_vcross(u1, u2, u3, v1, v2, v3, w1, w2, w3, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: v1, v2, v3
    real(kind=rp), dimension(n), intent(in) :: w1, w2, w3
    real(kind=rp), dimension(n), intent(out) :: u1, u2, u3
    integer :: i

    u1 = v2 * w3 - v3 * w2
    u2 = v3 * w1 - v1 * w3
    u3 = v1 * w2 - v2 * w1

  end subroutine host_vcross

  !> Compute a dot product \f$ dot = u \cdot v \f$ (2-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  subroutine host_vdot2(dot, u1, u2, v1, v2, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: u1, u2
    real(kind=rp), dimension(n), intent(in) :: v1, v2
    real(kind=rp), dimension(n), intent(out) :: dot

    dot = u1 * v1 + u2 * v2

  end subroutine host_vdot2

  !> Compute a dot product \f$ dot = u \cdot v \f$ (3-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  subroutine host_vdot3(dot, u1, u2, u3, v1, v2, v3, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: u1, u2, u3
    real(kind=rp), dimension(n), intent(in) :: v1, v2, v3
    real(kind=rp), dimension(n), intent(out) :: dot

    dot = u1 * v1 + u2 * v2 + u3 * v3

  end subroutine host_vdot3

  !> Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
  function host_vlsc3(u, v, w, n) result(s)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: u, v, w
    real(kind=rp) :: s

    s = sum(u * v * w)

  end function host_vlsc3

  !> Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
  function host_vlsc2(u, v, n) result(s)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: u, v
    real(kind=rp) :: s
    integer :: i

    s = sum(u * v)

  end function host_vlsc2

  !> Vector addition \f$ a = a + b \f$
  subroutine host_add2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b

    a = a + b

  end subroutine host_add2

  !> Vector addition \f$ a = b + c \f$
  subroutine host_add3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c

    a = b + c

  end subroutine host_add3

  !> Vector addition \f$ a = b + c + d\f$
  subroutine host_add4(a, b, c, d, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c, d

    a = b + c + d

  end subroutine host_add4

  !> Vector substraction \f$ a = a - b \f$
  subroutine host_sub2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b

    a = a - b

  end subroutine host_sub2

  !> Vector subtraction \f$ a = b - c \f$
  subroutine host_sub3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c

    a = b - c

  end subroutine host_sub3


  !> Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
  !! (multiplication on first argument)
  subroutine host_add2s1(a, b, c1, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: c1

    a = c1 * a + b

  end subroutine host_add2s1

  !> Vector addition with scalar multiplication \f$ a = a + c_1 b \f$
  !! (multiplication on second argument)
  subroutine host_add2s2(a, b, c1, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: c1

    a = a + c1 * b

  end subroutine host_add2s2

  !> Returns \f$ a = a + c1 * (b * b)\f$
  subroutine host_addsqr2s2(a, b, c1, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: c1

    a = a + c1 * (b * b)

  end subroutine host_addsqr2s2

  !> Multiplication by constant c \f$ a = c \cdot b \f$
  subroutine host_cmult2(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: c

    a = c * b

  end subroutine host_cmult2

  !> Vector division \f$ a = a / b \f$
  subroutine host_invcol2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b

    a = a / b

  end subroutine host_invcol2


  !> Vector multiplication \f$ a = a \cdot b \f$
  subroutine host_col2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b

    a = a * b

  end subroutine host_col2

  !> Vector multiplication with 3 vectors \f$ a = b \cdot c \f$
  subroutine host_col3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c

    a = b * c

  end subroutine host_col3

  !> Returns \f$ a = a - b*c \f$
  subroutine host_subcol3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c

    a = a - b * c

  end subroutine host_subcol3

  !> Returns \f$ a = c1 * b + c2 * c \f$
  subroutine host_add3s2(a, b, c, c1, c2 ,n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c
    real(kind=rp), intent(in) :: c1, c2

    a = c1 * b + c2 * c

  end subroutine host_add3s2


  !> Returns \f$ a = a - b*c*d \f$
  subroutine host_subcol4(a, b, c, d, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c, d

    a = a - b * c * d

  end subroutine host_subcol4

  !> Returns \f$ a = a + b*c \f$
  subroutine host_addcol3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c

    a = a + b * c

  end subroutine host_addcol3

  !> Returns \f$ a = a + b*c*d \f$
  subroutine host_addcol4(a, b, c, d, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c, d

    a = a + b * c * d

  end subroutine host_addcol4

  !> Returns \f$ a = b \dot c - d \cdot e \f$
  subroutine host_ascol5(a, b, c, d, e, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c, d, e

    a = b * c - d * e

  end subroutine host_ascol5

  !> Returns \f$ a = b \dot c1 (a - c2 \cdot c)\f$
  subroutine host_p_update(a, b, c, c1, c2, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c
    real(kind=rp), intent(in) :: c1, c2

    a = b + c1 * (a - c2 * c)

  end subroutine host_p_update

  !> Returns \f$ a = b \dot c1 (a - c2 \cdot c)\f$
  subroutine host_x_update(a, b, c, c1, c2, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c
    real(kind=rp), intent(in) :: c1, c2

    a = a + c1 * b + c2 * c

  end subroutine host_x_update

  !> Weighted inner product \f$ a^T b \f$
  function host_glsc2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: a, b
    real(kind=rp) :: host_glsc2, tmp
    integer :: ierr

    tmp = sum(a * b)

    call MPI_Allreduce(tmp, host_glsc2, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

  end function host_glsc2

  !> Weighted inner product \f$ a^T b c \f$
  function host_glsc3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: a, b, c
    real(kind=rp) :: host_glsc3, tmp
    integer :: ierr

    tmp = sum(a * b * c)

    call MPI_Allreduce(tmp, host_glsc3, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

  end function host_glsc3

  !> Weighted inner product \f$ a^T b c d \f$
  function host_glsc4(a, b, c, d, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: a, b, c, d
    real(kind=rp) :: host_glsc4, tmp
    integer :: ierr

    tmp = sum(a * b * c * d)

    call MPI_Allreduce(tmp, host_glsc4, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

  end function host_glsc4

  !> Heap Sort for double precision arrays
  !! @details Following p 231 Num. Rec., 1st Ed.
  !! @param[inout]  a   vector to be sorted
  !! @param[out]   ind  permutation array
  !! @param[in]   n   array size
  subroutine host_sortrp(a, ind, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: a(n)
    integer, intent(out) :: ind(n)
    real(kind=rp) :: aa
    integer :: j, ir, i, ii, l

    do j = 1, n
       ind(j) = j
    end do

    if (n .le. 1) return

    l = n / 2 + 1
    ir = n
    do while (.true.)
       if (l .gt. 1) then
          l = l - 1
          aa = a(l)
          ii = ind(l)
       else
          aa = a(ir)
          ii = ind(ir)
          a(ir) = a(1)
          ind(ir) = ind(1)
          ir = ir - 1
          if (ir .eq. 1) then
             a(1) = aa
             ind(1) = ii
             return
          end if
       end if
       i = l
       j = l+l
       do while (j .le. ir)
          if (j .lt. ir .and. a(j) .lt. a(j + 1)) j = j + 1
          if (aa .lt. a(j)) then
             a(i) = a(j)
             ind(i) = ind(j)
             i = j
             j = j + j
          else
             j = ir + 1
          end if
       end do
       a(i) = aa
       ind(i) = ii
    end do
  end subroutine host_sortrp

  !> Heap Sort for single integer arrays
  !! @details Following p 231 Num. Rec., 1st Ed.
  !! @param[inout]  a   vector to be sorted
  !! @param[out]   ind  permutation array
  !! @param[in]   n   array size
  subroutine host_sorti4(a, ind, n)
    integer, intent(in) :: n
    integer(i4), intent(inout) :: a(n)
    integer, intent(out) :: ind(n)
    integer(i4) :: aa
    integer :: j, ir, i, ii, l

    do j = 1, n
       ind(j) = j
    end do

    if (n .le. 1) return

    l = n / 2 + 1
    ir = n
    do while (.true.)
       if (l .gt. 1) then
          l = l - 1
          aa = a (l)
          ii = ind(l)
       else
          aa = a(ir)
          ii = ind(ir)
          a(ir) = a(1)
          ind(ir) = ind(1)
          ir = ir - 1
          if (ir .eq. 1) then
             a(1) = aa
             ind(1) = ii
             return
          end if
       end if
       i = l
       j = l + l
       do while (j .le. ir)
          if (j .lt. ir .and. a(j) .lt. a(j + 1)) j = j + 1
          if (aa .lt. a(j)) then
             a(i) = a(j)
             ind(i) = ind(j)
             i = j
             j = j + j
          else
             j = ir + 1
          end if
       end do
       a(i) = aa
       ind(i) = ii
    end do
  end subroutine host_sorti4

  !> host_host_sort double precision array acording to ind vector
  !! @param[inout]  b   vector to be reordered
  !! @param[in]   ind  permutation array
  !! @param[in]   n   array size
  subroutine host_swapdp(b, ind, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: b(n)
    integer, intent(in) :: ind(n)
    real(kind=rp) :: temp(n)
    integer :: i, jj

    temp = b
    do i = 1, n
       jj = ind(i)
       b(i) = temp(jj)
    end do
  end subroutine host_swapdp

  !> host_host_sort single integer array acording to ind vector
  !! @param[inout]  b   vector to be reordered
  !! @param[in]   ind  permutation array
  !! @param[in]   n   array size
  subroutine host_swapi4(b, ind, n)
    integer, intent(in) :: n
    integer(i4), intent(inout) :: b(n)
    integer, intent(in) :: ind(n)
    integer(i4) :: temp(n)
    integer :: i, jj

    temp = b
    do i = 1, n
       jj = ind(i)
       b(i) = temp(jj)
    end do
  end subroutine host_swapi4

  !> reorder double precision array - inverse of host_swap
  !! @param[inout]  b   vector to be reordered
  !! @param[in]   ind  permutation array
  !! @param[in]   n   array size
  subroutine host_reorddp(b, ind, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: b(n)
    integer, intent(in) :: ind(n)
    real(kind=rp) :: temp(n)
    integer :: i, jj

    temp = b
    do i = 1, n
       jj = ind(i)
       b(jj) = temp(i)
    end do
  end subroutine host_reorddp

  !> reorder single integer array - inverse of host_swap
  !! @param[inout]  b   vector to be reordered
  !! @param[in]   ind  permutation array
  !! @param[in]   n   array size
  subroutine host_reordi4(b, ind, n)
    integer, intent(in) :: n
    integer(i4), intent(inout) :: b(n)
    integer, intent(in) :: ind(n)
    integer(i4) :: temp(n)
    integer :: i, jj

    temp = b
    do i = 1, n
       jj = ind(i)
       b(jj) = temp(i)
    end do
  end subroutine host_reordi4

  !> Flip double precision vector b and ind
  !! @param[inout]  b   vector to be reordered
  !! @param[inout]  ind  permutation array
  !! @param[in]   n   array size
  subroutine host_flipvdp(b, ind, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: b(n)
    integer, intent(inout) :: ind(n)
    real(kind=rp) :: temp(n)
    integer :: tempind(n)
    integer :: i, jj

    do i = 1, n
       jj = n + 1 - i
       temp(jj) = b(i)
       tempind(jj) = ind(i)
    end do

    b = temp
    ind = tempind

  end subroutine host_flipvdp

  !> Flip single integer vector b and ind
  !! @param[inout]  b   vector to be reordered
  !! @param[inout]  ind  permutation array
  !! @param[in]   n   array size
  subroutine host_flipvi4(b, ind, n)
    integer, intent(in) :: n
    integer(i4), intent(inout) :: b(n)
    integer, intent(inout) :: ind(n)
    integer(i4) :: temp(n)
    integer :: tempind(n)
    integer :: i, jj

    do i = 1, n
       jj = n + 1 - i
       temp(jj) = b(i)
       tempind(jj) = ind(i)
    end do

    b = temp
    ind = tempind
  end subroutine host_flipvi4

end module host_math
