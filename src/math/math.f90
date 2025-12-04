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
module math
  use num_types, only : rp, dp, sp, qp, i4, xp
  use comm, only: NEKO_COMM, MPI_REAL_PRECISION, MPI_EXTRA_PRECISION
  use mpi_f08, only: MPI_MIN, MPI_MAX, MPI_SUM, MPI_IN_PLACE, MPI_INTEGER, &
       MPI_Allreduce
  implicit none
  private

  !> Machine epsilon \f$ \epsilon \f$
  real(kind=rp), public, parameter :: NEKO_EPS = epsilon(1.0_rp)

  !> \f$ ln(2) \f$
  real(kind=rp), public, parameter :: NEKO_M_LN2 = log(2.0_rp)

  !> \f$ \pi \f$
  real(kind=rp), public, parameter :: pi = 4._rp*atan(1._rp)

  interface abscmp
     module procedure sabscmp, dabscmp, qabscmp
  end interface abscmp

  interface sort
     module procedure sortrp, sorti4
  end interface sort

  interface swap
     module procedure swapdp, swapi4
  end interface swap

  interface reord
     module procedure reorddp, reordi4
  end interface reord

  interface flipv
     module procedure flipvdp, flipvi4
  end interface flipv

  interface relcmp
     module procedure srelcmp, drelcmp, qrelcmp
  end interface relcmp

  public :: abscmp, rzero, izero, row_zero, rone, copy, cmult, cadd, cfill, &
       glsum, glmax, glmin, chsign, vlmax, vlmin, invcol1, invcol3, invers2, &
       vcross, vdot2, vdot3, vlsc3, vlsc2, add2, add3, add4, sub2, sub3, &
       add2s1, add2s2, addsqr2s2, cmult2, invcol2, col2, col3, subcol3, &
       add3s2, add4s3, add5s4, subcol4, addcol3, addcol4, addcol3s2, ascol5, &
       p_update, x_update, glsc2, glsc3, glsc4, sort, masked_copy_0, &
       cfill_mask, relcmp, glimax, glimin, swap, reord, flipv, cadd2, &
       masked_gather_copy_0, absval, matinv3, matinv39, &
       pwmax2, pwmax3, cpwmax2, cpwmax3, pwmin2, pwmin3, cpwmin2, cpwmin3, &
       masked_scatter_copy_0, cdiv, cdiv2, glsubnorm, &
       masked_copy, masked_gather_copy, masked_scatter_copy, sabscmp, dabscmp

contains

  !> Return single precision absolute comparison \f$ | x - y | < \epsilon \f$
  pure function sabscmp(x, y, tol)
    real(kind=sp), intent(in) :: x
    real(kind=sp), intent(in) :: y
    real(kind=sp), intent(in), optional :: tol
    logical :: sabscmp

    if (present(tol)) then
       sabscmp = abs(x - y) .lt. tol
    else
       sabscmp = abs(x - y) .lt. NEKO_EPS
    end if

  end function sabscmp

  !> Return double precision absolute comparison \f$ | x - y | < \epsilon \f$
  pure function dabscmp(x, y, tol)
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(in) :: y
    real(kind=dp), intent(in), optional :: tol
    logical :: dabscmp

    if (present(tol)) then
       dabscmp = abs(x - y) .lt. tol
    else
       dabscmp = abs(x - y) .lt. NEKO_EPS
    end if

  end function dabscmp

  !> Return double precision absolute comparison \f$ | x - y | < \epsilon \f$
  pure function qabscmp(x, y, tol)
    real(kind=qp), intent(in) :: x
    real(kind=qp), intent(in) :: y
    real(kind=qp), intent(in), optional :: tol
    logical :: qabscmp

    if (present(tol)) then
       qabscmp = abs(x - y) .lt. tol
    else
       qabscmp = abs(x - y) .lt. NEKO_EPS
    end if

  end function qabscmp

  !> Return single precision relative comparison \f$ | x - y |<= \epsilon*|y| \f$
  pure function srelcmp(x, y, eps)
    real(kind=sp), intent(in) :: x
    real(kind=sp), intent(in) :: y
    real(kind=sp), intent(in), optional :: eps
    logical :: srelcmp
    if (present(eps)) then
       srelcmp = abs(x - y) .le. eps*abs(y)
    else
       srelcmp = abs(x - y) .le. NEKO_EPS*abs(y)
    end if

  end function srelcmp

  !> Return double precision relative comparison \f$ | x - y |/|y| < \epsilon \f$
  pure function drelcmp(x, y, eps)
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(in) :: y
    real(kind=dp), intent(in), optional :: eps
    logical :: drelcmp
    if (present(eps)) then
       drelcmp = abs(x - y) .le. eps*abs(y)
    else
       drelcmp = abs(x - y) .le. NEKO_EPS*abs(y)
    end if

  end function drelcmp


  !> Return quad precision relative comparison \f$ | x - y |/|y| < \epsilon \f$
  pure function qrelcmp(x, y, eps)
    real(kind=qp), intent(in) :: x
    real(kind=qp), intent(in) :: y
    real(kind=qp), intent(in), optional :: eps
    logical :: qrelcmp
    if (present(eps)) then
       qrelcmp = abs(x - y)/abs(y) .lt. eps
    else
       qrelcmp = abs(x - y)/abs(y) .lt. NEKO_EPS
    end if

  end function qrelcmp

  !> Zero a real vector
  subroutine rzero(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    integer :: i

    do i = 1, n
       a(i) = 0.0_rp
    end do
  end subroutine rzero

  !> Zero an integer vector
  subroutine izero(a, n)
    integer, intent(in) :: n
    integer, dimension(n), intent(inout) :: a
    integer :: i

    do i = 1, n
       a(i) = 0
    end do
  end subroutine izero

  !> Sets row e to 0 in matrix a
  subroutine row_zero(a, m, n, e)
    integer, intent(in) :: m, n, e
    real(kind=rp), intent(inout) :: a(m,n)
    integer :: j

    do j = 1,n
       a(e,j) = 0.0_rp
    end do
  end subroutine row_zero

  !> Set all elements to one
  subroutine rone(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    integer :: i

    do i = 1, n
       a(i) = 1.0_rp
    end do
  end subroutine rone

  !> Copy a vector \f$ a = b \f$
  subroutine copy(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(inout) :: a
    integer :: i

    do i = 1, n
       a(i) = b(i)
    end do

  end subroutine copy

  !> Copy a masked vector \f$ a(mask) = b(mask) \f$.
  !! @param a Destination array of size `n`.
  !! @param b Source array of size `n`.
  !! @param mask Mask array of length n_mask + 1, where `mask(0) = n_mask`
  !! the length of the mask array.
  !! @param n Size of the arrays `a` and `b`.
  !! @param n_mask Size of the mask array `mask`.
  subroutine masked_copy_0(a, b, mask, n, n_mask)
    integer, intent(in) :: n, n_mask
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(inout) :: a
    integer, dimension(0:n_mask) :: mask
    integer :: i, j

    do i = 1, n_mask
       j = mask(i)
       a(j) = b(j)
    end do

  end subroutine masked_copy_0

  !> Copy a masked vector \f$ a(mask) = b(mask) \f$.
  !! @param a Destination array of size `n`.
  !! @param b Source array of size `n`.
  !! @param mask Mask array of length n_mask + 1, where `mask(0) = n_mask`
  !! the length of the mask array.
  !! @param n Size of the arrays `a` and `b`.
  !! @param n_mask Size of the mask array `mask`.
  subroutine masked_copy(a, b, mask, n, n_mask)
    integer, intent(in) :: n, n_mask
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(inout) :: a
    integer, dimension(n_mask) :: mask
    integer :: i, j

    do i = 1, n_mask
       j = mask(i)
       a(j) = b(j)
    end do

  end subroutine masked_copy

  !> Gather a masked vector to reduced contigous vector
  !! \f$ a = b(mask) \f$.
  !! @param a Destination array of size `n_mask`.
  !! @param b Source array of size `n`.
  !! @param mask Mask array of length n_mask + 1, where `mask(0) = n_mask`
  !! the length of the mask array.
  !! @param n Size of the array `b`.
  !! @param n_mask Size of the mask array `mask` and `a`.
  !! @note Assumes `n .ge. n_mask`.
  subroutine masked_gather_copy_0(a, b, mask, n, n_mask)
    integer, intent(in) :: n, n_mask
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n_mask), intent(inout) :: a
    integer, dimension(0:n_mask) :: mask
    integer :: i, j

    do i = 1, n_mask
       j = mask(i)
       a(i) = b(j)
    end do

  end subroutine masked_gather_copy_0

  !> Gather a masked vector to reduced contigous vector
  !! \f$ a = b(mask) \f$.
  !! @param a Destination array of size `n_mask`.
  !! @param b Source array of size `n`.
  !! @param mask Mask array of length n_mask + 1, where `mask(0) = n_mask`
  !! the length of the mask array.
  !! @param n Size of the array `b`.
  !! @param n_mask Size of the mask array `mask` and `a`.
  !! @note Assumes `n .ge. n_mask`.
  subroutine masked_gather_copy(a, b, mask, n, n_mask)
    integer, intent(in) :: n, n_mask
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n_mask), intent(inout) :: a
    integer, dimension(n_mask) :: mask
    integer :: i, j

    do i = 1, n_mask
       j = mask(i)
       a(i) = b(j)
    end do

  end subroutine masked_gather_copy

  !> Scatter a contigous vector to masked positions in a target array
  !! \f$ a(mask) = b \f$.
  !! @param a Destination array of size `n`.
  !! @param b Source array of size `n_mask`.
  !! @param mask Mask array of length n_mask + 1, where `mask(0) = n_mask + 1`
  !! the length of the mask array.
  !! @param n Size of the array `mask`and array `b`.
  !! @param m Size of the mask array `a`.
  !! @note Assumes `n .ge. n_mask`.
  subroutine masked_scatter_copy_0(a, b, mask, n, n_mask)
    integer, intent(in) :: n, n_mask
    real(kind=rp), dimension(n_mask), intent(in) :: b
    real(kind=rp), dimension(n), intent(inout) :: a
    integer, dimension(0:n_mask) :: mask
    integer :: i, j

    do i = 1, n_mask
       j = mask(i)
       a(j) = b(i)
    end do

  end subroutine masked_scatter_copy_0

  !> Scatter a contigous vector to masked positions in a target array
  !! \f$ a(mask) = b \f$.
  !! @param a Destination array of size `n`.
  !! @param b Source array of size `n_mask`.
  !! @param mask Mask array of length n_mask + 1, where `mask(0) = n_mask + 1`
  !! the length of the mask array.
  !! @param n Size of the array `mask`and array `b`.
  !! @param m Size of the mask array `a`.
  !! @note Assumes `n .ge. n_mask`.
  subroutine masked_scatter_copy(a, b, mask, n, n_mask)
    integer, intent(in) :: n, n_mask
    real(kind=rp), dimension(n_mask), intent(in) :: b
    real(kind=rp), dimension(n), intent(inout) :: a
    integer, dimension(n_mask) :: mask
    integer :: i, j

    do i = 1, n_mask
       j = mask(i)
       a(j) = b(i)
    end do

  end subroutine masked_scatter_copy

  !> @brief Fill a constant to a masked vector.
  !! \f$ a_i = c, for i in mask \f$
  subroutine cfill_mask(a, c, n, mask, n_mask)
    integer, intent(in) :: n, n_mask
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer, dimension(n_mask), intent(in) :: mask
    integer :: i

    do i = 1, n_mask
       a(mask(i)) = c
    end do

  end subroutine cfill_mask

  !> Multiplication by constant c \f$ a = c \cdot a \f$
  subroutine cmult(a, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer :: i

    do i = 1, n
       a(i) = c * a(i)
    end do
  end subroutine cmult

  !> Multiplication by constant c \f$ a = c \cdot b \f$
  subroutine cmult2(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: c
    integer :: i

    do i = 1, n
       a(i) = c * b(i)
    end do

  end subroutine cmult2

  !> Division of constant c by elements of a \f$ a = c / a \f$
  subroutine cdiv(a, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer :: i

    do i = 1, n
       a(i) = c / a(i)
    end do
  end subroutine cdiv

  !> Division of constant c by elements of a \f$ a = c / b \f$
  subroutine cdiv2(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: c
    integer :: i

    do i = 1, n
       a(i) = c / b(i)
    end do
  end subroutine cdiv2

  !> Add a scalar to vector \f$ a_i = a_i + s \f$
  subroutine cadd(a, s, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), intent(in) :: s
    integer :: i

    do i = 1, n
       a(i) = a(i) + s
    end do
  end subroutine cadd

  !> Add a scalar to vector \f$ a_i = b_i + s \f$
  subroutine cadd2(a, b, s, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: s
    integer :: i

    do i = 1, n
       a(i) = b(i) + s
    end do
  end subroutine cadd2

  !> Set all elements to a constant c \f$ a = c \f$
  subroutine cfill(a, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer :: i

    do i = 1, n
       a(i) = c
    end do
  end subroutine cfill

  !> Sum a vector of length n
  function glsum(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n) :: a
    real(kind=rp) :: glsum
    real(kind=xp) :: tmp
    integer :: i, ierr
    tmp = 0.0_rp
    do i = 1, n
       tmp = tmp + a(i)
    end do
    call MPI_Allreduce(MPI_IN_PLACE, tmp, 1, &
         MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    glsum = tmp

  end function glsum

  !>Max of a vector of length n
  function glmax(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n) :: a
    real(kind=rp) :: tmp, glmax
    integer :: i, ierr

    tmp = -huge(0.0_rp)
    do i = 1, n
       tmp = max(tmp,a(i))
    end do
    call MPI_Allreduce(tmp, glmax, 1, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
  end function glmax

  !>Max of an integer vector of length n
  function glimax(a, n)
    integer, intent(in) :: n
    integer, dimension(n) :: a
    integer :: tmp, glimax
    integer :: i, ierr

    tmp = -huge(0)
    do i = 1, n
       tmp = max(tmp,a(i))
    end do
    call MPI_Allreduce(tmp, glimax, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)
  end function glimax

  !>Min of a vector of length n
  function glmin(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n) :: a
    real(kind=rp) :: tmp, glmin
    integer :: i, ierr

    tmp = huge(0.0_rp)
    do i = 1, n
       tmp = min(tmp,a(i))
    end do
    call MPI_Allreduce(tmp, glmin, 1, &
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
  end function glmin

  !>Min of an integer vector of length n
  function glimin(a, n)
    integer, intent(in) :: n
    integer, dimension(n) :: a
    integer :: tmp, glimin
    integer :: i, ierr

    tmp = huge(0)
    do i = 1, n
       tmp = min(tmp,a(i))
    end do
    call MPI_Allreduce(tmp, glimin, 1, &
         MPI_INTEGER, MPI_MIN, NEKO_COMM, ierr)
  end function glimin




  !> Change sign of vector \f$ a = -a \f$
  subroutine chsign(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    integer :: i

    do i = 1, n
       a(i) = -a(i)
    end do

  end subroutine chsign

  !> maximum value of a vector of length @a n
  function vlmax(vec,n) result(tmax)
    integer :: n, i
    real(kind=rp), intent(in) :: vec(n)
    real(kind=rp) :: tmax
    tmax = real(-99d20, rp)
    do i = 1, n
       tmax = max(tmax, vec(i))
    end do
  end function vlmax

  !> minimun value of a vector of length @a n
  function vlmin(vec,n) result(tmin)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: vec(n)
    real(kind=rp) :: tmin
    integer :: i
    tmin = real(99.0e20, rp)
    do i = 1, n
       tmin = min(tmin, vec(i))
    end do
  end function vlmin

  !> Invert a vector \f$ a = 1 / a \f$
  subroutine invcol1(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    integer :: i

    do i = 1, n
       a(i) = 1.0_xp / real(a(i),xp)
    end do

  end subroutine invcol1

  !> Invert a vector \f$ a = b / c \f$
  subroutine invcol3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b,c
    integer :: i

    do i = 1, n
       a(i) = real(b(i),xp) / c(i)
    end do

  end subroutine invcol3

  !> Compute inverted vector \f$ a = 1 / b \f$
  subroutine invers2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    integer :: i

    do i = 1, n
       a(i) = 1.0_xp / real(b(i),xp)
    end do

  end subroutine invers2

  !> Compute a cross product \f$ u = v \times w \f$
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  subroutine vcross(u1, u2, u3, v1, v2, v3, w1, w2, w3, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: v1, v2, v3
    real(kind=rp), dimension(n), intent(in) :: w1, w2, w3
    real(kind=rp), dimension(n), intent(out) :: u1, u2, u3
    integer :: i

    do i = 1, n
       u1(i) = v2(i)*w3(i) - v3(i)*w2(i)
       u2(i) = v3(i)*w1(i) - v1(i)*w3(i)
       u3(i) = v1(i)*w2(i) - v2(i)*w1(i)
    end do

  end subroutine vcross

  !> Compute a dot product \f$ dot = u \cdot v \f$ (2-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  subroutine vdot2(dot, u1, u2, v1, v2, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: u1, u2
    real(kind=rp), dimension(n), intent(in) :: v1, v2
    real(kind=rp), dimension(n), intent(out) :: dot
    integer :: i
    do i = 1, n
       dot(i) = u1(i)*v1(i) + u2(i)*v2(i)
    end do

  end subroutine vdot2

  !> Compute a dot product \f$ dot = u \cdot v \f$ (3-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  subroutine vdot3(dot, u1, u2, u3, v1, v2, v3, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: u1, u2, u3
    real(kind=rp), dimension(n), intent(in) :: v1, v2, v3
    real(kind=rp), dimension(n), intent(out) :: dot
    integer :: i

    do i = 1, n
       dot(i) = u1(i)*v1(i) + u2(i)*v2(i) + u3(i)*v3(i)
    end do

  end subroutine vdot3

  !> Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
  function vlsc3(u, v, w, n) result(s)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: u, v, w
    real(kind=rp) :: s
    integer :: i

    s = 0.0_rp
    do i = 1, n
       s = s + u(i)*v(i)*w(i)
    end do

  end function vlsc3

  !> Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
  function vlsc2(u, v, n) result(s)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: u, v
    real(kind=rp) :: s
    integer :: i

    s = 0.0_rp
    do i = 1, n
       s = s + u(i)*v(i)
    end do

  end function vlsc2

  !> Vector addition \f$ a = a + b \f$
  subroutine add2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    integer :: i

    do i = 1, n
       a(i) = a(i) + b(i)
    end do

  end subroutine add2

  !> Vector addition \f$ a = b + c \f$
  subroutine add3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    integer :: i

    do i = 1, n
       a(i) = b(i) + c(i)
    end do

  end subroutine add3

  !> Vector addition \f$ a = b + c + d\f$
  subroutine add4(a, b, c, d, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: a
    real(kind=rp), dimension(n), intent(in) :: d
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp), dimension(n), intent(in) :: b
    integer :: i

    do i = 1, n
       a(i) = b(i) + c(i) + d(i)
    end do

  end subroutine add4

  !> Vector substraction \f$ a = a - b \f$
  subroutine sub2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    integer :: i

    do i = 1, n
       a(i) = a(i) - b(i)
    end do

  end subroutine sub2

  !> Vector subtraction \f$ a = b - c \f$
  subroutine sub3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    integer :: i

    do i = 1, n
       a(i) = b(i) - c(i)
    end do

  end subroutine sub3


  !> Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
  !! (multiplication on first argument)
  subroutine add2s1(a, b, c1, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: c1
    integer :: i

    do i = 1, n
       a(i) = c1 * a(i) + b(i)
    end do

  end subroutine add2s1

  !> Vector addition with scalar multiplication  \f$ a = a + c_1 b \f$
  !! (multiplication on second argument)
  subroutine add2s2(a, b, c1, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: c1
    integer :: i

    do i = 1, n
       a(i) = a(i) + c1 * b(i)
    end do

  end subroutine add2s2

  !> Returns \f$ a = a + c1 * (b * b )\f$
  subroutine addsqr2s2(a, b, c1, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: c1
    integer :: i

    do i = 1,n
       a(i) = a(i) + c1 * ( b(i) * b(i) )
    end do

  end subroutine addsqr2s2

  !> Vector division \f$ a = a / b \f$
  subroutine invcol2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    integer :: i

    do i = 1, n
       a(i) = real(a(i),xp) /b(i)
    end do

  end subroutine invcol2


  !> Vector multiplication \f$ a = a \cdot b \f$
  subroutine col2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    integer :: i

    do i = 1, n
       a(i) = a(i) * b(i)
    end do

  end subroutine col2

  !> Vector multiplication with 3 vectors \f$ a =  b \cdot c \f$
  subroutine col3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    integer :: i

    do i = 1, n
       a(i) = b(i) * c(i)
    end do

  end subroutine col3

  !> Returns \f$ a = a - b*c \f$
  subroutine subcol3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    integer :: i

    do i = 1,n
       a(i) = a(i) - b(i) * c(i)
    end do

  end subroutine subcol3

  !> Returns \f$ a = c1 * b + c2 * c \f$
  subroutine add3s2(a, b, c, c1, c2 ,n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp), intent(in) :: c1, c2
    integer :: i

    do i = 1,n
       a(i) = c1 * b(i) + c2 * c(i)
    end do

  end subroutine add3s2

  !> Returns \f$ a = c1 * b + c2 * c + c3 * d \f$
  subroutine add4s3(a, b, c, d, c1, c2, c3, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp), dimension(n), intent(in) :: d
    real(kind=rp), intent(in) :: c1, c2, c3
    integer :: i

    do i = 1,n
       a(i) = c1 * b(i) + c2 * c(i) + c3 * d(i)
    end do

  end subroutine add4s3

  !> Returns \f$ a = a + c1 * b + c2 * c + c3 * d + c4 * e\f$
  subroutine add5s4(a, b, c, d, e, c1, c2, c3, c4, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp), dimension(n), intent(in) :: d
    real(kind=rp), dimension(n), intent(in) :: e
    real(kind=rp), intent(in) :: c1, c2, c3, c4
    integer :: i

    do i = 1,n
       a(i) = a(i) + c1 * b(i) + c2 * c(i) + c3 * d(i) + c4 * e(i)
    end do

  end subroutine add5s4

  !> Returns \f$ a = a - b*c*d \f$
  subroutine subcol4(a, b, c, d, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp), dimension(n), intent(in) :: d
    integer :: i

    do i = 1,n
       a(i) = a(i) - b(i) * c(i) * d(i)
    end do

  end subroutine subcol4

  !> Returns \f$ a = a + b*c \f$
  subroutine addcol3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    integer :: i

    do i = 1,n
       a(i) = a(i) + b(i) * c(i)
    end do

  end subroutine addcol3

  !> Returns \f$ a = a + b*c*d \f$
  subroutine addcol4(a, b, c, d, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp), dimension(n), intent(in) :: d
    integer :: i

    do i = 1,n
       a(i) = a(i) + b(i) * c(i) * d(i)
    end do

  end subroutine addcol4

  !> Returns \f$ a = a + s(b*c) \f$
  subroutine addcol3s2(a, b, c, s, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp), intent(in) :: s
    integer :: i

    do i = 1,n
       a(i) = a(i) + s * b(i) * c(i)
    end do

  end subroutine addcol3s2

  !> Returns \f$ a = b \dot c - d \cdot e \f$
  subroutine ascol5(a, b, c, d, e, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp), dimension(n), intent(in) :: d
    real(kind=rp), dimension(n), intent(in) :: e
    integer :: i

    do i = 1,n
       a(i) = b(i)*c(i)-d(i)*e(i)
    end do

  end subroutine ascol5

  !> Returns \f$ a = b \dot c1 ( a - c2 \cdot c )\f$
  subroutine p_update(a, b, c, c1, c2, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp), intent(in) :: c1, c2
    integer :: i

    do i = 1,n
       a(i) = b(i) + c1*(a(i)-c2*c(i))
    end do

  end subroutine p_update

  !> Returns \f$ a = b \dot c1 ( a - c2 \cdot c )\f$
  subroutine x_update(a, b, c, c1, c2, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp), intent(in) :: c1, c2
    integer :: i

    do i = 1,n
       a(i) = a(i) + c1*b(i)+c2*c(i)
    end do

  end subroutine x_update

  !> Weighted inner product \f$ a^T b \f$
  function glsc2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp) :: glsc2
    real(kind=xp) :: tmp
    integer :: i, ierr

    tmp = 0.0_xp
    do i = 1, n
       tmp = tmp + a(i) * b(i)
    end do

    call MPI_Allreduce(MPI_IN_PLACE, tmp, 1, &
         MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    glsc2 = tmp
  end function glsc2

  !> Weighted inner product \f$ a^T b c \f$
  function glsc3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp) :: glsc3
    real(kind=xp) :: tmp
    integer :: i, ierr

    tmp = 0.0_xp
    do i = 1, n
       tmp = tmp + a(i) * b(i) * c(i)
    end do

    call MPI_Allreduce(MPI_IN_PLACE, tmp, 1, &
         MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    glsc3 = tmp

  end function glsc3
  function glsc4(a, b, c, d, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp), dimension(n), intent(in) :: d
    real(kind=rp) :: glsc4
    real(kind=xp) :: tmp
    integer :: i, ierr

    tmp = 0.0_xp
    do i = 1, n
       tmp = tmp + a(i) * b(i) * c(i) * d(i)
    end do

    call MPI_Allreduce(MPI_IN_PLACE, tmp, 1, &
         MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    glsc4 = tmp

  end function glsc4

  !> Returns the norm of the difference of two vectors
  !! \f$ \sqrt{(a-b)^T (a-b)} \f$
  function glsubnorm(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp) :: glsubnorm
    real(kind=xp) :: tmp
    integer :: i, ierr

    tmp = 0.0_xp
    do i = 1, n
       tmp = tmp + (a(i) - b(i))**2
    end do

    call MPI_Allreduce(MPI_IN_PLACE, tmp, 1, &
         MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    glsubnorm = sqrt(tmp)

  end function glsubnorm

  !> Heap Sort for double precision arrays
  !! @details Following p 231 Num. Rec., 1st Ed.
  !! @param[inout]   a     vector to be sorted
  !! @param[out]     ind   permutation array
  !! @param[in]      n     array size
  subroutine sortrp(a, ind, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: a(n)
    integer, intent(out) :: ind(n)
    real(kind=rp) :: aa
    integer :: j, ir, i, ii, l

    do j = 1, n
       ind(j) = j
    end do

    if (n .le. 1) return


    l = n/2+1
    ir = n
    do while (.true.)
       if (l .gt. 1) then
          l = l-1
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
          if (j .lt. ir) then
             if ( a(j) .lt. a(j+1) ) j = j + 1
          end if
          if (aa .lt. a(j)) then
             a(i) = a(j)
             ind(i) = ind(j)
             i = j
             j = j+j
          else
             j = ir+1
          end if
       end do
       a(i) = aa
       ind(i) = ii
    end do
  end subroutine sortrp

  !> Heap Sort for single integer arrays
  !! @details Following p 231 Num. Rec., 1st Ed.
  !! @param[inout]   a     vector to be sorted
  !! @param[out]     ind   permutation array
  !! @param[in]      n     array size
  subroutine sorti4(a, ind, n)
    integer, intent(in) :: n
    integer(i4), intent(inout) :: a(n)
    integer, intent(out) :: ind(n)
    integer(i4) :: aa
    integer :: j, ir, i, ii, l

    do j = 1, n
       ind(j) = j
    end do

    if (n .le. 1) return

    l = n/2+1
    ir = n
    do while (.true.)
       if (l .gt. 1) then
          l = l - 1
          aa = a (l)
          ii = ind(l)
       else
          aa = a(ir)
          ii = ind(ir)
          a(ir) = a( 1)
          ind(ir) = ind( 1)
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
          if (j .lt. ir) then
             if ( a(j) .lt. a(j + 1) ) j = j + 1
          end if
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
  end subroutine sorti4

  !> sort double precision array acording to ind vector
  !! @param[inout]   b     vector to be reordered
  !! @param[in]      ind   permutation array
  !! @param[in]      n     array size
  subroutine swapdp(b, ind, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: b(n)
    integer, intent(in) :: ind(n)
    real(kind=rp) :: temp(n)
    integer :: i, jj

    do i = 1, n
       temp(i) = b(i)
    end do
    do i = 1, n
       jj = ind(i)
       b(i) = temp(jj)
    end do
  end subroutine swapdp

  !> sort single integer array acording to ind vector
  !! @param[inout]   b     vector to be reordered
  !! @param[in]      ind   permutation array
  !! @param[in]      n     array size
  subroutine swapi4(b, ind, n)
    integer, intent(in) :: n
    integer(i4), intent(inout) :: b(n)
    integer, intent(in) :: ind(n)
    integer(i4) :: temp(n)
    integer :: i, jj

    do i = 1, n
       temp(i) = b(i)
    end do
    do i = 1, n
       jj = ind(i)
       b(i) = temp(jj)
    end do
  end subroutine swapi4

  !> reorder double precision array - inverse of swap
  !! @param[inout]   b     vector to be reordered
  !! @param[in]      ind   permutation array
  !! @param[in]      n     array size
  subroutine reorddp(b, ind, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: b(n)
    integer, intent(in) :: ind(n)
    real(kind=rp) :: temp(n)
    integer :: i, jj

    do i = 1, n
       temp(i) = b(i)
    end do
    do i = 1, n
       jj = ind(i)
       b(jj) = temp(i)
    end do
  end subroutine reorddp

  !> reorder single integer array - inverse of swap
  !! @param[inout]   b     vector to be reordered
  !! @param[in]      ind   permutation array
  !! @param[in]      n     array size
  subroutine reordi4(b, ind, n)
    integer, intent(in) :: n
    integer(i4), intent(inout) :: b(n)
    integer, intent(in) :: ind(n)
    integer(i4) :: temp(n)
    integer :: i, jj

    do i = 1, n
       temp(i) = b(i)
    end do
    do i = 1, n
       jj = ind(i)
       b(jj) = temp(i)
    end do
  end subroutine reordi4

  !> Flip double precision vector b and ind
  !! @param[inout]   b     vector to be reordered
  !! @param[inout]   ind   permutation array
  !! @param[in]      n     array size
  subroutine flipvdp(b, ind, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: b(n)
    integer, intent(inout) :: ind(n)
    real(kind=rp) :: temp(n)
    integer :: tempind(n)
    integer :: i, jj

    do i = 1, n
       jj = n+1-i
       temp(jj) = b(i)
       tempind(jj) = ind(i)
    end do
    do i = 1,n
       b(i) = temp(i)
       ind(i) = tempind(i)
    end do
  end subroutine flipvdp

  !> Flip single integer vector b and ind
  !! @param[inout]   b     vector to be reordered
  !! @param[inout]   ind   permutation array
  !! @param[in]      n     array size
  subroutine flipvi4(b, ind, n)
    integer, intent(in) :: n
    integer(i4), intent(inout) :: b(n)
    integer, intent(inout) :: ind(n)
    integer(i4) :: temp(n)
    integer :: tempind(n)
    integer :: i, jj

    do i = 1, n
       jj = n+1-i
       temp(jj) = b(i)
       tempind(jj) = ind(i)
    end do
    do i = 1,n
       b(i) = temp(i)
       ind(i) = tempind(i)
    end do
  end subroutine flipvi4

  !> Take the absolute value of an array
  !! @param[inout]   a     vector to be manipulated
  !! @param[in]      n     array size
  subroutine absval(a, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    integer :: i
    do i = 1, n
       a(i) = abs(a(i))
    end do
  end subroutine absval

  ! ========================================================================== !
  ! Point-wise operations

  !> Point-wise maximum of two vectors \f$ a = \max(a, b) \f$
  subroutine pwmax2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    integer :: i

    do i = 1, n
       a(i) = max(a(i), b(i))
    end do
  end subroutine pwmax2

  !> Point-wise maximum of two vectors \f$ a = \max(b, c) \f$
  subroutine pwmax3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c
    integer :: i

    do i = 1, n
       a(i) = max(b(i), c(i))
    end do
  end subroutine pwmax3

  !> Point-wise maximum of scalar and vector \f$ a = \max(a, b) \f$
  subroutine cpwmax2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), intent(in) :: b
    integer :: i

    do i = 1, n
       a(i) = max(a(i), b)
    end do
  end subroutine cpwmax2

  !> Point-wise maximum of scalar and vector \f$ a = \max(b, c) \f$
  subroutine cpwmax3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: c
    integer :: i

    do i = 1, n
       a(i) = max(b(i), c)
    end do
  end subroutine cpwmax3

  !> Point-wise minimum of two vectors \f$ a = \min(a, b) \f$
  subroutine pwmin2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    integer :: i

    do i = 1, n
       a(i) = min(a(i), b(i))
    end do
  end subroutine pwmin2

  !> Point-wise minimum of two vectors \f$ a = \min(b, c) \f$
  subroutine pwmin3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b, c
    integer :: i

    do i = 1, n
       a(i) = min(b(i), c(i))
    end do
  end subroutine pwmin3

  !> Point-wise minimum of scalar and vector \f$ a = \min(a, b) \f$
  subroutine cpwmin2(a, b, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), intent(in) :: b
    integer :: i

    do i = 1, n
       a(i) = min(a(i), b)
    end do
  end subroutine cpwmin2

  !> Point-wise minimum of scalar and vector \f$ a = \min(b, c) \f$
  subroutine cpwmin3(a, b, c, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), intent(in) :: c
    integer :: i

    do i = 1, n
       a(i) = min(b(i), c)
    end do
  end subroutine cpwmin3

  ! M33INV and M44INV by David G. Simpson pure function version from
  ! https://fortranwiki.org/fortran/show/Matrix+inversion
  ! Invert 3x3 matrix
  function matinv39(a11, a12, a13, a21, a22, a23, a31, a32, a33) &
       result(B)
    real(kind=rp), intent(in) :: a11, a12, a13, a21, a22, a23, a31, a32, a33
    real(xp) :: A(3,3) !! Matrix
    real(rp) :: B(3,3) !! Inverse matrix
    A(1,1) = a11
    A(1,2) = a12
    A(1,3) = a13
    A(2,1) = a21
    A(2,2) = a22
    A(2,3) = a23
    A(3,1) = a31
    A(3,2) = a32
    A(3,3) = a33
    B = matinv3(A)
  end function matinv39

  !> Performs a direct calculation of the inverse of a 3×3 matrix.
  !! M33INV and M44INV by David G. Simpson pure function version from
  !! https://fortranwiki.org/fortran/show/Matrix+inversion
  !! Invert 3x3 matrix
  function matinv3(A) result(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(kind=xp), intent(in) :: A(3,3) !! Matrix
    real(kind=xp) :: B(3,3) !! Inverse matrix
    real(kind=xp) :: detinv

    ! Calculate the inverse determinant of the matrix
    ! first index x,y,z, second r, s, t
    detinv = 1.0_xp/real(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1),xp)
    ! Calculate the inverse of the matrix
    ! first index r, s, t, second x, y, z
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  end function matinv3

end module math
