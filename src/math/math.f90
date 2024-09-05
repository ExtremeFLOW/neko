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
  use num_types, only : rp, dp, sp, qp, i4
  use comm
  implicit none
  private

  public :: abscmp, rzero, izero, row_zero, rone, copy, cmult, cadd, cfill, &
       glsum, glmax, glmin, chsign, vlmax, vlmin, invcol1, invcol3, invers2, &
       vcross, vdot2, vdot3, vlsc3, vlsc2, add2, add3, add4, sub2, sub3, &
       add2s1, add2s2, addsqr2s2, cmult2, invcol2, col2, col3, subcol3, &
       add3s2, subcol4, addcol3, addcol4, ascol5, p_update, x_update, glsc2, &
       glsc3, glsc4, sort, masked_copy, cfill_mask, relcmp, glimax, glimin, &
       swap, reord, flipv, cadd2

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

  !> Return single precision absolute comparison \f$ | x - y | < \epsilon \f$
  interface sabscmp
     module function host_sabscmp(x, y) result(compare)
       real(kind=sp), intent(in) :: x
       real(kind=sp), intent(in) :: y
       logical :: compare
     end function host_sabscmp
  end interface sabscmp

  !> Return double precision absolute comparison \f$ | x - y | < \epsilon \f$
  interface dabscmp
     module function host_dabscmp(x, y) result(compare)
       real(kind=dp), intent(in) :: x
       real(kind=dp), intent(in) :: y
       logical :: compare
     end function host_dabscmp
  end interface dabscmp

  !> Return double precision absolute comparison \f$ | x - y | < \epsilon \f$
  interface qabscmp
     module function host_qabscmp(x, y) result(compare)
       real(kind=qp), intent(in) :: x
       real(kind=qp), intent(in) :: y
       logical :: compare
     end function host_qabscmp
  end interface qabscmp

  !> Return single precision relative comparison \f$ | x - y |<= \epsilon*|y| \f$
  interface srelcmp
     module function host_srelcmp(x, y, eps) result(compare)
       real(kind=sp), intent(in) :: x
       real(kind=sp), intent(in) :: y
       real(kind=sp), intent(in), optional :: eps
       logical :: compare
     end function host_srelcmp
  end interface srelcmp

  !> Return double precision relative comparison \f$ | x - y |/|y| < \epsilon \f$
  interface drelcmp
     module function host_drelcmp(x, y, eps) result(compare)
       real(kind=dp), intent(in) :: x
       real(kind=dp), intent(in) :: y
       real(kind=dp), intent(in), optional :: eps
       logical :: compare
     end function host_drelcmp
  end interface drelcmp


  !> Return quad precision relative comparison \f$ | x - y |/|y| < \epsilon \f$
  interface qrelcmp
     module function host_qrelcmp(x, y, eps) result(compare)
       real(kind=qp), intent(in) :: x
       real(kind=qp), intent(in) :: y
       real(kind=qp), intent(in), optional :: eps
       logical :: compare
     end function host_qrelcmp
  end interface qrelcmp

  !> Zero a real vector
  interface rzero
     module subroutine host_rzero(a, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
     end subroutine host_rzero
  end interface rzero

  !> Zero an integer vector
  interface izero
     module subroutine host_izero(a, n)
       integer, intent(in) :: n
       integer, dimension(n), intent(inout) :: a
     end subroutine host_izero
  end interface izero

  !> Sets row e to 0 in matrix a
  interface row_zero
     module subroutine host_row_zero(a, m, n, e)
       integer, intent(in) :: m, n, e
       real(kind=rp), intent(inout) :: a(m,n)
       integer :: j
     end subroutine host_row_zero
  end interface row_zero

  !> Set all elements to one
  interface rone
     module subroutine host_rone(a, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
     end subroutine host_rone
  end interface rone

  !> Copy a vector \f$ a = b \f$
  interface copy
     module subroutine host_copy(a, b, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(in) :: b
       real(kind=rp), dimension(n), intent(out) :: a
     end subroutine host_copy
  end interface copy

  !> Copy a masked vector \f$ a(mask) = b(mask) \f$.
  !! @param a Destination array of size `n`.
  !! @param b Source array of size `n`.
  !! @param mask Mask array of length m+1, where `mask(0) =m`
  !! the length of the mask array.
  !! @param n Size of the arrays `a` and `b`.
  !! @param m Size of the mask array `mask`.
  interface masked_copy
     module subroutine host_masked_copy(a, b, mask, n, m)
       integer, intent(in) :: n, m
       real(kind=rp), dimension(n), intent(in) :: b
       real(kind=rp), dimension(n), intent(out) :: a
       integer, dimension(0:m), intent(in) :: mask
     end subroutine host_masked_copy
  end interface masked_copy

  !> @brief Fill a constant to a masked vector.
  !! \f$ a_i = c, for i in mask \f$
  interface cfill_mask
     module subroutine host_cfill_mask(a, c, size, mask, mask_size)
       integer, intent(in) :: size, mask_size
       real(kind=rp), dimension(size), intent(inout) :: a
       real(kind=rp), intent(in) :: c
       integer, dimension(mask_size), intent(in) :: mask
     end subroutine host_cfill_mask
  end interface cfill_mask

  !> Multiplication by constant c \f$ a = c \cdot a \f$
  interface cmult
     module subroutine host_cmult(a, c, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), intent(in) :: c
     end subroutine host_cmult
  end interface cmult

  !> Add a scalar to vector \f$ a_i = a_i + s \f$
  interface cadd
     module subroutine host_cadd(a, s, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), intent(in) :: s
     end subroutine host_cadd
  end interface cadd

  !> Add a scalar to vector \f$ a_i = b_i + s \f$
  interface cadd2
     module subroutine host_cadd2(a, b, s, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b
       real(kind=rp), intent(in) :: s
       integer :: i
     end subroutine host_cadd2
  end interface cadd2

  !> Set all elements to a constant c \f$ a = c \f$
  interface cfill
     module subroutine host_cfill(a, c, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(out) :: a
       real(kind=rp), intent(in) :: c
     end subroutine host_cfill
  end interface cfill

  !> Sum a vector of length n
  interface glsum
     module function host_glsum(a, n) result(global_sum)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(in) :: a
       real(kind=rp) :: global_sum
     end function host_glsum
  end interface glsum

  !> Max of a vector of length n
  interface glmax
     module function host_glmax(a, n) result(glmax)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(in) :: a
       real(kind=rp) :: global_sum
       integer :: ierr
     end function host_glmax
  end interface glmax

  !> Max of an integer vector of length n
  interface glimax
     module function host_glimax(a, n) result(glmax)
       integer, intent(in) :: n
       integer, dimension(n), intent(in) :: a
       integer :: global_sum
       integer :: ierr
     end function host_glimax
  end interface glimax

  !> Min of a vector of length n
  interface glmin
     module function host_glmin(a, n) result(glmax)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(in) :: a
       real(kind=rp) :: global_sum
     end function host_glmin
  end interface glmin

  !> Min of an integer vector of length n
  interface glimin
     module function host_glimin(a, n) result(glmax)
       integer, intent(in) :: n
       integer, dimension(n), intent(in) :: a
       integer :: global_sum
     end function host_glimin
  end interface glimin




  !> Change sign of vector \f$ a = -a \f$
  interface chsign
     module subroutine host_chsign(a, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
     end subroutine host_chsign
  end interface chsign

  !> maximum value of a vector of length @a n
  interface vlmax
     module function host_vlmax(vec, n) result(tmax)
       integer :: n, i
       real(kind=rp), intent(in) :: vec(n)
       real(kind=rp) :: tmax
     end function host_vlmax
  end interface vlmax

  !> minimun value of a vector of length @a n
  interface vlmin
     module function host_vlmin(vec, n) result(tmin)
       integer, intent(in) :: n
       real(kind=rp), intent(in) :: vec(n)
       real(kind=rp) :: tmin
     end function host_vlmin
  end interface vlmin

  !> Invert a vector \f$ a = 1 / a \f$
  interface invcol1
     module subroutine host_invcol1(a, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
     end subroutine host_invcol1
  end interface invcol1

  !> Invert a vector \f$ a = b / c \f$
  interface invcol3
     module subroutine host_invcol3(a, b, c, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(out) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c
     end subroutine host_invcol3
  end interface invcol3

  !> Compute inverted vector \f$ a = 1 / b \f$
  interface invers2
     module subroutine host_invers2(a, b, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(out) :: a
       real(kind=rp), dimension(n), intent(in) :: b
     end subroutine host_invers2
  end interface invers2

  !> Compute a cross product \f$ u = v \times w \f$
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  interface vcross
     module subroutine host_vcross(u1, u2, u3, v1, v2, v3, w1, w2, w3, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(in) :: v1, v2, v3
       real(kind=rp), dimension(n), intent(in) :: w1, w2, w3
       real(kind=rp), dimension(n), intent(out) :: u1, u2, u3
     end subroutine host_vcross
  end interface vcross

  !> Compute a dot product \f$ dot = u \cdot v \f$ (2-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  interface vdot2
     module subroutine host_vdot2(dot, u1, u2, v1, v2, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(in) :: u1, u2
       real(kind=rp), dimension(n), intent(in) :: v1, v2
       real(kind=rp), dimension(n), intent(out) :: dot
     end subroutine host_vdot2
  end interface vdot2

  !> Compute a dot product \f$ dot = u \cdot v \f$ (3-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  interface vdot3
     module subroutine host_vdot3(dot, u1, u2, u3, v1, v2, v3, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(in) :: u1, u2, u3
       real(kind=rp), dimension(n), intent(in) :: v1, v2, v3
       real(kind=rp), dimension(n), intent(out) :: dot
     end subroutine host_vdot3
  end interface vdot3

  !> Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
  interface vlsc3
     module function host_vlsc3(u, v, w, n) result(s)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(in) :: u, v, w
       real(kind=rp) :: s
     end function host_vlsc3
  end interface vlsc3

  !> Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
  interface vlsc2
     module function host_vlsc2(u, v, n) result(s)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(in) :: u, v
       real(kind=rp) :: s
     end function host_vlsc2
  end interface vlsc2

  !> Vector addition \f$ a = a + b \f$
  interface add2
     module subroutine host_add2(a, b, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b
     end subroutine host_add2
  end interface add2

  !> Vector addition \f$ a = b + c \f$
  interface add3
     module subroutine host_add3(a, b, c, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(out) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c
     end subroutine host_add3
  end interface add3

  !> Vector addition \f$ a = b + c + d\f$
  interface add4
     module subroutine host_add4(a, b, c, d, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(out) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c, d
     end subroutine host_add4
  end interface add4

  !> Vector substraction \f$ a = a - b \f$
  interface sub2
     module subroutine host_sub2(a, b, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b
     end subroutine host_sub2
  end interface sub2

  !> Vector subtraction \f$ a = b - c \f$
  interface sub3
     module subroutine host_sub3(a, b, c, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(out) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c
     end subroutine host_sub3
  end interface sub3


  !> Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
  !! (multiplication on first argument)
  interface add2s1
     module subroutine host_add2s1(a, b, c1, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b
       real(kind=rp), intent(in) :: c1
     end subroutine host_add2s1
  end interface add2s1

  !> Vector addition with scalar multiplication \f$ a = a + c_1 b \f$
  !! (multiplication on second argument)
  interface add2s2
     module subroutine host_add2s2(a, b, c1, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b
       real(kind=rp), intent(in) :: c1
     end subroutine host_add2s2
  end interface add2s2

  !> Returns \f$ a = a + c1 * (b * b)\f$
  interface addsqr2s2
     module subroutine host_addsqr2s2(a, b, c1, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b
       real(kind=rp), intent(in) :: c1
     end subroutine host_addsqr2s2
  end interface addsqr2s2

  !> Multiplication by constant c \f$ a = c \cdot b \f$
  interface cmult2
     module subroutine host_cmult2(a, b, c, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(out) :: a
       real(kind=rp), dimension(n), intent(in) :: b
       real(kind=rp), intent(in) :: c
     end subroutine host_cmult2
  end interface cmult2

  !> Vector division \f$ a = a / b \f$
  interface invcol2
     module subroutine host_invcol2(a, b, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b
     end subroutine host_invcol2
  end interface invcol2


  !> Vector multiplication \f$ a = a \cdot b \f$
  interface col2
     module subroutine host_col2(a, b, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b
     end subroutine host_col2
  end interface col2

  !> Vector multiplication with 3 vectors \f$ a = b \cdot c \f$
  interface col3
     module subroutine host_col3(a, b, c, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(out) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c
     end subroutine host_col3
  end interface col3

  !> Returns \f$ a = a - b*c \f$
  interface subcol3
     module subroutine host_subcol3(a, b, c, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c
     end subroutine host_subcol3
  end interface subcol3

  !> Returns \f$ a = c1 * b + c2 * c \f$
  interface add3s2
     module subroutine host_add3s2(a, b, c, c1, c2 ,n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(out) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c
       real(kind=rp), intent(in) :: c1, c2
     end subroutine host_add3s2
  end interface add3s2


  !> Returns \f$ a = a - b*c*d \f$
  interface subcol4
     module subroutine host_subcol4(a, b, c, d, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c, d
     end subroutine host_subcol4
  end interface subcol4

  !> Returns \f$ a = a + b*c \f$
  interface addcol3
     module subroutine host_addcol3(a, b, c, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c
     end subroutine host_addcol3
  end interface addcol3

  !> Returns \f$ a = a + b*c*d \f$
  interface addcol4
     module subroutine host_addcol4(a, b, c, d, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c, d
     end subroutine host_addcol4
  end interface addcol4

  !> Returns \f$ a = b \dot c - d \cdot e \f$
  interface ascol5
     module subroutine host_ascol5(a, b, c, d, e, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(out) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c, d, e
     end subroutine host_ascol5
  end interface ascol5

  !> Returns \f$ a = b \dot c1 (a - c2 \cdot c)\f$
  interface p_update
     module subroutine host_p_update(a, b, c, c1, c2, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c
       real(kind=rp), intent(in) :: c1, c2
     end subroutine host_p_update
  end interface p_update

  !> Returns \f$ a = b \dot c1 (a - c2 \cdot c)\f$
  interface x_update
     module subroutine host_x_update(a, b, c, c1, c2, n)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: a
       real(kind=rp), dimension(n), intent(in) :: b, c
       real(kind=rp), intent(in) :: c1, c2
     end subroutine host_x_update
  end interface x_update

  !> Weighted inner product \f$ a^T b \f$
  interface glsc2
     module function host_glsc2(a, b, n) result(global_sum)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(in) :: a, b
       real(kind=rp) :: global_sum
     end function host_glsc2
  end interface glsc2

  !> Weighted inner product \f$ a^T b c \f$
  interface glsc3
     module function host_glsc3(a, b, c, n) result(global_sum)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(in) :: a, b, c
       real(kind=rp) :: global_sum
     end function host_glsc3
  end interface glsc3

  !> Weighted inner product \f$ a^T b c d \f$
  interface glsc4
     module function host_glsc4(a, b, c, d, n) result(global_sum)
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(in) :: a, b, c, d
       real(kind=rp) :: global_sum
     end function host_glsc4
  end interface glsc4

  !> Heap Sort for double precision arrays
  !! @details Following p 231 Num. Rec., 1st Ed.
  !! @param[inout]  a   vector to be sorted
  !! @param[out]   ind  permutation array
  !! @param[in]   n   array size
  interface sortrp
     module subroutine host_sortrp(a, ind, n)
       integer, intent(in) :: n
       real(kind=rp), intent(inout) :: a(n)
       integer, intent(out) :: ind(n)
     end subroutine host_sortrp
  end interface sortrp

  !> Heap Sort for single integer arrays
  !! @details Following p 231 Num. Rec., 1st Ed.
  !! @param[inout]  a   vector to be sorted
  !! @param[out]   ind  permutation array
  !! @param[in]   n   array size
  interface sorti4
     module subroutine host_sorti4(a, ind, n)
       integer, intent(in) :: n
       integer(i4), intent(inout) :: a(n)
       integer, intent(out) :: ind(n)
     end subroutine host_sorti4
  end interface sorti4

  !> sort double precision array acording to ind vector
  !! @param[inout]  b   vector to be reordered
  !! @param[in]   ind  permutation array
  !! @param[in]   n   array size
  interface swapdp
     module subroutine host_swapdp(b, ind, n)
       integer, intent(in) :: n
       real(kind=rp), intent(inout) :: b(n)
       integer, intent(in) :: ind(n)
     end subroutine host_swapdp
  end interface swapdp

  !> sort single integer array acording to ind vector
  !! @param[inout]  b   vector to be reordered
  !! @param[in]   ind  permutation array
  !! @param[in]   n   array size
  interface swapi4
     module subroutine host_swapi4(b, ind, n)
       integer, intent(in) :: n
       integer(i4), intent(inout) :: b(n)
       integer, intent(in) :: ind(n)
     end subroutine host_swapi4
  end interface swapi4

  !> reorder double precision array - inverse of swap
  !! @param[inout]  b   vector to be reordered
  !! @param[in]   ind  permutation array
  !! @param[in]   n   array size
  interface reorddp
     module subroutine host_reorddp(b, ind, n)
       integer, intent(in) :: n
       real(kind=rp), intent(inout) :: b(n)
       integer, intent(in) :: ind(n)
     end subroutine host_reorddp
  end interface reorddp

  !> reorder single integer array - inverse of swap
  !! @param[inout]  b   vector to be reordered
  !! @param[in]   ind  permutation array
  !! @param[in]   n   array size
  interface reordi4
     module subroutine host_reordi4(b, ind, n)
       integer, intent(in) :: n
       integer(i4), intent(inout) :: b(n)
       integer, intent(in) :: ind(n)
     end subroutine host_reordi4
  end interface reordi4

  !> Flip double precision vector b and ind
  !! @param[inout]  b   vector to be reordered
  !! @param[inout]  ind  permutation array
  !! @param[in]   n   array size
  interface flipvdp
     module subroutine host_flipvdp(b, ind, n)
       integer, intent(in) :: n
       real(kind=rp), intent(inout) :: b(n)
       integer, intent(inout) :: ind(n)
     end subroutine host_flipvdp
  end interface flipvdp

  !> Flip single integer vector b and ind
  !! @param[inout]  b   vector to be reordered
  !! @param[inout]  ind  permutation array
  !! @param[in]   n   array size
  interface flipvi4
     module subroutine host_flipvi4(b, ind, n)
       integer, intent(in) :: n
       integer(i4), intent(inout) :: b(n)
       integer, intent(inout) :: ind(n)
     end subroutine host_flipvi4
  end interface flipvi4

end module math
