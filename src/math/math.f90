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
  use num_types, only: rp, dp, sp, qp, i4
  use host_math, only: host_abscmp, host_relcmp, host_rzero, host_izero, host_row_zero, &
       host_rone, host_copy, host_cmult, host_cadd, host_cfill, host_glsum, host_glmax, &
       host_glmin, host_chsign, host_vlmax, host_vlmin, host_invcol1, host_invcol3, &
       host_invers2, host_vcross, host_vdot2, host_vdot3, host_vlsc3, host_vlsc2, &
       host_add2, host_add3, host_add4, host_sub2, host_sub3, host_add2s1, host_add2s2, &
       host_addsqr2s2, host_cmult2, host_invcol2, host_col2, host_col3, host_subcol3, &
       host_add3s2, host_subcol4, host_addcol3, host_addcol4, host_ascol5, host_p_update, &
       host_x_update, host_glsc2, host_glsc3, host_glsc4, host_sortrp, host_sorti4, &
       host_swapdp, host_swapi4, host_reorddp, host_reordi4, host_flipvdp, host_flipvi4, &
       host_masked_copy, host_cfill_mask, host_glimax, host_glimin, host_cadd2, &
       host_qabscmp, host_qrelcmp, host_sabscmp, host_dabscmp, host_srelcmp, host_drelcmp, &
       host_masked_red_copy
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

  interface sort
     module procedure host_sortrp, host_sorti4
  end interface sort

  interface swap
     module procedure host_swapdp, host_swapi4
  end interface swap

  interface reord
     module procedure host_reorddp, host_reordi4
  end interface reord

  interface flipv
     module procedure host_flipvdp, host_flipvi4
  end interface flipv

  !> Return absolute comparison \f$ | x - y | < \epsilon \f$
  interface abscmp
     module procedure host_sabscmp, host_dabscmp, host_qabscmp
  end interface abscmp

  !> Return relative comparison \f$ | x - y |<= \epsilon*|y| \f$
  interface relcmp
     module procedure host_srelcmp, host_drelcmp, host_qrelcmp
  end interface relcmp

  !> Zero a real vector
  interface rzero
     module procedure host_rzero
  end interface rzero

  !> Zero an integer vector
  interface izero
     module procedure host_izero
  end interface izero

  !> Sets row e to 0 in matrix a
  interface row_zero
     module procedure host_row_zero
  end interface row_zero

  !> Set all elements to one
  interface rone
     module procedure host_rone
  end interface rone

  !> Copy a vector \f$ a = b \f$
  interface copy
     module procedure host_copy
  end interface copy

  !> Copy a masked vector \f$ a(mask) = b(mask) \f$.
  !! @param a Destination array of size `n`.
  !! @param b Source array of size `n`.
  !! @param mask Mask array of length m+1, where `mask(0) =m`
  !! the length of the mask array.
  !! @param n Size of the arrays `a` and `b`.
  !! @param m Size of the mask array `mask`.
  interface masked_copy
     module procedure host_masked_copy
  end interface masked_copy

  !> Copy a masked vector to reduced contigous vector
  !! \f$ a = b(mask) \f$.
  !! @param a Destination array of size `m`.
  !! @param b Source array of size `n`.
  !! @param mask Mask array of length m+1, where `mask(0) =m`
  !! the length of the mask array.
  !! @param n Size of the array `b`.
  !! @param m Size of the mask array `mask` and `a`.
  interface masked_red_copy
     module procedure host_masked_red_copy
  end interface masked_red_copy


  !> @brief Fill a constant to a masked vector.
  !! \f$ a_i = c, for i in mask \f$
  interface cfill_mask
     module procedure host_cfill_mask
  end interface cfill_mask

  !> Multiplication by constant c \f$ a = c \cdot a \f$
  interface cmult
     module procedure host_cmult
  end interface cmult

  !> Add a scalar to vector \f$ a_i = a_i + s \f$
  interface cadd
     module procedure host_cadd
  end interface cadd

  !> Add a scalar to vector \f$ a_i = b_i + s \f$
  interface cadd2
     module procedure host_cadd2

  end interface cadd2

  !> Set all elements to a constant c \f$ a = c \f$
  interface cfill
     module procedure host_cfill
  end interface cfill

  !> Sum a vector of length n
  interface glsum
     module procedure host_glsum
  end interface glsum

  !> Max of a vector of length n
  interface glmax
     module procedure host_glmax

  end interface glmax

  !> Max of an integer vector of length n
  interface glimax
     module procedure host_glimax

  end interface glimax

  !> Min of a vector of length n
  interface glmin
     module procedure host_glmin
  end interface glmin

  !> Min of an integer vector of length n
  interface glimin
     module procedure host_glimin
  end interface glimin


  !> Change sign of vector \f$ a = -a \f$
  interface chsign
     module procedure host_chsign

  end interface chsign

  !> maximum value of a vector of length @a n
  interface vlmax
     module procedure host_vlmax
  end interface vlmax

  !> minimun value of a vector of length @a n
  interface vlmin
     module procedure host_vlmin
  end interface vlmin

  !> Invert a vector \f$ a = 1 / a \f$
  interface invcol1
     module procedure host_invcol1
  end interface invcol1

  !> Invert a vector \f$ a = b / c \f$
  interface invcol3
     module procedure host_invcol3
  end interface invcol3

  !> Compute inverted vector \f$ a = 1 / b \f$
  interface invers2
     module procedure host_invers2
  end interface invers2

  !> Compute a cross product \f$ u = v \times w \f$
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  interface vcross
     module procedure host_vcross
  end interface vcross

  !> Compute a dot product \f$ dot = u \cdot v \f$ (2-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  interface vdot2
     module procedure host_vdot2
  end interface vdot2

  !> Compute a dot product \f$ dot = u \cdot v \f$ (3-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  interface vdot3
     module procedure host_vdot3
  end interface vdot3

  !> Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
  interface vlsc3
     module procedure host_vlsc3
  end interface vlsc3

  !> Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
  interface vlsc2
     module procedure host_vlsc2
  end interface vlsc2

  !> Vector addition \f$ a = a + b \f$
  interface add2
     module procedure host_add2
  end interface add2

  !> Vector addition \f$ a = b + c \f$
  interface add3
     module procedure host_add3
  end interface add3

  !> Vector addition \f$ a = b + c + d\f$
  interface add4
     module procedure host_add4
  end interface add4

  !> Vector substraction \f$ a = a - b \f$
  interface sub2
     module procedure host_sub2
  end interface sub2

  !> Vector subtraction \f$ a = b - c \f$
  interface sub3
     module procedure host_sub3
  end interface sub3

  !> Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
  !! (multiplication on first argument)
  interface add2s1
     module procedure host_add2s1
  end interface add2s1

  !> Vector addition with scalar multiplication \f$ a = a + c_1 b \f$
  !! (multiplication on second argument)
  interface add2s2
     module procedure host_add2s2
  end interface add2s2

  !> Returns \f$ a = a + c1 * (b * b)\f$
  interface addsqr2s2
     module procedure host_addsqr2s2
  end interface addsqr2s2

  !> Multiplication by constant c \f$ a = c \cdot b \f$
  interface cmult2
     module procedure host_cmult2
  end interface cmult2

  !> Vector division \f$ a = a / b \f$
  interface invcol2
     module procedure host_invcol2
  end interface invcol2

  !> Vector multiplication \f$ a = a \cdot b \f$
  interface col2
     module procedure host_col2
  end interface col2

  !> Vector multiplication with 3 vectors \f$ a = b \cdot c \f$
  interface col3
     module procedure host_col3
  end interface col3

  !> Returns \f$ a = a - b*c \f$
  interface subcol3
     module procedure host_subcol3
  end interface subcol3

  !> Returns \f$ a = c1 * b + c2 * c \f$
  interface add3s2
     module procedure host_add3s2
  end interface add3s2

  !> Returns \f$ a = a - b*c*d \f$
  interface subcol4
     module procedure host_subcol4
  end interface subcol4

  !> Returns \f$ a = a + b*c \f$
  interface addcol3
     module procedure host_addcol3
  end interface addcol3

  !> Returns \f$ a = a + b*c*d \f$
  interface addcol4
     module procedure host_addcol4
  end interface addcol4

  !> Returns \f$ a = b \dot c - d \cdot e \f$
  interface ascol5
     module procedure host_ascol5
  end interface ascol5

  !> Returns \f$ a = b \dot c1 (a - c2 \cdot c)\f$
  interface p_update
     module procedure host_p_update
  end interface p_update

  !> Returns \f$ a = b \dot c1 (a - c2 \cdot c)\f$
  interface x_update
     module procedure host_x_update
  end interface x_update

  !> Weighted inner product \f$ a^T b \f$
  interface glsc2
     module procedure host_glsc2
  end interface glsc2

  !> Weighted inner product \f$ a^T b c \f$
  interface glsc3
     module procedure host_glsc3
  end interface glsc3

  !> Weighted inner product \f$ a^T b c d \f$
  interface glsc4
     module procedure host_glsc4
  end interface glsc4
end module math
