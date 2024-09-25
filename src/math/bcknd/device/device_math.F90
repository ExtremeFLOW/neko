! Copyright (c) 2021-2023, The Neko Authors
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
module device_math
  use comm
  use utils, only : neko_error
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_copy(a_d, b_d, n) &
          bind(c, name='hip_copy')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine hip_copy
  end interface

  interface
     subroutine hip_masked_copy(a_d, b_d, mask_d, n, m) &
          bind(c, name='hip_masked_copy')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d, b_d, mask_d
       integer(c_int) :: n, m
     end subroutine hip_masked_copy
  end interface
   
  interface
     subroutine hip_masked_red_copy(a_d, b_d, mask_d, n, m) &
          bind(c, name='hip_masked_red_copy')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d, b_d, mask_d
       integer(c_int) :: n, m
     end subroutine hip_masked_red_copy
  end interface
 
  interface
     subroutine hip_cfill_mask(a_d, c, size, mask_d, mask_size) &
          bind(c, name='hip_cfill_mask')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: size
       type(c_ptr), value :: mask_d
       integer(c_int) :: mask_size
     end subroutine hip_cfill_mask
  end interface

  interface
     subroutine hip_cmult(a_d, c, n) &
          bind(c, name='hip_cmult')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine hip_cmult
  end interface

  interface
     subroutine hip_cmult2(a_d, b_d, c, n) &
          bind(c, name='hip_cmult2')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine hip_cmult2
  end interface

  interface
     subroutine hip_cadd(a_d, c, n) &
          bind(c, name='hip_cadd')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine hip_cadd
  end interface

  interface
     subroutine hip_cadd2(a_d, b_d, c, n) &
          bind(c, name='hip_cadd2')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       type(c_ptr), value :: b_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine hip_cadd2
  end interface

  interface
     subroutine hip_cfill(a_d, c, n) &
          bind(c, name='hip_cfill')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine hip_cfill
  end interface

  interface
     subroutine hip_rzero(a_d, n) &
          bind(c, name='hip_rzero')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine hip_rzero
  end interface

  interface
     subroutine hip_add2(a_d, b_d, n) &
          bind(c, name='hip_add2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine hip_add2
  end interface

  interface
     subroutine hip_add4(a_d, b_d, c_d, d_d, n) &
          bind(c, name='hip_add4')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d, d_d
       integer(c_int) :: n
     end subroutine hip_add4
  end interface

  interface
     subroutine hip_add2s1(a_d, b_d, c1, n) &
          bind(c, name='hip_add2s1')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c1
       integer(c_int) :: n
     end subroutine hip_add2s1
  end interface

  interface
     subroutine hip_add2s2(a_d, b_d, c1, n) &
          bind(c, name='hip_add2s2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c1
       integer(c_int) :: n
     end subroutine hip_add2s2
  end interface

  interface
     subroutine hip_add2s2_many(y_d,x_d_d,a_d,j,n) &
          bind(c, name='hip_add2s2_many')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: y_d, x_d_d, a_d
       integer(c_int) :: j, n
     end subroutine hip_add2s2_many
  end interface

  interface
     subroutine hip_addsqr2s2(a_d, b_d, c1, n) &
          bind(c, name='hip_addsqr2s2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c1
       integer(c_int) :: n
     end subroutine hip_addsqr2s2
  end interface

  interface
     subroutine hip_add3s2(a_d, b_d, c_d, c1, c2, n) &
          bind(c, name='hip_add3s2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       real(c_rp) :: c1, c2
       integer(c_int) :: n
     end subroutine hip_add3s2
  end interface

  interface
     subroutine hip_invcol1(a_d, n) &
          bind(c, name='hip_invcol1')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine hip_invcol1
  end interface

  interface
     subroutine hip_invcol2(a_d, b_d, n) &
          bind(c, name='hip_invcol2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine hip_invcol2
  end interface

  interface
     subroutine hip_col2(a_d, b_d, n) &
          bind(c, name='hip_col2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine hip_col2
  end interface

  interface
     subroutine hip_col3(a_d, b_d, c_d, n) &
          bind(c, name='hip_col3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine hip_col3
  end interface

  interface
     subroutine hip_subcol3(a_d, b_d, c_d, n) &
          bind(c, name='hip_subcol3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine hip_subcol3
  end interface

  interface
     subroutine hip_sub2(a_d, b_d, n) &
          bind(c, name='hip_sub2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine hip_sub2
  end interface

  interface
     subroutine hip_sub3(a_d, b_d, c_d, n) &
          bind(c, name='hip_sub3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine hip_sub3
  end interface

  interface
     subroutine hip_add3(a_d, b_d, c_d, n) &
          bind(c, name='hip_add3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine hip_add3
  end interface

  interface
     subroutine hip_addcol3(a_d, b_d, c_d, n) &
          bind(c, name='hip_addcol3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine hip_addcol3
  end interface

  interface
     subroutine hip_addcol4(a_d, b_d, c_d, d_d, n) &
          bind(c, name='hip_addcol4')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d, d_d
       integer(c_int) :: n
     end subroutine hip_addcol4
  end interface

  interface
     subroutine hip_vdot3(dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, n) &
          bind(c, name='hip_vdot3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d
       integer(c_int) :: n
     end subroutine hip_vdot3
  end interface

  interface
     subroutine hip_vcross(u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, &
                           w1_d, w2_d, w3_d,  n) &
          bind(c, name='hip_vcross')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: u1_d, u2_d, u3_d
       type(c_ptr), value :: v1_d, v2_d, v3_d
       type(c_ptr), value :: w1_d, w2_d, w3_d
       integer(c_int) :: n
     end subroutine hip_vcross
  end interface

  interface
     real(c_rp) function hip_vlsc3(u_d, v_d, w_d, n) &
          bind(c, name='hip_vlsc3')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d
       integer(c_int) :: n
     end function hip_vlsc3
  end interface

  interface
     real(c_rp) function hip_glsc3(a_d, b_d, c_d, n) &
          bind(c, name='hip_glsc3')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end function hip_glsc3
  end interface

  interface
     subroutine hip_glsc3_many(h,w_d,v_d_d,mult_d,j,n) &
          bind(c, name='hip_glsc3_many')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: w_d, v_d_d, mult_d
       integer(c_int) :: j, n
       real(c_rp) :: h(j)
     end subroutine hip_glsc3_many
  end interface

  interface
     real(c_rp) function hip_glsc2(a_d, b_d, n) &
          bind(c, name='hip_glsc2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end function hip_glsc2
  end interface

  interface
     real(c_rp) function hip_glsum(a_d, n) &
          bind(c, name='hip_glsum')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end function hip_glsum
  end interface

  interface
     subroutine hip_absval(a_d, n) &
          bind(c, name = 'hip_absval')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine hip_absval
  end interface

#elif HAVE_CUDA
  interface
     subroutine cuda_copy(a_d, b_d, n) &
          bind(c, name='cuda_copy')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine cuda_copy
  end interface

  interface
     subroutine cuda_masked_copy(a_d, b_d, mask_d, n, m) &
          bind(c, name='cuda_masked_copy')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d, b_d, mask_d
       integer(c_int) :: n, m
     end subroutine cuda_masked_copy
  end interface

  interface
     subroutine cuda_masked_red_copy(a_d, b_d, mask_d, n, m) &
          bind(c, name='cuda_masked_red_copy')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d, b_d, mask_d
       integer(c_int) :: n, m
     end subroutine cuda_masked_red_copy
  end interface

  interface
     subroutine cuda_cfill_mask(a_d, c, size, mask_d, mask_size) &
          bind(c, name='cuda_cfill_mask')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: size
       type(c_ptr), value :: mask_d
       integer(c_int) :: mask_size
     end subroutine cuda_cfill_mask
  end interface

  interface
     subroutine cuda_cmult(a_d, c, n) &
          bind(c, name='cuda_cmult')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine cuda_cmult
  end interface

  interface
     subroutine cuda_cmult2(a_d, b_d, c, n) &
          bind(c, name='cuda_cmult2')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine cuda_cmult2
  end interface

  interface
     subroutine cuda_cadd(a_d, c, n) &
          bind(c, name='cuda_cadd')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine cuda_cadd
  end interface

  interface
     subroutine cuda_cadd2(a_d, b_d, c, n) &
          bind(c, name='cuda_cadd2')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       type(c_ptr), value :: b_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine cuda_cadd2
  end interface

  interface
     subroutine cuda_cfill(a_d, c, n) &
          bind(c, name='cuda_cfill')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine cuda_cfill
  end interface

  interface
     subroutine cuda_rzero(a_d, n) &
          bind(c, name='cuda_rzero')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine cuda_rzero
  end interface

  interface
     subroutine cuda_add2(a_d, b_d, n) &
          bind(c, name='cuda_add2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine cuda_add2
  end interface

  interface
     subroutine cuda_add4(a_d, b_d, c_d, d_d, n) &
          bind(c, name='cuda_add4')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d, d_d
       integer(c_int) :: n
     end subroutine cuda_add4
  end interface

  interface
     subroutine cuda_add2s1(a_d, b_d, c1, n) &
          bind(c, name='cuda_add2s1')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c1
       integer(c_int) :: n
     end subroutine cuda_add2s1
  end interface

  interface
     subroutine cuda_add2s2(a_d, b_d, c1, n) &
          bind(c, name='cuda_add2s2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c1
       integer(c_int) :: n
     end subroutine cuda_add2s2
  end interface

  interface
     subroutine cuda_addsqr2s2(a_d, b_d, c1, n) &
          bind(c, name='cuda_addsqr2s2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c1
       integer(c_int) :: n
     end subroutine cuda_addsqr2s2
  end interface

  interface
     subroutine cuda_add3s2(a_d, b_d, c_d, c1, c2, n) &
          bind(c, name='cuda_add3s2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       real(c_rp) :: c1, c2
       integer(c_int) :: n
     end subroutine cuda_add3s2
  end interface

  interface
     subroutine cuda_invcol1(a_d, n) &
          bind(c, name='cuda_invcol1')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine cuda_invcol1
  end interface

  interface
     subroutine cuda_invcol2(a_d, b_d, n) &
          bind(c, name='cuda_invcol2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine cuda_invcol2
  end interface

  interface
     subroutine cuda_col2(a_d, b_d, n) &
          bind(c, name='cuda_col2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine cuda_col2
  end interface

  interface
     subroutine cuda_col3(a_d, b_d, c_d, n) &
          bind(c, name='cuda_col3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine cuda_col3
  end interface

  interface
     subroutine cuda_subcol3(a_d, b_d, c_d, n) &
          bind(c, name='cuda_subcol3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine cuda_subcol3
  end interface

  interface
     subroutine cuda_sub2(a_d, b_d, n) &
          bind(c, name='cuda_sub2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine cuda_sub2
  end interface

  interface
     subroutine cuda_sub3(a_d, b_d, c_d, n) &
          bind(c, name='cuda_sub3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine cuda_sub3
  end interface

  interface
     subroutine cuda_add3(a_d, b_d, c_d, n) &
          bind(c, name='cuda_add3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine cuda_add3
  end interface

  interface
     subroutine cuda_addcol3(a_d, b_d, c_d, n) &
          bind(c, name='cuda_addcol3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine cuda_addcol3
  end interface

  interface
     subroutine cuda_addcol4(a_d, b_d, c_d, d_d, n) &
          bind(c, name='cuda_addcol4')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d, d_d
       integer(c_int) :: n
     end subroutine cuda_addcol4
  end interface

  interface
     subroutine cuda_vdot3(dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, n) &
          bind(c, name='cuda_vdot3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d
       integer(c_int) :: n
     end subroutine cuda_vdot3
  end interface

  interface
     subroutine cuda_vcross(u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, &
                           w1_d, w2_d, w3_d,  n) &
          bind(c, name='cuda_vcross')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: u1_d, u2_d, u3_d
       type(c_ptr), value :: v1_d, v2_d, v3_d
       type(c_ptr), value :: w1_d, w2_d, w3_d
       integer(c_int) :: n
     end subroutine cuda_vcross
  end interface

  interface
     real(c_rp) function cuda_vlsc3(u_d, v_d, w_d, n) &
          bind(c, name='cuda_vlsc3')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d
       integer(c_int) :: n
     end function cuda_vlsc3
  end interface

  interface
     subroutine cuda_add2s2_many(y_d,x_d_d,a_d,j,n) &
          bind(c, name='cuda_add2s2_many')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: y_d, x_d_d, a_d
       integer(c_int) :: j, n
     end subroutine cuda_add2s2_many
  end interface

  interface
     real(c_rp) function cuda_glsc3(a_d, b_d, c_d, n) &
          bind(c, name='cuda_glsc3')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end function cuda_glsc3
  end interface

  interface
     subroutine cuda_glsc3_many(h,w_d,v_d_d,mult_d,j,n) &
          bind(c, name='cuda_glsc3_many')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: w_d, v_d_d, mult_d
       integer(c_int) :: j, n
       real(c_rp) :: h(j)
     end subroutine cuda_glsc3_many
  end interface

  interface
     real(c_rp) function cuda_glsc2(a_d, b_d, n) &
          bind(c, name='cuda_glsc2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end function cuda_glsc2
  end interface

  interface
     real(c_rp) function cuda_glsum(a_d, n) &
          bind(c, name='cuda_glsum')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end function cuda_glsum
  end interface

  interface
     subroutine cuda_absval(a_d, n) &
          bind(c, name = 'cuda_absval')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine cuda_absval
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_copy(a_d, b_d, n) &
          bind(c, name='opencl_copy')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine opencl_copy
  end interface

  interface
     subroutine opencl_masked_copy(a_d, b_d, mask_d, n, m) &
          bind(c, name='opencl_masked_copy')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d, b_d, mask_d
       integer(c_int) :: n, m
     end subroutine opencl_masked_copy
  end interface

  interface
     subroutine opencl_cfill_mask(a_d, c, size, mask_d, mask_size) &
          bind(c, name='opencl_cfill_mask')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: size
       type(c_ptr), value :: mask_d
       integer(c_int) :: mask_size
     end subroutine opencl_cfill_mask
  end interface

  interface
     subroutine opencl_cmult(a_d, c, n) &
          bind(c, name='opencl_cmult')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine opencl_cmult
  end interface

  interface
     subroutine opencl_cmult2(a_d, b_d, c, n) &
          bind(c, name='opencl_cmult2')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine opencl_cmult2
  end interface

  interface
     subroutine opencl_cadd(a_d, c, n) &
          bind(c, name='opencl_cadd')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine opencl_cadd
  end interface

  interface
     subroutine opencl_cadd2(a_d, b_d, c, n) &
          bind(c, name='opencl_cadd2')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       type(c_ptr), value :: b_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine opencl_cadd2
  end interface

  interface
     subroutine opencl_cfill(a_d, c, n) &
          bind(c, name='opencl_cfill')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine opencl_cfill
  end interface

  interface
     subroutine opencl_rzero(a_d, n) &
          bind(c, name='opencl_rzero')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine opencl_rzero
  end interface

  interface
     subroutine opencl_rone(a_d, n) &
          bind(c, name='opencl_rone')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine opencl_rone
  end interface

  interface
     subroutine opencl_add2(a_d, b_d, n) &
          bind(c, name='opencl_add2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine opencl_add2
  end interface

  interface
     subroutine opencl_add4(a_d, b_d, c_d, d_d, n) &
          bind(c, name='opencl_add4')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d, d_d
       integer(c_int) :: n
     end subroutine opencl_add4
  end interface

  interface
     subroutine opencl_add2s1(a_d, b_d, c1, n) &
          bind(c, name='opencl_add2s1')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c1
       integer(c_int) :: n
     end subroutine opencl_add2s1
  end interface

  interface
     subroutine opencl_add2s2(a_d, b_d, c1, n) &
          bind(c, name='opencl_add2s2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c1
       integer(c_int) :: n
     end subroutine opencl_add2s2
  end interface

  interface
     subroutine opencl_add2s2_many(y_d,x_d_d,a_d,j,n) &
          bind(c, name='opencl_add2s2_many')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: y_d, x_d_d, a_d
       integer(c_int) :: j, n
     end subroutine opencl_add2s2_many
  end interface

  interface
     subroutine opencl_addsqr2s2(a_d, b_d, c1, n) &
          bind(c, name='opencl_addsqr2s2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c1
       integer(c_int) :: n
     end subroutine opencl_addsqr2s2
  end interface

  interface
     subroutine opencl_add3s2(a_d, b_d, c_d, c1, c2, n) &
          bind(c, name='opencl_add3s2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       real(c_rp) :: c1, c2
       integer(c_int) :: n
     end subroutine opencl_add3s2
  end interface

  interface
     subroutine opencl_invcol1(a_d, n) &
          bind(c, name='opencl_invcol1')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine opencl_invcol1
  end interface

  interface
     subroutine opencl_invcol2(a_d, b_d, n) &
          bind(c, name='opencl_invcol2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine opencl_invcol2
  end interface

  interface
     subroutine opencl_col2(a_d, b_d, n) &
          bind(c, name='opencl_col2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine opencl_col2
  end interface

  interface
     subroutine opencl_col3(a_d, b_d, c_d, n) &
          bind(c, name='opencl_col3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine opencl_col3
  end interface

  interface
     subroutine opencl_subcol3(a_d, b_d, c_d, n) &
          bind(c, name='opencl_subcol3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine opencl_subcol3
  end interface

  interface
     subroutine opencl_sub2(a_d, b_d, n) &
          bind(c, name='opencl_sub2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine opencl_sub2
  end interface

  interface
     subroutine opencl_sub3(a_d, b_d, c_d, n) &
          bind(c, name='opencl_sub3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine opencl_sub3
  end interface

  interface
     subroutine opencl_add3(a_d, b_d, c_d, n) &
          bind(c, name='opencl_add3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine opencl_add3
  end interface

  interface
     subroutine opencl_addcol3(a_d, b_d, c_d, n) &
          bind(c, name='opencl_addcol3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine opencl_addcol3
  end interface

  interface
     subroutine opencl_addcol4(a_d, b_d, c_d, d_d, n) &
          bind(c, name='opencl_addcol4')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d, d_d
       integer(c_int) :: n
     end subroutine opencl_addcol4
  end interface

  interface
     subroutine opencl_vdot3(dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, n) &
          bind(c, name='opencl_vdot3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d
       integer(c_int) :: n
     end subroutine opencl_vdot3
  end interface

  interface
     real(c_rp) function opencl_glsc3(a_d, b_d, c_d, n) &
          bind(c, name='opencl_glsc3')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end function opencl_glsc3
  end interface

  interface
     subroutine opencl_glsc3_many(h, w_d, v_d_d, mult_d, j, n) &
          bind(c, name='opencl_glsc3_many')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       integer(c_int) :: j, n
       type(c_ptr), value :: w_d, v_d_d, mult_d
       real(c_rp) :: h(j)
     end subroutine opencl_glsc3_many
  end interface

  interface
     real(c_rp) function opencl_glsc2(a_d, b_d, n) &
          bind(c, name='opencl_glsc2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end function opencl_glsc2
  end interface

  interface
     real(c_rp) function opencl_glsum(a_d, n) &
          bind(c, name='opencl_glsum')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end function opencl_glsum
  end interface
#endif

  public :: device_copy, device_rzero, device_rone, device_cmult, device_cmult2,&
       device_cadd, device_cadd2, device_cfill, device_add2, device_add3, &
       device_add4, device_add2s1, device_add2s2, device_addsqr2s2, &
       device_add3s2, device_invcol1, device_invcol2, device_col2, &
       device_col3, device_subcol3, device_sub2, device_sub3, device_addcol3, &
       device_addcol4, device_vdot3, device_vlsc3, device_glsc3, &
       device_glsc3_many, device_add2s2_many, device_glsc2, device_glsum, &
       device_masked_copy, device_cfill_mask, &
       device_masked_red_copy, device_vcross, &
       device_absval

contains

  !> Copy a vector \f$ a = b \f$
  subroutine device_copy(a_d, b_d, n)
    type(c_ptr) :: a_d, b_d
    integer :: n
#ifdef HAVE_HIP
    call hip_copy(a_d, b_d, n)
#elif HAVE_CUDA
    call cuda_copy(a_d, b_d, n)
#elif HAVE_OPENCL
    call opencl_copy(a_d, b_d, n)
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_copy

  !> Copy a masked vector \f$ a(mask) = b(mask) \f$.
  subroutine device_masked_copy(a_d, b_d, mask_d, n, m)
    type(c_ptr) :: a_d, b_d, mask_d
    integer :: n, m
#ifdef HAVE_HIP
    call hip_masked_copy(a_d, b_d, mask_d, n, m)
#elif HAVE_CUDA
    call cuda_masked_copy(a_d, b_d, mask_d, n, m)
#elif HAVE_OPENCL
    call opencl_masked_copy(a_d, b_d, mask_d, n, m)
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_masked_copy

  subroutine device_masked_red_copy(a_d, b_d, mask_d, n, m)
    type(c_ptr) :: a_d, b_d, mask_d
    integer :: n, m
#ifdef HAVE_HIP
    call hip_masked_red_copy(a_d, b_d, mask_d, n, m)
#elif HAVE_CUDA
    call cuda_masked_red_copy(a_d, b_d, mask_d, n, m)
#elif HAVE_OPENCL
    call neko_error('No OpenCL bcknd, masked red copy')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_masked_red_copy

  !> @brief Fill a constant to a masked vector.
  !! \f$ a_i = c, for i in mask \f$
  subroutine device_cfill_mask(a_d, c, size, mask_d, mask_size)
    type(c_ptr) :: a_d
    real(kind=rp), intent(in) :: c
    integer :: size
    type(c_ptr) :: mask_d
    integer :: mask_size
#ifdef HAVE_HIP
    call hip_cfill_mask(a_d, c, size, mask_d, mask_size)
#elif HAVE_CUDA
    call cuda_cfill_mask(a_d, c, size, mask_d, mask_size)
#elif HAVE_OPENCL
    call opencl_cfill_mask(a_d, c, size, mask_d, mask_size)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cfill_mask

  !> Zero a real vector
  subroutine device_rzero(a_d, n)
    type(c_ptr) :: a_d
    integer :: n
#ifdef HAVE_HIP
    call hip_rzero(a_d, n)
#elif HAVE_CUDA
    call cuda_rzero(a_d, n)
#elif HAVE_OPENCL
    call opencl_rzero(a_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_rzero

  !> Set all elements to one
  subroutine device_rone(a_d, n)
    type(c_ptr) :: a_d
    integer :: n
    real(kind=rp) :: one = 1.0_rp
#if defined(HAVE_HIP) || defined(HAVE_CUDA) || defined(HAVE_OPENCL)
    call device_cfill(a_d, one, n)
#elif HAVE_OPENCL
    call opencl_rone(a_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_rone

  !> Multiplication by constant c \f$ a = c \cdot a \f$
  subroutine device_cmult(a_d, c, n)
    type(c_ptr) :: a_d
    real(kind=rp), intent(in) :: c
    integer :: n
#ifdef HAVE_HIP
    call hip_cmult(a_d, c, n)
#elif HAVE_CUDA
    call cuda_cmult(a_d, c, n)
#elif HAVE_OPENCL
    call opencl_cmult(a_d, c, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cmult

  !> Multiplication by constant c \f$ a = c \cdot b \f$
  subroutine device_cmult2(a_d, b_d, c, n)
    type(c_ptr) :: a_d, b_d
    real(kind=rp), intent(in) :: c
    integer :: n
#ifdef HAVE_HIP
    call hip_cmult2(a_d, b_d, c, n)
#elif HAVE_CUDA
    call cuda_cmult2(a_d, b_d, c, n)
#elif HAVE_OPENCL
    call opencl_cmult2(a_d, b_d, c, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cmult2

  !> Add a scalar to vector \f$ a = a + s \f$
  subroutine device_cadd(a_d, c, n)
    type(c_ptr) :: a_d
    real(kind=rp), intent(in) :: c
    integer :: n
#ifdef HAVE_HIP
    call hip_cadd(a_d, c, n)
#elif HAVE_CUDA
    call cuda_cadd(a_d, c, n)
#elif HAVE_OPENCL
    call opencl_cadd(a_d, c, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cadd

  !> Add a scalar to vector \f$ a = b + s \f$
  subroutine device_cadd2(a_d, b_d, c, n)
    type(c_ptr) :: a_d
    type(c_ptr) :: b_d
    real(kind=rp), intent(in) :: c
    integer :: n
#ifdef HAVE_HIP
    call hip_cadd2(a_d, b_d, c, n)
#elif HAVE_CUDA
    call cuda_cadd2(a_d, b_d, c, n)
#elif HAVE_OPENCL
    call opencl_cadd2(a_d, b_d, c, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cadd2

  !> Set all elements to a constant c \f$ a = c \f$
  subroutine device_cfill(a_d, c, n)
    type(c_ptr) :: a_d
    real(kind=rp), intent(in) :: c
    integer :: n
#ifdef HAVE_HIP
    call hip_cfill(a_d, c, n)
#elif HAVE_CUDA
    call cuda_cfill(a_d, c, n)
#elif HAVE_OPENCL
    call opencl_cfill(a_d, c, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cfill

  !> Vector addition \f$ a = a + b \f$
  subroutine device_add2(a_d, b_d, n)
    type(c_ptr) :: a_d, b_d
    integer :: n
#ifdef HAVE_HIP
    call hip_add2(a_d, b_d, n)
#elif HAVE_CUDA
    call cuda_add2(a_d, b_d, n)
#elif HAVE_OPENCL
    call opencl_add2(a_d, b_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add2

  subroutine device_add4(a_d, b_d, c_d, d_d, n)
    type(c_ptr) :: a_d, b_d, c_d, d_d
    integer :: n
#ifdef HAVE_HIP
    call hip_add4(a_d, b_d, c_d, d_d, n)
#elif HAVE_CUDA
    call cuda_add4(a_d, b_d, c_d, d_d, n)
#elif HAVE_OPENCL
    call opencl_add4(a_d, b_d, c_d, d_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add4

  subroutine device_add2s1(a_d, b_d, c1, n)
    type(c_ptr) :: a_d, b_d
    real(kind=rp) :: c1
    integer :: n
#ifdef HAVE_HIP
    call hip_add2s1(a_d, b_d, c1, n)
#elif HAVE_CUDA
    call cuda_add2s1(a_d, b_d, c1, n)
#elif HAVE_OPENCL
    call opencl_add2s1(a_d, b_d, c1, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add2s1

  !> Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
  !! (multiplication on first argument)
  subroutine device_add2s2(a_d, b_d, c1, n)
    type(c_ptr) :: a_d, b_d
    real(kind=rp) :: c1
    integer :: n
#ifdef HAVE_HIP
    call hip_add2s2(a_d, b_d, c1, n)
#elif HAVE_CUDA
    call cuda_add2s2(a_d, b_d, c1, n)
#elif HAVE_OPENCL
    call opencl_add2s2(a_d, b_d, c1, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add2s2

  !> Returns \f$ a = a + c1 * (b * b )\f$
  subroutine device_addsqr2s2(a_d, b_d, c1, n)
    type(c_ptr) :: a_d, b_d
    real(kind=rp) :: c1
    integer :: n
#ifdef HAVE_HIP
    call hip_addsqr2s2(a_d, b_d, c1, n)
#elif HAVE_CUDA
    call cuda_addsqr2s2(a_d, b_d, c1, n)
#elif HAVE_OPENCL
    call opencl_addsqr2s2(a_d, b_d, c1, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_addsqr2s2

  !> Vector addition \f$ a = b + c \f$
  subroutine device_add3(a_d, b_d, c_d, n)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
#ifdef HAVE_HIP
    call hip_add3(a_d, b_d, c_d, n)
#elif HAVE_CUDA
    call cuda_add3(a_d, b_d, c_d, n)
#elif HAVE_OPENCL
    call opencl_add3(a_d, b_d, c_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add3

  !> Returns \f$ a = c1 * b + c2 * c \f$
  subroutine device_add3s2(a_d, b_d, c_d, c1, c2 ,n)
    type(c_ptr) :: a_d, b_d, c_d
    real(kind=rp) :: c1, c2
    integer :: n
#ifdef HAVE_HIP
    call hip_add3s2(a_d, b_d, c_d, c1, c2, n)
#elif HAVE_CUDA
    call cuda_add3s2(a_d, b_d, c_d, c1, c2, n)
#elif HAVE_OPENCL
    call opencl_add3s2(a_d, b_d, c_d, c1, c2, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add3s2

  !> Invert a vector \f$ a = 1 / a \f$
  subroutine device_invcol1(a_d, n)
    type(c_ptr) :: a_d
    integer :: n
#ifdef HAVE_HIP
    call hip_invcol1(a_d, n)
#elif HAVE_CUDA
    call cuda_invcol1(a_d, n)
#elif HAVE_OPENCL
    call opencl_invcol1(a_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_invcol1

  !> Vector division \f$ a = a / b \f$
  subroutine device_invcol2(a_d, b_d, n)
    type(c_ptr) :: a_d, b_d
    integer :: n
#ifdef HAVE_HIP
    call hip_invcol2(a_d, b_d, n)
#elif HAVE_CUDA
    call cuda_invcol2(a_d, b_d, n)
#elif HAVE_OPENCL
    call opencl_invcol2(a_d, b_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_invcol2

  !> Vector multiplication \f$ a = a \cdot b \f$
  subroutine device_col2(a_d, b_d, n)
    type(c_ptr) :: a_d, b_d
    integer :: n
#ifdef HAVE_HIP
    call hip_col2(a_d, b_d, n)
#elif HAVE_CUDA
    call cuda_col2(a_d, b_d, n)
#elif HAVE_OPENCL
    call opencl_col2(a_d, b_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_col2

  !> Vector multiplication with 3 vectors \f$ a =  b \cdot c \f$
  subroutine device_col3(a_d, b_d, c_d, n)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
#ifdef HAVE_HIP
    call hip_col3(a_d, b_d, c_d, n)
#elif HAVE_CUDA
    call cuda_col3(a_d, b_d, c_d, n)
#elif HAVE_OPENCL
    call opencl_col3(a_d, b_d, c_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_col3

  !> Returns \f$ a = a - b*c \f$
  subroutine device_subcol3(a_d, b_d, c_d, n)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
#ifdef HAVE_HIP
    call hip_subcol3(a_d, b_d, c_d, n)
#elif HAVE_CUDA
    call cuda_subcol3(a_d, b_d, c_d, n)
#elif HAVE_OPENCL
    call opencl_subcol3(a_d, b_d, c_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_subcol3

  !> Vector substraction \f$ a = a - b \f$
  subroutine device_sub2(a_d, b_d, n)
    type(c_ptr) :: a_d, b_d
    integer :: n
#ifdef HAVE_HIP
    call hip_sub2(a_d, b_d, n)
#elif HAVE_CUDA
    call cuda_sub2(a_d, b_d, n)
#elif HAVE_OPENCL
    call opencl_sub2(a_d, b_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_sub2

  !> Vector subtraction \f$ a = b - c \f$
  subroutine device_sub3(a_d, b_d, c_d, n)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
#ifdef HAVE_HIP
    call hip_sub3(a_d, b_d, c_d, n)
#elif HAVE_CUDA
    call cuda_sub3(a_d, b_d, c_d, n)
#elif HAVE_OPENCL
    call opencl_sub3(a_d, b_d, c_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_sub3

  !> Returns \f$ a = a + b*c \f$
  subroutine device_addcol3(a_d, b_d, c_d, n)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
#ifdef HAVE_HIP
    call hip_addcol3(a_d, b_d, c_d, n)
#elif HAVE_CUDA
    call cuda_addcol3(a_d, b_d, c_d, n)
#elif HAVE_OPENCL
    call opencl_addcol3(a_d, b_d, c_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_addcol3

  !> Returns \f$ a = a + b*c*d \f$
  subroutine device_addcol4(a_d, b_d, c_d, d_d, n)
    type(c_ptr) :: a_d, b_d, c_d, d_D
    integer :: n
#ifdef HAVE_HIP
    call hip_addcol4(a_d, b_d, c_d, d_d, n)
#elif HAVE_CUDA
    call cuda_addcol4(a_d, b_d, c_d, d_d, n)
#elif HAVE_OPENCL
    call opencl_addcol4(a_d, b_d, c_d, d_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_addcol4

  !> Compute a dot product \f$ dot = u \cdot v \f$ (3-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  subroutine device_vdot3(dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, n)
    type(c_ptr) :: dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d
    integer :: n
#ifdef HAVE_HIP
    call hip_vdot3(dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, n)
#elif HAVE_CUDA
    call cuda_vdot3(dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, n)
#elif HAVE_OPENCL
    call opencl_vdot3(dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_vdot3

  !> Compute a cross product \f$ u1, u2, u3 = v1,v2,v3 \cross w1,w2,w3 \f$ (3-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  subroutine device_vcross(u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, &
                           w1_d, w2_d, w3_d, n)
    type(c_ptr) :: u1_d, u2_d, u3_d
    type(c_ptr) :: v1_d, v2_d, v3_d
    type(c_ptr) :: w1_d, w2_d, w3_d
    integer :: n
#ifdef HAVE_HIP
    call hip_vcross(u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, & 
                    w1_d, w2_d, w3_d, n)
#elif HAVE_CUDA
    call cuda_vcross(u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, & 
                     w1_d, w2_d, w3_d, n)
#elif HAVE_OPENCL
    call neko_error("no opencl backedn vcross")
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_vcross


  !> Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
  function device_vlsc3(u_d, v_d, w_d, n) result(res)
    type(c_ptr) :: u_d, v_d, w_d
    integer :: n
    real(kind=rp) :: res
    res = 0.0_rp
#ifdef HAVE_HIP
    res = hip_vlsc3(u_d, v_d, w_d, n)
#elif HAVE_CUDA
    res = cuda_vlsc3(u_d, v_d, w_d, n)
#elif HAVE_OPENCL
    ! Same kernel as glsc3 (currently no device MPI for OpenCL)
    res = opencl_glsc3(u_d, v_d, w_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end function device_vlsc3

  !> Weighted inner product \f$ a^T b c \f$
  function device_glsc3(a_d, b_d, c_d, n) result(res)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n, ierr
    real(kind=rp) :: res
#ifdef HAVE_HIP
    res = hip_glsc3(a_d, b_d, c_d, n)
#elif HAVE_CUDA
    res = cuda_glsc3(a_d, b_d, c_d, n)
#elif HAVE_OPENCL
    res = opencl_glsc3(a_d, b_d, c_d, n)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif
  end function device_glsc3

  subroutine device_glsc3_many(h,w_d,v_d_d,mult_d,j,n)
    type(c_ptr), value :: w_d, v_d_d, mult_d
    integer(c_int) :: j, n
    real(c_rp) :: h(j)
    integer :: ierr
#ifdef HAVE_HIP
    call hip_glsc3_many(h,w_d,v_d_d,mult_d,j,n)
#elif HAVE_CUDA
    call cuda_glsc3_many(h,w_d,v_d_d,mult_d,j,n)
#elif HAVE_OPENCL
    call opencl_glsc3_many(h,w_d,v_d_d,mult_d,j,n)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, h, j, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif
  end subroutine device_glsc3_many

  subroutine device_add2s2_many(y_d,x_d_d,a_d,j,n)
    type(c_ptr), value :: y_d, x_d_d, a_d
    integer(c_int) :: j, n
#ifdef HAVE_HIP
    call hip_add2s2_many(y_d,x_d_d,a_d,j,n)
#elif HAVE_CUDA
    call cuda_add2s2_many(y_d,x_d_d,a_d,j,n)
#elif HAVE_OPENCL
    call opencl_add2s2_many(y_d,x_d_d,a_d,j,n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add2s2_many

  !> Weighted inner product \f$ a^T b \f$
  function device_glsc2(a_d, b_d, n) result(res)
    type(c_ptr) :: a_d, b_d
    integer :: n, ierr
    real(kind=rp) :: res
#ifdef HAVE_HIP
    res = hip_glsc2(a_d, b_d, n)
#elif HAVE_CUDA
    res = cuda_glsc2(a_d, b_d, n)
#elif HAVE_OPENCL
    res = opencl_glsc2(a_d, b_d, n)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif
  end function device_glsc2

  !> Sum a vector of length n
  function device_glsum(a_d, n) result(res)
    type(c_ptr) :: a_d
    integer :: n, ierr
    real(kind=rp) :: res
#ifdef HAVE_HIP
    res = hip_glsum(a_d, n)
#elif HAVE_CUDA
    res = cuda_glsum(a_d, n)
#elif HAVE_OPENCL
    res = opencl_glsum(a_d, n)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif
  end function device_glsum

  subroutine device_absval(a_d, n)
    integer, intent(in) :: n
    type(c_ptr) :: a_d
#ifdef HAVE_HIP
    call hip_absval(a_d, n)
#elif HAVE_CUDA
    call cuda_absval(a_d, n)
#elif HAVE_OPENCL
    call neko_error('OPENCL is not implemented for device_absval')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_absval


end module device_math
