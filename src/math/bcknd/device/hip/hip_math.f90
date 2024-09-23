! Copyright (c) 2024, The Neko Authors
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
module hip_math
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding
  implicit none
  public
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
          w1_d, w2_d, w3_d, n) &
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

end module hip_math
