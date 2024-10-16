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
module cuda_math
  use num_types, only: rp, c_rp
  implicit none
  public

  interface
     subroutine cuda_copy(a_d, b_d, n) &
          bind(c, name = 'cuda_copy')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine cuda_copy

     subroutine cuda_masked_copy(a_d, b_d, mask_d, n, m) &
          bind(c, name = 'cuda_masked_copy')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d, b_d, mask_d
       integer(c_int) :: n, m
     end subroutine cuda_masked_copy

     subroutine cuda_masked_red_copy(a_d, b_d, mask_d, n, m) &
          bind(c, name = 'cuda_masked_red_copy')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d, b_d, mask_d
       integer(c_int) :: n, m
     end subroutine cuda_masked_red_copy

     subroutine cuda_cfill_mask(a_d, c, size, mask_d, mask_size) &
          bind(c, name = 'cuda_cfill_mask')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: size
       type(c_ptr), value :: mask_d
       integer(c_int) :: mask_size
     end subroutine cuda_cfill_mask

     subroutine cuda_cmult(a_d, c, n) &
          bind(c, name = 'cuda_cmult')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine cuda_cmult

     subroutine cuda_cmult2(a_d, b_d, c, n) &
          bind(c, name = 'cuda_cmult2')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine cuda_cmult2

     subroutine cuda_cadd(a_d, c, n) &
          bind(c, name = 'cuda_cadd')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine cuda_cadd

     subroutine cuda_cadd2(a_d, b_d, c, n) &
          bind(c, name = 'cuda_cadd2')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d
       type(c_ptr), value :: b_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine cuda_cadd2

     subroutine cuda_cfill(a_d, c, n) &
          bind(c, name = 'cuda_cfill')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: n
     end subroutine cuda_cfill

     subroutine cuda_rzero(a_d, n) &
          bind(c, name = 'cuda_rzero')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine cuda_rzero

     subroutine cuda_add2(a_d, b_d, n) &
          bind(c, name = 'cuda_add2')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine cuda_add2

     subroutine cuda_add4(a_d, b_d, c_d, d_d, n) &
          bind(c, name = 'cuda_add4')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d, b_d, c_d, d_d
       integer(c_int) :: n
     end subroutine cuda_add4

     subroutine cuda_add2s1(a_d, b_d, c1, n) &
          bind(c, name = 'cuda_add2s1')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c1
       integer(c_int) :: n
     end subroutine cuda_add2s1

     subroutine cuda_add2s2(a_d, b_d, c1, n) &
          bind(c, name = 'cuda_add2s2')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c1
       integer(c_int) :: n
     end subroutine cuda_add2s2

     subroutine cuda_addsqr2s2(a_d, b_d, c1, n) &
          bind(c, name = 'cuda_addsqr2s2')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d, b_d
       real(c_rp) :: c1
       integer(c_int) :: n
     end subroutine cuda_addsqr2s2

     subroutine cuda_add3s2(a_d, b_d, c_d, c1, c2, n) &
          bind(c, name = 'cuda_add3s2')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d, b_d, c_d
       real(c_rp) :: c1, c2
       integer(c_int) :: n
     end subroutine cuda_add3s2

     subroutine cuda_invcol1(a_d, n) &
          bind(c, name = 'cuda_invcol1')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine cuda_invcol1

     subroutine cuda_invcol2(a_d, b_d, n) &
          bind(c, name = 'cuda_invcol2')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine cuda_invcol2

     subroutine cuda_col2(a_d, b_d, n) &
          bind(c, name = 'cuda_col2')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine cuda_col2

     subroutine cuda_col3(a_d, b_d, c_d, n) &
          bind(c, name = 'cuda_col3')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine cuda_col3

     subroutine cuda_subcol3(a_d, b_d, c_d, n) &
          bind(c, name = 'cuda_subcol3')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine cuda_subcol3

     subroutine cuda_sub2(a_d, b_d, n) &
          bind(c, name = 'cuda_sub2')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine cuda_sub2

     subroutine cuda_sub3(a_d, b_d, c_d, n) &
          bind(c, name = 'cuda_sub3')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine cuda_sub3

     subroutine cuda_add3(a_d, b_d, c_d, n) &
          bind(c, name = 'cuda_add3')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine cuda_add3

     subroutine cuda_addcol3(a_d, b_d, c_d, n) &
          bind(c, name = 'cuda_addcol3')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine cuda_addcol3

     subroutine cuda_addcol4(a_d, b_d, c_d, d_d, n) &
          bind(c, name = 'cuda_addcol4')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: a_d, b_d, c_d, d_d
       integer(c_int) :: n
     end subroutine cuda_addcol4

     subroutine cuda_vdot3(dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, n) &
          bind(c, name = 'cuda_vdot3')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d
       integer(c_int) :: n
     end subroutine cuda_vdot3

     subroutine cuda_vcross(u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, &
          w1_d, w2_d, w3_d, n) &
          bind(c, name = 'cuda_vcross')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: u1_d, u2_d, u3_d
       type(c_ptr), value :: v1_d, v2_d, v3_d
       type(c_ptr), value :: w1_d, w2_d, w3_d
       integer(c_int) :: n
     end subroutine cuda_vcross

     real(c_rp) function cuda_vlsc3(u_d, v_d, w_d, n) &
          bind(c, name = 'cuda_vlsc3')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d
       integer(c_int) :: n
     end function cuda_vlsc3

     subroutine cuda_add2s2_many(y_d, x_d_d, a_d, j, n) &
          bind(c, name = 'cuda_add2s2_many')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: y_d, x_d_d, a_d
       integer(c_int) :: j, n
     end subroutine cuda_add2s2_many

     real(c_rp) function cuda_glsc3(a_d, b_d, c_d, n) &
          bind(c, name = 'cuda_glsc3')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end function cuda_glsc3

     subroutine cuda_glsc3_many(h, w_d, v_d_d, mult_d, j, n) &
          bind(c, name = 'cuda_glsc3_many')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: w_d, v_d_d, mult_d
       integer(c_int) :: j, n
       real(c_rp) :: h(j)
     end subroutine cuda_glsc3_many

     real(c_rp) function cuda_glsc2(a_d, b_d, n) &
          bind(c, name = 'cuda_glsc2')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end function cuda_glsc2

     real(c_rp) function cuda_glsum(a_d, n) &
          bind(c, name = 'cuda_glsum')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end function cuda_glsum

     subroutine cuda_absval(a_d, n) &
          bind(c, name = 'cuda_absval')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine cuda_absval
  end interface
end module cuda_math
