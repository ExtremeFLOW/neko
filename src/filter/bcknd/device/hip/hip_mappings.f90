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

!> Hip interface binding for mappings
module hip_mappings
  use num_types, only: rp, c_rp
  implicit none
  private

  public :: hip_smooth_step, hip_step_function, hip_permeability

  ! Interfaces for the backend functions
  interface
     subroutine hip_smooth_step(x, edge0, edge1, n) &
          bind(c, name = "hip_smooth_step")
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: x
       real(c_rp) :: edge0, edge1
       integer(c_int) :: n
     end subroutine hip_smooth_step

     subroutine hip_step_function(x, edge, left, right, n) &
          bind(c, name = "hip_step_function")
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: x
       real(c_rp) :: edge, left, right
       integer(c_int) :: n
     end subroutine hip_step_function

     subroutine hip_permeability(x, k_0, k_1, q, n) &
          bind(c, name = "hip_permeability")
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: x
       real(c_rp) :: k_0, k_1, q
       integer(c_int) :: n
     end subroutine hip_permeability
  end interface

end module hip_mappings
