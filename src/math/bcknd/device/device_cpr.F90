! Copyright (c) 2021-2022, The Neko Authors
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
module device_cpr
  use num_types
  use utils
  use comm
  use, intrinsic :: iso_c_binding
  implicit none

  interface
     subroutine cuda_glsc3_elem(res_d,a_d,b_d,c_d,lx,nelv) &
       bind(c, name='cuda_glsc3_elem')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: res_d, a_d, b_d, c_d
       integer(c_int) :: lx, nelv
     end subroutine cuda_glsc3_elem
  end interface


contains

  subroutine device_glsc3_elem(res_d,a_d,b_d,c_d,lx,nelv)
    type(c_ptr), value :: res_d,a_d, b_d, c_d
    integer(c_int) :: lx, nelv
#ifdef HAVE_HIP
    call neko_error('No elem glsc3_elem in this device')
#elif HAVE_CUDA
    call cuda_glsc3_elem(res_d,a_d,b_d,c_d,lx,nelv)
#elif HAVE_OPENCL
    call neko_error('No elem glsc3_elem in this device')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_glsc3_elem
  
  
 
end module device_cpr
