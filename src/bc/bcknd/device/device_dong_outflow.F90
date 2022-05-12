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
module device_dong_outflow
  use num_types
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  interface
     subroutine hip_dong_outflow_apply_scalar(msk, x, normal_x, normal_y,&
                                              normal_z, u, v, w, uinf, delta, m) &
          bind(c, name='hip_dong_outflow_apply_scalar')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       integer(c_int) :: m
       real(kind=c_rp) :: uinf, delta 
       type(c_ptr), value :: msk, x, u, v, w, normal_x, normal_y, normal_z
     end subroutine hip_dong_outflow_apply_scalar
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_dong_outflow_apply_scalar(msk, x, normal_x, normal_y,&
                                              normal_z, u, v, w, uinf, delta, m) &
          bind(c, name='cuda_dong_outflow_apply_scalar')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       integer(c_int) :: m
       real(kind=c_rp) :: uinf, delta 
       type(c_ptr), value :: msk, x, u, v, w, normal_x, normal_y, normal_z
     end subroutine cuda_dong_outflow_apply_scalar
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_dong_outflow_apply_scalar(msk, x, normal_x, normal_y,&
                                              normal_z, u, v, w, uinf, delta, m) &
          bind(c, name='opencl_dong_outflow_apply_scalar')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       integer(c_int) :: m
       real(kind=c_rp) :: uinf, delta 
       type(c_ptr), value :: msk, x, u, v, w, normal_x, normal_y, normal_z
     end subroutine opencl_dong_outflow_apply_scalar
  end interface
#endif

contains

  subroutine device_dong_outflow_apply_scalar(msk, x, normal_x, normal_y,&
                                              normal_z, u, v, w, uinf, delta, m)
    integer(c_int) :: m
    real(kind=c_rp) :: uinf, delta 
    type(c_ptr) :: msk, x, u, v, w, normal_x, normal_y, normal_z

#ifdef HAVE_HIP
    call hip_dong_outflow_apply_scalar(msk, x, normal_x, normal_y,&
                                              normal_z, u, v, w, uinf, delta, m)
#elif HAVE_CUDA
    call cuda_dong_outflow_apply_scalar(msk, x, normal_x, normal_y,&
                                              normal_z, u, v, w, uinf, delta, m)
#elif HAVE_OPENCL
    !call opencl_dong_outflow_apply_scalar(msk, x, normal_x, normal_y,&
    !                                          normal_z, u, v, w, uinf, delta, m)
    call neko_error('No opencl device backend support for dong')
#else
    call neko_error('No device backend configured')
#endif
    
  end subroutine device_dong_outflow_apply_scalar
  
end module device_dong_outflow
