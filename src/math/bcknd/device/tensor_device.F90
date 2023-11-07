! Copyright (c) 2022, The Neko Authors
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
module tensor_device
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  implicit none
  private

  public :: tnsr3d_device, tnsr3d_el_list_device
  
#ifdef HAVE_HIP
   interface
     subroutine hip_tnsr3d_el_list(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, elements, n_points) &
          bind(c, name='hip_tnsr3d_el_list')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: v_d, u_d, A_d, Bt_d, Ct_d, elements
       integer(c_int) :: nu, nv, n_points
     end subroutine hip_tnsr3d_el_list
  end interface
  interface
     subroutine hip_tnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv) &
          bind(c, name='hip_tnsr3d')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: v_d, u_d, A_d, Bt_d, Ct_d
       integer(c_int) :: nu, nv, nelv
     end subroutine hip_tnsr3d
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_tnsr3d_el_list(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, elements, n_points) &
          bind(c, name='cuda_tnsr3d_el_list')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: v_d, u_d, A_d, Bt_d, Ct_d, elements
       integer(c_int) :: nu, nv, n_points
     end subroutine cuda_tnsr3d_el_list
  end interface
  interface
     subroutine cuda_tnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv) &
          bind(c, name='cuda_tnsr3d')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: v_d, u_d, A_d, Bt_d, Ct_d
       integer(c_int) :: nu, nv, nelv
     end subroutine cuda_tnsr3d
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_tnsr3d_el_list(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, elements, n_points) &
          bind(c, name='opencl_tnsr3d_el_list')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: v_d, u_d, A_d, Bt_d, Ct_d, elements
       integer(c_int) :: nu, nv, n_points
     end subroutine opencl_tnsr3d_el_list
  end interface
  interface
     subroutine opencl_tnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv) &
          bind(c, name='opencl_tnsr3d')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: v_d, u_d, A_d, Bt_d, Ct_d
       integer(c_int) :: nu, nv, nelv
     end subroutine opencl_tnsr3d
  end interface
#endif
contains

  subroutine tnsr3d_device(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv)
    type(c_ptr) :: v_d, u_d, A_d, Bt_d, Ct_d
    integer(c_int) :: nu, nv, nelv
#ifdef HAVE_HIP
    call hip_tnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv)
#elif HAVE_CUDA
    call cuda_tnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv)
#elif HAVE_OPENCL
    call opencl_tnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine tnsr3d_device

  subroutine tnsr3d_el_list_device(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, elements, n_points)
    type(c_ptr) :: v_d, u_d, A_d, Bt_d, Ct_d, elements
    integer(c_int) :: nu, nv, n_points
#ifdef HAVE_HIP
    call hip_tnsr3d_el_list(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, elements, n_points)
#elif HAVE_CUDA
    call cuda_tnsr3d_el_list( v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, elements, n_points)
#elif HAVE_OPENCL
    call opencl_tnsr3d_el_list( v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, elements, n_points)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine tnsr3d_el_list_device


end module tensor_device
