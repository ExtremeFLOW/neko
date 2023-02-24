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
     subroutine cuda_lcsc3(res_d,a_d,b_d,c_d,lx,nelv) &
       bind(c, name='cuda_lcsc3')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: res_d, a_d, b_d, c_d
       integer(c_int) :: lx, nelv
     end subroutine cuda_lcsc3
  end interface

  interface
     subroutine cuda_lcsum(res_d,a_d,lx,nelv) &
       bind(c, name='cuda_lcsum')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: res_d, a_d
       integer(c_int) :: lx, nelv
     end subroutine cuda_lcsum
  end interface

  interface
     subroutine cuda_lcsort_abs(a_sort_d,key_sort_d,a_d,key_d,lx,nelv) &
       bind(c, name='cuda_lcsort_abs')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_sort_d, key_sort_d, a_d, key_d
       integer(c_int) :: lx, nelv
     end subroutine cuda_lcsort_abs
  end interface

  interface
     subroutine cuda_lcsort_bykey(a_sort_d,key_sort_d,a_d,key_d,lx,nelv) &
       bind(c, name='cuda_lcsort_bykey')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_sort_d, key_sort_d, a_d, key_d
       integer(c_int) :: lx, nelv
     end subroutine cuda_lcsort_bykey
  end interface

  interface
     subroutine cuda_lctnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d,bp_d,bp_key_d, nelv) &
          bind(c, name='cuda_lctnsr3d')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: v_d, u_d, A_d, Bt_d, Ct_d, bp_d, bp_key_d
       integer(c_int) :: nu, nv, nelv
     end subroutine cuda_lctnsr3d
  end interface

contains

  subroutine device_lcsc3(res_d,a_d,b_d,c_d,lx,nelv)
    type(c_ptr), value :: res_d,a_d, b_d, c_d
    integer(c_int) :: lx, nelv
#ifdef HAVE_HIP
    call neko_error('No elem glsc3_elem in this device')
#elif HAVE_CUDA
    call cuda_lcsc3(res_d,a_d,b_d,c_d,lx,nelv)
#elif HAVE_OPENCL
    call neko_error('No elem glsc3_elem in this device')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_lcsc3
   
  subroutine device_lcsum(res_d,a_d,lx,nelv)
    type(c_ptr), value :: res_d,a_d
    integer(c_int) :: lx, nelv
#ifdef HAVE_HIP
    call neko_error('No elem glsc3_elem in this device')
#elif HAVE_CUDA
    call cuda_lcsum(res_d,a_d,lx,nelv)
#elif HAVE_OPENCL
    call neko_error('No elem glsc3_elem in this device')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_lcsum
 
  subroutine device_lcsort_abs(a_sort_d,key_sort_d,a_d,key_d,lx,nelv)
    type(c_ptr), value :: a_sort_d,key_sort_d,a_d,key_d
    integer(c_int) :: lx, nelv
#ifdef HAVE_HIP
    call neko_error('No elem glsc3_elem in this device')
#elif HAVE_CUDA
    call cuda_lcsort_abs(a_sort_d,key_sort_d,a_d,key_d,lx,nelv)
#elif HAVE_OPENCL
    call neko_error('No elem glsc3_elem in this device')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_lcsort_abs

  subroutine device_lcsort_bykey(a_sort_d,key_sort_d,a_d,key_d,lx,nelv)
    type(c_ptr), value :: a_sort_d,key_sort_d,a_d,key_d
    integer(c_int) :: lx, nelv
#ifdef HAVE_HIP
    call neko_error('No elem glsc3_elem in this device')
#elif HAVE_CUDA
    call cuda_lcsort_bykey(a_sort_d,key_sort_d,a_d,key_d,lx,nelv)
#elif HAVE_OPENCL
    call neko_error('No elem glsc3_elem in this device')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_lcsort_bykey

  subroutine device_lctnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d,bp_d,bp_key_d, nelv)
    type(c_ptr) :: v_d, u_d, A_d, Bt_d, Ct_d,bp_d,bp_key_d
    integer(c_int) :: nu, nv, nelv
#ifdef HAVE_HIP
    call neko_error('No device backend configured')
#elif HAVE_CUDA
    call cuda_lctnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d,bp_d,bp_key_d, nelv)
#elif HAVE_OPENCL
    call neko_error('No device backend configured')
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_lctnsr3d

end module device_cpr
