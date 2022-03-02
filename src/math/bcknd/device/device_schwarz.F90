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
module device_schwarz
  use num_types
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  interface
     subroutine hip_schwarz_extrude(arr1_d,l1,f1,arr2_d,l2,f2,nx, nelv) &
          bind(c, name='hip_schwarz_extrude')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: arr1_d, arr2_d
       integer(c_int) :: l1, l2, nx, nelv
       real(c_rp) :: f1, f2
     end subroutine hip_schwarz_extrude
     subroutine hip_schwarz_toext3d(a_d,b_d,nx, nelv) &
          bind(c, name='hip_schwarz_toext3d')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d 
       integer(c_int) :: nx, nelv
     end subroutine hip_schwarz_toext3d
     subroutine hip_schwarz_toreg3d(b_d,a_d,nx, nelv) &
          bind(c, name='hip_schwarz_toreg3d')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d 
       integer(c_int) :: nx, nelv
     end subroutine hip_schwarz_toreg3d
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_schwarz_extrude(arr1_d,l1,f1,arr2_d,l2,f2,nx, nelv) &
          bind(c, name='cuda_schwarz_extrude')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: arr1_d, arr2_d
       integer(c_int) :: l1, l2, nx, nelv
       real(c_rp) :: f1, f2
     end subroutine cuda_schwarz_extrude
     subroutine cuda_schwarz_toext3d(a_d,b_d,nx, nelv) &
          bind(c, name='cuda_schwarz_toext3d')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d 
       integer(c_int) :: nx, nelv
     end subroutine cuda_schwarz_toext3d
     subroutine cuda_schwarz_toreg3d(b_d,a_d,nx, nelv) &
          bind(c, name='cuda_schwarz_toreg3d')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d 
       integer(c_int) :: nx, nelv
     end subroutine cuda_schwarz_toreg3d
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_schwarz_extrude(arr1_d,l1,f1,arr2_d,l2,f2,nx, nelv) &
          bind(c, name='opencl_schwarz_extrude')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: arr1_d, arr2_d
       integer(c_int) :: l1, l2, nx, nelv
       real(c_rp) :: f1, f2
     end subroutine opencl_schwarz_extrude
     subroutine opencl_schwarz_toext3d(a_d,b_d,nx, nelv) &
          bind(c, name='opencl_schwarz_toext3d')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d 
       integer(c_int) :: nx, nelv
     end subroutine opencl_schwarz_toext3d
     subroutine opencl_schwarz_toreg3d(b_d,a_d,nx, nelv) &
          bind(c, name='opencl_schwarz_toreg3d')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d 
       integer(c_int) :: nx, nelv
     end subroutine opencl_schwarz_toreg3d
  end interface
#endif
contains
  subroutine device_schwarz_extrude(arr1_d,l1,f1,arr2_d,l2,f2,nx,ny,nz, nelv)
    integer, intent(in) :: l1,l2,nx,ny,nz, nelv
    type(c_ptr), intent(inout) :: arr1_d,arr2_d
    real(kind=rp), intent(in) :: f1,f2
#ifdef HAVE_HIP
    call hip_schwarz_extrude(arr1_d,l1,f1,arr2_d,l2,f2,nx,nelv)
#elif HAVE_CUDA
    call cuda_schwarz_extrude(arr1_d,l1,f1,arr2_d,l2,f2,nx,nelv)
#elif HAVE_OPENCL
    call opencl_schwarz_extrude(arr1_d,l1,f1,arr2_d,l2,f2,nx,nelv)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_schwarz_extrude

  subroutine device_schwarz_toext3d(a_d,b_d,nx, nelv)
    integer, intent(in) :: nx, nelv
    type(c_ptr) :: a_d, b_d
#ifdef HAVE_HIP
    call hip_schwarz_toext3d(a_d,b_d,nx,nelv)
#elif HAVE_CUDA
    call cuda_schwarz_toext3d(a_d,b_d,nx,nelv)
#elif HAVE_OPENCL
    call opencl_schwarz_toext3d(a_d,b_d,nx,nelv)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_schwarz_toext3d

  subroutine device_schwarz_toreg3d(b_d,a_d,nx, nelv)
    integer, intent(in) :: nx, nelv
    type(c_ptr) :: a_d, b_d
#ifdef HAVE_HIP
    call hip_schwarz_toreg3d(b_d,a_d,nx,nelv)
#elif HAVE_CUDA
    call cuda_schwarz_toreg3d(b_d,a_d,nx,nelv)
#elif HAVE_OPENCL
    call opencl_schwarz_toreg3d(b_d,a_d,nx,nelv)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_schwarz_toreg3d
end module device_schwarz
