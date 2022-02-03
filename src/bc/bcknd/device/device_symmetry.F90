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
module device_symmetry
  use num_types
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  interface
     subroutine hip_symmetry_apply_vector(xmsk, ymsk, zmsk, x, y, z, m, n, l) &
          bind(c, name='hip_symmetry_apply_vector')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, l
       type(c_ptr), value :: xmsk, ymsk, zmsk, x, y, z
     end subroutine hip_symmetry_apply_vector
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_symmetry_apply_vector(xmsk, ymsk, zmsk, x, y, z, m, n, l) &
          bind(c, name='cuda_symmetry_apply_vector')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, l
       type(c_ptr), value :: xmsk, ymsk, zmsk, x, y, z
     end subroutine cuda_symmetry_apply_vector
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_symmetry_apply_vector(xmsk, ymsk, zmsk, x, y, z, m, n, l) &
          bind(c, name='opencl_symmetry_apply_vector')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, l
       type(c_ptr), value :: xmsk, ymsk, zmsk, x, y, z
     end subroutine opencl_symmetry_apply_vector
  end interface
#endif
  
contains

  subroutine device_symmetry_apply_vector(xmsk, ymsk, zmsk, x, y, z, m, n, l)
    integer, intent(in) :: m, n, l
    type(c_ptr) :: xmsk, ymsk, zmsk, x, y, z

#ifdef HAVE_HIP
    call hip_symmetry_apply_vector(xmsk, ymsk, zmsk, x, y, z, m, n, l)
#elif HAVE_CUDA
    call cuda_symmetry_apply_vector(xmsk, ymsk, zmsk, x, y, z, m, n, l)
#elif HAVE_OPENCL
    call opencl_symmetry_apply_vector(xmsk, ymsk, zmsk, x, y, z, m, n, l)
#else
    call neko_error('No device backend configured')
#endif
    
  end subroutine device_symmetry_apply_vector
  
end module device_symmetry
