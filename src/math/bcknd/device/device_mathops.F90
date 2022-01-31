! Copyright (c) 2021, The Neko Authors
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
module device_mathops
  use num_types
  use utils
  use, intrinsic :: iso_c_binding
  implicit none
  
#ifdef HAVE_HIP
  interface
     subroutine hip_opchsign(a1_d, a2_d, a3_d, gdim, n) &
          bind(c, name='hip_opchsign')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a1_d, a2_d, a3_d
       integer(c_int) :: gdim, n
     end subroutine hip_opchsign
  end interface

  interface
     subroutine hip_opcolv(a1_d, a2_d, a3_d, c_d, gdim, n) &
          bind(c, name='hip_opcolv')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a1_d, a2_d, a3_d, c_d
       integer(c_int) :: gdim, n
     end subroutine hip_opcolv
  end interface

  interface
     subroutine hip_opcolv3c(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, d, gdim, n) &
          bind(c, name='hip_opcolv3c')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d
       real(c_rp) :: d
       integer(c_int) :: gdim, n
     end subroutine hip_opcolv3c
  end interface

  interface
     subroutine hip_opadd2cm(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c, gdim, n) &
          bind(c, name='hip_opadd2cm')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d
       real(c_rp) :: c
       integer(c_int) :: gdim, n
     end subroutine hip_opadd2cm
  end interface

  interface
     subroutine hip_opadd2col(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, gdim, n) &
          bind(c, name='hip_opadd2col')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d
       integer(c_int) :: gdim, n
     end subroutine hip_opadd2col
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_opchsign(a1_d, a2_d, a3_d, gdim, n) &
          bind(c, name='cuda_opchsign')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a1_d, a2_d, a3_d
       integer(c_int) :: gdim, n
     end subroutine cuda_opchsign
  end interface

  interface
     subroutine cuda_opcolv(a1_d, a2_d, a3_d, c_d, gdim, n) &
          bind(c, name='cuda_opcolv')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a1_d, a2_d, a3_d, c_d
       integer(c_int) :: gdim, n
     end subroutine cuda_opcolv
  end interface

  interface
     subroutine cuda_opcolv3c(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, d, gdim, n) &
          bind(c, name='cuda_opcolv3c')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d
       real(c_rp) :: d
       integer(c_int) :: gdim, n
     end subroutine cuda_opcolv3c
  end interface

  interface
     subroutine cuda_opadd2cm(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c, gdim, n) &
          bind(c, name='cuda_opadd2cm')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d
       real(c_rp) :: c
       integer(c_int) :: gdim, n
     end subroutine cuda_opadd2cm
  end interface

  interface
     subroutine cuda_opadd2col(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, gdim, n) &
          bind(c, name='cuda_opadd2col')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d
       integer(c_int) :: gdim, n
     end subroutine cuda_opadd2col
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_opchsign(a1_d, a2_d, a3_d, gdim, n) &
          bind(c, name='opencl_opchsign')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a1_d, a2_d, a3_d
       integer(c_int) :: gdim, n
     end subroutine opencl_opchsign
  end interface

  interface
     subroutine opencl_opcolv(a1_d, a2_d, a3_d, c_d, gdim, n) &
          bind(c, name='opencl_opcolv')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a1_d, a2_d, a3_d, c_d
       integer(c_int) :: gdim, n
     end subroutine opencl_opcolv
  end interface

  interface
     subroutine opencl_opcolv3c(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, d, gdim, n) &
          bind(c, name='opencl_opcolv3c')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d
       real(c_rp) :: d
       integer(c_int) :: gdim, n
     end subroutine opencl_opcolv3c
  end interface

  interface
     subroutine opencl_opadd2cm(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c, gdim, n) &
          bind(c, name='opencl_opadd2cm')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d
       real(c_rp) :: c
       integer(c_int) :: gdim, n
     end subroutine opencl_opadd2cm
  end interface

  interface
     subroutine opencl_opadd2col(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, gdim, n) &
          bind(c, name='opencl_opadd2col')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d
       integer(c_int) :: gdim, n
     end subroutine opencl_opadd2col
  end interface
#endif

contains

  !> \f$ a = -a \f$
  subroutine device_opchsign(a1_d, a2_d, a3_d, gdim, n)
    type(c_ptr) :: a1_d, a2_d, a3_d
    integer :: n, gdim
#ifdef HAVE_HIP
    call hip_opchsign(a1_d, a2_d, a3_d, gdim, n)
#elif HAVE_CUDA
    call cuda_opchsign(a1_d, a2_d, a3_d, gdim, n)
#elif HAVE_OPENCL
    call opencl_opchsign(a1_d, a2_d, a3_d, gdim, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_opchsign

  !> \f$ a = a * c \f$
  subroutine device_opcolv(a1_d, a2_d, a3_d, c_d, gdim, n)
    type(c_ptr) :: a1_d, a2_d, a3_d, c_d
    integer :: n, gdim
#ifdef HAVE_HIP
    call hip_opcolv(a1_d, a2_d, a3_d, c_d, gdim, n)
#elif HAVE_CUDA
    call cuda_opcolv(a1_d, a2_d, a3_d, c_d, gdim, n)
#elif HAVE_OPENCL
    call opencl_opcolv(a1_d, a2_d, a3_d, c_d, gdim, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_opcolv

  !> \f$ a(i) = b(i) * c(i) * d \f$ 
  subroutine device_opcolv3c(a1_d, a2_d, a3_d, &
                             b1_d, b2_d, b3_d, c_d, d, n, gdim)
    type(c_ptr) :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d
    real(kind=rp) :: d
    integer :: n, gdim
#ifdef HAVE_HIP
    call hip_opcolv3c(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, d, gdim, n)
#elif HAVE_CUDA
    call cuda_opcolv3c(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, d, gdim, n)
#elif HAVE_OPENCL
    call opencl_opcolv3c(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, d, gdim, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_opcolv3c

  !> \f$ a(i) = a + b(i) * c \f$ 
  subroutine device_opadd2cm (a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c, n, gdim)
    type(c_ptr) :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d
    real(kind=rp) :: c
    integer :: n, gdim
#ifdef HAVE_HIP
    call hip_opadd2cm(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c, gdim, n)
#elif HAVE_CUDA
    call cuda_opadd2cm(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c, gdim, n)
#elif HAVE_OPENCL
    call opencl_opadd2cm(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c, gdim, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_opadd2cm

  !> \f$ a(i) = a + b(i) * c(i) \f$
  subroutine device_opadd2col (a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, n, gdim)
    type(c_ptr) :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d
    integer :: n, gdim
#ifdef HAVE_HIP
    call hip_opadd2col(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, gdim, n)
#elif HAVE_CUDA
    call cuda_opadd2col(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, gdim, n)
#elif HAVE_OPENCL
    call opencl_opadd2col(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, gdim, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_opadd2col

end module device_mathops
