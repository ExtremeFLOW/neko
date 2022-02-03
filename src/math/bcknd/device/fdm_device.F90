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
module fdm_device
  use num_types
  use utils
  use device
  use, intrinsic :: iso_c_binding
  implicit none
#ifdef HAVE_HIP
  interface
     subroutine hip_fdm_do_fast(e_d, r_d, s_d, d_d, nl, nelv) &
          bind(c, name='hip_fdm_do_fast')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: e_d, r_d, s_d, d_d
       integer(c_int) :: nl, nelv
     end subroutine hip_fdm_do_fast
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_fdm_do_fast(e_d, r_d, s_d, d_d, nl, nelv) &
          bind(c, name='cuda_fdm_do_fast')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: e_d, r_d, s_d, d_d
       integer(c_int) :: nl, nelv
     end subroutine cuda_fdm_do_fast
  end interface
#endif
contains

  subroutine fdm_do_fast_device(e, r, s, d, nl, ldim, nelv)
    integer, intent(in) :: nl, nelv, ldim
    real(kind=rp), intent(inout) :: e(nl**ldim, nelv)
    real(kind=rp), intent(inout) :: r(nl**ldim, nelv)
    real(kind=rp), intent(inout) :: s(nl*nl,2,ldim, nelv)
    real(kind=rp), intent(inout) :: d(nl**ldim, nelv)    
    integer ::  ie, nn, i
    type(c_ptr) :: e_d, r_d, s_d, d_d

    e_d = device_get_ptr(e, nelv)
    r_d = device_get_ptr(r, nelv)
    s_d = device_get_ptr(s, nelv)
    d_d = device_get_ptr(d, nelv)
    if (ldim .ne. 3) call neko_error('fdm dim not supported')

#ifdef HAVE_HIP
    call hip_fdm_do_fast(e_d, r_d, s_d, d_d, nl, nelv)
#elif HAVE_CUDA
    call cuda_fdm_do_fast(e_d, r_d, s_d, d_d, nl, nelv)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine fdm_do_fast_device

end module fdm_device
