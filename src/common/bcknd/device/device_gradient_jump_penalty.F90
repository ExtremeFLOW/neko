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
module device_gradient_jump_penalty
  use num_types, only : c_rp, rp
  use utils, only : neko_error
  use device, only : glb_cmd_queue
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_pick_facet_value_hex(b_d, a_d, nx, nelv) &
          bind(c, name = 'hip_pick_facet_value_hex')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: nx, nelv
     end subroutine hip_pick_facet_value_hex
  end interface
  interface
     subroutine hip_gradient_jump_penalty_finalize(penalty_d, &
                                           penalty_facet_d, &
                                           dphidxi_d, &
                                           nx, nelv) &
          bind(c, name = 'hip_gradient_jump_penalty_finalize')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: penalty_d, &
                             penalty_facet_d, dphidxi_d
       integer(c_int) :: nx, nelv
     end subroutine hip_gradient_jump_penalty_finalize
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_pick_facet_value_hex(b_d, a_d, nx, nelv) &
          bind(c, name = 'cuda_pick_facet_value_hex')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: nx, nelv
     end subroutine cuda_pick_facet_value_hex
  end interface
  interface
     subroutine cuda_gradient_jump_penalty_finalize(penalty_d, &
                                           penalty_facet_d, &
                                           dphidxi_d, &
                                           nx, nelv) &
          bind(c, name = 'cuda_gradient_jump_penalty_finalize')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: penalty_d, &
                             penalty_facet_d, dphidxi_d
       integer(c_int) :: nx, nelv
     end subroutine cuda_gradient_jump_penalty_finalize
  end interface
#elif HAVE_OPENCL

#endif

  public :: device_pick_facet_value_hex, &
            device_gradient_jump_penalty_finalize

contains

  subroutine device_pick_facet_value_hex(b_d, a_d, nx, nelv)
    integer, intent(in) :: nx, nelv
    type(c_ptr) :: a_d, b_d
#ifdef HAVE_HIP
    call hip_pick_facet_value_hex(b_d, a_d, nx, nelv)
#elif HAVE_CUDA
    call cuda_pick_facet_value_hex(b_d, a_d, nx, nelv)
#elif HAVE_OPENCL
    call neko_error('OPENCL is not implemented for gradient jump penalty')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_pick_facet_value_hex

  subroutine device_gradient_jump_penalty_finalize(penalty_d, &
                                           penalty_facet_d, &
                                           dphidxi_d, &
                                           nx, nelv)
    integer, intent(in) :: nx, nelv
    type(c_ptr) :: penalty_d, penalty_facet_d, dphidxi_d
#ifdef HAVE_HIP
    call hip_gradient_jump_penalty_finalize(penalty_d, &
                                           penalty_facet_d, &
                                           dphidxi_d, &
                                           nx, nelv)
#elif HAVE_CUDA
    call cuda_gradient_jump_penalty_finalize(penalty_d, &
                                           penalty_facet_d, &
                                           dphidxi_d, &
                                           nx, nelv)
#elif HAVE_OPENCL
    call neko_error('OPENCL is not implemented for gradient jump penalty')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_gradient_jump_penalty_finalize

end module device_gradient_jump_penalty
