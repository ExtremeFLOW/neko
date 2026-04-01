! Copyright (c) 2026, The Neko Authors
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
module device_coupled_vector_bc_resolver
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_null_ptr
  implicit none
  private

#ifdef HAVE_CUDA
  interface
     subroutine cuda_coupled_vector_bc_resolver_apply(mixed_msk, x, y, z, &
          constraint_n, constraint_t1, constraint_t2, n, t1, t2, m, strm) &
          bind(c, name='cuda_coupled_vector_bc_resolver_apply')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: mixed_msk, x, y, z
       type(c_ptr), value :: constraint_n, constraint_t1, constraint_t2
       type(c_ptr), value :: n, t1, t2, strm
     end subroutine cuda_coupled_vector_bc_resolver_apply
  end interface
#endif

  public :: device_coupled_vector_bc_resolver_apply

contains

  subroutine device_coupled_vector_bc_resolver_apply(mixed_msk, x, y, z, &
       constraint_n, constraint_t1, constraint_t2, n, t1, t2, m, strm)
    integer, intent(in) :: m
    type(c_ptr), intent(in) :: mixed_msk, x, y, z
    type(c_ptr), intent(in) :: constraint_n, constraint_t1, constraint_t2
    type(c_ptr), intent(in) :: n, t1, t2
    type(c_ptr), intent(in), optional :: strm
    type(c_ptr) :: strm_

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = c_null_ptr
    end if

#ifdef HAVE_CUDA
    call cuda_coupled_vector_bc_resolver_apply(mixed_msk, x, y, z, &
         constraint_n, constraint_t1, constraint_t2, n, t1, t2, m, strm_)
#else
    call neko_error('CUDA backend not configured for coupled vector BC apply')
#endif

  end subroutine device_coupled_vector_bc_resolver_apply

end module device_coupled_vector_bc_resolver
