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
!> Device wrappers for basis-aware mixed boundary-condition constraint kernels.
!! These routines operate on a resolved mixed-node mask together with a local
!! orthonormal basis stored as nodewise normal and tangent vectors. They are
!! used to impose homogeneous or prescribed constraints on vector fields in the
!! local `(n, t1, t2)` basis on device backends.
module device_constrain_mixed_bc
  use num_types, only : c_rp, rp
  use utils, only : neko_error
  use device, only : glb_cmd_queue
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  implicit none
  private

#ifdef HAVE_CUDA
  interface
     subroutine cuda_constrain_mixed_bc_zero(mixed_msk, x, y, z, &
          constraint_n, constraint_t1, constraint_t2, n, t1, t2, m, strm) &
          bind(c, name='cuda_constrain_mixed_bc_zero')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       implicit none
       integer(c_int) :: constraint_n, constraint_t1, constraint_t2
       integer(c_int) :: m
       type(c_ptr), value :: mixed_msk, x, y, z, n, t1, t2, strm
     end subroutine cuda_constrain_mixed_bc_zero

     subroutine cuda_constrain_mixed_bc_set(mixed_msk, x, y, z, &
          constraint_n, constraint_t1, constraint_t2, n, t1, t2, &
          values_n, values_t1, values_t2, m, strm) &
          bind(c, name='cuda_constrain_mixed_bc_set')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       implicit none
       integer(c_int) :: constraint_n, constraint_t1, constraint_t2
       integer(c_int) :: m
       type(c_ptr), value :: mixed_msk, x, y, z, n, t1, t2
       type(c_ptr), value :: values_n, values_t1, values_t2, strm
     end subroutine cuda_constrain_mixed_bc_set

     subroutine cuda_constrain_mixed_bc_set_const(mixed_msk, x, y, z, &
          constraint_n, constraint_t1, constraint_t2, n, t1, t2, &
          value_n, value_t1, value_t2, m, strm) &
          bind(c, name='cuda_constrain_mixed_bc_set_const')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       implicit none
       real(c_rp) :: value_n, value_t1, value_t2
       integer(c_int) :: constraint_n, constraint_t1, constraint_t2
       integer(c_int) :: m
       type(c_ptr), value :: mixed_msk, x, y, z, n, t1, t2, strm
     end subroutine cuda_constrain_mixed_bc_set_const
  end interface
#endif

  public :: device_constrain_mixed_bc_zero
  public :: device_constrain_mixed_bc_set
  public :: device_constrain_mixed_bc_set_const

contains

  !> Constrain mixed-boundary projections by zeroing the selected local
  !! components on the device.
  !! @param[in] mixed_msk Device pointer to the resolved mixed-node mask.
  !! @param[inout] x Device pointer to the x-component field.
  !! @param[inout] y Device pointer to the y-component field.
  !! @param[inout] z Device pointer to the z-component field.
  !! @param[in] constraint_n Nonzero if the normal component is constrained.
  !! @param[in] constraint_t1 Nonzero if the first tangent component is
  !! constrained.
  !! @param[in] constraint_t2 Nonzero if the second tangent component is
  !! constrained.
  !! @param[in] n Device pointer to the nodewise normal vectors.
  !! @param[in] t1 Device pointer to the nodewise first tangent vectors.
  !! @param[in] t2 Device pointer to the nodewise second tangent vectors.
  !! @param[in] m Number of resolved mixed nodes.
  !! @param[in] strm Optional backend stream/queue to launch on.
  subroutine device_constrain_mixed_bc_zero(mixed_msk, x, y, z, &
       constraint_n, constraint_t1, constraint_t2, n, t1, t2, m, strm)
    integer, intent(in) :: constraint_n, constraint_t1, constraint_t2
    integer, intent(in) :: m
    type(c_ptr), intent(in) :: mixed_msk, x, y, z
    type(c_ptr), intent(in) :: n, t1, t2
    type(c_ptr), intent(in), optional :: strm
    type(c_ptr) :: strm_

    if (m .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#ifdef HAVE_CUDA
    call cuda_constrain_mixed_bc_zero(mixed_msk, x, y, z, &
         constraint_n, constraint_t1, constraint_t2, n, t1, t2, m, strm_)
#else
    call neko_error('CUDA backend not configured for mixed BC constraint')
#endif

  end subroutine device_constrain_mixed_bc_zero

  !> Constrain mixed-boundary projections by assigning prescribed local
  !! components derived from nodewise target vectors on the device.
  !! @param[in] mixed_msk Device pointer to the resolved mixed-node mask.
  !! @param[inout] x Device pointer to the x-component field.
  !! @param[inout] y Device pointer to the y-component field.
  !! @param[inout] z Device pointer to the z-component field.
  !! @param[in] constraint_n Nonzero if the normal component is constrained.
  !! @param[in] constraint_t1 Nonzero if the first tangent component is
  !! constrained.
  !! @param[in] constraint_t2 Nonzero if the second tangent component is
  !! constrained.
  !! @param[in] n Device pointer to the nodewise normal vectors.
  !! @param[in] t1 Device pointer to the nodewise first tangent vectors.
  !! @param[in] t2 Device pointer to the nodewise second tangent vectors.
  !! @param[in] values_n Device pointer to the x-components of the target
  !! vectors.
  !! @param[in] values_t1 Device pointer to the y-components of the target
  !! vectors.
  !! @param[in] values_t2 Device pointer to the z-components of the target
  !! vectors.
  !! @param[in] m Number of resolved mixed nodes.
  !! @param[in] strm Optional backend stream/queue to launch on.
  subroutine device_constrain_mixed_bc_set(mixed_msk, x, y, z, &
       constraint_n, constraint_t1, constraint_t2, n, t1, t2, &
       values_n, values_t1, values_t2, m, strm)
    integer, intent(in) :: constraint_n, constraint_t1, constraint_t2
    integer, intent(in) :: m
    type(c_ptr), intent(in) :: mixed_msk, x, y, z
    type(c_ptr), intent(in) :: n, t1, t2
    type(c_ptr), intent(in) :: values_n, values_t1, values_t2
    type(c_ptr), intent(in), optional :: strm
    type(c_ptr) :: strm_

    if (m .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#ifdef HAVE_CUDA
    call cuda_constrain_mixed_bc_set(mixed_msk, x, y, z, &
         constraint_n, constraint_t1, constraint_t2, n, t1, t2, &
         values_n, values_t1, values_t2, m, strm_)
#else
    call neko_error('CUDA backend not configured for mixed BC constraint')
#endif

  end subroutine device_constrain_mixed_bc_set

  !> Constrain mixed-boundary projections by assigning prescribed local
  !! components derived from a constant target vector on the device.
  !! @param[in] mixed_msk Device pointer to the resolved mixed-node mask.
  !! @param[inout] x Device pointer to the x-component field.
  !! @param[inout] y Device pointer to the y-component field.
  !! @param[inout] z Device pointer to the z-component field.
  !! @param[in] constraint_n Nonzero if the normal component is constrained.
  !! @param[in] constraint_t1 Nonzero if the first tangent component is
  !! constrained.
  !! @param[in] constraint_t2 Nonzero if the second tangent component is
  !! constrained.
  !! @param[in] n Device pointer to the nodewise normal vectors.
  !! @param[in] t1 Device pointer to the nodewise first tangent vectors.
  !! @param[in] t2 Device pointer to the nodewise second tangent vectors.
  !! @param[in] value_n x-component of the constant target vector.
  !! @param[in] value_t1 y-component of the constant target vector.
  !! @param[in] value_t2 z-component of the constant target vector.
  !! @param[in] m Number of resolved mixed nodes.
  !! @param[in] strm Optional backend stream/queue to launch on.
  subroutine device_constrain_mixed_bc_set_const(mixed_msk, x, y, z, &
       constraint_n, constraint_t1, constraint_t2, n, t1, t2, &
       value_n, value_t1, value_t2, m, strm)
    integer, intent(in) :: constraint_n, constraint_t1, constraint_t2
    integer, intent(in) :: m
    real(kind=rp), intent(in) :: value_n, value_t1, value_t2
    type(c_ptr), intent(in) :: mixed_msk, x, y, z
    type(c_ptr), intent(in) :: n, t1, t2
    type(c_ptr), intent(in), optional :: strm
    type(c_ptr) :: strm_

    if (m .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#ifdef HAVE_CUDA
    call cuda_constrain_mixed_bc_set_const(mixed_msk, x, y, z, &
         constraint_n, constraint_t1, constraint_t2, n, t1, t2, &
         value_n, value_t1, value_t2, m, strm_)
#else
    call neko_error('CUDA backend not configured for mixed BC constraint')
#endif

  end subroutine device_constrain_mixed_bc_set_const

end module device_constrain_mixed_bc
