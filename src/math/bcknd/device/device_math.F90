! Copyright (c) 2021-2025, The Neko Authors
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
module device_math
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
  use num_types, only : rp, c_rp
  use utils, only : neko_error
  use comm, only : NEKO_COMM, pe_size, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_SUM, MPI_IN_PLACE, MPI_Allreduce
  use device, only : glb_cmd_queue
  ! ========================================================================== !
  ! Device math interfaces

  use hip_math
  use cuda_math
  use opencl_math

  implicit none
  private

  interface device_cadd
     module procedure device_radd, device_iadd
  end interface device_cadd

  public :: device_copy, device_rzero, device_rone, device_cmult, &
       device_cmult2, device_cadd, device_cadd2, device_cfill, device_add2, &
       device_add3, device_add4, device_add2s1, device_add2s2, &
       device_addsqr2s2, device_add3s2, device_add4s3, device_add5s4, &
       device_invcol1, device_invcol2, device_col2, device_col3, &
       device_subcol3, device_sub2, device_sub3, device_addcol3, &
       device_addcol4, device_addcol3s2, device_vdot3, device_vlsc3, &
       device_glsc3, device_glsc3_many, device_add2s2_many, device_glsc2, &
       device_glsum, device_masked_copy_0, device_cfill_mask, &
       device_vcross, device_absval, device_masked_atomic_reduction_0, &
       device_masked_gather_copy_0, device_masked_scatter_copy_0, &
       device_invcol3, device_cdiv, device_cdiv2, device_glsubnorm, &
       device_pwmax2, device_pwmax3, device_cpwmax2, device_cpwmax3, &
       device_pwmin2, device_pwmin3, device_cpwmin2, device_cpwmin3

contains

  !> Copy a vector \f$ a = b \f$
  subroutine device_copy(a_d, b_d, n, strm)
    type(c_ptr) :: a_d, b_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_copy(a_d, b_d, n, strm_)
#elif HAVE_CUDA
    call cuda_copy(a_d, b_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_copy(a_d, b_d, n, strm_)
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_copy

  !> Copy a masked vector \f$ a(mask) = b(mask) \f$.
  subroutine device_masked_copy_0(a_d, b_d, mask_d, n, n_mask, strm)
    type(c_ptr) :: a_d, b_d, mask_d
    integer :: n, n_mask
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1 .or. n_mask .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_masked_copy(a_d, b_d, mask_d, n, n_mask, strm_)
#elif HAVE_CUDA
    call cuda_masked_copy(a_d, b_d, mask_d, n, n_mask, strm_)
#elif HAVE_OPENCL
    call opencl_masked_copy(a_d, b_d, mask_d, n, n_mask, strm_)
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_masked_copy_0

  !> Gather a masked vector \f$ a(i) = b(mask(i)) \f$.
  subroutine device_masked_gather_copy_0(a_d, b_d, mask_d, n, n_mask, strm)
    type(c_ptr) :: a_d, b_d, mask_d
    integer :: n, n_mask
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1 .or. n_mask .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_masked_gather_copy(a_d, b_d, mask_d, n, n_mask, strm_)
#elif HAVE_CUDA
    call cuda_masked_gather_copy(a_d, b_d, mask_d, n, n_mask, strm_)
#elif HAVE_OPENCL
    call opencl_masked_gather_copy(a_d, b_d, mask_d, n, n_mask, strm_)
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_masked_gather_copy_0

  !> Scatter a masked vector \f$ a((mask(i)) = b(i) \f$.
  subroutine device_masked_scatter_copy_0(a_d, b_d, mask_d, n, n_mask, strm)
    type(c_ptr) :: a_d, b_d, mask_d
    integer :: n, n_mask
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1 .or. n_mask .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_masked_scatter_copy(a_d, b_d, mask_d, n, n_mask, strm_)
#elif HAVE_CUDA
    call cuda_masked_scatter_copy(a_d, b_d, mask_d, n, n_mask, strm_)
#elif HAVE_OPENCL
    call opencl_masked_scatter_copy(a_d, b_d, mask_d, n, n_mask, strm_)
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_masked_scatter_copy_0

  subroutine device_masked_atomic_reduction_0(a_d, b_d, mask_d, n, n_mask, strm)
    type(c_ptr) :: a_d, b_d, mask_d
    integer :: n, n_mask
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1 .or. n_mask .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_masked_atomic_reduction(a_d, b_d, mask_d, n, n_mask, strm_)
#elif HAVE_CUDA
    call cuda_masked_atomic_reduction(a_d, b_d, mask_d, n, n_mask, strm_)
#elif HAVE_OPENCL
    call neko_error('No OpenCL bcknd, masked atomic reduction')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_masked_atomic_reduction_0

  !> @brief Fill a constant to a masked vector.
  !! \f$ a_i = c, for i in mask \f$
  subroutine device_cfill_mask(a_d, c, n, mask_d, n_mask, strm)
    type(c_ptr) :: a_d
    real(kind=rp), intent(in) :: c
    integer :: n
    type(c_ptr) :: mask_d
    integer :: n_mask
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1 .or. n_mask .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_cfill_mask(a_d, c, n, mask_d, n_mask, strm_)
#elif HAVE_CUDA
    call cuda_cfill_mask(a_d, c, n, mask_d, n_mask, strm_)
#elif HAVE_OPENCL
    call opencl_cfill_mask(a_d, c, n, mask_d, n_mask, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cfill_mask

  !> Zero a real vector
  subroutine device_rzero(a_d, n, strm)
    type(c_ptr) :: a_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_rzero(a_d, n, strm_)
#elif HAVE_CUDA
    call cuda_rzero(a_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_rzero(a_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_rzero

  !> Set all elements to one
  subroutine device_rone(a_d, n, strm)
    type(c_ptr) :: a_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_
    real(kind=rp), parameter :: one = 1.0_rp

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP || HAVE_CUDA || HAVE_OPENCL
    call device_cfill(a_d, one, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_rone

  !> Multiplication by constant c \f$ a = c \cdot a \f$
  subroutine device_cmult(a_d, c, n, strm)
    type(c_ptr) :: a_d
    real(kind=rp), intent(in) :: c
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_cmult(a_d, c, n, strm_)
#elif HAVE_CUDA
    call cuda_cmult(a_d, c, n, strm_)
#elif HAVE_OPENCL
    call opencl_cmult(a_d, c, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cmult

  !> Multiplication by constant c \f$ a = c \cdot b \f$
  subroutine device_cmult2(a_d, b_d, c, n, strm)
    type(c_ptr) :: a_d, b_d
    real(kind=rp), intent(in) :: c
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_cmult2(a_d, b_d, c, n, strm_)
#elif HAVE_CUDA
    call cuda_cmult2(a_d, b_d, c, n, strm_)
#elif HAVE_OPENCL
    call opencl_cmult2(a_d, b_d, c, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cmult2

  !> Division of constant c by array \f$ a = c / a \f$
  subroutine device_cdiv(a_d, c, n, strm)
    type(c_ptr) :: a_d
    real(kind=rp), intent(in) :: c
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_cdiv(a_d, c, n, strm_)
#elif HAVE_CUDA
    call cuda_cdiv(a_d, c, n, strm_)
#elif HAVE_OPENCL
    call opencl_cdiv(a_d, c, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cdiv

  !> Division of constant c by array \f$ a = c / b \f$
  subroutine device_cdiv2(a_d, b_d, c, n, strm)
    type(c_ptr) :: a_d, b_d
    real(kind=rp), intent(in) :: c
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_cdiv2(a_d, b_d, c, n, strm_)
#elif HAVE_CUDA
    call cuda_cdiv2(a_d, b_d, c, n, strm_)
#elif HAVE_OPENCL
    call opencl_cdiv2(a_d, b_d, c, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cdiv2

  !> Add a scalar to vector \f$ a = a + s \f$
  subroutine device_radd(a_d, c, n, strm)
    type(c_ptr) :: a_d
    real(kind=rp), intent(in) :: c
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_radd(a_d, c, n, strm_)
#elif HAVE_CUDA
    call cuda_radd(a_d, c, n, strm_)
#elif HAVE_OPENCL
    call opencl_radd(a_d, c, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_radd

  !> Add a scalar to vector \f$ a = b + s \f$
  subroutine device_cadd2(a_d, b_d, c, n, strm)
    type(c_ptr) :: a_d
    type(c_ptr) :: b_d
    real(kind=rp), intent(in) :: c
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_cadd2(a_d, b_d, c, n, strm_)
#elif HAVE_CUDA
    call cuda_cadd2(a_d, b_d, c, n, strm_)
#elif HAVE_OPENCL
    call opencl_cadd2(a_d, b_d, c, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cadd2

  !> Set all elements to a constant c \f$ a = c \f$
  subroutine device_cfill(a_d, c, n, strm)
    type(c_ptr) :: a_d
    real(kind=rp), intent(in) :: c
    integer :: n
    type(c_ptr), optional ::strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_cfill(a_d, c, n, strm_)
#elif HAVE_CUDA
    call cuda_cfill(a_d, c, n, strm_)
#elif HAVE_OPENCL
    call opencl_cfill(a_d, c, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cfill

  !> Vector addition \f$ a = a + b \f$
  subroutine device_add2(a_d, b_d, n, strm)
    type(c_ptr) :: a_d, b_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_add2(a_d, b_d, n, strm_)
#elif HAVE_CUDA
    call cuda_add2(a_d, b_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_add2(a_d, b_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add2

  subroutine device_add4(a_d, b_d, c_d, d_d, n, strm)
    type(c_ptr) :: a_d, b_d, c_d, d_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_add4(a_d, b_d, c_d, d_d, n, strm_)
#elif HAVE_CUDA
    call cuda_add4(a_d, b_d, c_d, d_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_add4(a_d, b_d, c_d, d_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add4

  subroutine device_add2s1(a_d, b_d, c1, n, strm)
    type(c_ptr) :: a_d, b_d
    real(kind=rp) :: c1
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_add2s1(a_d, b_d, c1, n, strm_)
#elif HAVE_CUDA
    call cuda_add2s1(a_d, b_d, c1, n, strm_)
#elif HAVE_OPENCL
    call opencl_add2s1(a_d, b_d, c1, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add2s1

  !> Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
  !! (multiplication on first argument)
  subroutine device_add2s2(a_d, b_d, c1, n, strm)
    type(c_ptr) :: a_d, b_d
    real(kind=rp) :: c1
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_add2s2(a_d, b_d, c1, n, strm_)
#elif HAVE_CUDA
    call cuda_add2s2(a_d, b_d, c1, n, strm_)
#elif HAVE_OPENCL
    call opencl_add2s2(a_d, b_d, c1, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add2s2

  !> Returns \f$ a = a + c1 * (b * b )\f$
  subroutine device_addsqr2s2(a_d, b_d, c1, n, strm)
    type(c_ptr) :: a_d, b_d
    real(kind=rp) :: c1
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_addsqr2s2(a_d, b_d, c1, n, strm_)
#elif HAVE_CUDA
    call cuda_addsqr2s2(a_d, b_d, c1, n, strm_)
#elif HAVE_OPENCL
    call opencl_addsqr2s2(a_d, b_d, c1, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_addsqr2s2

  !> Vector addition \f$ a = b + c \f$
  subroutine device_add3(a_d, b_d, c_d, n, strm)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_add3(a_d, b_d, c_d, n, strm_)
#elif HAVE_CUDA
    call cuda_add3(a_d, b_d, c_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_add3(a_d, b_d, c_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add3

  !> Returns \f$ a = c1 * b + c2 * c \f$
  subroutine device_add3s2(a_d, b_d, c_d, c1, c2 , n, strm)
    type(c_ptr) :: a_d, b_d, c_d
    real(kind=rp) :: c1, c2
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_add3s2(a_d, b_d, c_d, c1, c2, n, strm_)
#elif HAVE_CUDA
    call cuda_add3s2(a_d, b_d, c_d, c1, c2, n, strm_)
#elif HAVE_OPENCL
    call opencl_add3s2(a_d, b_d, c_d, c1, c2, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add3s2

  !> Returns \f$ a = c1 * b + c2 * c + c3 * d\f$
  subroutine device_add4s3(a_d, b_d, c_d, d_d, c1, c2 , c3, n, strm)
    type(c_ptr) :: a_d, b_d, c_d, d_d
    real(kind=rp) :: c1, c2, c3
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_add4s3(a_d, b_d, c_d, d_d, c1, c2, c3, n, strm_)
#elif HAVE_CUDA
    call cuda_add4s3(a_d, b_d, c_d, d_d, c1, c2, c3, n, strm_)
#elif HAVE_OPENCL
    call opencl_add4s3(a_d, b_d, c_d, d_d, c1, c2, c3, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add4s3

  !> Returns \f$ a = a + c1 * b + c2 * c + c3 * d + c4 * e\f$
  subroutine device_add5s4(a_d, b_d, c_d, d_d, e_d, c1, c2 , c3, c4, n, strm)
    type(c_ptr) :: a_d, b_d, c_d, d_d, e_d
    real(kind=rp) :: c1, c2, c3, c4
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_add5s4(a_d, b_d, c_d, d_d, e_d, c1, c2, c3, c4, n, strm_)
#elif HAVE_CUDA
    call cuda_add5s4(a_d, b_d, c_d, d_d, e_d, c1, c2, c3, c4, n, strm_)
#elif HAVE_OPENCL
    call opencl_add5s4(a_d, b_d, c_d, d_d, e_d, c1, c2, c3, c4, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add5s4

  !> Invert a vector \f$ a = 1 / a \f$
  subroutine device_invcol1(a_d, n, strm)
    type(c_ptr) :: a_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_invcol1(a_d, n, strm_)
#elif HAVE_CUDA
    call cuda_invcol1(a_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_invcol1(a_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_invcol1

  !> Vector division \f$ a = a / b \f$
  subroutine device_invcol2(a_d, b_d, n, strm)
    type(c_ptr) :: a_d, b_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_invcol2(a_d, b_d, n, strm_)
#elif HAVE_CUDA
    call cuda_invcol2(a_d, b_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_invcol2(a_d, b_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_invcol2

  !> Vector division \f$ a = b / c \f$
  subroutine device_invcol3(a_d, b_d, c_d, n, strm)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#ifdef HAVE_HIP
    call hip_invcol3(a_d, b_d, c_d, n, strm_)
#elif HAVE_CUDA
    call cuda_invcol3(a_d, b_d, c_d, n, strm_)
#elif HAVE_OPENCL
    ! call opencl_invcol3(a_d, b_d, c_d, n)
    call neko_error('opencl_invcol3 not implemented')
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_invcol3

  !> Vector multiplication \f$ a = a \cdot b \f$
  subroutine device_col2(a_d, b_d, n, strm)
    type(c_ptr) :: a_d, b_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

    if (n .lt. 1) return
#if HAVE_HIP
    call hip_col2(a_d, b_d, n, strm_)
#elif HAVE_CUDA
    call cuda_col2(a_d, b_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_col2(a_d, b_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_col2

  !> Vector multiplication with 3 vectors \f$ a =  b \cdot c \f$
  subroutine device_col3(a_d, b_d, c_d, n, strm)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_col3(a_d, b_d, c_d, n, strm_)
#elif HAVE_CUDA
    call cuda_col3(a_d, b_d, c_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_col3(a_d, b_d, c_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_col3

  !> Returns \f$ a = a - b*c \f$
  subroutine device_subcol3(a_d, b_d, c_d, n, strm)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_subcol3(a_d, b_d, c_d, n, strm_)
#elif HAVE_CUDA
    call cuda_subcol3(a_d, b_d, c_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_subcol3(a_d, b_d, c_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_subcol3

  !> Vector substraction \f$ a = a - b \f$
  subroutine device_sub2(a_d, b_d, n, strm)
    type(c_ptr) :: a_d, b_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_sub2(a_d, b_d, n, strm_)
#elif HAVE_CUDA
    call cuda_sub2(a_d, b_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_sub2(a_d, b_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_sub2

  !> Vector subtraction \f$ a = b - c \f$
  subroutine device_sub3(a_d, b_d, c_d, n, strm)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_sub3(a_d, b_d, c_d, n, strm_)
#elif HAVE_CUDA
    call cuda_sub3(a_d, b_d, c_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_sub3(a_d, b_d, c_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_sub3

  !> Returns \f$ a = a + b*c \f$
  subroutine device_addcol3(a_d, b_d, c_d, n, strm)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_addcol3(a_d, b_d, c_d, n, strm_)
#elif HAVE_CUDA
    call cuda_addcol3(a_d, b_d, c_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_addcol3(a_d, b_d, c_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_addcol3

  !> Returns \f$ a = a + b*c*d \f$
  subroutine device_addcol4(a_d, b_d, c_d, d_d, n, strm)
    type(c_ptr) :: a_d, b_d, c_d, d_D
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_addcol4(a_d, b_d, c_d, d_d, n, strm_)
#elif HAVE_CUDA
    call cuda_addcol4(a_d, b_d, c_d, d_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_addcol4(a_d, b_d, c_d, d_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_addcol4

  !> Returns \f$ a = a + s(b*c) \f$
  subroutine device_addcol3s2(a_d, b_d, c_d, s, n, strm)
    type(c_ptr) :: a_d, b_d, c_d
    real(kind=rp) :: s
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_addcol3s2(a_d, b_d, c_d, s, n, strm_)
#elif HAVE_CUDA
    call cuda_addcol3s2(a_d, b_d, c_d, s, n, strm_)
#elif HAVE_OPENCL
    call opencl_addcol3s2(a_d, b_d, c_d, s, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_addcol3s2

  !> Compute a dot product \f$ dot = u \cdot v \f$ (3-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  subroutine device_vdot3(dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, n, strm)
    type(c_ptr) :: dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_vdot3(dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, n, strm_)
#elif HAVE_CUDA
    call cuda_vdot3(dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_vdot3(dot_d, u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_vdot3

  !> Compute a cross product \f$ u = v \times w \f$ (3-d version)
  !! assuming vector components \f$ u = (u_1, u_2, u_3) \f$ etc.
  subroutine device_vcross(u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, &
       w1_d, w2_d, w3_d, n, strm)
    type(c_ptr) :: u1_d, u2_d, u3_d
    type(c_ptr) :: v1_d, v2_d, v3_d
    type(c_ptr) :: w1_d, w2_d, w3_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_vcross(u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, &
         w1_d, w2_d, w3_d, n, strm_)
#elif HAVE_CUDA
    call cuda_vcross(u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, &
         w1_d, w2_d, w3_d, n, strm_)
#elif HAVE_OPENCL
    call opencl_vcross(u1_d, u2_d, u3_d, v1_d, v2_d, v3_d, &
         w1_d, w2_d, w3_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_vcross


  !> Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
  function device_vlsc3(u_d, v_d, w_d, n, strm) result(res)
    type(c_ptr) :: u_d, v_d, w_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_
    real(kind=rp) :: res

    res = 0.0_rp

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    res = hip_vlsc3(u_d, v_d, w_d, n, strm_)
#elif HAVE_CUDA
    res = cuda_vlsc3(u_d, v_d, w_d, n, strm_)
#elif HAVE_OPENCL
    ! Same kernel as glsc3 (currently no device MPI for OpenCL)
    res = opencl_glsc3(u_d, v_d, w_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end function device_vlsc3

  !> Weighted inner product \f$ a^T b c \f$
  function device_glsc3(a_d, b_d, c_d, n, strm) result(res)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n, ierr
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_
    real(kind=rp) :: res

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

    res = 0.0_rp
#if HAVE_HIP
    res = hip_glsc3(a_d, b_d, c_d, n, strm_)
#elif HAVE_CUDA
    res = cuda_glsc3(a_d, b_d, c_d, n, strm_)
#elif HAVE_OPENCL
    res = opencl_glsc3(a_d, b_d, c_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif
  end function device_glsc3

  subroutine device_glsc3_many(h, w_d, v_d_d, mult_d, j, n, strm)
    type(c_ptr), value :: w_d, v_d_d, mult_d
    integer(c_int) :: j, n
    real(c_rp) :: h(j)
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_
    integer :: ierr

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_glsc3_many(h, w_d, v_d_d, mult_d, j, n, strm_)
#elif HAVE_CUDA
    call cuda_glsc3_many(h, w_d, v_d_d, mult_d, j, n, strm_)
#elif HAVE_OPENCL
    call opencl_glsc3_many(h, w_d, v_d_d, mult_d, j, n, strm_)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, h, j, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif
  end subroutine device_glsc3_many

  subroutine device_add2s2_many(y_d, x_d_d, a_d, j, n, strm)
    type(c_ptr), value :: y_d, x_d_d, a_d
    integer(c_int) :: j, n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_add2s2_many(y_d, x_d_d, a_d, j, n, strm_)
#elif HAVE_CUDA
    call cuda_add2s2_many(y_d, x_d_d, a_d, j, n, strm_)
#elif HAVE_OPENCL
    call opencl_add2s2_many(y_d, x_d_d, a_d, j, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add2s2_many

  !> Weighted inner product \f$ a^T b \f$
  function device_glsc2(a_d, b_d, n, strm) result(res)
    type(c_ptr) :: a_d, b_d
    integer :: n, ierr
    real(kind=rp) :: res
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

    res = 0.0_rp
#if HAVE_HIP
    res = hip_glsc2(a_d, b_d, n, strm_)
#elif HAVE_CUDA
    res = cuda_glsc2(a_d, b_d, n, strm_)
#elif HAVE_OPENCL
    res = opencl_glsc2(a_d, b_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif
  end function device_glsc2

  !> Returns the norm of the difference of two vectors
  !! \f$ \sqrt{(a-b)^T (a-b)} \f$
  function device_glsubnorm(a_d, b_d, n, strm) result(res)
    type(c_ptr), intent(in) :: a_d, b_d
    integer, intent(in) :: n
    integer :: ierr
    real(kind=rp) :: res
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

    res = 0.0_rp
#if HAVE_HIP
    res = hip_glsubnorm2(a_d, b_d, n, strm_)
#elif HAVE_CUDA
    res = cuda_glsubnorm2(a_d, b_d, n, strm_)
#elif HAVE_OPENCL
    res = opencl_glsubnorm2(a_d, b_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif

    res = sqrt(res)
  end function device_glsubnorm

  !> Sum a vector of length n
  function device_glsum(a_d, n, strm) result(res)
    type(c_ptr) :: a_d
    integer :: n, ierr
    real(kind=rp) :: res
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

    res = 0.0_rp
#if HAVE_HIP
    res = hip_glsum(a_d, n, strm_)
#elif HAVE_CUDA
    res = cuda_glsum(a_d, n, strm_)
#elif HAVE_OPENCL
    res = opencl_glsum(a_d, n, strm_)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif
  end function device_glsum

  subroutine device_absval(a_d, n, strm)
    integer, intent(in) :: n
    type(c_ptr) :: a_d
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#ifdef HAVE_HIP
    call hip_absval(a_d, n, strm_)
#elif HAVE_CUDA
    call cuda_absval(a_d, n, strm_)
#elif HAVE_OPENCL
    call neko_error('OPENCL is not implemented for device_absval')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_absval

  ! ========================================================================== !
  ! Device point-wise max

  !> Compute the point-wise maximum of two vectors
  !! \f$ a_i = \max(a_i, b_i) \f$
  subroutine device_pwmax2(a_d, b_d, n, strm)
    type(c_ptr) :: a_d, b_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_pwmax_vec2(a_d, b_d, n, strm_)
#elif HAVE_CUDA
    call cuda_pwmax_vec2(a_d, b_d, n, strm_)
#elif HAVE_OPENCL
    call neko_error('No OpenCL backend for device_pwmax2')
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_pwmax2

  !> Compute the point-wise maximum of two vectors
  !! \f$ a_i = \max(b_i, c_i) \f$
  subroutine device_pwmax3(a_d, b_d, c_d, n, strm)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_pwmax_vec3(a_d, b_d, c_d, n, strm_)
#elif HAVE_CUDA
    call cuda_pwmax_vec3(a_d, b_d, c_d, n, strm_)
#elif HAVE_OPENCL
    call neko_error('No OpenCL backend for device_pwmax3')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_pwmax3

  !> Compute the point-wise maximum of a vector and a scalar
  !! \f$ a_i = \max(a_i, c) \f$
  subroutine device_cpwmax2(a_d, c, n, strm)
    type(c_ptr) :: a_d
    real(kind=rp), intent(in) :: c
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_pwmax_sca2(a_d, c, n, strm_)
#elif HAVE_CUDA
    call cuda_pwmax_sca2(a_d, c, n, strm_)
#elif HAVE_OPENCL
    call neko_error('No OpenCL backend for device_cpwmax2')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_cpwmax2

  !> Compute the point-wise maximum of a vector and a scalar
  !! \f$ a_i = \max(b_i, c) \f$
  subroutine device_cpwmax3(a_d, b_d, c, n, strm)
    type(c_ptr) :: a_d, b_d
    real(kind=rp), intent(in) :: c
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_pwmax_sca3(a_d, b_d, c, n, strm_)
#elif HAVE_CUDA
    call cuda_pwmax_sca3(a_d, b_d, c, n, strm_)
#elif HAVE_OPENCL
    call neko_error('No OpenCL backend for device_cpwmax3')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_cpwmax3

  ! ========================================================================== !
  ! Device point-wise min

  !> Compute the point-wise minimum of two vectors
  !! \f$ a_i = \min(a_i, b_i) \f$
  subroutine device_pwmin2(a_d, b_d, n, strm)
    type(c_ptr) :: a_d, b_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_pwmin_vec2(a_d, b_d, n, strm_)
#elif HAVE_CUDA
    call cuda_pwmin_vec2(a_d, b_d, n, strm_)
#elif HAVE_OPENCL
    call neko_error('No OpenCL backend for device_pwmin2')
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_pwmin2

  !> Compute the point-wise minimum of two vectors
  !! \f$ a_i = \min(b_i, c_i) \f$
  subroutine device_pwmin3(a_d, b_d, c_d, n, strm)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_pwmin_vec3(a_d, b_d, c_d, n, strm_)
#elif HAVE_CUDA
    call cuda_pwmin_vec3(a_d, b_d, c_d, n, strm_)
#elif HAVE_OPENCL
    call neko_error('No OpenCL backend for device_pwmin3')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_pwmin3

  !> Compute the point-wise minimum of a vector and a scalar
  !! \f$ a_i = \min(a_i, c) \f$
  subroutine device_cpwmin2(a_d, c, n, strm)
    type(c_ptr) :: a_d
    real(kind=rp), intent(in) :: c
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_pwmin_sca2(a_d, c, n, strm_)
#elif HAVE_CUDA
    call cuda_pwmin_sca2(a_d, c, n, strm_)
#elif HAVE_OPENCL
    call neko_error('No OpenCL backend for device_cpwmin2')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_cpwmin2

  !> Compute the point-wise minimum of a vector and a scalar
  !! \f$ a_i = \min(b_i, c) \f$
  subroutine device_cpwmin3(a_d, b_d, c, n, strm)
    type(c_ptr) :: a_d, b_d
    real(kind=rp), intent(in) :: c
    integer :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_pwmin_sca3(a_d, b_d, c, n, strm_)
#elif HAVE_CUDA
    call cuda_pwmin_sca3(a_d, b_d, c, n, strm_)
#elif HAVE_OPENCL
    call neko_error('No OpenCL backend for device_cpwmin3')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_cpwmin3

  ! ========================================================================== !
  ! Integer operations

  !> Add an integer scalar to vector \f$ a = a + s \f$
  subroutine device_iadd(a_d, c, n, strm)
    type(c_ptr), intent(inout) :: a_d
    integer, intent(in) :: c
    integer, intent(in) :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_
    if (n .lt. 1) return

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

#if HAVE_HIP
    call hip_iadd(a_d, c, n, strm_)
#elif HAVE_CUDA
    call cuda_iadd(a_d, c, n, strm_)
#elif HAVE_OPENCL
    call opencl_iadd(a_d, c, n, strm_)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_iadd


end module device_math
