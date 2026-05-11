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
!> Residual weighting and normalization policies for KSP solvers.
module ksp_residual
  use num_types, only : rp
  use field, only : field_t
  use coefs, only : coef_t
  use utils, only : neko_error, neko_type_error
  implicit none
  private

  character(len=16), parameter :: KSP_RESIDUAL_WEIGHT_TYPES(2) = &
       [character(len=16) :: 'unit', 'coef_mult']
  character(len=16), parameter :: KSP_RESIDUAL_NORMALIZATION_TYPES(4) = &
       [character(len=16) :: 'none', 'volume', 'initial', 'rhs']

  !> KSP residual policy.
  !!
  !! The policy does not own any field-sized state. Solvers are expected to
  !! provide scratch storage for the weight field and use it in their own norm
  !! reductions.
  type, public :: ksp_residual_t
     character(len=:), allocatable :: weight_type
     character(len=:), allocatable :: normalization_type
     procedure(ksp_residual_weight_intrf), nopass, pointer :: &
          compute_weight => null()
     procedure(ksp_residual_normalization_intrf), nopass, pointer :: &
          compute_normalization => null()
   contains
     procedure, pass(this) :: init => ksp_residual_init
     procedure, pass(this) :: free => ksp_residual_free
  end type ksp_residual_t

  abstract interface
     !> Fill a caller-owned field with per-DOF norm weights.
     subroutine ksp_residual_weight_intrf(weight, coef, n)
       import field_t
       import coef_t
       implicit none
       type(field_t), intent(inout) :: weight
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: n
     end subroutine ksp_residual_weight_intrf

     !> Compute a scalar residual normalization factor.
     function ksp_residual_normalization_intrf(coef, reference_norm) &
          result(scale)
       import coef_t
       import rp
       implicit none
       type(coef_t), intent(in) :: coef
       real(kind=rp), optional, intent(in) :: reference_norm
       real(kind=rp) :: scale
     end function ksp_residual_normalization_intrf
  end interface

contains

  !> Initialize a residual policy from weight and normalization names.
  subroutine ksp_residual_init(this, weight_type, normalization_type)
    class(ksp_residual_t), intent(inout) :: this
    character(len=*), intent(in) :: weight_type
    character(len=*), intent(in) :: normalization_type

    call this%free()

    select case (trim(weight_type))
    case ('unit')
       this%compute_weight => ksp_residual_weight_unit
    case ('coef_mult')
       this%compute_weight => ksp_residual_weight_coef_mult
    case default
       call neko_type_error('KSP residual weight', weight_type, &
            KSP_RESIDUAL_WEIGHT_TYPES)
    end select

    select case (trim(normalization_type))
    case ('none')
       this%compute_normalization => ksp_residual_normalization_none
    case ('volume')
       this%compute_normalization => ksp_residual_normalization_volume
    case ('initial')
       this%compute_normalization => ksp_residual_normalization_reference
    case ('rhs')
       this%compute_normalization => ksp_residual_normalization_reference
    case default
       call neko_type_error('KSP residual normalization', &
            normalization_type, KSP_RESIDUAL_NORMALIZATION_TYPES)
    end select

    this%weight_type = trim(weight_type)
    this%normalization_type = trim(normalization_type)

  end subroutine ksp_residual_init

  !> Release a residual policy.
  subroutine ksp_residual_free(this)
    class(ksp_residual_t), intent(inout) :: this

    this%compute_weight => null()
    this%compute_normalization => null()

    if (allocated(this%weight_type)) then
       deallocate(this%weight_type)
    end if

    if (allocated(this%normalization_type)) then
       deallocate(this%normalization_type)
    end if

  end subroutine ksp_residual_free

  !> Unit weights for plain Euclidean-like residual norms.
  subroutine ksp_residual_weight_unit(weight, coef, n)
    type(field_t), intent(inout) :: weight
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n

    if (.not. allocated(weight%x)) then
      call neko_error('KSP residual weight field is not allocated')
    end if

    weight%x = 1.0_rp

  end subroutine ksp_residual_weight_unit

  !> Use the SEM multiplicity weight from the coefficient structure.
  subroutine ksp_residual_weight_coef_mult(weight, coef, n)
    type(field_t), intent(inout) :: weight
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n

    if (.not. allocated(weight%x)) then
      call neko_error('KSP residual weight field is not allocated')
    end if

    weight%x = coef%mult

  end subroutine ksp_residual_weight_coef_mult

  !> No normalization.
  function ksp_residual_normalization_none(coef, reference_norm) result(scale)
    type(coef_t), intent(in) :: coef
    real(kind=rp), optional, intent(in) :: reference_norm
    real(kind=rp) :: scale

    scale = 1.0_rp

  end function ksp_residual_normalization_none

  !> Volume-based normalization, e.g. sqrt(r^TWr / V).
  function ksp_residual_normalization_volume(coef, reference_norm) &
       result(scale)
    type(coef_t), intent(in) :: coef
    real(kind=rp), optional, intent(in) :: reference_norm
    real(kind=rp) :: scale

    scale = sqrt(coef%volume)

  end function ksp_residual_normalization_volume

  !> Reference-based normalization using an externally supplied scalar.
  function ksp_residual_normalization_reference(coef, reference_norm) &
       result(scale)
    type(coef_t), intent(in) :: coef
    real(kind=rp), optional, intent(in) :: reference_norm
    real(kind=rp) :: scale

    if (.not. present(reference_norm)) then
       call neko_error('KSP residual reference normalization requires a ' // &
            'reference norm')
    end if

    scale = max(reference_norm, tiny(1.0_rp))

  end function ksp_residual_normalization_reference

end module ksp_residual
