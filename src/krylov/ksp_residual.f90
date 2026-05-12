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
  use ax_product, only : ax_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc_list, only : bc_list_t
  use scratch_registry, only : neko_scratch_registry
  use math, only : cfill, col2, copy, glsc3, glsum
  use utils, only : neko_error, neko_type_error
  implicit none
  private

  character(len=16), parameter :: KSP_RESIDUAL_WEIGHT_TYPES(3) = &
       [character(len=16) :: 'none', 'volume', 'inverse_volume']
  character(len=16), parameter :: KSP_RESIDUAL_NORMALIZATION_TYPES(4) = &
       [character(len=16) :: 'none', 'volume', 'initial', 'equation_scale']

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

     !> Compute a scalar residual normalization factor from the solver state.
      function ksp_residual_normalization_intrf(Ax, ref, rhs, coef, gs_h, blst, &
           weight, n) result(scale)
       import ax_t
       import field_t
       import coef_t
       import gs_t
       import bc_list_t
       import rp
       implicit none
       class(ax_t), intent(in) :: Ax
        type(field_t), intent(in) :: ref
       integer, intent(in) :: n
       real(kind=rp), intent(in) :: rhs(n)
       type(coef_t), intent(in) :: coef
       type(gs_t), intent(inout) :: gs_h
       type(bc_list_t), intent(inout) :: blst
       type(field_t), intent(in) :: weight
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
    case ('none')
       this%compute_weight => ksp_residual_weight_none
    case ('volume')
       this%compute_weight => ksp_residual_weight_volume
    case ('inverse_volume')
       this%compute_weight => ksp_residual_weight_inverse_volume
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
       this%compute_normalization => ksp_residual_normalization_initial
    case ('equation_scale')
       this%compute_normalization => ksp_residual_normalization_equation_scale
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

  !> Default assembled-CG residual weight based on multiplicity.
  subroutine ksp_residual_weight_none(weight, coef, n)
    type(field_t), intent(inout) :: weight
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n

    if (.not. allocated(weight%x)) then
      call neko_error('KSP residual weight field is not allocated')
    end if

    call copy(weight%x, coef%mult, n)

  end subroutine ksp_residual_weight_none

  !> Use the SEM mass/volume matrix from the coefficient structure.
  subroutine ksp_residual_weight_volume(weight, coef, n)
    type(field_t), intent(inout) :: weight
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n

    if (.not. allocated(weight%x)) then
      call neko_error('KSP residual weight field is not allocated')
    end if

    call copy(weight%x, coef%B, n)
    call col2(weight%x, coef%mult, n)

  end subroutine ksp_residual_weight_volume

  !> Use the SEM mass/volume matrix from the coefficient structure.
  subroutine ksp_residual_weight_inverse_volume(weight, coef, n)
    type(field_t), intent(inout) :: weight
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n

    if (.not. allocated(weight%x)) then
      call neko_error('KSP residual weight field is not allocated')
    end if

    call copy(weight%x, coef%Binv, n)
    call col2(weight%x, coef%mult, n)

  end subroutine ksp_residual_weight_inverse_volume

  !> No normalization.
  function ksp_residual_normalization_none(Ax, ref, rhs, coef, gs_h, blst, &
       weight, n) result(scale)
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(in) :: ref
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: rhs(n)
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs_h
    type(bc_list_t), intent(inout) :: blst
    type(field_t), intent(in) :: weight
    real(kind=rp) :: scale

    scale = 1.0_rp

  end function ksp_residual_normalization_none

  !> Volume-based normalization, e.g. sqrt(r^TWr / V).
  function ksp_residual_normalization_volume(Ax, ref, rhs, coef, gs_h, blst, &
       weight, n) &
       result(scale)
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(in) :: ref
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: rhs(n)
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs_h
    type(bc_list_t), intent(inout) :: blst
    type(field_t), intent(in) :: weight
    real(kind=rp) :: scale

    scale = sqrt(coef%volume)

  end function ksp_residual_normalization_volume

  !> Normalize using an equation-scale weighted L2 scaling.
  !!
  !! Similar structurally to what OpenFOAM does, but the
  !! resulting scale uses the same weighted discrete L2 measure as the solver
  !! residual itself. For correction solves with rhs = b - A(ref), the
  !! full-system right-hand side b is reconstructed internally.
  function ksp_residual_normalization_equation_scale(Ax, ref, rhs, coef, gs_h, &
       blst, weight, n) result(scale)
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(in) :: ref
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: rhs(n)
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs_h
    type(bc_list_t), intent(inout) :: blst
    type(field_t), intent(in) :: weight
    real(kind=rp) :: scale
    real(kind=rp) :: mean_x, a_ref, a_mean
    integer :: ax_idx, ref_idx, const_idx, i
    type(field_t), pointer :: ax_x, ax_ref, x_const

    call neko_scratch_registry%request_field(ax_x, ax_idx, .false.)
    call neko_scratch_registry%request_field(ax_ref, ref_idx, .false.)
    call neko_scratch_registry%request_field(x_const, const_idx, .false.)

    mean_x = glsum(ref%x, n) / real(n, rp)

    call Ax%compute(ax_x%x, ref%x, coef, ref%msh, ref%Xh)
    call gs_h%op(ax_x%x, n, GS_OP_ADD)
    call blst%apply(ax_x)

    call cfill(x_const%x, mean_x, n)
    call Ax%compute(ax_ref%x, x_const%x, coef, ref%msh, ref%Xh)
    call gs_h%op(ax_ref%x, n, GS_OP_ADD)
    call blst%apply(ax_ref)

    do i = 1, n
       a_ref = ax_x%x(i, 1, 1, 1)
       a_mean = ax_ref%x(i, 1, 1, 1)
       ! This is now Ax - A x_const
       ax_x%x(i, 1, 1, 1) = a_ref - a_mean
       ! This is now b - A x_const
       ax_ref%x(i, 1, 1, 1) = rhs(i) + a_ref - a_mean
    end do

    ! Sum the norms of the two contributions to get the final scale.
    scale = sqrt(glsc3(ax_x%x, ax_x%x, weight%x, n) + &
         glsc3(ax_ref%x, ax_ref%x, weight%x, n))
    scale = max(scale, tiny(1.0_rp))

    call neko_scratch_registry%relinquish_field(const_idx)
    call neko_scratch_registry%relinquish_field(ref_idx)
    call neko_scratch_registry%relinquish_field(ax_idx)

  end function ksp_residual_normalization_equation_scale

  !> Normalize by the weighted norm of the initial residual.
  !! Since we solve for dx with guess 0, this is just the norm of the RHS,
  !! i.e. b - A(ref). So this scales the initial residual to 1,
  !! effectively making the absolute tolerance a relative one.
  function ksp_residual_normalization_initial(Ax, ref, rhs, coef, gs_h, blst, &
       weight, n) result(scale)
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(in) :: ref
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: rhs(n)
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs_h
    type(bc_list_t), intent(inout) :: blst
    type(field_t), intent(in) :: weight
    real(kind=rp) :: scale

    scale = sqrt(glsc3(rhs, rhs, weight%x, n))
    scale = max(scale, tiny(1.0_rp))

  end function ksp_residual_normalization_initial

end module ksp_residual
