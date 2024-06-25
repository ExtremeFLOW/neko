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
!
!> Implements `gradient_jump_penalty_t`.
module gradient_jump_penalty
  use num_types, only : rp
  use math
  use field, only : field_t
  use dofmap , only : dofmap_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use coefs, only : coef_t

  implicit none
  private

  !> Implements the gradient jump penalty.
  !! @note Reference DOI: 10.1016/j.cma.2021.114200
  type, public :: gradient_jump_penalty_t
     !> Coefficient of the penalty.
     real :: tau
     !> Polynomial order
     integer :: p
     !> Penalty terms
     real(kind=rp), allocatable :: penalty(:, :, :, :)
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()


  contains
     !> Constructor.
     procedure, pass(this) :: init => gradient_jump_penalty_init
     !> Destructor
     procedure, pass(this) :: free => gradient_jump_penalty_free
    !  !> Compute gradient jump penalty term.
    !  procedure, pass(this) :: compute => gradient_jump_penalty_compute

  end type gradient_jump_penalty_t

contains
  !> Constructor.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  subroutine gradient_jump_penalty_init(this, dofmap, coef)
    implicit none
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), target, intent(in) :: coef
    
    call this%free()

    this%p = dofmap%xh%lx - 1
    if (this%p .gt. 1) then
       this%tau = 0.8_rp * (this%p + 1) ** (-4)
    else
       this%tau = 0.02_rp
    end if
    this%coef => coef

    allocate(this%penalty(this%p + 1, this%p + 1 , this%p + 1 , this%coef%msh%nelv))


  end subroutine gradient_jump_penalty_init

  !> Destructor for the gradient_jump_penalty_t class.
  subroutine gradient_jump_penalty_free(this)
    implicit none
    class(gradient_jump_penalty_t), intent(inout) :: this
    
    if (allocated(this%penalty)) then
       deallocate(this%penalty)
    end if
    nullify(this%coef)

  end subroutine gradient_jump_penalty_free

!   !> Compute eddy viscosity.
!   !! @param t The time value.
!   !! @param tstep The current time-step.
!   subroutine gradient_jump_penalty_compute(this, t, tstep)
!     class(grad_jump_penalty_t), intent(inout) :: this
!     real(kind=rp), intent(in) :: t
!     integer, intent(in) :: tstep

!     call set_ds_filt(this%test_filter)

!     if (NEKO_BCKND_DEVICE .eq. 1) then
!         call neko_error("Dynamic Smagorinsky model not implemented on accelarators.")
!     else
!         call dynamic_smagorinsky_compute_cpu(t, tstep, this%coef, this%nut, &
!                                 this%delta, this%c_dyn, this%test_filter, &
!                                 this%mij, this%lij, this%num, this%den)
!     end if

!   end subroutine gradient_jump_penalty_compute

end module gradient_jump_penalty