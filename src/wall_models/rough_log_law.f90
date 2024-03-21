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
!> Implements `rough_log_law_t`.
module rough_log_law
  use num_types, only : rp
  use json_module, only : json_file
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use wall_model, only : wall_model_t
  implicit none
  private

  !> Wall model based on the log-law for a rough wall.
  type, public, extends(wall_model_t) :: rough_log_law_t

    !> The von Karman coefficient.
    real(kind=rp) :: kappa = 0.41
   contains
     !> Constructor.
     procedure, pass(this) :: init => rough_log_law_init
     !> Destructor.
     procedure, pass(this) :: free => rough_log_law_free
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => rough_log_law_compute
  end type rough_log_law_t

contains
  !> Constructor
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  !! @param nu_name The name of the turbulent viscosity field.
  !! @param json A dictionary with parameters.
  subroutine rough_log_law_init(this, dofmap, coef, msk, facet, json)
    class(rough_log_law_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    type(json_file), intent(inout) :: json

    call this%init_base(dofmap, coef, msk, facet)

    ! Make RTS
    this%kappa = 0.41

  end subroutine rough_log_law_init

  !> Destructor for the rough_log_law_t (base) class.
  subroutine rough_log_law_free(this)
    class(rough_log_law_t), intent(inout) :: this

    call this%free_base()

  end subroutine rough_log_law_free

  !> Compute the wall shear stress.
  subroutine rough_log_law_compute(this, t, tstep)
    class(rough_log_law_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    write(*,*) "Computing stress"

    this%tau_x = 1.0_rp
    this%tau_y = 1.0_rp
    this%tau_z = 1.0_rp
  end subroutine rough_log_law_compute


end module rough_log_law