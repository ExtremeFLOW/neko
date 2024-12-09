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
!> Implements `adv_dummy_t`
module adv_dummy
  use advection, only: advection_t
  use num_types, only: rp
  use space, only: space_t
  use field, only: field_t
  use coefs, only: coef_t
  implicit none
  private

  !> A zero-valued advection that can be used to kill the advection term.
  type, public, extends(advection_t) :: adv_dummy_t
   contains
     !> Constructor.
     procedure, pass(this) :: init => init_adv_dummy
     !> Destructor.
     procedure, pass(this) :: free => free_adv_dummy
     !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
     !! the RHS.
     procedure, pass(this) :: compute => compute_adv_dummy
     !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to
     !! the RHS.
     procedure, pass(this) :: compute_scalar => &
       compute_scalar_adv_dummy
  end type adv_dummy_t

contains

  !> Constructor.
  !! @param coef The coefficients of the (space, mesh) pair.
  subroutine init_adv_dummy(this, coef)
    class(adv_dummy_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

  end subroutine init_adv_dummy

  !> Destructor.
  subroutine free_adv_dummy(this)
    class(adv_dummy_t), intent(inout) :: this

  end subroutine free_adv_dummy

  !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$ to the
  !! RHS.
  !! @details Here, does nothing.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param fx The x component of source term.
  !! @param fy The y component of source term.
  !! @param fz The z component of source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step, not required for this method.
  subroutine compute_adv_dummy(this, vx, vy, vz, fx, fy, fz, Xh, &
                                          coef, n, dt)
    class(adv_dummy_t), intent(inout) :: this
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: fx, fy, fz
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt

  end subroutine compute_adv_dummy

  !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to the
  !! RHS.
  !! @details Here, does nothing.
  !! @param this The object.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param s The scalar.
  !! @param fs The source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step, not required for this method.
  subroutine compute_scalar_adv_dummy(this, vx, vy, vz, s, fs, Xh, &
                                                 coef, n, dt)
    class(adv_dummy_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: s
    type(field_t), intent(inout) :: fs
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt

  end subroutine compute_scalar_adv_dummy

end module adv_dummy
