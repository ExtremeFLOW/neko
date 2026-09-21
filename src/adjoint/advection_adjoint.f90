!> @file advection_adjoint.f90
!! @copyright
!! Copyright (c) 2024-2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> Subroutines to add advection terms to the RHS of a transport equation.
module advection_adjoint
  use num_types, only: rp
  use space, only: space_t
  use field, only: field_t
  use coefs, only: coef_t
  use json_module, only : json_file
  use field_series, only: field_series_t
  use time_scheme_controller, only: time_scheme_controller_t
  implicit none
  private

  !> Base abstract type for computing the advection operator
  type, abstract :: advection_adjoint_t
   contains
     procedure(compute_adv_lin), pass(this), deferred :: compute_linear
     procedure(compute_adv_lin), pass(this), deferred :: compute_adjoint
     procedure(compute_scalar_adv_lin), pass(this), deferred :: &
          compute_adjoint_scalar
     procedure(advection_adjoint_free), pass(this), deferred :: free
  end type advection_adjoint_t

  interface
     !> A factory for \ref advection_adjoint_t descendants. Both creates and
     !! initializes the object.
     !! @param object The object allocated by the factory.
     !! @param json The parameter file.
     !! @param coef The coefficients of the (space, mesh) pair.
     !! @param ulag, vlag, wlag The lagged velocity fields.
     !! @param dtlag The lagged time steps.
     !! @param tlag The lagged times.
     !! @param time_scheme The bdf-ext time scheme used in the method.
     !! @param use_dummy If true, a dummy zero-valued advection type is
     !! @param slag The lagged scalar field.
     !! allocated. This can be used to kill the advection term.
     !! @note The factory both allocates and initializes `object`.
     module subroutine advection_adjoint_factory(object, json, coef, &
          ulag, vlag, wlag, dtlag, tlag, time_scheme, use_dummy, slag)
       class(advection_adjoint_t), allocatable, intent(inout) :: object
       type(json_file), intent(inout) :: json
       type(coef_t), intent(inout), target :: coef
       type(field_series_t), intent(in), target :: ulag, vlag, wlag
       real(kind=rp), intent(in), target :: dtlag(10)
       real(kind=rp), intent(in), target :: tlag(10)
       type(time_scheme_controller_t), intent(in), target :: time_scheme
       logical, optional, intent(in) :: use_dummy
       type(field_series_t), target, optional, intent(in) :: slag
     end subroutine advection_adjoint_factory
  end interface

  ! ========================================================================== !
  ! Linearized advection operator interface

  abstract interface
     !> Add advection operator to the right-hand-side for a fluld.
     !! @param this The object.
     !! @param vx The x component of velocity.
     !! @param vy The y component of velocity.
     !! @param vz The z component of velocity.
     !! @param fx The x component of source term.
     !! @param fy The y component of source term.
     !! @param fz The z component of source term.
     !! @param Xh The function space.
     !! @param coef The coefficients of the (Xh, mesh) pair.
     !! @param n Typically the size of the mesh.
     subroutine compute_adv_lin(this, vx, vy, vz, vxb, vyb, vzb, fx, fy, fz, &
          Xh, coef, n)
       import :: advection_adjoint_t
       import :: coef_t
       import :: space_t
       import :: field_t
       import :: rp
       class(advection_adjoint_t), intent(inout) :: this
       type(space_t), intent(inout) :: Xh
       type(coef_t), intent(inout) :: coef
       type(field_t), intent(inout) :: vx, vy, vz
       type(field_t), intent(inout) :: vxb, vyb, vzb
       type(field_t), intent(inout) :: fx, fy, fz
       integer, intent(in) :: n
     end subroutine compute_adv_lin
  end interface

  abstract interface
     !> Add advection operator to the right-hand-side for a scalar.
     !! @param this The object.
     !! @param vxb The x component of velocity.
     !! @param vyb The y component of velocity.
     !! @param vzb The z component of velocity.
     !! @param s The scalar.
     !! @param fs The source term.
     !! @param Xh The function space.
     !! @param coef The coefficients of the (Xh, mesh) pair.
     !! @param n Typically the size of the mesh.
     subroutine compute_scalar_adv_lin(this, vxb, vyb, vzb, s, fs, Xh, coef, &
          n, dt)
       import :: advection_adjoint_t
       import :: coef_t
       import :: space_t
       import :: field_t
       import :: rp
       class(advection_adjoint_t), intent(inout) :: this
       type(field_t), intent(inout) :: vxb, vyb, vzb
       type(field_t), intent(inout) :: s
       type(field_t), intent(inout) :: fs
       type(space_t), intent(inout) :: Xh
       type(coef_t), intent(inout) :: coef
       real(kind=rp), intent(in), optional :: dt
       integer, intent(in) :: n
     end subroutine compute_scalar_adv_lin
  end interface

  abstract interface
     !> Destructor
     subroutine advection_adjoint_free(this)
       import :: advection_adjoint_t
       class(advection_adjoint_t), intent(inout) :: this
     end subroutine advection_adjoint_free
  end interface

  public :: advection_adjoint_t, advection_adjoint_factory
end module advection_adjoint
