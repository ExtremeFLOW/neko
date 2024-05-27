! Copyright (c) 2021-2024, The Neko Authors
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
!> Subroutines to add advection terms to the RHS of a transport equation.
module advection
  use num_types, only : rp
  use space, only : space_t
  use field, only : field_t
  use coefs, only : coef_t
  implicit none
  private

  !> Base abstract type for computing the advection operator
  type, public, abstract :: advection_t
   contains
     procedure(compute_adv), pass(this), deferred :: compute
     procedure(compute_scalar_adv), pass(this), deferred :: compute_scalar
     procedure(advection_free), pass(this), deferred :: free
  end type advection_t

  !> Base abstract type for computing the advection operator
  type, public, abstract :: advection_lin_t
   contains
     procedure(compute_adv_lin), pass(this), deferred :: compute_linear
     procedure(compute_adv_lin), pass(this), deferred :: compute_adjoint
! TODO
!     procedure(compute_scalar_adv_lin), pass(this), deferred :: compute_scalar
     procedure(advection_lin_free), pass(this), deferred :: free
  end type advection_lin_t

  ! ========================================================================== !
  ! Advection operator interface

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
     subroutine compute_adv(this, vx, vy, vz, fx, fy, fz, Xh, coef, n)
       import :: advection_t
       import :: coef_t
       import :: space_t
       import :: field_t
       import :: rp
       class(advection_t), intent(inout) :: this
       type(space_t), intent(inout) :: Xh
       type(coef_t), intent(inout) :: coef
       type(field_t), intent(inout) :: vx, vy, vz
       type(field_t), intent(inout) :: fx, fy, fz
       integer, intent(in) :: n
     end subroutine compute_adv
  end interface

  abstract interface
     !> Add advection operator to the right-hand-side for a scalar.
     !! @param this The object.
     !! @param vx The x component of velocity.
     !! @param vy The y component of velocity.
     !! @param vz The z component of velocity.
     !! @param s The scalar.
     !! @param fs The source term.
     !! @param Xh The function space.
     !! @param coef The coefficients of the (Xh, mesh) pair.
     !! @param n Typically the size of the mesh.
     subroutine compute_scalar_adv(this, vx, vy, vz, s, fs, Xh, coef, n)
       import :: advection_t
       import :: coef_t
       import :: space_t
       import :: field_t
       import :: rp
       class(advection_t), intent(inout) :: this
       type(field_t), intent(inout) :: vx, vy, vz
       type(field_t), intent(inout) :: s
       type(field_t), intent(inout) :: fs
       type(space_t), intent(inout) :: Xh
       type(coef_t), intent(inout) :: coef
       integer, intent(in) :: n
     end subroutine compute_scalar_adv
  end interface

  abstract interface
     !> Destructor
     subroutine advection_free(this)
       import :: advection_t
       class(advection_t), intent(inout) :: this
     end subroutine advection_free
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
       import :: advection_lin_t
       import :: coef_t
       import :: space_t
       import :: field_t
       import :: rp
       class(advection_lin_t), intent(inout) :: this
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
     !! @param vx The x component of velocity.
     !! @param vy The y component of velocity.
     !! @param vz The z component of velocity.
     !! @param s The scalar.
     !! @param fs The source term.
     !! @param Xh The function space.
     !! @param coef The coefficients of the (Xh, mesh) pair.
     !! @param n Typically the size of the mesh.
     subroutine compute_scalar_adv_lin(this, vx, vy, vz, s, fs, Xh, coef, n)
       import :: advection_lin_t
       import :: coef_t
       import :: space_t
       import :: field_t
       import :: rp
       class(advection_lin_t), intent(inout) :: this
       type(field_t), intent(inout) :: vx, vy, vz
       type(field_t), intent(inout) :: s
       type(field_t), intent(inout) :: fs
       type(space_t), intent(inout) :: Xh
       type(coef_t), intent(inout) :: coef
       integer, intent(in) :: n
     end subroutine compute_scalar_adv_lin
  end interface

  abstract interface
     !> Destructor
     subroutine advection_lin_free(this)
       import :: advection_lin_t
       class(advection_lin_t), intent(inout) :: this
     end subroutine advection_lin_free
  end interface

end module advection
