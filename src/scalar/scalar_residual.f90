! Copyright (c) 2022-2024, The Neko Authors
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
!> Defines the residual for the scalar transport equation
module scalar_residual
  use gather_scatter, only : gs_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use source_scalar, only : source_scalar_t
  use facet_normal, only : facet_normal_t
  use space, only : space_t
  use mesh, only : mesh_t
  use num_types, only : rp
  implicit none
  private

  !> Abstract type to compute scalar residual
  type, public, abstract :: scalar_residual_t
   contains
     procedure(scalar_residual_interface), nopass, deferred :: compute
  end type scalar_residual_t

  abstract interface
     !> Interface for computing the residual of a scalar transport equation.
     !! @param Ax The Helmholtz operator.
     !! @param s The values of the scalar.
     !! @param s_res The values of the scalar residual.
     !! @param f_xH The right hand side.
     !! @param c_xH The SEM coefficients.
     !! @param msh The mesh.
     !! @param Xh The SEM function space.
     !! @param lambda The thermal conductivity.
     !! @param rhocp The density multiplied by the specific heat capacity.
     !! @param bd The coefficeints from the BDF differencing scheme.
     !! @param dt The timestep.
     !! @param n The total number of degrees of freedom.
     subroutine scalar_residual_interface(Ax, s, s_res, f_Xh, c_Xh, msh, Xh, &
                                          lambda, rhocp, bd, dt, n)
       import field_t
       import Ax_t
       import gs_t
       import facet_normal_t
       import source_scalar_t
       import space_t
       import coef_t
       import mesh_t
       import rp
       class(ax_t), intent(in) :: Ax
       type(mesh_t), intent(inout) :: msh
       type(space_t), intent(inout) :: Xh
       type(field_t), intent(inout) :: s
       type(field_t), intent(inout) :: s_res
       type(field_t), intent(inout) :: f_Xh
       type(coef_t), intent(inout) :: c_Xh
       type(field_t), intent(in) :: lambda
       real(kind=rp), intent(in) :: rhocp
       real(kind=rp), intent(in) :: bd
       real(kind=rp), intent(in) :: dt
       integer, intent(in) :: n
     end subroutine scalar_residual_interface
  end interface

  interface
     !> Factory for the scalar advection-diffusion residual.
     !! @details Only selects the compute backend.
     !! @param object The object to be allocated by the factory.
     module subroutine scalar_residual_factory(object)
       class(scalar_residual_t), allocatable, intent(inout) :: object
     end subroutine scalar_residual_factory
  end interface

  public :: scalar_residual_factory
  
end module scalar_residual
