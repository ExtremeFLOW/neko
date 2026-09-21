!> @file adjoint_pnpn_res.f90
!! @copyright
!! Copyright (c) 2024-2026, The Neko-TOP Authors
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
module adjoint_pnpn_residual
  use gather_scatter, only : gs_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use facet_normal, only : facet_normal_t
  use space, only : space_t
  use mesh, only : mesh_t
  use num_types, only : rp
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Abstract type to compute pressure residual
  type, public, abstract :: adjoint_pnpn_prs_res_t
   contains
     procedure(adjoint_prs_res), nopass, deferred :: compute
  end type adjoint_pnpn_prs_res_t

  !> Abstract type to compute velocity residual
  type, public, abstract :: adjoint_pnpn_vel_res_t
   contains
     procedure(adjoint_vel_res), nopass, deferred :: compute
  end type adjoint_pnpn_vel_res_t

  abstract interface
     !> Compute adjoint pressure residual.
     !! @param p Adjoint pressure field.
     !! @param p_res Pressure residual output.
     !! @param u Adjoint velocity x-component.
     !! @param v Adjoint velocity y-component.
     !! @param w Adjoint velocity z-component.
     !! @param f_x Explicit forcing x-component.
     !! @param f_y Explicit forcing y-component.
     !! @param f_z Explicit forcing z-component.
     !! @param c_xh Coefficients on the pressure space.
     !! @param gs_Xh Gather-scatter operator on the pressure space.
     !! @param bc_prs_surface Pressure boundary surface normals.
     !! @param bc_sym_surface Symmetry boundary surface normals.
     !! @param Ax Helmholtz operator.
     !! @param bd BDF coefficient for the current step.
     !! @param dt Time-step size.
     !! @param mu Dynamic viscosity field.
     !! @param rho Density field.
     !! @param event Backend event handle (optional for device backends).
     subroutine adjoint_prs_res(p, p_res, u, v, w, f_x, f_y, f_z, c_xh,&
          gs_Xh, bc_prs_surface, bc_sym_surface, Ax, bd, dt, mu, rho, event)
       import field_t
       import Ax_t
       import gs_t
       import facet_normal_t
       import coef_t
       import rp
       import c_ptr
       type(field_t), intent(inout) :: p !< Adjoint pressure field.
       type(field_t), intent(inout) :: u !< Adjoint velocity x-component.
       type(field_t), intent(inout) :: v !< Adjoint velocity y-component.
       type(field_t), intent(inout) :: w !< Adjoint velocity z-component.
       type(field_t), intent(inout) :: p_res !< Pressure residual output.
       type(field_t), intent(in) :: f_x !< Explicit forcing x-component.
       type(field_t), intent(in) :: f_y !< Explicit forcing y-component.
       type(field_t), intent(in) :: f_z !< Explicit forcing z-component.
       type(coef_t), intent(inout) :: c_Xh
       !< Coefficients on the pressure space.
       type(gs_t), intent(inout) :: gs_Xh
       !< Gather-scatter operator on the pressure space.
       type(facet_normal_t), intent(in) :: bc_prs_surface
       !< Pressure boundary surface normals.
       type(facet_normal_t), intent(in) :: bc_sym_surface
       !< Symmetry boundary surface normals.
       class(Ax_t), intent(inout) :: Ax !< Helmholtz operator.
       real(kind=rp), intent(in) :: bd !< BDF coefficient for the current step.
       real(kind=rp), intent(in) :: dt !< Time-step size.
       type(field_t), intent(in) :: mu !< Dynamic viscosity field.
       type(field_t), intent(in) :: rho !< Density field.
       type(c_ptr), intent(inout) :: event !< Backend event handle.
     end subroutine adjoint_prs_res
  end interface

  abstract interface
     !> Compute adjoint velocity residual.
     !! @param Ax Helmholtz operator.
     !! @param u Adjoint velocity x-component.
     !! @param v Adjoint velocity y-component.
     !! @param w Adjoint velocity z-component.
     !! @param u_res Residual for adjoint velocity x-component.
     !! @param v_res Residual for adjoint velocity y-component.
     !! @param w_res Residual for adjoint velocity z-component.
     !! @param p Adjoint pressure field.
     !! @param f_x Explicit forcing x-component.
     !! @param f_y Explicit forcing y-component.
     !! @param f_z Explicit forcing z-component.
     !! @param c_Xh Coefficients on the velocity space.
     !! @param msh Mesh object.
     !! @param Xh Velocity space.
     !! @param mu Dynamic viscosity field.
     !! @param rho Density field.
     !! @param bd BDF coefficient for the current step.
     !! @param dt Time-step size.
     !! @param n Number of degrees of freedom.
     subroutine adjoint_vel_res(Ax, u, v, w, u_res, v_res, w_res, &
          p, f_x, f_y, f_z, c_Xh, msh, Xh, mu, rho, bd, dt, n)
       import field_t
       import Ax_t
       import space_t
       import coef_t
       import mesh_t
       import rp
       class(ax_t), intent(in) :: Ax !< Helmholtz operator.
       type(mesh_t), intent(inout) :: msh !< Mesh object.
       type(space_t), intent(inout) :: Xh !< Velocity space.
       type(field_t), intent(inout) :: p !< Adjoint pressure field.
       type(field_t), intent(inout) :: u !< Adjoint velocity x-component.
       type(field_t), intent(inout) :: v !< Adjoint velocity y-component.
       type(field_t), intent(inout) :: w !< Adjoint velocity z-component.
       type(field_t), intent(inout) :: u_res
       !< Residual for adjoint velocity x-component.
       type(field_t), intent(inout) :: v_res
       !< Residual for adjoint velocity y-component.
       type(field_t), intent(inout) :: w_res
       !< Residual for adjoint velocity z-component.
       type(field_t), intent(in) :: f_x !< Explicit forcing x-component.
       type(field_t), intent(in) :: f_y !< Explicit forcing y-component.
       type(field_t), intent(in) :: f_z !< Explicit forcing z-component.
       type(coef_t), intent(inout) :: c_Xh
       !< Coefficients on the velocity space.
       type(field_t), intent(in) :: mu !< Dynamic viscosity field.
       type(field_t), intent(in) :: rho !< Density field.
       real(kind=rp), intent(in) :: bd !< BDF coefficient for the current step.
       real(kind=rp), intent(in) :: dt !< Time-step size.
       integer, intent(in) :: n !< Number of degrees of freedom.
     end subroutine adjoint_vel_res

  end interface

  interface

     !> Factory for the pressure residual computation routine for the PnPn fluid
     !! scheme with the constant-viscosity stress formulation.
     !! @details Only selects the compute backend.
     !! @param object The object to be allocated by the factory.
     module subroutine adjoint_pnpn_prs_res_factory(object)
       class(adjoint_pnpn_prs_res_t), allocatable, intent(inout) :: object
     end subroutine adjoint_pnpn_prs_res_factory

     !> Factory for the velocity residual computation routine for the PnPn fluid
     !! scheme with the constant-viscosity stress formulation.
     !! @details Only selects the compute backend.
     !! @param object The object to be allocated by the factory.
     module subroutine adjoint_pnpn_vel_res_factory(object)
       class(adjoint_pnpn_vel_res_t), allocatable, intent(inout) :: object
     end subroutine adjoint_pnpn_vel_res_factory

  end interface

  public :: adjoint_pnpn_prs_res_factory, adjoint_pnpn_vel_res_factory

end module adjoint_pnpn_residual
