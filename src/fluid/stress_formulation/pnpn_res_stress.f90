! Copyright (c) 2022, The Neko Authors
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
!> Defines Pressure and velocity residuals in the Pn-Pn formulation an full
!! viscous stress term formulation.
module pnpn_residual_stress
  use gather_scatter, only : gs_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use facet_normal, only : facet_normal_t
  use space, only : space_t
  use mesh, only : mesh_t
  use num_types, only : rp
  implicit none
  private

  !> Abstract type to compute the pressure residual for the PnPn fluid with
  !! full viscous stress formulation.
  !! @details Descendants correspond to implementations for different compute
  !! backends.
  type, public, abstract :: pnpn_prs_res_stress_t
   contains
     procedure(prs_res_stress), nopass, deferred :: compute
  end type pnpn_prs_res_stress_t

  !> Abstract type to compute the velocity residual for the PnPn fluid with
  !! full viscous stress formulation.
  !! @details Descendants correspond to implementations for different compute
  !! backends.
  type, public, abstract :: pnpn_vel_res_stress_t
   contains
     procedure(vel_res_stress), nopass, deferred :: compute
  end type pnpn_vel_res_stress_t

  abstract interface
     !> Compute the residual of the pressure equation.
     !! @param p The pressure field.
     !! @param p_res The output residual field.
     !! @param u The x component of velocity.
     !! @param v The y component of velocity.
     !! @param w The z component of velocity.
     !! @param u_e The x component of the explicitly time-extrapolated velocity.
     !! @param v_e The y component of the explicitly time-extrapolated velocity.
     !! @param w_e The z component of the explicitly time-extrapolated velocity.
     !! @param f_x The x component of the right-hand side.
     !! @param f_y The y component of the right-hand side.
     !! @param f_z The z component of the right-hand side.
     !! @param c_Xh The SEM coefficients.
     !! @param gs_Xh The gather-scatter.
     !! @param bcs_prs_surface Pressure boundary condition.
     !! @param bcs_sym_surface Symmetry boundary conditions.
     !! @param Ax The implicit laplacian operator kernel.
     !! @param bd The first coefficient of the BDF scheme.
     !! @param dt The time step.
     !! @param mu The dynamic viscosity.
     !! @param rho The density.
     subroutine prs_res_stress(p, p_res, u, v, w, u_e, v_e, w_e, f_x, f_y, f_z,&
       c_Xh, gs_Xh, bc_prs_surface, bc_sym_surface, Ax, bd, dt, mu, rho)
       import field_t
       import Ax_t
       import gs_t
       import facet_normal_t
       import coef_t
       import rp
       type(field_t), intent(inout) :: p, u, v, w
       type(field_t), intent(inout) :: u_e, v_e, w_e !< time-extrapolated velocity
       type(field_t), intent(inout) :: p_res
       !> Momentum source terms
       type(field_t), intent(inout) :: f_x, f_y, f_z
       type(coef_t), intent(inout) :: c_Xh
       type(gs_t), intent(inout) :: gs_Xh
       type(facet_normal_t), intent(inout) :: bc_prs_surface
       type(facet_normal_t), intent(inout) :: bc_sym_surface
       class(Ax_t), intent(inout) :: Ax
       real(kind=rp), intent(inout) :: bd
       real(kind=rp), intent(in) :: dt
       type(field_t), intent(in) :: mu
       type(field_t), intent(in) :: rho
     end subroutine prs_res_stress
  end interface

  abstract interface
     !> Compute the residual of the velocity equation.
     !! @param Ax The implicit viscous stress operator kernel.
     !! @param u The x component of velocity.
     !! @param v The y component of velocity.
     !! @param w The z component of velocity.
     !! @param u_res The x component of the output residual field.
     !! @param v_res The y component of the output residual field.
     !! @param w_res The z component of the output residual field.
     !! @param f_x The x component of the right-hand side.
     !! @param f_y The y component of the right-hand side.
     !! @param f_z The z component of the right-hand side.
     !! @param c_Xh The SEM coefficients.
     !! @param msh The mesh.
     !! @param Xh The SEM space.
     !! @param mu The dynamic viscosity.
     !! @param rho The density.
     !! @param bd The first coefficient of the BDF scheme.
     !! @param dt The time step.
     subroutine vel_res_stress(Ax, u, v, w, u_res, v_res, w_res, &
          p, f_x, f_y, f_z, c_Xh, msh, Xh, mu, rho, bd, dt, n)
       import field_t
       import Ax_t
       import gs_t
       import facet_normal_t
       import space_t
       import coef_t
       import mesh_t
       import rp
       class(ax_t), intent(in) :: Ax
       type(mesh_t), intent(inout) :: msh
       type(space_t), intent(inout) :: Xh
       type(field_t), intent(inout) :: p, u, v, w
       type(field_t), intent(inout) :: u_res, v_res, w_res
       type(field_t), intent(inout) :: f_x, f_y, f_z
       type(coef_t), intent(inout) :: c_Xh
       type(field_t), intent(in) :: mu
       type(field_t), intent(in) :: rho
       real(kind=rp), intent(in) :: bd
       real(kind=rp), intent(in) :: dt
       integer, intent(in) :: n
     end subroutine vel_res_stress

  end interface

end module pnpn_residual_stress
