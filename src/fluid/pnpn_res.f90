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
!> Defines Pressure and velocity residuals in the Pn-Pn formulation
module pnpn_residual
  use gather_scatter, only : gs_t  
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use facet_normal, only : facet_normal_t
  use space, only : space_t
  use mesh, only : mesh_t
  use num_types, only : rp
  use scratch_registry, only : scratch_registry_t
  implicit none
  private
  
  !> Abstract type to compute pressure residual
  type, public, abstract :: pnpn_prs_res_t
   contains
     procedure(prs_res), nopass, deferred :: compute     
  end type pnpn_prs_res_t

  !> Abstract type to compute velocity residual
  type, public, abstract :: pnpn_vel_res_t
   contains
     procedure(vel_res), nopass, deferred :: compute
  end type pnpn_vel_res_t
    
  abstract interface
     subroutine prs_res(p, p_res, u, v, w, u_e, v_e, w_e, f_x, f_y, f_z, c_xh,&
          gs_Xh, bc_prs_surface, bc_sym_surface, Ax, bd, dt, mu, rho)
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
       real(kind=rp), intent(in) :: mu
       real(kind=rp), intent(in) :: rho
     end subroutine prs_res
  end interface

  abstract interface
     subroutine vel_res(Ax, u, v, w, u_res, v_res, w_res, &
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
       real(kind=rp), intent(in) :: mu
       real(kind=rp), intent(in) :: rho
       real(kind=rp), intent(in) :: bd
       real(kind=rp), intent(in) :: dt
       integer, intent(in) :: n
     end subroutine vel_res

  end interface
 
end module pnpn_residual
