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
!> Pressure and velocity residuals for the low-Mach Pn-Pn formulation.
!!
!! Mirrors pnpn_residual but accepts an additional Q_T field (the thermal
!! divergence source derived from the temperature equation). Subsequent
!! commits will populate Q_T's contribution to the pressure equation and
!! introduce the variable-density/viscosity viscous terms. In this initial
!! commit the implementations ignore Q_T and reproduce the pnpn residuals
!! exactly, so enabling the new path does not change numerics.
module lowmach_residual
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

  !> Abstract type to compute the low-Mach pressure residual.
  type, public, abstract :: lowmach_prs_res_t
   contains
     procedure(lm_prs_res), nopass, deferred :: compute
  end type lowmach_prs_res_t

  !> Abstract type to compute the low-Mach velocity residual.
  type, public, abstract :: lowmach_vel_res_t
   contains
     procedure(lm_vel_res), nopass, deferred :: compute
  end type lowmach_vel_res_t

  abstract interface
     subroutine lm_prs_res(p, p_res, u, v, w, u_e, v_e, w_e, f_x, f_y, f_z, &
          c_Xh, gs_Xh, bc_prs_surface, bc_sym_surface, Ax, bd, dt, mu, rho, &
          Q_T, event)
       import field_t
       import Ax_t
       import gs_t
       import facet_normal_t
       import coef_t
       import rp
       import c_ptr
       type(field_t), intent(inout) :: p, u, v, w
       type(field_t), intent(in) :: u_e, v_e, w_e
       type(field_t), intent(inout) :: p_res
       type(field_t), intent(in) :: f_x, f_y, f_z
       type(coef_t), intent(inout) :: c_Xh
       type(gs_t), intent(inout) :: gs_Xh
       type(facet_normal_t), intent(in) :: bc_prs_surface
       type(facet_normal_t), intent(in) :: bc_sym_surface
       class(ax_t), intent(inout) :: Ax
       real(kind=rp), intent(in) :: bd
       real(kind=rp), intent(in) :: dt
       type(field_t), intent(in) :: mu
       type(field_t), intent(in) :: rho
       type(field_t), intent(in) :: Q_T
       type(c_ptr), intent(inout) :: event
     end subroutine lm_prs_res
  end interface

  abstract interface
     subroutine lm_vel_res(Ax, u, v, w, u_res, v_res, w_res, &
          p, f_x, f_y, f_z, c_Xh, msh, Xh, mu, rho, Q_T, bd, dt, n)
       import field_t
       import Ax_t
       import space_t
       import coef_t
       import mesh_t
       import rp
       class(ax_t), intent(in) :: Ax
       type(mesh_t), intent(inout) :: msh
       type(space_t), intent(inout) :: Xh
       type(field_t), intent(inout) :: p, u, v, w
       type(field_t), intent(inout) :: u_res, v_res, w_res
       type(field_t), intent(in) :: f_x, f_y, f_z
       type(coef_t), intent(inout) :: c_Xh
       type(field_t), intent(in) :: mu
       type(field_t), intent(in) :: rho
       type(field_t), intent(in) :: Q_T
       real(kind=rp), intent(in) :: bd
       real(kind=rp), intent(in) :: dt
       integer, intent(in) :: n
     end subroutine lm_vel_res
  end interface

  interface
     !> Factory for the low-Mach pressure residual. Selects the backend.
     module subroutine lowmach_prs_res_factory(object)
       class(lowmach_prs_res_t), allocatable, intent(inout) :: object
     end subroutine lowmach_prs_res_factory

     !> Factory for the low-Mach velocity residual. Selects the backend.
     module subroutine lowmach_vel_res_factory(object)
       class(lowmach_vel_res_t), allocatable, intent(inout) :: object
     end subroutine lowmach_vel_res_factory
  end interface

  public :: lowmach_prs_res_factory, lowmach_vel_res_factory

end module lowmach_residual
