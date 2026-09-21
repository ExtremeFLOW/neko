!> @file adjoint_pnpn_res_device.F90
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
!
!> Residuals in the Pn-Pn formulation (device backend)
module adjoint_pnpn_res_device
  use gather_scatter, only : gs_t, GS_OP_ADD
  use operators, only : cdtp
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use facet_normal, only : facet_normal_t
  use adjoint_pnpn_residual, only : adjoint_pnpn_prs_res_t, &
       adjoint_pnpn_vel_res_t
  use scratch_registry, only: neko_scratch_registry
  use mesh, only : mesh_t
  use num_types, only : rp
  use space, only : space_t
  use device_math, only : device_cfill, device_cmult, device_cmult2, &
       device_col2, device_add2
  use device, only : device_event_sync
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  type, public, extends(adjoint_pnpn_prs_res_t) :: adjoint_pnpn_prs_res_device_t
   contains
     procedure, nopass :: compute => adjoint_pnpn_prs_res_device_compute
  end type adjoint_pnpn_prs_res_device_t

  type, public, extends(adjoint_pnpn_vel_res_t) :: adjoint_pnpn_vel_res_device_t
   contains
     procedure, nopass :: compute => adjoint_pnpn_vel_res_device_compute
  end type adjoint_pnpn_vel_res_device_t

contains

  !> Compute adjoint pressure residual (device backend).
  !! @param p Adjoint pressure field.
  !! @param p_res Pressure residual output.
  !! @param u Adjoint velocity x-component.
  !! @param v Adjoint velocity y-component.
  !! @param w Adjoint velocity z-component.
  !! @param f_x Explicit forcing x-component.
  !! @param f_y Explicit forcing y-component.
  !! @param f_z Explicit forcing z-component.
  !! @param c_Xh Coefficients on the pressure space.
  !! @param gs_Xh Gather-scatter operator on the pressure space.
  !! @param bc_prs_surface Pressure boundary surface normals.
  !! @param bc_sym_surface Symmetry boundary surface normals.
  !! @param Ax Helmholtz operator.
  !! @param bd BDF coefficient for the current step.
  !! @param dt Time-step size.
  !! @param mu Dynamic viscosity field (assumed constant).
  !! @param rho Density field (assumed constant).
  !! @param event Backend event handle for gather-scatter synchronization.
  subroutine adjoint_pnpn_prs_res_device_compute(p, p_res, u, v, w, &
       f_x, f_y, f_z, c_Xh, gs_Xh, bc_prs_surface, bc_sym_surface, &
       Ax, bd, dt, mu, rho, event)
    type(field_t), intent(inout) :: p, u, v, w
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
    type(c_ptr), intent(inout) :: event
    real(kind=rp) :: rho_val, inv_rho
    integer :: n
    type(field_t), pointer :: ta1, ta2, ta3, wa1, wa2, wa3
    integer :: temp_indices(6)

    call neko_scratch_registry%request_field(ta1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(ta2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ta3, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(wa1, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(wa2, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(wa3, temp_indices(6), .false.)

    n = u%dof%size()

    ! We assume the material properties are constant
    rho_val = rho%x(1,1,1,1)
    inv_rho = 1.0_rp / rho_val

    c_Xh%ifh2 = .false.
    call device_cfill(c_Xh%h1_d, inv_rho, n)
    call device_cfill(c_Xh%h2_d, 0.0_rp, n)

    call device_cmult2(ta1%x_d, f_x%x_d, inv_rho, n)
    call device_cmult2(ta2%x_d, f_y%x_d, inv_rho, n)
    call device_cmult2(ta3%x_d, f_z%x_d, inv_rho, n)

    call device_col2(ta1%x_d, c_Xh%B_d, n)
    call device_col2(ta2%x_d, c_Xh%B_d, n)
    call device_col2(ta3%x_d, c_Xh%B_d, n)

    call gs_Xh%op(ta1, GS_OP_ADD, event)
    call gs_Xh%op(ta2, GS_OP_ADD, event)
    call gs_Xh%op(ta3, GS_OP_ADD, event)
    call device_event_sync(event)

    call device_col2(ta1%x_d, c_Xh%Binv_d, n)
    call device_col2(ta2%x_d, c_Xh%Binv_d, n)
    call device_col2(ta3%x_d, c_Xh%Binv_d, n)

    call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call cdtp(wa3%x, ta3%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)

    call Ax%compute(p_res%x, p%x, c_Xh, p%msh, p%Xh)

    call device_cmult(p_res%x_d, -1.0_rp, n)
    call device_add2(p_res%x_d, wa1%x_d, n)
    call device_add2(p_res%x_d, wa2%x_d, n)
    call device_add2(p_res%x_d, wa3%x_d, n)

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine adjoint_pnpn_prs_res_device_compute

  !> Compute adjoint velocity residual (device backend).
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
  !! @param mu Dynamic viscosity field (assumed constant).
  !! @param rho Density field (assumed constant).
  !! @param bd BDF coefficient for the current step.
  !! @param dt Time-step size.
  !! @param n Number of degrees of freedom.
  subroutine adjoint_pnpn_vel_res_device_compute(Ax, u, v, w, u_res, &
       v_res, w_res, p, f_x, f_y, f_z, c_Xh, msh, Xh, mu, rho, bd, dt, n)
    class(ax_t), intent(in) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(inout) :: u_res, v_res, w_res
    type(field_t), intent(in) :: f_x, f_y, f_z
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: mu
    type(field_t), intent(in) :: rho
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n
    real(kind=rp) :: mu_val, rho_val

    ! We assume the material properties are constant
    mu_val = mu%x(1,1,1,1)
    rho_val = rho%x(1,1,1,1)

    call device_cfill(c_Xh%h1_d, mu_val, n)
    call device_cfill(c_Xh%h2_d, rho_val * (bd / dt), n)
    c_Xh%ifh2 = .true.

    call Ax%compute_vector(u_res%x, v_res%x, w_res%x, u%x, v%x, w%x, c_Xh, &
         msh, Xh)

    call device_cmult(u_res%x_d, -1.0_rp, n)
    call device_add2(u_res%x_d, f_x%x_d, n)

    call device_cmult(v_res%x_d, -1.0_rp, n)
    call device_add2(v_res%x_d, f_y%x_d, n)

    call device_cmult(w_res%x_d, -1.0_rp, n)
    call device_add2(w_res%x_d, f_z%x_d, n)
  end subroutine adjoint_pnpn_vel_res_device_compute

end module adjoint_pnpn_res_device
