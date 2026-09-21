!> @file adjoint_pnpn_res_cpu.f90
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
!> Residuals in the Pn-Pn formulation (CPU version)
module adjoint_pnpn_res_cpu
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
  use math, only : cfill, col3, cmult, col2, chsign, add2
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  type, public, extends(adjoint_pnpn_prs_res_t) :: adjoint_pnpn_prs_res_cpu_t
   contains
     procedure, nopass :: compute => adjoint_pnpn_prs_res_cpu_compute
  end type adjoint_pnpn_prs_res_cpu_t

  type, public, extends(adjoint_pnpn_vel_res_t) :: adjoint_pnpn_vel_res_cpu_t
   contains
     procedure, nopass :: compute => adjoint_pnpn_vel_res_cpu_compute
  end type adjoint_pnpn_vel_res_cpu_t

contains

  !> Compute adjoint pressure residual (CPU backend).
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
  !! @param event Backend event handle (unused on CPU).
  subroutine adjoint_pnpn_prs_res_cpu_compute(p, p_res, u, v, w, &
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
    real(kind=rp) :: rho_val
    integer :: n
    type(field_t), pointer :: ta1, ta2, ta3, wa1, wa2, wa3
    integer :: temp_indices(6)

    call neko_scratch_registry%request_field(ta1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(ta2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ta3, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(wa1, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(wa2, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(wa3, temp_indices(6), .false.)

    n = c_Xh%dof%size()

    ! We assume the material properties are constant
    rho_val = rho%x(1,1,1,1)
    call cfill(c_Xh%h1, 1.0_rp / rho_val, n)
    call cfill(c_Xh%h2, 0.0_rp, n)
    c_Xh%ifh2 = .false.

    call col3(ta1%x, f_x%x, c_Xh%B, n)
    call col3(ta2%x, f_y%x, c_Xh%B, n)
    call col3(ta3%x, f_z%x, c_Xh%B, n)
    call cmult(ta1%x, 1.0_rp / rho_val, n)
    call cmult(ta2%x, 1.0_rp / rho_val, n)
    call cmult(ta3%x, 1.0_rp / rho_val, n)

    call gs_Xh%op(ta1, GS_OP_ADD)
    call gs_Xh%op(ta2, GS_OP_ADD)
    call gs_Xh%op(ta3, GS_OP_ADD)

    call col2(ta1%x, c_Xh%Binv, n)
    call col2(ta2%x, c_Xh%Binv, n)
    call col2(ta3%x, c_Xh%Binv, n)

    call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call cdtp(wa3%x, ta3%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)

    ! maybe I'm missing a mu??

    call Ax%compute(p_res%x, p%x, c_Xh, p%msh, p%Xh)

    call chsign(p_res%x, n)
    call add2(p_res%x, wa1%x, n)
    call add2(p_res%x, wa2%x, n)
    call add2(p_res%x, wa3%x, n)

    ! This is commented out because I don't think it applies any more... but
    ! sym BCs haven't been tested thoroughly.
    !  !
    !  ! Surface velocity terms
    !  !
    !  do concurrent (i = 1:n)
    !     wa1%x(i,1,1,1) = 0.0_rp
    !     wa2%x(i,1,1,1) = 0.0_rp
    !     wa3%x(i,1,1,1) = 0.0_rp
    !  end do

    !  call bc_sym_surface%apply_surfvec(wa1%x, wa2%x, wa3%x, ta1%x, ta2%x, &
    !                                    ta3%x,&
    !                                    n)

    !  dtbd = bd / dt
    !  do concurrent (i = 1:n)
    !     ta1%x(i,1,1,1) = 0.0_rp
    !     ta2%x(i,1,1,1) = 0.0_rp
    !     ta3%x(i,1,1,1) = 0.0_rp
    !  end do

    !  call bc_prs_surface%apply_surfvec(ta1%x, ta2%x, ta3%x, u%x, v%x, w%x, n)

    !  do concurrent (i = 1:n)
    !     p_res%x(i,1,1,1) = p_res%x(i,1,1,1) &
    !          - (dtbd * (ta1%x(i,1,1,1) + ta2%x(i,1,1,1) + ta3%x(i,1,1,1))) &
    !          - (wa1%x(i,1,1,1) + wa2%x(i,1,1,1) + wa3%x(i,1,1,1))
    !  end do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine adjoint_pnpn_prs_res_cpu_compute

  !> Compute adjoint velocity residual (CPU backend).
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
  subroutine adjoint_pnpn_vel_res_cpu_compute(Ax, u, v, w, u_res, &
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
    real(kind=rp) :: rho_val, mu_val
    ! We assume the material properties are constant
    rho_val = rho%x(1,1,1,1)
    mu_val = mu%x(1,1,1,1)

    call cfill(c_Xh%h1, mu_val, n)
    call cfill(c_Xh%h2, rho_val * bd / dt, n)
    c_Xh%ifh2 = .true.

    call Ax%compute(u_res%x, u%x, c_Xh, msh, Xh)
    call Ax%compute(v_res%x, v%x, c_Xh, msh, Xh)
    call Ax%compute(w_res%x, w%x, c_Xh, msh, Xh)

    call chsign(u_res%x, n)
    call chsign(v_res%x, n)
    call chsign(w_res%x, n)
    call add2(u_res%x, f_x%x, n)
    call add2(v_res%x, f_y%x, n)
    call add2(w_res%x, f_z%x, n)

  end subroutine adjoint_pnpn_vel_res_cpu_compute

end module adjoint_pnpn_res_cpu
