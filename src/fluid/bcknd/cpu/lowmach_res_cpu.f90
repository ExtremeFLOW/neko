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
!> Low-Mach residuals, CPU backend.
!!
!! Initial implementation: the pressure and velocity residual bodies are
!! copied from pnpn_res_cpu and the Q_T argument is accepted but ignored.
!! When enabled via the scheme, this reproduces the constant-density pnpn
!! result bit-for-bit, providing a stable foundation onto which the Q_T
!! source, the nabla-Q_T correction, and variable-property stresses will
!! be layered in subsequent commits.
module lowmach_res_cpu
  use gather_scatter, only : gs_t, GS_OP_ADD
  use operators, only : opgrad, curl, cdtp, rotate_cyc
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use facet_normal, only : facet_normal_t
  use lowmach_residual, only : lowmach_prs_res_t, lowmach_vel_res_t
  use scratch_registry, only : neko_scratch_registry
  use mesh, only : mesh_t
  use num_types, only : rp
  use space, only : space_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  type, public, extends(lowmach_prs_res_t) :: lowmach_prs_res_cpu_t
   contains
     procedure, nopass :: compute => lowmach_prs_res_cpu_compute
  end type lowmach_prs_res_cpu_t

  type, public, extends(lowmach_vel_res_t) :: lowmach_vel_res_cpu_t
   contains
     procedure, nopass :: compute => lowmach_vel_res_cpu_compute
  end type lowmach_vel_res_cpu_t

contains

  subroutine lowmach_prs_res_cpu_compute(p, p_res, u, v, w, u_e, v_e, w_e, &
       f_x, f_y, f_z, c_Xh, gs_Xh, bc_prs_surface, bc_sym_surface, Ax, bd, &
       dt, mu, rho, Q_T, event)
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
    real(kind=rp) :: dtbd, rho_val, mu_val
    integer :: n
    integer :: i
    type(field_t), pointer :: ta1, ta2, ta3, wa1, wa2, wa3, work1, work2
    integer :: temp_indices(8)

    call neko_scratch_registry%request_field(ta1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(ta2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ta3, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(wa1, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(wa2, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(wa3, temp_indices(6), .false.)
    call neko_scratch_registry%request_field(work1, temp_indices(7), .false.)
    call neko_scratch_registry%request_field(work2, temp_indices(8), .false.)

    n = c_Xh%dof%size()

    ! Constant-property approximation, identical to pnpn_res_cpu.
    rho_val = rho%x(1,1,1,1)
    mu_val = mu%x(1,1,1,1)
    do i = 1, n
       c_Xh%h1(i,1,1,1) = 1.0_rp / rho_val
       c_Xh%h2(i,1,1,1) = 0.0_rp
    end do
    c_Xh%ifh2 = .false.

    call curl(ta1, ta2, ta3, u_e, v_e, w_e, work1, work2, c_Xh)
    call curl(wa1, wa2, wa3, ta1, ta2, ta3, work1, work2, c_Xh)

    do concurrent (i = 1:n)
       ta1%x(i,1,1,1) = f_x%x(i,1,1,1) / rho_val &
            - ((wa1%x(i,1,1,1) * (mu_val / rho_val)) * c_Xh%B(i,1,1,1))
       ta2%x(i,1,1,1) = f_y%x(i,1,1,1) / rho_val &
            - ((wa2%x(i,1,1,1) * (mu_val / rho_val)) * c_Xh%B(i,1,1,1))
       ta3%x(i,1,1,1) = f_z%x(i,1,1,1) / rho_val &
            - ((wa3%x(i,1,1,1) * (mu_val / rho_val)) * c_Xh%B(i,1,1,1))
    end do

    call rotate_cyc(ta1%x, ta2%x, ta3%x, 1, c_Xh)
    call gs_Xh%op(ta1, GS_OP_ADD)
    call gs_Xh%op(ta2, GS_OP_ADD)
    call gs_Xh%op(ta3, GS_OP_ADD)
    call rotate_cyc(ta1%x, ta2%x, ta3%x, 0, c_Xh)

    do concurrent (i = 1:n)
       ta1%x(i,1,1,1) = ta1%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
       ta2%x(i,1,1,1) = ta2%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
       ta3%x(i,1,1,1) = ta3%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
    end do

    call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call cdtp(wa3%x, ta3%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)

    call Ax%compute(p_res%x, p%x, c_Xh, p%msh, p%Xh)

    dtbd = bd / dt

    do concurrent (i = 1:n)
       p_res%x(i,1,1,1) = (-p_res%x(i,1,1,1)) &
            + wa1%x(i,1,1,1) + wa2%x(i,1,1,1) + wa3%x(i,1,1,1)
    end do

    ! Low-Mach thermal divergence source:
    !   p_res += (bd/dt) * Q_T * B
    ! Follows the convention of Nek5000 plan4.f:
    !   call admcol3(respr, QTL, bm1, dtbd, ntot1)
    do concurrent (i = 1:n)
       p_res%x(i,1,1,1) = p_res%x(i,1,1,1) &
            + dtbd * Q_T%x(i,1,1,1) * c_Xh%B(i,1,1,1)
    end do

    do concurrent (i = 1:n)
       wa1%x(i,1,1,1) = 0.0_rp
       wa2%x(i,1,1,1) = 0.0_rp
       wa3%x(i,1,1,1) = 0.0_rp
    end do

    call bc_sym_surface%apply_surfvec(wa1%x, wa2%x, wa3%x, &
         ta1%x, ta2%x, ta3%x, n)
    do concurrent (i = 1:n)
       ta1%x(i,1,1,1) = 0.0_rp
       ta2%x(i,1,1,1) = 0.0_rp
       ta3%x(i,1,1,1) = 0.0_rp
    end do

    call bc_prs_surface%apply_surfvec(ta1%x, ta2%x, ta3%x, u%x, v%x, w%x, n)

    do concurrent (i = 1:n)
       p_res%x(i,1,1,1) = p_res%x(i,1,1,1) &
            - (dtbd * (ta1%x(i,1,1,1) + ta2%x(i,1,1,1) + ta3%x(i,1,1,1))) &
            - (wa1%x(i,1,1,1) + wa2%x(i,1,1,1) + wa3%x(i,1,1,1))
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine lowmach_prs_res_cpu_compute

  subroutine lowmach_vel_res_cpu_compute(Ax, u, v, w, u_res, v_res, w_res, &
       p, f_x, f_y, f_z, c_Xh, msh, Xh, mu, rho, Q_T, bd, dt, n)
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
    real(kind=rp) :: rho_val, mu_val
    integer :: temp_indices(3)
    type(field_t), pointer :: ta1, ta2, ta3
    integer :: i

    associate (unused_Q_T => Q_T); end associate

    rho_val = rho%x(1,1,1,1)
    mu_val = mu%x(1,1,1,1)

    do concurrent (i = 1:n)
       c_Xh%h1(i,1,1,1) = mu_val
       c_Xh%h2(i,1,1,1) = rho_val * bd / dt
    end do
    c_Xh%ifh2 = .true.

    call Ax%compute(u_res%x, u%x, c_Xh, msh, Xh)
    call Ax%compute(v_res%x, v%x, c_Xh, msh, Xh)
    call Ax%compute(w_res%x, w%x, c_Xh, msh, Xh)

    call neko_scratch_registry%request_field(ta1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(ta2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ta3, temp_indices(3), .false.)

    call opgrad(ta1%x, ta2%x, ta3%x, p%x, c_Xh)

    do concurrent (i = 1:n)
       u_res%x(i,1,1,1) = (-u_res%x(i,1,1,1)) - ta1%x(i,1,1,1) + f_x%x(i,1,1,1)
       v_res%x(i,1,1,1) = (-v_res%x(i,1,1,1)) - ta2%x(i,1,1,1) + f_y%x(i,1,1,1)
       w_res%x(i,1,1,1) = (-w_res%x(i,1,1,1)) - ta3%x(i,1,1,1) + f_z%x(i,1,1,1)
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine lowmach_vel_res_cpu_compute

end module lowmach_res_cpu
