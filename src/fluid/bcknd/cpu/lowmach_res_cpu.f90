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
!! Pressure residual: variable-density Poisson operator (h1 = 1/rho). The
!! viscous force uses the rotational stress form mu(curl curl u) - S^T(2 grad mu)
!! plus the low-Mach correction -(4/3) mu grad(Q_T) + (2/3) Q_T grad(mu), and the
!! continuity constraint div u = Q_T enters as the source (bd/dt) Q_T B.
!!
!! Velocity residual (momentum equation): the deviatoric viscous stress
!!   tau = mu (grad u + grad u^T) - (2/3) mu (div u) I
!! is split into
!!   * div[ mu (grad u + grad u^T) ]  -- evaluated by the EXISTING coupled
!!     full-stress operator (ax_helm_full) via Ax%compute_vector, with mu in
!!     coef%h1; no extra arguments are threaded into ax_helm, and
!!   * the dilatation -(2/3) grad(mu Q_T), assembled explicitly as a body
!!     force (Q_T is known) with the generic weak-gradient operator opgrad.
!! Using the low-Mach scheme therefore requires full_stress_formulation = true
!! so that Ax_vel is the coupled stress operator.
module lowmach_res_cpu
  use gather_scatter, only : gs_t, GS_OP_ADD
  use operators, only : opgrad, curl, cdtp, rotate_cyc, dudxyz, strain_rate
  use math, only : copy, cmult2, rzero, vdot3, cmult, sub2, col2, invers2
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
    real(kind=rp), parameter :: c43 = 4.0_rp / 3.0_rp
    real(kind=rp), parameter :: c13 = 1.0_rp / 3.0_rp
    real(kind=rp) :: dtbd
    integer :: n, nelv, lxyz
    integer :: i, e
    type(field_t), pointer :: ta1, ta2, ta3, wa1, wa2, wa3, work1, work2, work3
    type(field_t), pointer :: s11, s22, s33, s12, s13, s23
    integer :: temp_indices(15)

    call neko_scratch_registry%request_field(ta1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(ta2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ta3, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(wa1, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(wa2, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(wa3, temp_indices(6), .false.)
    call neko_scratch_registry%request_field(work1, temp_indices(7), .false.)
    call neko_scratch_registry%request_field(work2, temp_indices(8), .false.)
    call neko_scratch_registry%request_field(work3, temp_indices(9), .false.)
    call neko_scratch_registry%request_field(s11, temp_indices(10), .false.)
    call neko_scratch_registry%request_field(s22, temp_indices(11), .false.)
    call neko_scratch_registry%request_field(s33, temp_indices(12), .false.)
    call neko_scratch_registry%request_field(s12, temp_indices(13), .false.)
    call neko_scratch_registry%request_field(s13, temp_indices(14), .false.)
    call neko_scratch_registry%request_field(s23, temp_indices(15), .false.)

    n = c_Xh%dof%size()
    lxyz = c_Xh%Xh%lxyz
    nelv = c_Xh%msh%nelv

    ! Variable-density pressure-Poisson coefficient h1 = 1/rho.
    call invers2(c_Xh%h1, rho%x, n)
    call rzero(c_Xh%h2, n)
    c_Xh%ifh2 = .false.

    ! Deformation stress in rotational form: mu * (double curl of u_e), minus
    ! the variable-viscosity term S^T (2 grad mu). Together these give
    ! -div[ mu (grad u + grad u^T) ] for a divergence-free field; the
    ! divergence (low-Mach) correction is added separately below.
    call curl(ta1, ta2, ta3, u_e, v_e, w_e, work1, work2, c_Xh)
    call curl(wa1, wa2, wa3, ta1, ta2, ta3, work1, work2, c_Xh)
    call col2(wa1%x, mu%x, n)
    call col2(wa2%x, mu%x, n)
    call col2(wa3%x, mu%x, n)

    call strain_rate(s11%x, s22%x, s33%x, s12%x, s13%x, s23%x, &
         u_e, v_e, w_e, c_Xh)

    ! ta = 2 grad(mu)  (kept for both S^T grad mu and the low-Mach term).
    call dudxyz(ta1%x, mu%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call dudxyz(ta2%x, mu%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call dudxyz(ta3%x, mu%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
    call cmult(ta1%x, 2.0_rp, n)
    call cmult(ta2%x, 2.0_rp, n)
    call cmult(ta3%x, 2.0_rp, n)

    do e = 1, nelv
       call vdot3(work1%x(:,:,:,e), ta1%x(:,:,:,e), ta2%x(:,:,:,e), &
            ta3%x(:,:,:,e), s11%x(:,:,:,e), s12%x(:,:,:,e), s13%x(:,:,:,e), lxyz)
       call vdot3(work2%x(:,:,:,e), ta1%x(:,:,:,e), ta2%x(:,:,:,e), &
            ta3%x(:,:,:,e), s12%x(:,:,:,e), s22%x(:,:,:,e), s23%x(:,:,:,e), lxyz)
       call vdot3(work3%x(:,:,:,e), ta1%x(:,:,:,e), ta2%x(:,:,:,e), &
            ta3%x(:,:,:,e), s13%x(:,:,:,e), s23%x(:,:,:,e), s33%x(:,:,:,e), lxyz)
    end do

    ! wa = mu (curl curl u) - S^T (2 grad mu) = -div[mu(grad u + grad u^T)].
    call sub2(wa1%x, work1%x, n)
    call sub2(wa2%x, work2%x, n)
    call sub2(wa3%x, work3%x, n)

    ! Low-Mach correction to -div(tau) (div u = Q_T):
    !   wa += -(4/3) mu grad(Q_T) + (2/3) Q_T grad(mu).
    ! grad(mu) is still held as ta/2 (ta = 2 grad mu); grad(Q_T) via dudxyz.
    call dudxyz(work1%x, Q_T%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call dudxyz(work2%x, Q_T%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call dudxyz(work3%x, Q_T%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
    do concurrent (i = 1:n)
       wa1%x(i,1,1,1) = wa1%x(i,1,1,1) - c43 * mu%x(i,1,1,1) * work1%x(i,1,1,1) &
            + c13 * Q_T%x(i,1,1,1) * ta1%x(i,1,1,1)
       wa2%x(i,1,1,1) = wa2%x(i,1,1,1) - c43 * mu%x(i,1,1,1) * work2%x(i,1,1,1) &
            + c13 * Q_T%x(i,1,1,1) * ta2%x(i,1,1,1)
       wa3%x(i,1,1,1) = wa3%x(i,1,1,1) - c43 * mu%x(i,1,1,1) * work3%x(i,1,1,1) &
            + c13 * Q_T%x(i,1,1,1) * ta3%x(i,1,1,1)
    end do

    ! Momentum acceleration (f + div tau)/rho = (f - wa)/rho.
    do concurrent (i = 1:n)
       ta1%x(i,1,1,1) = f_x%x(i,1,1,1) / rho%x(i,1,1,1) &
            - (wa1%x(i,1,1,1) / rho%x(i,1,1,1)) * c_Xh%B(i,1,1,1)
       ta2%x(i,1,1,1) = f_y%x(i,1,1,1) / rho%x(i,1,1,1) &
            - (wa2%x(i,1,1,1) / rho%x(i,1,1,1)) * c_Xh%B(i,1,1,1)
       ta3%x(i,1,1,1) = f_z%x(i,1,1,1) / rho%x(i,1,1,1) &
            - (wa3%x(i,1,1,1) / rho%x(i,1,1,1)) * c_Xh%B(i,1,1,1)
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

    ! Low-Mach thermal divergence source (div u = Q_T):
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
    integer :: temp_indices(7)
    type(field_t), pointer :: ta1, ta2, ta3, gx, gy, gz, muqt
    integer :: i
    real(kind=rp), parameter :: c_dil = 2.0_rp / 3.0_rp

    ! Helmholtz coefficients for the implicit viscous/unsteady operator:
    !   h1 = mu  (viscous),  h2 = rho * bd/dt  (unsteady mass term).
    call copy(c_Xh%h1, mu%x, n)
    call cmult2(c_Xh%h2, rho%x, bd / dt, n)
    c_Xh%ifh2 = .true.

    ! Implicit deviatoric stress: weak div[ mu (grad u + grad u^T) ].
    ! This is the EXISTING coupled full-stress operator (ax_helm_full); mu is
    ! carried in c_Xh%h1, so ax_helm's interface is untouched. Requires Ax_vel
    ! to be the full-stress operator (full_stress_formulation = true).
    call Ax%compute_vector(u_res%x, v_res%x, w_res%x, u%x, v%x, w%x, &
         c_Xh, msh, Xh)

    call neko_scratch_registry%request_field(ta1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(ta2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ta3, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(gx, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(gy, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(gz, temp_indices(6), .false.)
    call neko_scratch_registry%request_field(muqt, temp_indices(7), .false.)

    ! Pressure gradient (weak).
    call opgrad(ta1%x, ta2%x, ta3%x, p%x, c_Xh)

    ! Low-Mach dilatation part of the deviatoric stress:
    !   div(tau)_dil = -(2/3) grad( mu * Q_T ),   since div u = Q_T.
    ! Q_T is known explicitly this step, so it is assembled as a body force
    ! with the generic weak-gradient operator -- no change to Ax.
    do concurrent (i = 1:n)
       muqt%x(i,1,1,1) = mu%x(i,1,1,1) * Q_T%x(i,1,1,1)
    end do
    call opgrad(gx%x, gy%x, gz%x, muqt%x, c_Xh)

    ! Momentum residual:  res = f - grad p + div(tau) - (rho bd/dt) M u.
    ! -Ax(u) supplies +div[mu(grad u+grad u^T)] - (rho bd/dt) M u (weak);
    ! the dilatation -(2/3) grad(mu Q_T) is added here in the same weak sense.
    do concurrent (i = 1:n)
       u_res%x(i,1,1,1) = (-u_res%x(i,1,1,1)) - ta1%x(i,1,1,1) &
            + f_x%x(i,1,1,1) - c_dil * gx%x(i,1,1,1)
       v_res%x(i,1,1,1) = (-v_res%x(i,1,1,1)) - ta2%x(i,1,1,1) &
            + f_y%x(i,1,1,1) - c_dil * gy%x(i,1,1,1)
       w_res%x(i,1,1,1) = (-w_res%x(i,1,1,1)) - ta3%x(i,1,1,1) &
            + f_z%x(i,1,1,1) - c_dil * gz%x(i,1,1,1)
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine lowmach_vel_res_cpu_compute

end module lowmach_res_cpu
