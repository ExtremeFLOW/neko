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
!> Leray projection of a velocity field onto the discretely divergence-free
!! subspace.
!!
!! @note Not to be confused with the module `projection`, which implements the
!! Krylov-subspace projection used to generate good initial guesses for the
!! linear solvers. This module implements the Helmholtz--Leray decomposition.
module div_free_projection
  use num_types, only : rp, i8
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use ax_product, only : ax_t
  use krylov, only : ksp_t, ksp_monitor_t
  use precon, only : pc_t
  use scalar_bc_projector, only : scalar_bc_projector_t
  use facet_normal, only : facet_normal_t
  use operators, only : cdtp, opgrad, ortho, div, rotate_cyc
  use opr_device, only : device_ortho
  use neko_config, only : NEKO_BCKND_DEVICE
  use math, only : add2, sub2, col2, copy, rzero, glsc3, glsum, cfill
  use device_math, only : device_add2, device_sub2, device_col2, &
       device_copy, device_rzero, device_glsc3, device_glsum, device_cfill
  use scratch_registry, only : neko_scratch_registry
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_warning
  implicit none
  private

  public :: project_div_free, div_norm

contains

  !> Project a velocity field onto the space of (weakly) divergence-free
  !! fields, using the same discrete operators as the Pn-Pn pressure step.
  !!
  !! @details
  !! The Helmholtz--Leray decomposition of the velocity field reads
  !! \f$ \mathbf{u} = \mathbf{u}_{df} + \nabla \phi \f$, where
  !! \f$ \mathbf{u}_{df} \f$ is divergence free. The potential \f$ \phi \f$
  !! solves the Poisson problem
  !! \f{eqnarray*}{
  !!    \nabla^2 \phi &=& \nabla \cdot \mathbf{u} \quad &\text{in } \Omega, \\
  !!    \partial_n \phi &=& 0 \quad &\text{on } \Gamma_D, \\
  !!    \phi &=& 0 \quad &\text{on } \Gamma_{out},
  !! \f}
  !! with \f$ \Gamma_D \f$ the boundaries where the velocity is prescribed
  !! (walls, inflow) and \f$ \Gamma_{out} \f$ the boundaries where the
  !! pressure is prescribed (outflow). The homogeneous Neumann condition
  !! guarantees that the correction has no normal component on \f$ \Gamma_D
  !! \f$, so a prescribed inflow profile is left untouched. The tangential
  !! velocity on \f$ \Gamma_D \f$ is in general modified by the projection.
  !!
  !! Discretely, the weak form
  !! \f$ (\nabla \phi, \nabla \psi) = (\mathbf{u}, \nabla \psi) -
  !! \langle \mathbf{u}\cdot\mathbf{n}, \psi \rangle_{\Gamma_D} \f$
  !! is assembled with `cdtp` and the facet-normal surface operator, i.e. with
  !! exactly the operators that build the Pn-Pn pressure residual. The
  !! resulting field is therefore divergence free in the same discrete sense
  !! that the pressure step of the time loop enforces, which is the strongest
  !! statement available for an equal-order (Pn-Pn) discretisation.
  !!
  !! @param u The x component of the velocity, modified in place.
  !! @param v The y component of the velocity, modified in place.
  !! @param w The z component of the velocity, modified in place.
  !! @param phi Work field holding the projection potential on exit.
  !! @param rhs Work field for the right-hand side of the Poisson problem.
  !! @param coef The SEM coefficients.
  !! @param gs The gather-scatter handle.
  !! @param Ax The Helmholtz operator used for the pressure.
  !! @param ksp The Krylov solver used for the pressure.
  !! @param pc The preconditioner used for the pressure.
  !! @param bc_prs_projector Projector zeroing the strong pressure boundaries.
  !! @param bc_prs_surface Facet-normal operator on the velocity Dirichlet
  !! boundaries, as used by the pressure residual.
  !! @param prs_dirichlet Whether any strong pressure boundary is present.
  !! @param glb_n_points The global number of non-unique GLL points.
  subroutine project_div_free(u, v, w, phi, rhs, coef, gs, Ax, ksp, pc, &
       bc_prs_projector, bc_prs_surface, prs_dirichlet, glb_n_points)
    type(field_t), intent(inout) :: u, v, w
    type(field_t), intent(inout) :: phi, rhs
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    class(ax_t), intent(in) :: Ax
    class(ksp_t), intent(inout) :: ksp
    class(pc_t), intent(inout) :: pc
    class(scalar_bc_projector_t), intent(inout) :: bc_prs_projector
    type(facet_normal_t), intent(in) :: bc_prs_surface
    logical, intent(in) :: prs_dirichlet
    integer(kind=i8), intent(in) :: glb_n_points
    type(field_t), pointer :: ta1, ta2, ta3
    type(ksp_monitor_t) :: ksp_result
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: div_before, div_after, net_flux
    integer :: temp_indices(3)
    integer :: n

    n = coef%dof%size()

    call neko_log%section('Divergence-free projection')

    div_before = div_norm(u, v, w, coef, gs)

    call neko_scratch_registry%request_field(ta1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(ta2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ta3, temp_indices(3), .false.)

    !
    ! The Poisson operator is the plain (unscaled) stiffness matrix. The
    ! density drops out of the projection, so there is no need to weight it.
    !
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(coef%h1_d, 1.0_rp, n)
       call device_rzero(coef%h2_d, n)
    else
       call cfill(coef%h1, 1.0_rp, n)
       call rzero(coef%h2, n)
    end if
    coef%ifh2 = .false.

    !
    ! Volume contribution, (u, grad(psi)), assembled exactly as the
    ! divergence term of the pressure residual.
    !
    call cdtp(rhs%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call cdtp(ta1%x, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call cdtp(ta2%x, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(rhs%x_d, ta1%x_d, n)
       call device_add2(rhs%x_d, ta2%x_d, n)
    else
       call add2(rhs%x, ta1%x, n)
       call add2(rhs%x, ta2%x, n)
    end if

    !
    ! Surface contribution, -<u.n, psi> on the boundaries where the velocity
    ! is prescribed. Without it the projection would drive the normal velocity
    ! to zero on those boundaries, wiping out a prescribed inflow.
    !
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(ta1%x_d, n)
       call device_rzero(ta2%x_d, n)
       call device_rzero(ta3%x_d, n)
       call bc_prs_surface%apply_surfvec_dev(ta1%x_d, ta2%x_d, ta3%x_d, &
            u%x_d, v%x_d, w%x_d)
       call device_sub2(rhs%x_d, ta1%x_d, n)
       call device_sub2(rhs%x_d, ta2%x_d, n)
       call device_sub2(rhs%x_d, ta3%x_d, n)
    else
       call rzero(ta1%x, n)
       call rzero(ta2%x, n)
       call rzero(ta3%x, n)
       call bc_prs_surface%apply_surfvec(ta1%x, ta2%x, ta3%x, &
            u%x, v%x, w%x, n)
       call sub2(rhs%x, ta1%x, n)
       call sub2(rhs%x, ta2%x, n)
       call sub2(rhs%x, ta3%x, n)
    end if

    !
    ! With no strong pressure boundary the problem is a pure Neumann one. It
    ! is solvable only if the net flux through the boundary vanishes, so warn
    ! when it does not and remove the resulting constant from the residual.
    !
    if (.not. prs_dirichlet) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          net_flux = device_glsum(rhs%x_d, n)
       else
          net_flux = glsum(rhs%x, n)
       end if

       write (log_buf, '(A,ES13.6)') 'Net boundary flux :', -net_flux
       call neko_log%message(log_buf)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_ortho(rhs%x_d, glb_n_points, n)
       else
          call ortho(rhs%x, glb_n_points, n)
       end if
    end if

    call gs%op(rhs, GS_OP_ADD)
    call bc_prs_projector%apply(rhs%x, n)

    !
    ! Solve for the projection potential.
    !
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(phi%x_d, n)
    else
       call rzero(phi%x, n)
    end if

    call pc%update()
    ksp_result = ksp%solve(Ax, phi, rhs%x, n, coef, bc_prs_projector, gs)

    if (.not. prs_dirichlet) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_ortho(phi%x_d, glb_n_points, n)
       else
          call ortho(phi%x, glb_n_points, n)
       end if
    end if

    write (log_buf, '(A,I6,2(A,ES13.6))') 'Poisson solve     : iters ', &
         ksp_result%iter, ', res ', ksp_result%res_start, ' -> ', &
         ksp_result%res_final
    call neko_log%message(log_buf)

    if (.not. ksp_result%converged) then
       call neko_warning('The divergence-free projection did not converge, ' &
            // 'the initial condition may still be divergent.')
    end if

    !
    ! Subtract the gradient of the potential. `opgrad` returns the weak
    ! gradient, so it is assembled and scaled by the inverse mass matrix to
    ! recover a continuous, point-wise gradient, mirroring how the velocity is
    ! corrected by the pressure gradient in the time loop.
    !
    call opgrad(ta1%x, ta2%x, ta3%x, phi%x, coef)

    call rotate_cyc(ta1%x, ta2%x, ta3%x, 1, coef)
    call gs%op(ta1%x, ta2%x, ta3%x, n, GS_OP_ADD)
    call rotate_cyc(ta1%x, ta2%x, ta3%x, 0, coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(ta1%x_d, coef%Binv_d, n)
       call device_col2(ta2%x_d, coef%Binv_d, n)
       call device_col2(ta3%x_d, coef%Binv_d, n)
       call device_sub2(u%x_d, ta1%x_d, n)
       call device_sub2(v%x_d, ta2%x_d, n)
       call device_sub2(w%x_d, ta3%x_d, n)
    else
       call col2(ta1%x, coef%Binv, n)
       call col2(ta2%x, coef%Binv, n)
       call col2(ta3%x, coef%Binv, n)
       call sub2(u%x, ta1%x, n)
       call sub2(v%x, ta2%x, n)
       call sub2(w%x, ta3%x, n)
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

    div_after = div_norm(u, v, w, coef, gs)

    write (log_buf, '(A,ES13.6)') 'div(u) L2, before :', div_before
    call neko_log%message(log_buf)
    write (log_buf, '(A,ES13.6)') 'div(u) L2, after  :', div_after
    call neko_log%message(log_buf)

    call neko_log%end_section()

  end subroutine project_div_free

  !> Compute the \f$ L^2 \f$ norm of the point-wise divergence of a velocity
  !! field, \f$ \sqrt{\int_\Omega (\nabla \cdot \mathbf{u})^2 \, dV} \f$.
  !! @param u The x component of the velocity.
  !! @param v The y component of the velocity.
  !! @param w The z component of the velocity.
  !! @param coef The SEM coefficients.
  !! @param gs The gather-scatter handle.
  function div_norm(u, v, w, coef, gs) result(norm)
    type(field_t), intent(in) :: u, v, w
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    real(kind=rp) :: norm
    type(field_t), pointer :: d
    integer :: temp_index, n

    n = coef%dof%size()

    call neko_scratch_registry%request_field(d, temp_index, .false.)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call div(d%x_d, u%x_d, v%x_d, w%x_d, coef)
       call gs%op(d, GS_OP_ADD)
       call device_col2(d%x_d, coef%mult_d, n)
       norm = device_glsc3(d%x_d, coef%B_d, d%x_d, n)
    else
       call div(d%x, u%x, v%x, w%x, coef)
       call gs%op(d, GS_OP_ADD)
       call col2(d%x, coef%mult, n)
       norm = glsc3(d%x, coef%B, d%x, n)
    end if

    norm = sqrt(max(norm, 0.0_rp))

    call neko_scratch_registry%relinquish_field(temp_index)

  end function div_norm

end module div_free_projection
