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
!> Leray projection of a velocity field onto the divergence-free subspace.
!! Unrelated to the module `projection`, which projects Krylov solutions.
module div_free_projection
  use num_types, only : rp, i8
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use ax_product, only : ax_t
  use krylov, only : ksp_t, ksp_monitor_t
  use precon, only : pc_t
  use scalar_bc_projector, only : scalar_bc_projector_t
  use vector_bc_projector, only : vector_bc_projector_t
  use facet_normal, only : facet_normal_t
  use operators, only : cdtp, opgrad, ortho, div, rotate_cyc
  use opr_device, only : device_ortho
  use neko_config, only : NEKO_BCKND_DEVICE
  use math, only : add2, sub2, col2, cmult, rzero, glsc3, glsum, cfill, &
       NEKO_EPS
  use device_math, only : device_add2, device_sub2, device_col2, &
       device_cmult, device_rzero, device_glsc3, device_glsum, device_cfill
  use scratch_registry, only : neko_scratch_registry
  use logger, only : neko_log, LOG_SIZE
  implicit none
  private

  public :: project_div_free, div_norm

contains

  !> Project a velocity field onto the divergence-free subspace.
  !!
  !! @details Solves \f$ \nabla^2 \phi = \nabla \cdot \mathbf{u} \f$ with
  !! \f$ \partial_n \phi = 0 \f$ on the strong velocity boundaries
  !! \f$ \Gamma_u \f$ and \f$ \phi = 0 \f$ on the strong pressure
  !! boundaries, then sets \f$ \mathbf{u} \leftarrow \mathbf{u} - \nabla
  !! \phi \f$. The weak form is \f$ (\nabla \phi, \nabla \psi) =
  !! (\mathbf{u}, \nabla \psi) - \langle \mathbf{u} \cdot \mathbf{n}, \psi
  !! \rangle_{\Gamma_u} \f$. It is assembled with the operators, boundary
  !! conditions and \f$ 1/\rho \f$ scaling of the pnpn pressure step. The
  !! correction is masked on \f$ \Gamma_u \f$ as the momentum residual is, so
  !! the velocity conditions imposed there are kept. Its tangential part is
  !! lost there, which leaves a one-point layer of divergence along
  !! \f$ \Gamma_u \f$, as in the time loop.
  !!
  !! @param u,v,w The velocity, modified in place.
  !! @param coef SEM coefficients.
  !! @param gs Gather-scatter handle.
  !! @param Ax Operator of the pressure.
  !! @param ksp Krylov solver of the pressure.
  !! @param pc Preconditioner of the pressure.
  !! @param bc_prs_projector Projector onto the strong pressure boundaries.
  !! @param bc_vel_projector Projector onto the strong velocity boundaries.
  !! @param bc_prs_surface Facet normals of the strong velocity boundaries.
  !! @param prs_dirichlet Whether any strong pressure boundary exists.
  !! @param glb_n_points Global number of (non-unique) GLL points.
  !! @param rho The (constant) density.
  !! @param rel_tol Residual reduction to solve to, relative to the initial
  !! residual. The solver's own tolerance is absolute and tuned for the
  !! pressure increment of a time step.
  !! @param max_iter Iteration cap for the solve.
  subroutine project_div_free(u, v, w, coef, gs, Ax, ksp, pc, &
       bc_prs_projector, bc_vel_projector, bc_prs_surface, prs_dirichlet, &
       glb_n_points, rho, rel_tol, max_iter)
    type(field_t), intent(inout) :: u, v, w
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    class(ax_t), intent(in) :: Ax
    class(ksp_t), intent(inout) :: ksp
    class(pc_t), intent(inout) :: pc
    class(scalar_bc_projector_t), intent(inout) :: bc_prs_projector
    class(vector_bc_projector_t), intent(inout) :: bc_vel_projector
    type(facet_normal_t), intent(in) :: bc_prs_surface
    logical, intent(in) :: prs_dirichlet
    integer(kind=i8), intent(in) :: glb_n_points
    real(kind=rp), intent(in) :: rho
    real(kind=rp), intent(in) :: rel_tol
    integer, intent(in) :: max_iter
    type(field_t), pointer :: ta1, ta2, ta3, phi, rhs
    type(ksp_monitor_t) :: ksp_result
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: div_before, div_after, net_flux, rhs_norm, rhs_scale
    real(kind=rp) :: abs_tol, res_target, res_final
    integer :: temp_indices(5), ksp_max_iter
    integer :: n

    n = coef%dof%size()

    div_before = div_norm(u, v, w, coef, gs)

    call neko_scratch_registry%request_field(ta1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(ta2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ta3, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(phi, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(rhs, temp_indices(5), .false.)

    ! The operator of the pressure step, 1/rho times the stiffness matrix. The
    ! scaling makes no difference to the projection, but the preconditioner
    ! may estimate spectral bounds on its first use and keep them.
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(coef%h1_d, 1.0_rp / rho, n)
       call device_rzero(coef%h2_d, n)
    else
       call cfill(coef%h1, 1.0_rp / rho, n)
       call rzero(coef%h2, n)
    end if
    coef%ifh2 = .false.

    ! Volume term (u, grad psi), as in the pressure residual. The size of its
    ! parts tells a right-hand side that is pure round-off from a real one.
    call cdtp(rhs%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call cdtp(ta1%x, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call cdtp(ta2%x, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       rhs_scale = sqrt(device_glsc3(rhs%x_d, coef%mult_d, rhs%x_d, n)) &
            + sqrt(device_glsc3(ta1%x_d, coef%mult_d, ta1%x_d, n)) &
            + sqrt(device_glsc3(ta2%x_d, coef%mult_d, ta2%x_d, n))
       call device_add2(rhs%x_d, ta1%x_d, n)
       call device_add2(rhs%x_d, ta2%x_d, n)
    else
       rhs_scale = sqrt(glsc3(rhs%x, coef%mult, rhs%x, n)) &
            + sqrt(glsc3(ta1%x, coef%mult, ta1%x, n)) &
            + sqrt(glsc3(ta2%x, coef%mult, ta2%x, n))
       call add2(rhs%x, ta1%x, n)
       call add2(rhs%x, ta2%x, n)
    end if

    ! Surface term -<u.n, psi> on the strong velocity boundaries. Without it
    ! the projection would enforce zero flux through them, e.g. an inflow.
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

    ! Pure Neumann problem: the net flux has to vanish for it to be solvable.
    ! Report it, and remove the constant mode, which leaves a nonzero flux
    ! behind as uniform divergence.
    if (.not. prs_dirichlet) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          net_flux = device_glsum(rhs%x_d, n)
       else
          net_flux = glsum(rhs%x, n)
       end if

       write (log_buf, '(A,ES13.6)') 'Net boundary flux :', -net_flux
       call neko_log%message(log_buf)

       if (abs(net_flux) .gt. 100.0_rp * NEKO_EPS * rhs_scale) then
          call neko_log%warning('The net flux through the velocity ' // &
               'boundaries does not vanish, its share of the divergence ' // &
               'cannot be removed')
       end if

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_ortho(rhs%x_d, glb_n_points, n)
       else
          call ortho(rhs%x, glb_n_points, n)
       end if
    end if

    call gs%op(rhs, GS_OP_ADD)
    call bc_prs_projector%apply(rhs%x, n)

    ! Solve to a residual reduction of rel_tol. The initial residual is
    ! measured the way the solvers measure theirs, and the solver's absolute
    ! tolerance and iteration cap are set for the duration of this solve. A
    ! right-hand side at the round-off of the terms it is summed from means
    ! the field is already divergence free, and there is nothing to solve for.
    if (NEKO_BCKND_DEVICE .eq. 1) then
       rhs_norm = device_glsc3(rhs%x_d, coef%mult_d, rhs%x_d, n)
       call device_rzero(phi%x_d, n)
    else
       rhs_norm = glsc3(rhs%x, coef%mult, rhs%x, n)
       call rzero(phi%x, n)
    end if
    rhs_norm = sqrt(max(rhs_norm, 0.0_rp) / coef%volume)
    rhs_scale = rhs_scale / sqrt(coef%volume)

    if (rhs_norm .le. 100.0_rp * NEKO_EPS * rhs_scale) then
       write (log_buf, '(A,ES10.3,A)') '||div u||_L2      : ', div_before, &
            ' -> unchanged'
       call neko_log%message(log_buf)
       call neko_log%message('Already divergence free to working precision')
       call neko_scratch_registry%relinquish_field(temp_indices)
       return
    end if

    res_target = rel_tol * rhs_norm
    abs_tol = ksp%abs_tol
    ksp_max_iter = ksp%max_iter
    ksp%abs_tol = res_target
    ksp%max_iter = max_iter
    call pc%update()
    ksp_result = ksp%solve(Ax, phi, rhs%x, n, coef, bc_prs_projector, gs, &
         niter = max_iter)
    ksp%abs_tol = abs_tol
    ksp%max_iter = ksp_max_iter

    if (.not. prs_dirichlet) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_ortho(phi%x_d, glb_n_points, n)
       else
          call ortho(phi%x, glb_n_points, n)
       end if
    end if

    ! The residual of the solve, measured as the solvers measure theirs. Not
    ! every solver reports its own: cheby runs a fixed number of iterations
    ! without computing it, and gmres only estimates it.
    call Ax%compute(ta1%x, phi%x, coef, u%msh, u%Xh)
    call gs%op(ta1, GS_OP_ADD)
    call bc_prs_projector%apply(ta1%x, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_sub2(ta1%x_d, rhs%x_d, n)
       res_final = device_glsc3(ta1%x_d, coef%mult_d, ta1%x_d, n)
    else
       call sub2(ta1%x, rhs%x, n)
       res_final = glsc3(ta1%x, coef%mult, ta1%x, n)
    end if
    res_final = sqrt(max(res_final, 0.0_rp) / coef%volume)

    write (log_buf, '(A,I5,A,ES9.2,A,ES9.2,A,ES8.1,A)') &
         'Poisson solve: ', ksp_result%iter, ' iters, res ', &
         rhs_norm, ' -> ', res_final, ' (target ', res_target, ')'
    call neko_log%message(log_buf)

    if (res_final .gt. res_target) then
       call neko_log%warning('The divergence-free projection did not reach ' &
            // 'its tolerance, raise divergence_free_max_iterations or use ' &
            // 'a stronger pressure preconditioner')
    end if

    ! u <- u - grad(phi) / rho. opgrad is the weak gradient: assemble it, mask
    ! it on the strong velocity boundaries as the momentum residual is masked,
    ! and scale it by the inverse mass matrix.
    call opgrad(ta1%x, ta2%x, ta3%x, phi%x, coef)

    call rotate_cyc(ta1, ta2, ta3, 1, coef)
    call gs%op(ta1%x, ta2%x, ta3%x, n, GS_OP_ADD)
    call rotate_cyc(ta1, ta2, ta3, 0, coef)
    call bc_vel_projector%apply(ta1%x, ta2%x, ta3%x, n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(ta1%x_d, coef%Binv_d, n)
       call device_col2(ta2%x_d, coef%Binv_d, n)
       call device_col2(ta3%x_d, coef%Binv_d, n)
       call device_cmult(ta1%x_d, 1.0_rp / rho, n)
       call device_cmult(ta2%x_d, 1.0_rp / rho, n)
       call device_cmult(ta3%x_d, 1.0_rp / rho, n)
       call device_sub2(u%x_d, ta1%x_d, n)
       call device_sub2(v%x_d, ta2%x_d, n)
       call device_sub2(w%x_d, ta3%x_d, n)
    else
       call col2(ta1%x, coef%Binv, n)
       call col2(ta2%x, coef%Binv, n)
       call col2(ta3%x, coef%Binv, n)
       call cmult(ta1%x, 1.0_rp / rho, n)
       call cmult(ta2%x, 1.0_rp / rho, n)
       call cmult(ta3%x, 1.0_rp / rho, n)
       call sub2(u%x, ta1%x, n)
       call sub2(v%x, ta2%x, n)
       call sub2(w%x, ta3%x, n)
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

    div_after = div_norm(u, v, w, coef, gs)

    write (log_buf, '(A,2(ES10.3,A))') '||div u||_L2      : ', &
         div_before, ' -> ', div_after, ''
    call neko_log%message(log_buf)

  end subroutine project_div_free

  !> The \f$ L^2 \f$ norm \f$ \sqrt{\int_\Omega (\nabla \cdot \mathbf{u})^2
  !! \, dV} \f$ of the point-wise divergence, made continuous by averaging
  !! across elements before the mass-weighted integral.
  !! @param u,v,w The velocity.
  !! @param coef SEM coefficients.
  !! @param gs Gather-scatter handle.
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
