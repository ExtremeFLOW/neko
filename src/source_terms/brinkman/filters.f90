! Copyright (c) 2024, The Neko Authors
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
!> A module containing filter functions and subroutines. These functions
!! are used to modify fields in a way that is useful for various
!! simulations.
module filters
  use field, only: field_t
  use num_types, only: rp, sp
  implicit none

  private
  public :: smooth_step_field, permeability_field, step_function_field, PDE_filter

contains

  !> @brief Apply a smooth step function to a field.
  !! @details The smooth step function is defined as:
  !! \f[
  !! t = (x - edge0) / (edge1 - edge0)
  !!  f(t) = \begin{cases}
  !!            t^3 (t (6x - 15) + 10), & t \in [0, 1] \\
  !!              0, & t \leq 0 \\
  !!              1, & t \geq 1 \\
  !!          \end{cases}
  !! \f]
  !! @note The step can be inverted by swapping edge0 and edge1.
  !!
  !! @param[in,out] F Field to be modified.
  !! @param[in] edge0 Edge giving output 0.
  !! @param[in] edge1 Edge giving output 1.
  subroutine smooth_step_field(F, edge0, edge1)
    use filters_cpu, only: smooth_step_cpu

    type(field_t), intent(inout) :: F
    real(kind=rp), intent(in) :: edge0, edge1

    F%x = smooth_step_cpu(F%x, edge0, edge1)
  end subroutine smooth_step_field

  !> @brief Apply a permeability function to a field.
  !! @details The permeability function is defined as:
  !! \f[ k(x) = k_0 + (k_1 - k_0) x \frac{q + 1}{q + x}} \f]
  !! @param[in,out] F Field to be modified.
  !! @param[in] k_0 Permeability at x=0.
  !! @param[in] k_1 Permeability at x=1.
  !! @param[in] q Penalty factor.
  subroutine permeability_field(F_out, x, k_0, k_1, q)
    use filters_cpu, only: permeability_cpu

    type(field_t), intent(inout) :: F_out
    type(field_t), intent(in) :: x
    real(kind=rp), intent(in) :: k_0, k_1
    real(kind=rp), intent(in) :: q

    F_out%x = permeability_cpu(x%x, k_0, k_1, q)
  end subroutine permeability_field

  !> @brief Apply a step function to a field.
  !! @param[in,out] F Field to be modified.
  !! @param[in] x0 Position of the step.
  !! @param[in] value0 Value of the field before the step.
  !! @param[in] value1 Value of the field after the step.
  subroutine step_function_field(F, x0, value0, value1)
    use filters_cpu, only: step_function_cpu

    type(field_t), intent(inout) :: F
    real(kind=rp), intent(in) :: x0, value0, value1

    F%x = step_function_cpu(F%x, x0, value0, value1)
  end subroutine step_function_field

  !> @brief Apply a PDE based filter, see B.S. Lazarov & O. Sigmund 2010,
  ! of the form $\f -r^2 \nabla^2 \tilde{\rho} + \tilde{\rho} = \rho \f$
  !! @param[out] F_out, filtered field $ \f \tilde{\rho}  \f$.
  !! @param[in] F_in, unfiltered field $\f \rho \f$.
  !! @param[in] filter radius $\f r \f$.
  subroutine PDE_filter(F_out, F_in, r, coef)
  ! HARRY
  ! this is annoying... hopefully it will be resolved when I make my PDE filter class
  use coefs, only: coef_t
    use ax_helm_fctry, only: ax_helm_factory
    use ax_product, only : ax_t
    use krylov, only : ksp_t, ksp_monitor_t
    use krylov_fctry, only : krylov_solver_factory
    use precon, only : pc_t
    use bc
    !use bc, only : bc_list_t, bc_list_init, bc_list_apply_scalar
    use neumann, only : neumann_t
    use profiler, only : profiler_start_region, profiler_end_region
    use gather_scatter, only : gs_t, GS_OP_ADD
    use pnpn_residual, only: pnpn_prs_res_t
    use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
    use fld_file_output
    use field_registry, only : neko_field_registry
    use logger, only : neko_log, LOG_SIZE
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    implicit none
    type(field_t), intent(in) ::  F_in
    type(field_t), intent(inout) ::  F_out
    real(kind=rp), intent(in) :: r
    type(coef_t), intent(inout) :: coef

	 ! for now, I'm going allocate everything, then relinguish it.
	 ! if I can get this to work, I'll wrap it up more cleanly into a class with an init etc
    ! Harry

	 ! private
    !> Filter residual
    ! class(pnpn_prs_res_t), allocatable :: filt_res
    !> Ax
     class(ax_t), allocatable :: Ax
    !> Solver results monitors ( filter )
    type(ksp_monitor_t) :: ksp_results(1)
    !> Krylov solver for the filter
    class(ksp_t), allocatable  :: ksp_filt     
    !> Filter Preconditioner
    class(pc_t), allocatable :: pc_filt        
    !> They will all be Neumann conditions.
     type(neumann_t) :: filter_bcs
    !> Filter boundary conditions
    type(bc_list_t) :: bclst_filt              

    character(len=LOG_SIZE) :: log_buf

    ! I want to look at it, delete later
    type(fld_file_output_t) :: fout
    type(field_t), pointer :: fu,fv,fw




    real(kind=rp), dimension(coef%dof%size()) ::  RHS

    real(kind=rp) :: abstol_filt
    integer :: ksp_n, ksp_max_iter, n, i
    character(len=:), allocatable :: ksp_solver, precon_type_filt

	call fout%init(sp, "yofam", 2)

	 ! -----------this will eventually become an init----------!

	 abstol_filt = 0.0000000001_rp
	 ksp_max_iter = 200
	 ksp_solver = "gmres"
	 ! I did some testing...
	 precon_type_filt = "ident"

    n = F_in%dof%size()


	 ! HARRY
	 ! if I'm being honest... I don't think this is doing anything...
	 ! but also, a Neumann condition, in my mind, is a "do nothing" condition
	 !
	 ! so maybe it's ok if I've fucked this up.

	 ! init filter BCs (all Neumann)
	 ! Create list with just Neumann bcs
    call bc_list_init(bclst_filt)
    ! add all the neumann BCs
    call filter_bcs%init_base(coef)
    !call filter_bcs%init_neumann(0.0_rp)
    call filter_bcs%finalize()
    call bc_list_add(bclst_filt, filter_bcs)

	 ! Setup backend dependent Ax routines
    call ax_helm_factory(Ax, .false.)

	 ! set up krylov solver
	 call krylov_solver_factory(ksp_filt, n, ksp_solver, ksp_max_iter, abstol_filt)
	 ! set up preconditioner
	 print *, "setting up precon"
	 call filter_precon_factory(pc_filt, ksp_filt, &
            coef, coef%dof, coef%gs_h, bclst_filt, precon_type_filt)
	 print *, "finished"




	 ! -----------this will eventually become an init----------!

    ! set up Helmholtz operators and RHS
    do i = 1, n
    	 ! note, h1 is already negative in its definition
       coef%h1(i,1,1,1) = r**2
       coef%h2(i,1,1,1) = 1.0_rp
       RHS(i) = F_in%x(i,1,1,1)*coef%B(i,1,1,1)
    end do
    coef%ifh2 = .true.
    ! maybe one ax compute?
    call Ax%compute(F_out%x, RHS, coef, coef%msh, coef%Xh)




	 ! set BCs
    call bc_list_apply_scalar(bclst_filt, RHS, n)

    ! Solve Helmholtz equation
    call profiler_start_region('filter solve')
      ksp_results(1) = &
         ksp_filt%solve(Ax, F_out, RHS, n, coef, bclst_filt, coef%gs_h)

    call profiler_end_region

    ! update preconditioner (needed?)

    call pc_filt%update()

    ! write it all out
    call neko_log%message('Filter')

    write(log_buf, '(A,A,A)') 'Iterations:   ',&
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') ksp_results%iter, &
         ksp_results%res_start, ksp_results%res_final
    call neko_log%message(log_buf)


    ! Check for divergence
    !   if (ieee_is_nan(ksp_results%res_final)) then
    !      call neko_log%error("Filter solver diverged")
    !      stop
    !   end if



    ! take a look

    !call fout%fields%append(F_in)
    !call fout%fields%append(F_out)
    fu => neko_field_registry%get_field('unfiltered_brinkman_indicator')
    fv => neko_field_registry%get_field('brinkman_indicator')
    fout%fields%items(1)%ptr => fu 
    fout%fields%items(2)%ptr => fv 
    call fout%sample(1.0_rp)





  end subroutine PDE_filter


    !> Initialize a Krylov preconditioner
    ! copied from scalar scheme...
    ! I can't imagine why anything should be different here
  subroutine filter_precon_factory(pc, ksp, coef, dof, gs, bclst, pctype)
  ! HARRY
  ! this is annoying... hopefully it will be resolved when I make my PDE filter class
  use coefs, only: coef_t
    use ax_helm_fctry, only: ax_helm_factory
    use ax_product, only : ax_t
    use krylov, only : ksp_t, ksp_monitor_t
    use krylov_fctry, only : krylov_solver_factory
    use precon, only : pc_t
    use bc
    !use bc, only : bc_list_t, bc_list_init, bc_list_apply_scalar
    use neumann, only : neumann_t
    use profiler, only : profiler_start_region, profiler_end_region
    use gather_scatter, only : gs_t, GS_OP_ADD
    use pnpn_residual, only: pnpn_prs_res_t
    use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
    use fld_file_output
    use field_registry, only : neko_field_registry
    use logger, only : neko_log, LOG_SIZE
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    ! extra ones
      use dofmap, only :  dofmap_t
        use jacobi, only : jacobi_t
  use device_jacobi, only : device_jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use hsmg, only : hsmg_t
    use precon_fctry, only : precon_factory, pc_t, precon_destroy
      use utils, only : neko_error

    implicit none
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(inout) :: coef
    type(dofmap_t), target, intent(inout) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=*) :: pctype

    call precon_factory(pc, pctype)

    select type(pcp => pc)
    type is(jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (sx_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (device_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is(hsmg_t)
       if (len_trim(pctype) .gt. 4) then
          if (index(pctype, '+') .eq. 5) then
             call pcp%init(dof%msh, dof%Xh, coef, dof, gs, &
                  bclst, trim(pctype(6:)))
          else
             call neko_error('Unknown coarse grid solver')
          end if
       else
          call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst)
       end if
    end select

    call ksp%set_pc(pc)

  end subroutine filter_precon_factory


end module filters
