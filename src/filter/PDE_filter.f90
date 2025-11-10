! Copyright (c) 2023, The Neko Authors
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
!
!> A PDE based filter

module PDE_filter
  use num_types, only: rp
  use json_module, only: json_file
  use json_utils, only: json_get_or_default, json_get
  use field_registry, only: neko_field_registry
  use field, only: field_t
  use coefs, only: coef_t
  use ax_product, only: ax_t, ax_helm_factory
  use krylov, only: ksp_t, ksp_monitor_t, krylov_solver_factory
  use precon, only: pc_t, precon_factory, precon_destroy
  use bc_list, only : bc_list_t
  use neumann, only: neumann_t
  use profiler, only: profiler_start_region, profiler_end_region
  use gather_scatter, only: gs_t, GS_OP_ADD
  use pnpn_residual, only: pnpn_prs_res_t
  use mesh, only: mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use field_registry, only: neko_field_registry
  use filter, only: filter_t
  use scratch_registry, only: neko_scratch_registry
  use field_math, only: field_copy, field_add3
  use coefs, only: coef_t
  use logger, only: neko_log, LOG_SIZE
  use neko_config, only: NEKO_BCKND_DEVICE
  use dofmap, only: dofmap_t
  use jacobi, only: jacobi_t
  use device_jacobi, only: device_jacobi_t
  use sx_jacobi, only: sx_jacobi_t
  use utils, only: neko_error
  use device_math, only: device_cfill, device_subcol3, device_cmult
  implicit none
  private

  !> A PDE based filter mapping $\rho \mapsto \tilde{\rho}$,
  !! see Lazarov & O. Sigmund 2010,
  !! by solving an equation
  !! of the form $\f -r^2 \nabla^2 \tilde{\rho} + \tilde{\rho} = \rho \f$
  type, public, extends(filter_t) :: PDE_filter_t

     !> Ax
     class(ax_t), allocatable :: Ax
     !> Solver results monitors ( filter )
     type(ksp_monitor_t) :: ksp_results(1)
     !> Krylov solver for the filter
     class(ksp_t), allocatable :: ksp_filt
     !> Filter Preconditioner
     class(pc_t), allocatable :: pc_filt
     !> Filter boundary conditions (they will all be Neumann, so empty)
     type(bc_list_t) :: bclst_filt

     ! Inputs from the user
     !> filter radius
     real(kind=rp) :: r
     !> tolerance for PDE filter
     real(kind=rp) :: abstol_filt
     !> max iterations for PDE filter
     integer :: ksp_max_iter
     !> method for solving PDE
     character(len=:), allocatable :: ksp_solver
     ! > preconditioner type
     character(len=:), allocatable :: precon_type_filt
     integer :: ksp_n, n, i



   contains
     !> Constructor from json.
     procedure, pass(this) :: init => PDE_filter_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_components => &
          PDE_filter_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => PDE_filter_free
     !> Apply the filter
     procedure, pass(this) :: apply => PDE_filter_apply
  end type PDE_filter_t

contains

  !> Constructor from json.
  subroutine PDE_filter_init_from_json(this, json, coef)
    class(PDE_filter_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in) :: coef

    ! user parameters
    call json_get(json, "filter.radius", this%r)

    call json_get_or_default(json, "filter.tolerance", this%abstol_filt, &
         1.0e-10_rp)

    call json_get_or_default(json, "filter.max_iter", this%ksp_max_iter, 200)

    call json_get_or_default(json, "filter.solver", this%ksp_solver, 'cg')

    call json_get_or_default(json, "filter.preconditioner", &
         this%precon_type_filt, 'jacobi')

    call this%init_base(json, coef)
    call PDE_filter_init_from_components(this, coef)

  end subroutine PDE_filter_init_from_json

  !> Actual constructor.
  subroutine PDE_filter_init_from_components(this, coef)
    class(PDE_filter_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer :: n

    n = this%coef%dof%size()

    ! init the bc list (all Neuman BCs, will remain empty)
    call this%bclst_filt%init()

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%Ax, full_formulation = .false.)

    ! set up krylov solver
    call krylov_solver_factory(this%ksp_filt, n, this%ksp_solver, &
         this%ksp_max_iter, this%abstol_filt)

    ! set up preconditioner
    call filter_precon_factory(this%pc_filt, this%ksp_filt, &
         this%coef, this%coef%dof, &
         this%coef%gs_h, &
         this%bclst_filt, this%precon_type_filt)

  end subroutine PDE_filter_init_from_components

  !> Destructor.
  subroutine PDE_filter_free(this)
    class(PDE_filter_t), intent(inout) :: this

    if (allocated(this%Ax)) then
       deallocate(this%Ax)
    end if

    if (allocated(this%ksp_filt)) then
       call this%ksp_filt%free()
       deallocate(this%ksp_filt)
    end if

    if (allocated(this%pc_filt)) then
       call precon_destroy(this%pc_filt)
       deallocate(this%pc_filt)
    end if

    if (allocated(this%ksp_solver)) then
       deallocate(this%ksp_solver)
    end if

    if (allocated(this%precon_type_filt)) then
       deallocate(this%precon_type_filt)
    end if

    call this%bclst_filt%free()

    call this%free_base()

  end subroutine PDE_filter_free

  !> Apply the filter
  !! @param F_out filtered field
  !! @param F_in unfiltered field
  subroutine PDE_filter_apply(this, F_out, F_in)
    class(PDE_filter_t), intent(inout) :: this
    type(field_t), intent(in) :: F_in
    type(field_t), intent(inout) :: F_out
    integer :: n, i
    ! type(field_t), pointer :: RHS
    type(field_t) :: RHS, d_F_out
    character(len=LOG_SIZE) :: log_buf
    ! integer :: temp_indices(1)

    n = this%coef%dof%size()
    ! TODO
    ! This is a bit awkward, because the init for the source terms occurs
    ! before the init of the scratch registry.
    ! So we can't use the scratch registry here.
    ! call neko_scratch_registry%request_field(RHS, temp_indices(1))
    call RHS%init(this%coef%dof)
    call d_F_out%init(this%coef%dof)

    ! in a similar fasion to pressure/velocity, we will solve for d_F_out.

    ! to improve convergence, we use F_in as an initial guess for F_out.
    ! so F_out = F_in + d_F_in.

    ! Defining the operator A = -r^2 \nabla^2 + I
    ! the system changes from:
    ! A (F_out) = F_in
    ! to
    ! A (d_F_out) = F_in - A(F_in)

    ! set up Helmholtz operators and RHS
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! TODO
       ! I think this is correct but I've never tested it
       call device_cfill(this%coef%h1_d, this%r**2, n)
       call device_cfill(this%coef%h2_d, 1.0_rp, n)
    else
       do i = 1, n
          ! h1 is already negative in its definition
          this%coef%h1(i,1,1,1) = this%r**2
          ! ax_helm includes the mass matrix in h2
          this%coef%h2(i,1,1,1) = 1.0_rp
       end do
    end if
    this%coef%ifh2 = .true.

    ! compute the A(F_in) component of the RHS
    ! (note, to be safe with the inout intent we first copy F_in to the
    !  temporary d_F_out)
    call field_copy(d_F_out, F_in)
    call this%Ax%compute(RHS%x, d_F_out%x, this%coef, this%coef%msh, &
         this%coef%Xh)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_subcol3(RHS%x_d, F_in%x_d, this%coef%B_d, n)
       call device_cmult(RHS%x_d, -1.0_rp, n)
    else
       do i = 1, n
          ! mass matrix should be included here
          RHS%x(i,1,1,1) = F_in%x(i,1,1,1) * this%coef%B(i,1,1,1) &
               - RHS%x(i,1,1,1)
       end do
    end if

    ! gather scatter
    call this%coef%gs_h%op(RHS, GS_OP_ADD)

    ! set BCs
    call this%bclst_filt%apply_scalar(RHS%x, n)

    ! Solve Helmholtz equation
    call profiler_start_region('filter solve')
    this%ksp_results(1) = &
         this%ksp_filt%solve(this%Ax, d_F_out, RHS%x, n, this%coef, &
         this%bclst_filt, this%coef%gs_h)

    call profiler_end_region

    ! add result
    call field_add3(F_out, F_in, d_F_out)
    ! update preconditioner (needed?)
    call this%pc_filt%update()

    ! write it all out
    call neko_log%message('Filter')

    write(log_buf, '(A,A,A)') 'Iterations:   ',&
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') this%ksp_results%iter, &
         this%ksp_results%res_start, this%ksp_results%res_final
    call neko_log%message(log_buf)

    !call neko_scratch_registry%relinquish_field(temp_indices)
    call RHS%free()
    call d_F_out%free()

  end subroutine PDE_filter_apply

  !> Initialize a Krylov preconditioner
  subroutine filter_precon_factory(pc, ksp, coef, dof, gs, bclst, &
       pctype)
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(in) :: coef
    type(dofmap_t), target, intent(in) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=*) :: pctype

    call precon_factory(pc, pctype)

    select type (pcp => pc)
    type is (jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (sx_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (device_jacobi_t)
       call pcp%init(coef, dof, gs)
    end select

    call ksp%set_pc(pc)

  end subroutine filter_precon_factory

end module PDE_filter
