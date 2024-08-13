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
! these were copied for lambda2
  use num_types, only : rp
  use json_module, only : json_file
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : lambda2op

! these are the ones I want
    use json_utils, only: json_get, json_get_or_default, json_extract_item
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
    use filter, only: filter_t
    use scratch_registry, only: neko_scratch_registry


    implicit none
  private

  type, public, extends(filter_t) :: PDE_filter_t
    ! these are all the things I want
    ! probably don't need to be public yeh?


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




	 ! these guys come in as inputs from the user
    !> filter radius
    real(kind=rp)  :: r
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
     procedure, pass(this) :: init_from_attributes => &
          PDE_filter_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => PDE_filter_free
     !> Compute the lambda2 field
     procedure, pass(this) :: apply => PDE_filter_apply
  end type PDE_filter_t

contains

  !> Constructor from json.
  subroutine PDE_filter_init_from_json(this, json, coef)
    class(PDE_filter_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout) :: coef
    character(len=20) :: fields(1)
    real(kind=rp) :: tmp_real


	! I'll do the json stuff later...
	 this%r = 0.01 
	 this%abstol_filt = 0.0000000001_rp
    this%ksp_max_iter = 200
    this%ksp_solver = "gmres"
    this%precon_type_filt = "ident"

    ! JSON stuff
    call json_get_or_default(json, 'filter.radius', this%r, 1.0_rp)

    call this%init_base(json, coef)
    call PDE_filter_init_from_attributes(this, coef)
   
  end subroutine PDE_filter_init_from_json

  !> Actual constructor.
  subroutine PDE_filter_init_from_attributes(this, coef)
    class(PDE_filter_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    integer :: n


    ! HARRY
    ! if I'm being honest... I don't think this is doing anything...
    ! but also, a Neumann condition, in my mind, is a "do nothing" condition
    !
    ! so maybe it's ok if I've fucked this up.

    n = this%coef%dof%size()

    ! init filter BCs (all Neumann)
    ! Create list with just Neumann bcs
    call bc_list_init(this%bclst_filt)
    ! add all the neumann BCs
    call this%filter_bcs%init_base(this%coef)
    call this%filter_bcs%init_neumann(0.0_rp)
    call this%filter_bcs%finalize()
    call bc_list_add(this%bclst_filt, this%filter_bcs)

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%Ax)

    ! set up krylov solver
    call krylov_solver_factory(this%ksp_filt, n, this%ksp_solver, this%ksp_max_iter, this%abstol_filt)
    ! set up preconditioner
    call filter_precon_factory(this%pc_filt, this%ksp_filt, &
            this%coef, this%coef%dof, this%coef%gs_h, this%bclst_filt, this%precon_type_filt)

  end subroutine PDE_filter_init_from_attributes

  !> Destructor.
  subroutine PDE_filter_free(this)
    class(PDE_filter_t), intent(inout) :: this
    call this%free_base()
  end subroutine PDE_filter_free

  !> Apply the filter
  !! @param F_out filtered field
  !! @param F_in unfiltered field
  subroutine PDE_filter_apply(this, F_out, F_in)
    class(PDE_filter_t), intent(inout) :: this
    type(field_t), intent(in) ::  F_in
    type(field_t), intent(inout) ::  F_out
    integer :: n, i
    type(field_t), pointer ::  RHS 
    integer :: temp_indices
    type(field_t), pointer :: ta1, ta2, ta3
    character(len=LOG_SIZE) :: log_buf

    n = this%coef%dof%size()

    
    ! you should be using the scratch registry!!!
    ! but for some reason i segfault here if I try
    !call neko_scratch_registry%request_field(RHS, temp_indices(1))
    call neko_field_registry%add_field(this%coef%dof, 'RHS')
    RHS => neko_field_registry%get_field('RHS')


    ! set up Helmholtz operators and RHS
    do i = 1, n
       ! h1 is already negative in its definition
       this%coef%h1(i,1,1,1) = this%r**2
       ! ax_helm includes the mass matrix in h2
       this%coef%h2(i,1,1,1) = 1.0_rp
       ! mass matrix should be included here
       RHS%x(i,1,1,1) = F_in%x(i,1,1,1)*this%coef%B(i,1,1,1)
    end do
    this%coef%ifh2 = .true.

    ! set BCs
    call bc_list_apply_scalar(this%bclst_filt, RHS%x, n)

	 ! gather scatter
    call this%coef%gs_h%op(RHS, GS_OP_ADD)

    ! Solve Helmholtz equation
    call profiler_start_region('filter solve')
      this%ksp_results(1) = &
         this%ksp_filt%solve(this%Ax, F_out, RHS%x, n, this%coef, this%bclst_filt, this%coef%gs_h)

    call profiler_end_region

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


	 ! YOU SHOULD RELINGUISH IF YOU GET THIS TO WORK!
    !call neko_scratch_registry%relinquish_field(temp_indices)



  end subroutine PDE_filter_apply

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


end module PDE_filter
