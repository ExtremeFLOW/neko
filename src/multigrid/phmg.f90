! Copyright (c) 2024-2026, The Neko Authors
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
!> Hybrid ph-multigrid preconditioner
module phmg
  use num_types, only : rp
  use precon, only : pc_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use space, only : space_t, GLL
  use dofmap, only : dofmap_t
  use field, only : field_t
  use coefs, only : coef_t, COEF_FULL, COEF_OPERATOR
  use mesh, only : mesh_t
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use scalar_bc_projector, only : scalar_bc_projector_t
  use dirichlet, only : dirichlet_t
  use utils, only : neko_error, neko_warning
  use cheby, only : cheby_t
  use cheby_device, only : cheby_device_t
  use jacobi, only : jacobi_t
  use device_jacobi, only : device_jacobi_t
  use schwarz, only : schwarz_t
  use ax_product, only : ax_t, ax_helm_allocator
  use tree_amg_multigrid, only : tamg_solver_t
  use interpolation, only : interpolator_t
  use json_module, only : json_file
  use json_utils, only : json_get_or_default, json_get
  use math, only : copy, col2, add2, add2s2, add2s1
  use device, only : device_get_ptr, device_stream_wait_event, glb_cmd_queue, &
       glb_cmd_event
  use device_math, only : device_rzero, device_copy, device_add2, &
       device_add2s2, device_invcol2, device_glsc2, device_col2, device_add2s1
  use neko_config, only : NEKO_BCKND_DEVICE
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER, &
       krylov_solver_factory
  use profiler, only : profiler_start_region, profiler_end_region
  use logger, only : neko_log, LOG_SIZE
  use, intrinsic :: iso_c_binding
  implicit none
  private


  type, private :: phmg_lvl_t
     integer :: lvl = -1
     integer :: smoother_itrs = 10
     type(space_t), pointer :: Xh
     type(dofmap_t), pointer :: dm_Xh
     type(gs_t), pointer :: gs_h
     type(schwarz_t) :: schwarz
     type(cheby_t) :: cheby
     type(cheby_device_t) :: cheby_device
     type(jacobi_t) :: jacobi
     type(device_jacobi_t) :: device_jacobi
     type(coef_t), pointer :: coef
     type(scalar_bc_projector_t) :: bc_projector
     type(dirichlet_t) :: bc
     type(field_t) :: r, w, z
  end type phmg_lvl_t

  type, public :: phmg_hrchy_t
     type(phmg_lvl_t), allocatable :: lvl(:)
  end type phmg_hrchy_t


  type, public, extends(pc_t) :: phmg_t
     type(tamg_solver_t) :: amg_solver
     integer :: nlvls
     type(phmg_hrchy_t) :: phmg_hrchy
     class(ax_t), allocatable :: ax
     type(interpolator_t), allocatable :: intrp(:)
     !> Coordinate interpolators.
     type(interpolator_t), allocatable :: crd_intrp(:)
     type(mesh_t), pointer :: msh
     !> Fine-level metrics_version at the last refresh, to detect a change.
     integer :: last_metrics_version = -1
     !> When .false., update() does nothing and coarse levels stay frozen
     !> at their initial geometry. When .true., update() interpolates
     !> the fine-level geometry to the coarse levels.
     logical :: update_enabled = .false.
     !> Whether a refresh also re-estimates the Chebyshev eigenvalues.
     logical :: refresh_eigs = .true.
     !> Re-estimate eigenvalues only every N-th refresh.
     integer :: refresh_eigs_frequency = 20
     !> Restart eigenvalue estimation from the previous eigenvector.
     logical :: eigs_warm_start = .true.
     !> Power iterations used for warm-started re-estimations
     integer :: power_its_refresh = 20
     !> Refreshes done so far, for refresh_eigs_frequency
     integer :: n_refresh = 0
     !> Which accelerator the Chebyshev smoother uses
     character(len=:), allocatable :: cheby_acc
   contains
     procedure, pass(this) :: init => phmg_init
     procedure, pass(this) :: init_from_components => &
          phmg_init_from_components
     procedure, pass(this) :: free => phmg_free
     procedure, pass(this) :: solve => phmg_solve
     procedure, pass(this) :: update => phmg_update
     procedure, private, pass(this) :: mg_cycle => phmg_mg_cycle
  end type phmg_t

contains

  subroutine phmg_init(this, coef, bclst, phmg_params)
    class(phmg_t), intent(inout), target :: this
    type(coef_t), intent(in), target :: coef
    type(bc_list_t), intent(inout), target :: bclst
    type(json_file), intent(inout) :: phmg_params
    integer :: crs_tamg_lvls, crs_tamg_itrs, crs_tamg_cheby_degree
    integer :: smoother_itrs
    character(len=:), allocatable :: cheby_acc
    integer, allocatable :: pcrs_sched(:)
    logical :: update_enabled

    call json_get_or_default(phmg_params, 'smoother_iterations', &
         smoother_itrs, 3)

    call json_get_or_default(phmg_params, 'smoother_cheby_acc', &
         cheby_acc, "jacobi")

    call json_get_or_default(phmg_params, 'coarse_grid.levels', &
         crs_tamg_lvls, 3)

    call json_get_or_default(phmg_params, 'coarse_grid.iterations', &
         crs_tamg_itrs, 1)

    call json_get_or_default(phmg_params, 'coarse_grid.cheby_degree', &
         crs_tamg_cheby_degree, 4)

    if (phmg_params%valid_path('pcoarsening_schedule')) then
       call json_get(phmg_params, 'pcoarsening_schedule', pcrs_sched)
    else
       allocate(pcrs_sched(2))
       pcrs_sched(1) = 3
       pcrs_sched(2) = 1
    end if

    call json_get_or_default(phmg_params, 'update.enabled', &
         update_enabled, .false.)

    ! Control the eigenvalue re-estimation; geometry always refreshes.
    call json_get_or_default(phmg_params, 'update.eigs.enabled', &
         this%refresh_eigs, .true.)
    call json_get_or_default(phmg_params, 'update.eigs.frequency', &
         this%refresh_eigs_frequency, 20)
    call json_get_or_default(phmg_params, 'update.eigs.warm_start', &
         this%eigs_warm_start, .true.)
    call json_get_or_default(phmg_params, 'update.eigs.warm_start_iterations', &
         this%power_its_refresh, 20)

    call this%init_from_components(coef, bclst, smoother_itrs, &
         cheby_acc, crs_tamg_lvls, crs_tamg_itrs, crs_tamg_cheby_degree,&
         pcrs_sched, update_enabled)

  end subroutine phmg_init

  subroutine phmg_init_from_components(this, coef, bclst, smoother_itrs, &
       cheby_acc, crs_tamg_lvls, crs_tamg_itrs, crs_tamg_cheby_degree, &
       pcrs_sched, update_enabled)
    class(phmg_t), intent(inout), target :: this
    type(coef_t), intent(in), target :: coef
    type(bc_list_t), intent(inout), target :: bclst
    integer, intent(in) :: smoother_itrs
    character(len=:), allocatable :: cheby_acc
    integer, intent(in) :: crs_tamg_lvls, crs_tamg_itrs
    integer, intent(in) :: crs_tamg_cheby_degree
    integer, intent(in), allocatable :: pcrs_sched(:)
    logical, intent(in), optional :: update_enabled
    integer :: lx_crs, lx_mid
    integer, allocatable :: lx_lvls(:)
    integer :: n, i, j, st
    class(bc_t), pointer :: bc_j
    logical :: use_jacobi, use_cheby
    use_jacobi = .true.
    use_cheby = .true.

    this%msh => coef%msh
    this%cheby_acc = trim(cheby_acc)

    if (present(update_enabled)) then
       this%update_enabled = update_enabled
    else
       this%update_enabled = .false.
    end if

    this%nlvls = size(pcrs_sched) + 1
    allocate(lx_lvls(0:this%nlvls - 1))
    lx_lvls(1:) = pcrs_sched + 1

    allocate(this%phmg_hrchy%lvl(0:this%nlvls - 1))

    this%phmg_hrchy%lvl(0)%lvl = 0
    this%phmg_hrchy%lvl(0)%smoother_itrs = smoother_itrs
    this%phmg_hrchy%lvl(0)%Xh => coef%Xh
    this%phmg_hrchy%lvl(0)%coef => coef
    this%phmg_hrchy%lvl(0)%dm_Xh => coef%dof
    this%phmg_hrchy%lvl(0)%gs_h => coef%gs_h

    do i = 1, this%nlvls - 1
       allocate(this%phmg_hrchy%lvl(i)%Xh)
       allocate(this%phmg_hrchy%lvl(i)%dm_Xh)
       allocate(this%phmg_hrchy%lvl(i)%gs_h)
       allocate(this%phmg_hrchy%lvl(i)%coef)

       this%phmg_hrchy%lvl(i)%lvl = i
       this%phmg_hrchy%lvl(i)%smoother_itrs = smoother_itrs
       call this%phmg_hrchy%lvl(i)%Xh%init(GLL, lx_lvls(i), lx_lvls(i), &
            lx_lvls(i))
       call this%phmg_hrchy%lvl(i)%dm_Xh%init(coef%msh, &
            this%phmg_hrchy%lvl(i)%Xh)
       call this%phmg_hrchy%lvl(i)%gs_h%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       ! A coarse level reads G_ij, h1, h2, B and mult and nothing else, so
       ! retaining the rest costs memory for no purpose. The exception is a
       ! moving mesh, where update() has to rebuild the geometry from the
       ! derivative arrays.
       if (this%update_enabled) then
          call this%phmg_hrchy%lvl(i)%coef%init( &
               this%phmg_hrchy%lvl(i)%gs_h, COEF_FULL)
       else
          call this%phmg_hrchy%lvl(i)%coef%init( &
               this%phmg_hrchy%lvl(i)%gs_h, COEF_OPERATOR)
       end if
    end do

    do i = 0, this%nlvls - 1
       call this%phmg_hrchy%lvl(i)%r%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       call this%phmg_hrchy%lvl(i)%w%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       call this%phmg_hrchy%lvl(i)%z%init(this%phmg_hrchy%lvl(i)%dm_Xh)

       this%phmg_hrchy%lvl(i)%coef%ifh2 = coef%ifh2
       call copy(this%phmg_hrchy%lvl(i)%coef%h1, coef%h1, &
            this%phmg_hrchy%lvl(i)%dm_Xh%size())

       call this%phmg_hrchy%lvl(i)%bc%init_base(this%phmg_hrchy%lvl(i)%coef)
       if (bclst%size() .gt. 0 ) then
          do j = 1, bclst%size()
             bc_j => bclst%get(j)
             call this%phmg_hrchy%lvl(i)%bc%mark_facets(bc_j%marked_facet)
          end do
       end if
       call this%phmg_hrchy%lvl(i)%bc%finalize()
       call this%phmg_hrchy%lvl(i)%bc%set_g(0.0_rp)
       call this%phmg_hrchy%lvl(i)%bc_projector%mark(this%phmg_hrchy%lvl(i)%bc)

       !> Initialize Smoothers
       if (trim(cheby_acc) .eq. "schwarz") then
          call this%phmg_hrchy%lvl(i)%schwarz%init( &
               this%phmg_hrchy%lvl(i)%Xh, &
               this%phmg_hrchy%lvl(i)%dm_Xh, &
               this%phmg_hrchy%lvl(i)%gs_h, &
               this%phmg_hrchy%lvl(i)%bc_projector, &
               coef%msh)
       end if

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call this%phmg_hrchy%lvl(i)%device_jacobi%init(&
               this%phmg_hrchy%lvl(i)%coef, &
               this%phmg_hrchy%lvl(i)%dm_Xh, &
               this%phmg_hrchy%lvl(i)%gs_h)
       else
          call this%phmg_hrchy%lvl(i)%jacobi%init(&
               this%phmg_hrchy%lvl(i)%coef, &
               this%phmg_hrchy%lvl(i)%dm_Xh, &
               this%phmg_hrchy%lvl(i)%gs_h)
       end if

       if (NEKO_BCKND_DEVICE .eq. 1) then
          if (trim(cheby_acc) .eq. "jacobi") then
             call this%phmg_hrchy%lvl(i)%cheby_device%init( &
                  this%phmg_hrchy%lvl(i)%dm_Xh%size(), smoother_itrs, &
                  this%phmg_hrchy%lvl(i)%device_jacobi)
             st = 1
          else
             call this%phmg_hrchy%lvl(i)%cheby_device%init( &
                  this%phmg_hrchy%lvl(i)%dm_Xh%size(), smoother_itrs)
             st = 0
             if (trim(cheby_acc) .eq. "schwarz") then
                this%phmg_hrchy%lvl(i)%cheby_device%schwarz => &
                     this%phmg_hrchy%lvl(i)%schwarz
                st = 2
             end if
          end if
       else
          if (trim(cheby_acc) .eq. "jacobi") then
             call this%phmg_hrchy%lvl(i)%cheby%init( &
                  this%phmg_hrchy%lvl(i)%dm_Xh%size(), smoother_itrs, &
                  this%phmg_hrchy%lvl(i)%jacobi)
             st = 1
          else
             call this%phmg_hrchy%lvl(i)%cheby%init( &
                  this%phmg_hrchy%lvl(i)%dm_Xh%size(), smoother_itrs)
             st = 0
             if (trim(cheby_acc) .eq. "schwarz") then
                this%phmg_hrchy%lvl(i)%cheby%schwarz => &
                     this%phmg_hrchy%lvl(i)%schwarz
                st = 2
             end if
          end if
       end if

    end do

    call print_phmg_info(this%nlvls, st, this%phmg_hrchy)

    ! Create backend specific Ax operator
    call ax_helm_allocator(this%ax, type_name = "standard")

    ! Interpolator Fine + mg levels
    allocate(this%intrp(this%nlvls - 1))
    do i = 1, this%nlvls -1
       call this%intrp(i)%init(this%phmg_hrchy%lvl(i-1)%Xh, &
            this%phmg_hrchy%lvl(i)%Xh)
    end do

    ! Coarse space first. Each level maps from level 0.
    if (this%update_enabled) then
       allocate(this%crd_intrp(this%nlvls - 1))
       do i = 1, this%nlvls - 1
          call this%crd_intrp(i)%init(this%phmg_hrchy%lvl(i)%Xh, &
               this%phmg_hrchy%lvl(0)%Xh)
       end do
    end if

    call this%amg_solver%init(this%ax, this%phmg_hrchy%lvl(this%nlvls -1)%Xh, &
         this%phmg_hrchy%lvl(this%nlvls -1)%coef, this%msh, &
         this%phmg_hrchy%lvl(this%nlvls-1)%gs_h, crs_tamg_lvls, &
         this%phmg_hrchy%lvl(this%nlvls -1)%bc_projector, &
         crs_tamg_itrs, crs_tamg_cheby_degree)

    ! update() only refreshes when `lvl(0)%coef%metrics_version` changes.
    this%last_metrics_version = this%phmg_hrchy%lvl(0)%coef%metrics_version

    ! Hand the eigenvalue policy to every smoother.
    if (this%update_enabled) then
       do i = 0, this%nlvls - 1
          this%phmg_hrchy%lvl(i)%cheby%warm_start_eigs = this%eigs_warm_start
          this%phmg_hrchy%lvl(i)%cheby%power_its_refresh = &
               this%power_its_refresh
          this%phmg_hrchy%lvl(i)%cheby_device%warm_start_eigs = &
               this%eigs_warm_start
          this%phmg_hrchy%lvl(i)%cheby_device%power_its_refresh = &
               this%power_its_refresh
       end do
       call this%amg_solver%set_eig_refresh(this%eigs_warm_start, &
            this%power_its_refresh)

       if (trim(cheby_acc) .eq. "schwarz") then
          call neko_warning("PHMG: the Schwarz smoother is not refreshed " // &
               "when the mesh changes. Its local solves stay at the " // &
               "initial geometry.")
       end if
    end if

  end subroutine phmg_init_from_components

  subroutine phmg_free(this)
    class(phmg_t), intent(inout) :: this
    integer :: i

    call this%amg_solver%free()

    if (allocated(this%intrp)) then
       do i = 1, size(this%intrp)
          call this%intrp(i)%free()
       end do
       deallocate(this%intrp)
    end if

    if (allocated(this%crd_intrp)) then
       do i = 1, size(this%crd_intrp)
          call this%crd_intrp(i)%free()
       end do
       deallocate(this%crd_intrp)
    end if

    if (allocated(this%ax)) then
       deallocate(this%ax)
    end if

    if (allocated(this%phmg_hrchy%lvl)) then
       do i = lbound(this%phmg_hrchy%lvl, 1), ubound(this%phmg_hrchy%lvl, 1)
          call this%phmg_hrchy%lvl(i)%r%free()
          call this%phmg_hrchy%lvl(i)%w%free()
          call this%phmg_hrchy%lvl(i)%z%free()

          call this%phmg_hrchy%lvl(i)%cheby%free()
          call this%phmg_hrchy%lvl(i)%cheby_device%free()
          call this%phmg_hrchy%lvl(i)%jacobi%free()
          call this%phmg_hrchy%lvl(i)%device_jacobi%free()

          if (allocated(this%phmg_hrchy%lvl(i)%schwarz%work1)) then
             call this%phmg_hrchy%lvl(i)%schwarz%free()
          end if

          call this%phmg_hrchy%lvl(i)%bc_projector%free()
          call this%phmg_hrchy%lvl(i)%bc%free()

          ! Level 0 borrows Xh, dm_Xh, gs_h and coef from the caller,
          ! all other levels own them
          if (this%phmg_hrchy%lvl(i)%lvl .gt. 0) then
             call this%phmg_hrchy%lvl(i)%coef%free()
             call this%phmg_hrchy%lvl(i)%gs_h%free()
             call this%phmg_hrchy%lvl(i)%dm_Xh%free()
             call this%phmg_hrchy%lvl(i)%Xh%free()
             deallocate(this%phmg_hrchy%lvl(i)%coef)
             deallocate(this%phmg_hrchy%lvl(i)%gs_h)
             deallocate(this%phmg_hrchy%lvl(i)%dm_Xh)
             deallocate(this%phmg_hrchy%lvl(i)%Xh)
          end if

          nullify(this%phmg_hrchy%lvl(i)%coef)
          nullify(this%phmg_hrchy%lvl(i)%gs_h)
          nullify(this%phmg_hrchy%lvl(i)%dm_Xh)
          nullify(this%phmg_hrchy%lvl(i)%Xh)
       end do
       deallocate(this%phmg_hrchy%lvl)
    end if

    nullify(this%msh)

  end subroutine phmg_free

  subroutine phmg_solve(this, z, r, n)
    class(phmg_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(c_ptr) :: z_d, r_d
    type(ksp_monitor_t) :: ksp_results
    integer :: i

    call profiler_start_region('PHMG_solve', 8)
    associate( mglvl => this%phmg_hrchy%lvl)
      if (NEKO_BCKND_DEVICE .eq. 1) then
         z_d = device_get_ptr(z)
         r_d = device_get_ptr(r)
         !We should not work with the input
         call device_copy(mglvl(0)%r%x_d, r_d, n)

         call device_rzero(mglvl(0)%z%x_d, n)
         call device_rzero(mglvl(0)%w%x_d, n)

         call this%mg_cycle()

         call device_copy(z_d, mglvl(0)%z%x_d, n)
      else
         !OCL NORECURRENCE, NOVREC, NOALIAS
         !DIR$ CONCURRENT
         !DIR$ IVDEP
         !GCC$ ivdep
         !$omp parallel do
         do i = 1, n
            !We should not work with the input
            mglvl(0)%r%x(i,1,1,1) = r(i)

            mglvl(0)%z%x(i,1,1,1) = 0.0_rp
            mglvl(0)%w%x(i,1,1,1) = 0.0_rp
         end do
         !$omp end parallel do

         call this%mg_cycle()

         !OCL NORECURRENCE, NOVREC, NOALIAS
         !DIR$ CONCURRENT
         !DIR$ IVDEP
         !GCC$ ivdep
         !$omp parallel do
         do i = 1, n
            z(i) = mglvl(0)%z%x(i,1,1,1)
         end do
         !$omp end parallel do
      end if
    end associate
    call profiler_end_region('PHMG_solve', 8)

  end subroutine phmg_solve

  !> Bring the preconditioner back in sync after the mesh has changed.
  subroutine phmg_update(this)
    class(phmg_t), intent(inout) :: this
    integer :: fine_version
    logical :: do_eigs

    if (.not. this%update_enabled) return

    fine_version = this%phmg_hrchy%lvl(0)%coef%metrics_version

    ! Mesh has not changed since the last refresh.
    if (fine_version .eq. this%last_metrics_version) return

    call phmg_update_coarse_geometry(this)

    call phmg_update_smoother_acc(this)

    do_eigs = (this%refresh_eigs) .and. &
         (this%refresh_eigs_frequency .gt. 0) .and. &
         (mod(this%n_refresh, max(this%refresh_eigs_frequency, 1)) &
         .eq. 0)
    if (do_eigs) call phmg_update_smoother_eigs(this)

    this%n_refresh = this%n_refresh + 1
    this%last_metrics_version = fine_version

  end subroutine phmg_update


  !> Sample the new fine coordinates on the coarse levels and rebuild their
  !! metrics. Level 0 shares the caller's dofmap and coef, so it is skipped.
  subroutine phmg_update_coarse_geometry(this)
    class(phmg_t), intent(inout) :: this
    integer :: i

    if (.not. allocated(this%crd_intrp)) then
       call neko_error("PHMG: update requested but coordinate " // &
            "interpolators were not initialized.")
    end if

    call profiler_start_region('PHMG_update_geometry')
    associate (mg => this%phmg_hrchy%lvl, nelv => this%msh%nelv)
      do i = 1, this%nlvls - 1
         call this%crd_intrp(i)%map(mg(i)%dm_Xh%x, mg(0)%dm_Xh%x, &
              nelv, mg(i)%Xh)
         call this%crd_intrp(i)%map(mg(i)%dm_Xh%y, mg(0)%dm_Xh%y, &
              nelv, mg(i)%Xh)
         call this%crd_intrp(i)%map(mg(i)%dm_Xh%z, mg(0)%dm_Xh%z, &
              nelv, mg(i)%Xh)

         call mg(i)%coef%recompute_metrics()
      end do
    end associate
    call profiler_end_region('PHMG_update_geometry')

  end subroutine phmg_update_coarse_geometry


  !> Rebuild the accelerator the Chebyshev smoother uses, which depends on
  !! the geometry.
  subroutine phmg_update_smoother_acc(this)
    class(phmg_t), intent(inout) :: this
    integer :: i

    call profiler_start_region('PHMG_update_smoother_acc')
    associate (mg => this%phmg_hrchy%lvl)
      select case (trim(this%cheby_acc))
      case ("jacobi")
         do i = 0, this%nlvls - 1
            if (NEKO_BCKND_DEVICE .eq. 1) then
               call mg(i)%device_jacobi%update()
            else
               call mg(i)%jacobi%update()
            end if
         end do
      case ("schwarz")
         ! Not refreshed. The local solves are built by the FDM from the
         ! element coordinates, and fdm_t has no update path.
      case default
         ! No accelerator to refresh.
      end select
    end associate
    call profiler_end_region('PHMG_update_smoother_acc')

  end subroutine phmg_update_smoother_acc


  !> Mark every smoother to re-estimate its eigenvalues on the next solve.
  subroutine phmg_update_smoother_eigs(this)
    class(phmg_t), intent(inout) :: this
    integer :: i

    associate (mg => this%phmg_hrchy%lvl)
      do i = 0, this%nlvls - 1
         if (NEKO_BCKND_DEVICE .eq. 1) then
            mg(i)%cheby_device%recompute_eigs = .true.
         else
            mg(i)%cheby%recompute_eigs = .true.
         end if
      end do
    end associate

    call this%amg_solver%invalidate_eigs()

  end subroutine phmg_update_smoother_eigs


  subroutine phmg_mg_cycle(this)
    class(phmg_t), intent(inout) :: this
    type(ksp_monitor_t) :: ksp_results
    character(len=2) :: lvl_name
    integer :: lvl, i

    associate(mg => this%phmg_hrchy%lvl, intrp => this%intrp, &
         msh => this%msh, Ax => this%Ax)
      do lvl = 0, this%nlvls-2
         write(lvl_name, '(I0)') lvl
         call profiler_start_region( "PHMG_level_" // trim(lvl_name))
         associate(z => mg(lvl)%z, r => mg(lvl)%r, w => mg(lvl)%w)
           !------------!
           !   SMOOTH   !
           !------------!
           if (NEKO_BCKND_DEVICE .eq. 1) then
              mg(lvl)%cheby_device%zero_initial_guess = .true.
              ksp_results = mg(lvl)%cheby_device%solve(Ax, z, &
                   r%x, mg(lvl)%dm_Xh%size(), &
                   mg(lvl)%coef, mg(lvl)%bc_projector, &
                   mg(lvl)%gs_h, niter = mg(lvl)%smoother_itrs)
           else
              mg(lvl)%cheby%zero_initial_guess = .true.
              ksp_results = mg(lvl)%cheby%solve(Ax, z, &
                   r%x, mg(lvl)%dm_Xh%size(), &
                   mg(lvl)%coef, mg(lvl)%bc_projector, &
                   mg(lvl)%gs_h, niter = mg(lvl)%smoother_itrs)
           end if

           !------------!
           !  Residual  !
           !------------!
           call Ax%compute(w%x, z%x, mg(lvl)%coef, msh, mg(lvl)%Xh)
           call mg(lvl)%gs_h%op(w%x, mg(lvl)%dm_Xh%size(), GS_OP_ADD, &
                glb_cmd_event)
           call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)
           call mg(lvl)%bc_projector%apply(w%x, mg(lvl)%dm_Xh%size())

           if (NEKO_BCKND_DEVICE .eq. 1) then
              call device_add2s1(w%x_d, r%x_d, -1.0_rp, mg(lvl)%dm_Xh%size())
           else
              !OCL NORECURRENCE, NOVREC, NOALIAS
              !DIR$ CONCURRENT
              !DIR$ IVDEP
              !GCC$ ivdep
              !$omp parallel do
              do i = 1, mg(lvl)%dm_Xh%size()
                 w%x(i,1,1,1) = r%x(i,1,1,1) - w%x(i,1,1,1)
              end do
              !$omp end parallel do
           end if

           !------------!
           !  Restrict  !
           !------------!
           if (NEKO_BCKND_DEVICE .eq. 1) then
              call device_col2(w%x_d, mg(lvl)%coef%mult_d, mg(lvl)%dm_Xh%size())
           else
              call col2(w%x, mg(lvl)%coef%mult, mg(lvl)%dm_Xh%size())
           end if

           call intrp(lvl+1)%map(mg(lvl+1)%r%x, w%x, msh%nelv, mg(lvl+1)%Xh)

           call mg(lvl+1)%gs_h%op(mg(lvl+1)%r%x, mg(lvl+1)%dm_Xh%size(), &
                GS_OP_ADD, glb_cmd_event)
           call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)

           call mg(lvl+1)%bc_projector%apply( &
                mg(lvl+1)%r%x, &
                mg(lvl+1)%dm_Xh%size())

           if (NEKO_BCKND_DEVICE .eq. 1) then
              call device_rzero(mg(lvl+1)%z%x_d, mg(lvl+1)%dm_Xh%size())
           else
              !OCL NORECURRENCE, NOVREC, NOALIAS
              !DIR$ CONCURRENT
              !DIR$ IVDEP
              !GCC$ ivdep
              !$omp parallel do
              do i = 1, mg(lvl+1)%dm_Xh%size()
                 mg(lvl+1)%z%x(i,1,1,1) = 0.0_rp
              end do
              !$omp end parallel do
           end if
         end associate
         call profiler_end_region( "PHMG_level_" // trim(lvl_name))
      end do

      call profiler_start_region( 'PHMG_coarse-solve' )
      !------------!
      !   SOLVE    !
      !------------!
      call this%amg_solver%solve(mg(this%nlvls-1)%z%x, &
           mg(this%nlvls-1)%r%x, &
           mg(this%nlvls-1)%dm_Xh%size())
      call profiler_end_region( 'PHMG_coarse-solve' )

      do lvl = (this%nlvls-2), 0, -1
         write(lvl_name, '(I0)') lvl
         call profiler_start_region( "PHMG_level_" // trim(lvl_name))
         associate(z => mg(lvl)%z, r => mg(lvl)%r, w => mg(lvl)%w)
           !------------!
           !  Project   !
           !------------!
           call intrp(lvl+1)%map(w%x, mg(lvl+1)%z%x, msh%nelv, mg(lvl)%Xh)

           call mg(lvl)%gs_h%op(w%x, mg(lvl)%dm_Xh%size(), GS_OP_ADD, &
                glb_cmd_event)
           call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)

           if (NEKO_BCKND_DEVICE .eq. 1) then
              call device_col2(w%x_d, mg(lvl)%coef%mult_d, mg(lvl)%dm_Xh%size())
           else
              call col2(w%x, mg(lvl)%coef%mult, mg(lvl)%dm_Xh%size())
           end if

           !------------!
           !  Correct   !
           !------------!
           if (NEKO_BCKND_DEVICE .eq. 1) then
              call device_add2(z%x_d, w%x_d, mg(lvl)%dm_Xh%size())
           else
              !OCL NORECURRENCE, NOVREC, NOALIAS
              !DIR$ CONCURRENT
              !DIR$ IVDEP
              !GCC$ ivdep
              !$omp parallel do
              do i = 1, mg(lvl)%dm_Xh%size()
                 z%x(i,1,1,1) = z%x(i,1,1,1) + w%x(i,1,1,1)
              end do
              !$omp end parallel do
           end if

           !------------!
           !   SMOOTH   !
           !------------!
           if (NEKO_BCKND_DEVICE .eq. 1) then
              ksp_results = mg(lvl)%cheby_device%solve(Ax, z, &
                   r%x, mg(lvl)%dm_Xh%size(), &
                   mg(lvl)%coef, mg(lvl)%bc_projector, &
                   mg(lvl)%gs_h, niter = mg(lvl)%smoother_itrs)
           else
              ksp_results = mg(lvl)%cheby%solve(Ax, z, &
                   r%x, mg(lvl)%dm_Xh%size(), &
                   mg(lvl)%coef, mg(lvl)%bc_projector, &
                   mg(lvl)%gs_h, niter = mg(lvl)%smoother_itrs)
           end if
         end associate
         call profiler_end_region( "PHMG_level_" // trim(lvl_name))
      end do
    end associate

  end subroutine phmg_mg_cycle

  !> Wraps jacobi solve as a residual update relaxation method
  !! @param z solution vector (inout)
  !! @param r rhs vector
  !! @param w scratch work space vector
  !! @param mg phmg level from hierarchy
  !! @param msh mesh
  !! @param Ax matrix vector class object
  !! @param n vector length
  !! @param lvl not used
  subroutine phmg_jacobi_smoother(z, r, w, mg, msh, Ax, n, lvl)
    type(phmg_lvl_t) :: mg
    class(ax_t), intent(inout) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(field_t), intent(inout) :: z, r, w
    integer, intent(in) :: n, lvl
    integer :: i, j, iblk, ni, niblk

    ni = mg%smoother_itrs
    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, ni
          call Ax%compute(w%x, z%x, mg%coef, msh, mg%Xh)
          call mg%gs_h%op(w%x, n, GS_OP_ADD, glb_cmd_event)
          call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)
          call mg%bc_projector%apply(w%x, n)
          call device_add2s1(w%x_d, r%x_d, -1.0_rp, n)

          call mg%device_jacobi%solve(w%x, w%x, n)

          call device_add2s2(z%x_d, w%x_d, 0.6_rp, n)
       end do
    else
       do i = 1, ni
          call Ax%compute(w%x, z%x, mg%coef, msh, mg%Xh)
          call mg%gs_h%op(w%x, n, GS_OP_ADD)
          call mg%bc_projector%apply(w%x, n)
          call add2s1(w%x, r%x, -1.0_rp, n)

          call mg%jacobi%solve(w%x, w%x, n)

          call add2s2(z%x, w%x, 0.6_rp, n)
       end do
    end if
  end subroutine phmg_jacobi_smoother


  subroutine phmg_resid_monitor(z, r, w, mg, msh, Ax, lvl, typ)
    integer :: lvl, typ
    type(phmg_lvl_t) :: mg
    class(ax_t), intent(inout) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(field_t) :: z, r, w
    real(kind=rp) :: val
    character(len=LOG_SIZE) :: log_buf
    call Ax%compute(w%x, z%x, mg%coef, msh, mg%Xh)
    call mg%gs_h%op(w%x, mg%dm_Xh%size(), GS_OP_ADD)
    call mg%bc_projector%apply(w%x, mg%dm_Xh%size())
    call device_add2s1(w%x_d, r%x_d, -1.0_rp, mg%dm_Xh%size())
    val = device_glsc2(w%x_d, w%x_d, mg%dm_Xh%size())
    if (typ .eq. 1) then
       write(log_buf, '(A15,I4,F12.6)') 'PRESMOO - PRE', lvl, val
    else if (typ .eq. 2) then
       write(log_buf, '(A15,I4,F12.6)') 'PRESMOO -POST', lvl, val
    else if (typ .eq. 3) then
       write(log_buf, '(A15,I4,F12.6)') 'POSTSMOO- PRE', lvl, val
    else if (typ .eq. 4) then
       write(log_buf, '(A15,I4,F12.6)') 'POSTSMOO-POST', lvl, val
    else if (typ .eq. 5) then
       write(log_buf, '(A15,I4,F12.6)') 'TAMG - PRE', lvl, val
    else if (typ .eq. 6) then
       write(log_buf, '(A15,I4,F12.6)') 'TAMG -POST', lvl, val
    else
       write(log_buf, '(A15,I4,F12.6)') 'RESID', lvl, val
    end if
    call neko_log%message(log_buf)
  end subroutine phmg_resid_monitor

  subroutine print_phmg_info(nlvls, smoo_type, phmg)
    integer, intent(in) :: nlvls
    integer, intent(in) :: smoo_type
    type(phmg_hrchy_t) :: phmg
    integer :: i, clvl
    character(len=LOG_SIZE) :: log_buf, smoo_name

    call neko_log%section('PHMG')

    if (smoo_type .eq. 1) then
       write(smoo_name, '(A16)') 'CHEBY-acc JACOBI'
    else if (smoo_type .eq. 2) then
       write(smoo_name, '(A17)') 'CHEBY-acc SCHWARZ'
    else
       write(smoo_name, '(A5)') 'CHEBY'
    end if

    write(log_buf, '(A28,I2,A8)') &
         'Creating PHMG hierarchy with', &
         nlvls, 'levels.'
    call neko_log%message(log_buf)

    clvl = nlvls - 1
    do i = 0, nlvls-1
       write(log_buf, '(A8,I2,A8,I2)') &
            '-- level', i, '-- lx:', phmg%lvl(i)%Xh%lx
       call neko_log%message(log_buf)

       if (i .eq. clvl) then
          write(log_buf, '(A19,A20)') &
               'Solve:', 'tAMG'
          call neko_log%message(log_buf)
       else
          write(log_buf, '(A22,A20)') &
               'Smoother:', &
               trim(smoo_name)
          call neko_log%message(log_buf)

          write(log_buf, '(A28,I2)') &
               'Smoother Iters:', &
               phmg%lvl(i)%smoother_itrs
          call neko_log%message(log_buf)
       end if
    end do

    call neko_log%end_section()

  end subroutine print_phmg_info

end module phmg
