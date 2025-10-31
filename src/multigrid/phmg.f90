! Copyright (c) 2024-2025, The Neko Authors
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
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use dirichlet, only : dirichlet_t
  use utils, only : neko_error
  use cheby, only : cheby_t
  use cheby_device, only : cheby_device_t
  use jacobi, only : jacobi_t
  use device_jacobi, only : device_jacobi_t
  use schwarz, only : schwarz_t
  use ax_product, only : ax_t, ax_helm_factory
  use tree_amg_multigrid, only : tamg_solver_t
  use interpolation, only : interpolator_t
  use json_module, only : json_file
  use json_utils, only : json_get_or_default, json_get
  use math, only : copy, col2, add2, sub3, add2s2
  use device, only : device_get_ptr, device_stream_wait_event, glb_cmd_queue, &
       glb_cmd_event
  use device_math, only : device_rzero, device_copy, device_add2, device_sub3,&
       device_add2s2, device_invcol2, device_glsc2, device_col2
  use profiler, only : profiler_start_region, profiler_end_region
  use neko_config, only: NEKO_BCKND_DEVICE
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER, &
       krylov_solver_factory
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
     type(bc_list_t) :: bclst
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
     type(mesh_t), pointer :: msh
   contains
     procedure, pass(this) :: init => phmg_init
     procedure, pass(this) :: init_from_components => &
          phmg_init_from_components
     procedure, pass(this) :: free => phmg_free
     procedure, pass(this) :: solve => phmg_solve
     procedure, pass(this) :: update => phmg_update
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

    call json_get_or_default(phmg_params, 'smoother_iterations', &
         smoother_itrs, 10)

    call json_get_or_default(phmg_params, 'smoother_cheby_acc', &
         cheby_acc, "jacobi")

    call json_get_or_default(phmg_params, 'coarse_grid.levels', &
         crs_tamg_lvls, 3)

    call json_get_or_default(phmg_params, 'coarse_grid.iterations', &
         crs_tamg_itrs, 1)

    call json_get_or_default(phmg_params, 'coarse_grid.cheby_degree', &
         crs_tamg_cheby_degree, 5)

    if (phmg_params%valid_path('pcoarsening_schedule')) then
       call json_get(phmg_params, 'pcoarsening_schedule', pcrs_sched)
    else
       allocate(pcrs_sched(2))
       pcrs_sched(1) = 3
       pcrs_sched(2) = 1
    end if


    call this%init_from_components(coef, bclst, smoother_itrs, &
         cheby_acc, crs_tamg_lvls, crs_tamg_itrs, crs_tamg_cheby_degree,&
         pcrs_sched)

  end subroutine phmg_init

  subroutine phmg_init_from_components(this, coef, bclst, smoother_itrs, &
       cheby_acc, crs_tamg_lvls, crs_tamg_itrs, crs_tamg_cheby_degree, &
       pcrs_sched)
    class(phmg_t), intent(inout), target :: this
    type(coef_t), intent(in), target :: coef
    type(bc_list_t), intent(inout), target :: bclst
    integer, intent(in) :: smoother_itrs
    character(len=:), allocatable :: cheby_acc
    integer, intent(in) :: crs_tamg_lvls, crs_tamg_itrs
    integer, intent(in) :: crs_tamg_cheby_degree
    integer, intent(in), allocatable :: pcrs_sched(:)
    integer :: lx_crs, lx_mid
    integer, allocatable :: lx_lvls(:)
    integer :: n, i, j, st
    class(bc_t), pointer :: bc_j
    logical :: use_jacobi, use_cheby
    use_jacobi = .true.
    use_cheby = .true.

    this%msh => coef%msh

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
       call this%phmg_hrchy%lvl(i)%coef%init(this%phmg_hrchy%lvl(i)%gs_h)
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
       call this%phmg_hrchy%lvl(i)%bclst%init()
       call this%phmg_hrchy%lvl(i)%bclst%append(this%phmg_hrchy%lvl(i)%bc)

       !> Initialize Smoothers
       if (trim(cheby_acc) .eq. "schwarz") then
          call this%phmg_hrchy%lvl(i)%schwarz%init( &
               this%phmg_hrchy%lvl(i)%Xh, &
               this%phmg_hrchy%lvl(i)%dm_Xh, &
               this%phmg_hrchy%lvl(i)%gs_h, &
               this%phmg_hrchy%lvl(i)%bclst, &
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
    call ax_helm_factory(this%ax, full_formulation = .false.)

    ! Interpolator Fine + mg levels
    allocate(this%intrp(this%nlvls - 1))
    do i = 1, this%nlvls -1
       call this%intrp(i)%init(this%phmg_hrchy%lvl(i-1)%Xh, &
            this%phmg_hrchy%lvl(i)%Xh)
    end do

    call this%amg_solver%init(this%ax, this%phmg_hrchy%lvl(this%nlvls -1)%Xh, &
         this%phmg_hrchy%lvl(this%nlvls -1)%coef, this%msh, &
         this%phmg_hrchy%lvl(this%nlvls-1)%gs_h, crs_tamg_lvls, &
         this%phmg_hrchy%lvl(this%nlvls -1)%bclst, &
         crs_tamg_itrs, crs_tamg_cheby_degree)

  end subroutine phmg_init_from_components

  subroutine phmg_free(this)
    class(phmg_t), intent(inout) :: this
  end subroutine phmg_free

  subroutine phmg_solve(this, z, r, n)
    class(phmg_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(c_ptr) :: z_d, r_d
    type(ksp_monitor_t) :: ksp_results


    associate( mglvl => this%phmg_hrchy%lvl)
      if (NEKO_BCKND_DEVICE .eq. 1) then
         z_d = device_get_ptr(z)
         r_d = device_get_ptr(r)
         !We should not work with the input
         call device_copy(mglvl(0)%r%x_d, r_d, n)
         call device_rzero(mglvl(0)%z%x_d, n)
         call device_rzero(mglvl(0)%w%x_d, n)
         call phmg_mg_cycle(mglvl(0)%z, mglvl(0)%r, mglvl(0)%w, 0, &
              this%nlvls -1, mglvl, this%intrp, this%msh, this%Ax, &
              this%amg_solver)

         call device_copy(z_d, mglvl(0)%z%x_d, n)
      else
         !We should not work with the input
         call copy(mglvl(0)%r%x, r, n)

         mglvl(0)%z%x = 0.0_rp
         mglvl(0)%w%x = 0.0_rp

         call phmg_mg_cycle(mglvl(0)%z, mglvl(0)%r, mglvl(0)%w, 0, &
              this%nlvls -1, mglvl, this%intrp, this%msh, this%Ax, &
              this%amg_solver)

         call copy(z, mglvl(0)%z%x, n)
      end if
    end associate

  end subroutine phmg_solve

  subroutine phmg_update(this)
    class(phmg_t), intent(inout) :: this
  end subroutine phmg_update


  recursive subroutine phmg_mg_cycle(z, r, w, lvl, clvl, &
       mg, intrp, msh, Ax, amg_solver)
    type(ksp_monitor_t) :: ksp_results
    integer :: lvl, clvl
    type(phmg_lvl_t) :: mg(0:clvl)
    type(interpolator_t) :: intrp(1:clvl)
    type(tamg_solver_t), intent(inout) :: amg_solver
    class(ax_t), intent(inout) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(field_t) :: z, r, w
    integer :: i
    logical :: use_jacobi
    real(kind=rp) :: val

    use_jacobi = .false.
    call profiler_start_region('PHMG_cycle', 8)
    !>----------<!
    !> SMOOTH   <!
    !>----------<!
    call profiler_start_region('PHMG_PreSmooth', 9)
    if (use_jacobi) then
       call phmg_jacobi_smoother(z, r, w, mg(lvl), msh, Ax, &
            mg(lvl)%dm_Xh%size(), lvl)
    else
       if (NEKO_BCKND_DEVICE .eq. 1) then
          mg(lvl)%cheby_device%zero_initial_guess = .true.
          ksp_results = mg(lvl)%cheby_device%solve(Ax, z, &
               r%x, mg(lvl)%dm_Xh%size(), &
               mg(lvl)%coef, mg(lvl)%bclst, &
               mg(lvl)%gs_h, niter = mg(lvl)%smoother_itrs)
       else
          mg(lvl)%cheby%zero_initial_guess = .true.
          ksp_results = mg(lvl)%cheby%solve(Ax, z, &
               r%x, mg(lvl)%dm_Xh%size(), &
               mg(lvl)%coef, mg(lvl)%bclst, &
               mg(lvl)%gs_h, niter = mg(lvl)%smoother_itrs)
       end if
    end if
    call profiler_end_region('PHMG_PreSmooth', 9)

    !>----------<!
    !> Residual <!
    !>----------<!
    call Ax%compute(w%x, z%x, mg(lvl)%coef, msh, mg(lvl)%Xh)
    call mg(lvl)%gs_h%op(w%x, mg(lvl)%dm_Xh%size(), GS_OP_ADD, glb_cmd_event)
    call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)
    call mg(lvl)%bclst%apply_scalar(w%x, mg(lvl)%dm_Xh%size())

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_sub3(w%x_d, r%x_d, w%x_d, mg(lvl)%dm_Xh%size())
    else
       w%x = r%x - w%x
    end if

    !>----------<!
    !> Restrict <!
    !>----------<!
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(w%x_d, mg(lvl)%coef%mult_d, mg(lvl)%dm_Xh%size())
    else
       call col2(w%x, mg(lvl)%coef%mult, mg(lvl)%dm_Xh%size())
    end if

    call profiler_start_region('PHMG_map_to_coarse', 9)
    call intrp(lvl+1)%map(mg(lvl+1)%r%x, w%x, msh%nelv, mg(lvl+1)%Xh)
    call profiler_end_region('PHMG_map_to_coarse', 9)

    call mg(lvl+1)%gs_h%op(mg(lvl+1)%r%x, mg(lvl+1)%dm_Xh%size(), &
         GS_OP_ADD, glb_cmd_event)
    call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)

    call mg(lvl+1)%bclst%apply_scalar( &
         mg(lvl+1)%r%x, &
         mg(lvl+1)%dm_Xh%size())
    !>----------<!
    !> SOLVE    <!
    !>----------<!
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(mg(lvl+1)%z%x_d, mg(lvl+1)%dm_Xh%size())
    else
       mg(lvl+1)%z%x = 0.0_rp
    end if
    if (lvl+1 .eq. clvl) then
       call profiler_start_region('PHMG_tAMG_coarse_grid', 9)
       call amg_solver%solve(mg(lvl+1)%z%x, &
            mg(lvl+1)%r%x, &
            mg(lvl+1)%dm_Xh%size())
       call profiler_end_region('PHMG_tAMG_coarse_grid', 9)

    else
       call phmg_mg_cycle(mg(lvl+1)%z, mg(lvl+1)%r, mg(lvl+1)%w, lvl+1, &
            clvl, mg, intrp, msh, Ax, amg_solver)
    end if

    !>----------<!
    !> Project  <!
    !>----------<!
    call profiler_start_region('PHMG_map_to_fine', 9)
    call intrp(lvl+1)%map(w%x, mg(lvl+1)%z%x, msh%nelv, mg(lvl)%Xh)
    call profiler_end_region('PHMG_map_to_fine', 9)

    call mg(lvl)%gs_h%op(w%x, mg(lvl)%dm_Xh%size(), GS_OP_ADD, glb_cmd_event)
    call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(w%x_d, mg(lvl)%coef%mult_d, mg(lvl)%dm_Xh%size())
    else
       call col2(w%x, mg(lvl)%coef%mult, mg(lvl)%dm_Xh%size())
    end if

    !>----------<!
    !> Correct  <!
    !>----------<!
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(z%x_d, w%x_d, mg(lvl)%dm_Xh%size())
    else
       z%x = z%x + w%x
    end if

    !>----------<!
    !> SMOOTH   <!
    !>----------<!
    call profiler_start_region('PHMG_PostSmooth', 9)
    if (use_jacobi) then
       call phmg_jacobi_smoother(z, r, w, mg(lvl), msh, Ax, &
            mg(lvl)%dm_Xh%size(), lvl)
    else
       if (NEKO_BCKND_DEVICE .eq. 1) then
          ksp_results = mg(lvl)%cheby_device%solve(Ax, z, &
               r%x, mg(lvl)%dm_Xh%size(), &
               mg(lvl)%coef, mg(lvl)%bclst, &
               mg(lvl)%gs_h, niter = mg(lvl)%smoother_itrs)
       else
          ksp_results = mg(lvl)%cheby%solve(Ax, z, &
               r%x, mg(lvl)%dm_Xh%size(), &
               mg(lvl)%coef, mg(lvl)%bclst, &
               mg(lvl)%gs_h, niter = mg(lvl)%smoother_itrs)
       end if
    end if
    call profiler_end_region('PHMG_PostSmooth', 9)

    call profiler_end_region('PHMG_cycle', 8)
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
    integer :: i, iblk, ni, niblk

    ni = mg%smoother_itrs
    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, ni
          call Ax%compute(w%x, z%x, mg%coef, msh, mg%Xh)
          call mg%gs_h%op(w%x, n, GS_OP_ADD, glb_cmd_event)
          call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)
          call mg%bclst%apply_scalar(w%x, n)
          call device_sub3(w%x_d, r%x_d, w%x_d, n)

          call mg%device_jacobi%solve(w%x, w%x, n)

          call device_add2s2(z%x_d, w%x_d, 0.6_rp, n)
       end do
    else
       do i = 1, ni
          call Ax%compute(w%x, z%x, mg%coef, msh, mg%Xh)
          call mg%gs_h%op(w%x, n, GS_OP_ADD)
          call mg%bclst%apply_scalar(w%x, n)
          call sub3(w%x, r%x, w%x, n)

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
    call mg%bclst%apply_scalar(w%x, mg%dm_Xh%size())
    call device_sub3(w%x_d, r%x_d, w%x_d, mg%dm_Xh%size())
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
