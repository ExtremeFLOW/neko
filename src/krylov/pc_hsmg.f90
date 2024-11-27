! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC.
!
! The UChicago Argonne, LLC as Operator of Argonne National
! Laboratory holds copyright in the Software. The copyright holder
! reserves all rights except those expressly granted to licensees,
! and U.S. Government license rights.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the disclaimer below.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the disclaimer (as noted below)
! in the documentation and/or other materials provided with the
! distribution.
!
! 3. Neither the name of ANL nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
! UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF
! ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Additional BSD Notice
! ---------------------
! 1. This notice is required to be provided under our contract with
! the U.S. Department of Energy (DOE). This work was produced at
! Argonne National Laboratory under Contract
! No. DE-AC02-06CH11357 with the DOE.
!
! 2. Neither the United States Government nor UCHICAGO ARGONNE,
! LLC nor any of their employees, makes any warranty,
! express or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.
!
! 3. Also, reference herein to any specific commercial products, process,
! or services by trade name, trademark, manufacturer or otherwise does
! not necessarily constitute or imply its endorsement, recommendation,
! or favoring by the United States Government or UCHICAGO ARGONNE LLC.
! The views and opinions of authors expressed
! herein do not necessarily state or reflect those of the United States
! Government or UCHICAGO ARGONNE, LLC, and shall
! not be used for advertising or product endorsement purposes.
!
!> Krylov preconditioner
module hsmg
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use math, only : copy, col2, add2
  use utils, only : neko_error
  use precon, only : pc_t, precon_factory, precon_destroy
  use ax_product, only : ax_t, ax_helm_factory
  use gather_scatter, only : gs_t, GS_OP_ADD
  use interpolation, only : interpolator_t
  use bc, only : bc_t, bc_list_apply_scalar, bc_list_t, bc_list_add, &
       bc_list_init
  use dirichlet, only : dirichlet_t
  use schwarz, only : schwarz_t
  use jacobi, only : jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use device_jacobi, only : device_jacobi_t
  use device
  use device_math, only : device_copy, device_col2, device_add2
  use profiler, only : profiler_start_region, profiler_end_region
  use space, only : space_t, GLL
  use dofmap, only : dofmap_t
  use field, only : field_t
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER, &
       krylov_solver_factory, krylov_solver_destroy
  !$ use omp_lib
  implicit none
  private

  !Struct to arrange our multigridlevels
  type, private :: multigrid_t
     type(dofmap_t), pointer :: dof
     type(gs_t), pointer  :: gs_h
     type(space_t), pointer :: Xh
     type(coef_t), pointer :: coef
     type(bc_list_t), pointer :: bclst
     type(schwarz_t), pointer :: schwarz
     type(field_t), pointer :: e
  end type multigrid_t

  type, public, extends(pc_t) :: hsmg_t
     type(mesh_t), pointer :: msh
     integer :: nlvls !< Number of levels in the multigrid
     type(multigrid_t), allocatable :: grids(:) !< array for multigrids
     type(gs_t) :: gs_crs, gs_mg !< gather scatter for lower levels
     type(space_t) :: Xh_crs, Xh_mg !< spaces for lower levels
     type(dofmap_t) :: dm_crs, dm_mg
     type(coef_t) :: c_crs, c_mg
     type(dirichlet_t) :: bc_crs, bc_mg, bc_reg
     type(bc_list_t) :: bclst_crs, bclst_mg, bclst_reg
     type(schwarz_t) :: schwarz, schwarz_mg, schwarz_crs !< Schwarz decompostions
     type(field_t) :: e, e_mg, e_crs !< Solve fields
     type(field_t) :: wf !< Work fields
     class(ksp_t), allocatable :: crs_solver !< Solver for course problem
     integer :: niter = 10 !< Number of iter of crs sovlve
     class(pc_t), allocatable :: pc_crs !< Some basic precon for crs
     class(ax_t), allocatable :: ax !< Matrix for crs solve
     real(kind=rp), allocatable :: r(:)!< Residual work array
     type(interpolator_t) :: interp_fine_mid
     type(interpolator_t) :: interp_mid_crs
     real(kind=rp), allocatable :: w(:) !< work array
     type(c_ptr) :: w_d = C_NULL_PTR
     type(c_ptr) :: r_d = C_NULL_PTR
     type(c_ptr) :: hsmg_event
     type(c_ptr) :: gs_event
   contains
     procedure, pass(this) :: init => hsmg_init
     procedure, pass(this) :: free => hsmg_free
     procedure, pass(this) :: solve => hsmg_solve
     procedure, pass(this) :: update => hsmg_set_h
  end type hsmg_t

contains

  !> @note I do not think we actually use the same grids as they do in the original!
  subroutine hsmg_init(this, msh, Xh, coef, dof, gs_h, bclst, crs_pctype)
    class(hsmg_t), intent(inout), target :: this
    type(mesh_t), intent(inout), target :: msh
    type(space_t), intent(inout), target :: Xh
    type(coef_t), intent(inout), target :: coef
    type(dofmap_t), intent(inout), target :: dof
    type(gs_t), intent(inout), target :: gs_h
    type(bc_list_t), intent(inout), target :: bclst
    character(len=*), optional :: crs_pctype
    integer :: n, i
    integer :: lx_crs, lx_mid

    call this%free()
    this%nlvls = 3
    lx_crs = 2
    if (Xh%lx .lt. 5) then
       lx_mid = max(Xh%lx-1,3)

       if(Xh%lx .le. 2) then
          call neko_error('Polynomial order < 2 not supported for hsmg precon')
       end if

    else
       lx_mid = 4
    end if
    this%msh => msh
    allocate(this%grids(this%nlvls))
    allocate(this%w(dof%size()))
    allocate(this%r(dof%size()))


    ! Compute all elements as if they are deformed
    call msh%all_deformed()

    n = dof%size()
    call this%e%init(dof, 'work array')
    call this%wf%init(dof, 'work 2')

    call this%Xh_crs%init(GLL, lx_crs, lx_crs, lx_crs)
    call this%dm_crs%init(msh, this%Xh_crs)
    call this%gs_crs%init(this%dm_crs)
    call this%e_crs%init(this%dm_crs, 'work crs')
    call this%c_crs%init(this%gs_crs)

    call this%Xh_mg%init(GLL, lx_mid, lx_mid, lx_mid)
    call this%dm_mg%init(msh, this%Xh_mg)
    call this%gs_mg%init(this%dm_mg)
    call this%e_mg%init(this%dm_mg, 'work midl')
    call this%c_mg%init(this%gs_mg)

    ! Create backend specific Ax operator
    call ax_helm_factory(this%ax, full_formulation = .false.)

    ! Create a backend specific preconditioner
    call precon_factory(this%pc_crs, 'jacobi')

    ! Create a backend specific krylov solver
    if (present(crs_pctype)) then
       call krylov_solver_factory(this%crs_solver, &
            this%dm_crs%size(), trim(crs_pctype), KSP_MAX_ITER, M = this%pc_crs)
    else
       call krylov_solver_factory(this%crs_solver, &
            this%dm_crs%size(), 'cg', KSP_MAX_ITER, M = this%pc_crs)
    end if

    call this%bc_crs%init_base(this%c_crs)
    call this%bc_mg%init_base(this%c_mg)
    call this%bc_reg%init_base(coef)
    if (bclst%n .gt. 0) then
       do i = 1, bclst%n
          call this%bc_reg%mark_facets(bclst%bc(i)%bcp%marked_facet)
          call this%bc_crs%mark_facets(bclst%bc(i)%bcp%marked_facet)
          call this%bc_mg%mark_facets(bclst%bc(i)%bcp%marked_facet)
       end do
    end if
    call this%bc_reg%finalize()
    call this%bc_reg%set_g(real(0d0, rp))
    call bc_list_init(this%bclst_reg)
    call bc_list_add(this%bclst_reg, this%bc_reg)

    call this%bc_crs%finalize()
    call this%bc_crs%set_g(real(0d0, rp))
    call bc_list_init(this%bclst_crs)
    call bc_list_add(this%bclst_crs, this%bc_crs)


    call this%bc_mg%finalize()
    call this%bc_mg%set_g(0.0_rp)
    call bc_list_init(this%bclst_mg)
    call bc_list_add(this%bclst_mg, this%bc_mg)

    call this%schwarz%init(Xh, dof, gs_h, this%bclst_reg, msh)
    call this%schwarz_mg%init(this%Xh_mg, this%dm_mg, this%gs_mg,&
                              this%bclst_mg, msh)

    call this%interp_fine_mid%init(Xh, this%Xh_mg)
    call this%interp_mid_crs%init(this%Xh_mg, this%Xh_crs)

    call hsmg_fill_grid(dof, gs_h, Xh, coef, this%bclst_reg, this%schwarz, &
                        this%e, this%grids, 3)
    call hsmg_fill_grid(this%dm_mg, this%gs_mg, this%Xh_mg, this%c_mg, &
                        this%bclst_mg, this%schwarz_mg, this%e_mg, &
                        this%grids, 2)
    call hsmg_fill_grid(this%dm_crs, this%gs_crs, this%Xh_crs, &
                        this%c_crs, this%bclst_crs, this%schwarz_crs, &
                        this%e_crs, this%grids, 1)

    call hsmg_set_h(this)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%w, this%w_d, n)
       call device_map(this%r, this%r_d, n)
    end if

    select type(pc => this%pc_crs)
    type is (jacobi_t)
       call pc%init(this%c_crs, this%dm_crs, this%gs_crs)
    type is (sx_jacobi_t)
       call pc%init(this%c_crs, this%dm_crs, this%gs_crs)
    type is (device_jacobi_t)
       call pc%init(this%c_crs, this%dm_crs, this%gs_crs)
    end select
    
    call device_event_create(this%hsmg_event, 2)
    call device_event_create(this%gs_event, 2)
  end subroutine hsmg_init

  subroutine hsmg_set_h(this)
    class(hsmg_t), intent(inout) :: this
!    integer :: i
    !Yeah I dont really know what to do here. For incompressible flow not much happens
    this%grids(1)%coef%ifh2 = .false.
    call copy(this%grids(1)%coef%h1, this%grids(3)%coef%h1, &
         this%grids(1)%dof%size())
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%grids(1)%coef%h1_d, this%grids(3)%coef%h1_d, &
            this%grids(1)%dof%size())
    end if
  end subroutine hsmg_set_h


  subroutine hsmg_fill_grid(dof, gs_h, Xh, coef, bclst, schwarz, e, grids, l)
    type(dofmap_t), target, intent(in) :: dof
    type(gs_t), target, intent(in) :: gs_h
    type(space_t), target, intent(in) :: Xh
    type(coef_t), target, intent(in) :: coef
    type(bc_list_t), target, intent(in) :: bclst
    type(schwarz_t), target, intent(in) :: schwarz
    type(field_t), target, intent(in) :: e
    integer, intent(in) :: l
    type(multigrid_t), intent(inout), dimension(l) :: grids


    grids(l)%dof => dof
    grids(l)%gs_h => gs_h
    grids(l)%Xh => Xh
    grids(l)%coef => coef
    grids(l)%bclst => bclst
    grids(l)%schwarz => schwarz
    grids(l)%e => e

  end subroutine hsmg_fill_grid

  subroutine hsmg_free(this)
    class(hsmg_t), intent(inout) :: this

    if (allocated(this%ax)) then
       deallocate(this%ax)
    end if

    if (allocated(this%grids)) then
       deallocate(this%grids)
    end if

    if (allocated(this%w)) then
       deallocate(this%w)
    end if

    if (allocated(this%r)) then
       deallocate(this%r)
    end if

    call this%schwarz%free()
    call this%schwarz_mg%free()

    call this%c_crs%free()
    call this%c_mg%free()
    call this%e%free()
    call this%e_mg%free()
    call this%e_crs%free()

    call this%gs_crs%free()
    call this%gs_mg%free()
    call this%interp_mid_crs%free()
    call this%interp_fine_mid%free()

    if (allocated(this%crs_solver)) then
       call krylov_solver_destroy(this%crs_solver)
       deallocate(this%crs_solver)
    end if

    if (allocated(this%pc_crs)) then
       call precon_destroy(this%pc_crs)
    end if

  end subroutine hsmg_free

  !> The h1mg preconditioner from Nek5000.
  subroutine hsmg_solve(this, z, r, n)
    integer, intent(in) :: n
    class(hsmg_t), intent(inout) :: this
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(c_ptr) :: z_d, r_d
    type(ksp_monitor_t) :: crs_info
    integer :: thrdid, nthrds

    call profiler_start_region('HSMG_solve', 8)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       z_d = device_get_ptr(z)
       r_d = device_get_ptr(r)
       !We should not work with the input
       call device_copy(this%r_d, r_d, n)
       call bc_list_apply_scalar(this%bclst_reg, r, n)

       !OVERLAPPING Schwarz exchange and solve
       !! DOWNWARD Leg of V-cycle, we are pretty hardcoded here but w/e
       call device_col2(this%r_d, this%grids(3)%coef%mult_d, &
                        this%grids(3)%dof%size())
       !Restrict to middle level
       call this%interp_fine_mid%map(this%e%x, this%r, &
                                     this%msh%nelv, this%grids(2)%Xh)
       call this%grids(2)%gs_h%op(this%e%x, &
                  this%grids(2)%dof%size(), GS_OP_ADD, this%gs_event)
       call device_event_sync(this%gs_event)
       call device_copy(this%r_d, r_d, n)
       call bc_list_apply_scalar(this%bclst_reg, r, n)
       call device_copy(this%w_d, this%e%x_d, this%grids(2)%dof%size())
       call bc_list_apply_scalar(this%bclst_mg, this%w, &
                                 this%grids(2)%dof%size())
       !OVERLAPPING Schwarz exchange and solve
       call device_col2(this%w_d, this%grids(2)%coef%mult_d, &
                        this%grids(2)%dof%size())
       !restrict residual to crs
       call this%interp_mid_crs%map(this%wf%x, this%w, this%msh%nelv, &
                                    this%grids(1)%Xh)
       !Crs solve
       call device_copy(this%w_d, this%e%x_d, this%grids(2)%dof%size())
       call bc_list_apply_scalar(this%bclst_mg, this%w, &
                                 this%grids(2)%dof%size())

       !$omp parallel private(thrdid, nthrds)

       thrdid = 0
       nthrds = 1
       !$ thrdid = omp_get_thread_num()
       !$ nthrds = omp_get_num_threads()

       if (thrdid .eq. 0) then
          call profiler_start_region('HSMG_schwarz', 9)
          call this%grids(3)%schwarz%compute(z, this%r)
          call this%grids(2)%schwarz%compute(this%grids(2)%e%x, this%w)
          call profiler_end_region('HSMG_schwarz', 9)
       end if
       if (nthrds .eq. 1 .or. thrdid .eq. 1) then
          call profiler_start_region('HSMG_coarse_grid', 10)
          call this%grids(1)%gs_h%op(this%wf%x, &
               this%grids(1)%dof%size(), GS_OP_ADD, this%gs_event)
          call device_event_sync(this%gs_event)
          call bc_list_apply_scalar(this%grids(1)%bclst, this%wf%x, &
                                    this%grids(1)%dof%size())
          call profiler_start_region('HSMG_coarse_solve', 11)
          crs_info = this%crs_solver%solve(this%Ax, this%grids(1)%e, &
                                       this%wf%x, &
                                       this%grids(1)%dof%size(), &
                                       this%grids(1)%coef, &
                                       this%grids(1)%bclst, &
                                       this%grids(1)%gs_h, this%niter)
          call profiler_end_region('HSMG_coarse_solve', 11)
          call bc_list_apply_scalar(this%grids(1)%bclst, this%grids(1)%e%x,&
                                    this%grids(1)%dof%size())
          call profiler_end_region('HSMG_coarse_grid', 10)
       end if
       !$omp end parallel

       call this%interp_mid_crs%map(this%w, this%grids(1)%e%x, &
                                    this%msh%nelv, this%grids(2)%Xh)
       call device_add2(this%grids(2)%e%x_d, this%w_d, this%grids(2)%dof%size())

       call this%interp_fine_mid%map(this%w, this%grids(2)%e%x, &
                                     this%msh%nelv, this%grids(3)%Xh)
       call device_add2(z_d, this%w_d, this%grids(3)%dof%size())
       call this%grids(3)%gs_h%op(z, this%grids(3)%dof%size(), &
                                     GS_OP_ADD, this%gs_event)
       call device_event_sync(this%gs_event)
       call device_col2(z_d, this%grids(3)%coef%mult_d, &
                        this%grids(3)%dof%size())
    else
       !We should not work with the input
       call copy(this%r, r, n)

       !OVERLAPPING Schwarz exchange and solve
       call this%grids(3)%schwarz%compute(z, this%r)
       ! DOWNWARD Leg of V-cycle, we are pretty hardcoded here but w/e
       call col2(this%r, this%grids(3)%coef%mult, &
                 this%grids(3)%dof%size())
       !Restrict to middle level
       call this%interp_fine_mid%map(this%w, this%r, &
                                     this%msh%nelv, this%grids(2)%Xh)
       call this%grids(2)%gs_h%op(this%w, this%grids(2)%dof%size(), GS_OP_ADD)
       !OVERLAPPING Schwarz exchange and solve
       call this%grids(2)%schwarz%compute(this%grids(2)%e%x, this%w)
       call col2(this%w, this%grids(2)%coef%mult, this%grids(2)%dof%size())
       !restrict residual to crs
       call this%interp_mid_crs%map(this%r, this%w, &
            this%msh%nelv, this%grids(1)%Xh)
       !Crs solve

       call this%grids(1)%gs_h%op(this%r, this%grids(1)%dof%size(), GS_OP_ADD)
       call bc_list_apply_scalar(this%grids(1)%bclst, this%r, &
                                 this%grids(1)%dof%size())
       call profiler_start_region('HSMG_coarse-solve', 11)
       crs_info = this%crs_solver%solve(this%Ax, this%grids(1)%e, this%r, &
                                    this%grids(1)%dof%size(), &
                                    this%grids(1)%coef, &
                                    this%grids(1)%bclst, &
                                    this%grids(1)%gs_h, this%niter)
       call profiler_end_region('HSMG_coarse-solve', 11)
       call bc_list_apply_scalar(this%grids(1)%bclst, this%grids(1)%e%x,&
                                 this%grids(1)%dof%size())


       call this%interp_mid_crs%map(this%w, this%grids(1)%e%x, &
                                    this%msh%nelv, this%grids(2)%Xh)
       call add2(this%grids(2)%e%x, this%w, this%grids(2)%dof%size())

       call this%interp_fine_mid%map(this%w, this%grids(2)%e%x, &
                                     this%msh%nelv, this%grids(3)%Xh)
       call add2(z, this%w, this%grids(3)%dof%size())
       call this%grids(3)%gs_h%op(z, this%grids(3)%dof%size(), GS_OP_ADD)
       call col2(z, this%grids(3)%coef%mult, this%grids(3)%dof%size())

    end if
    call profiler_end_region('HSMG_solve', 8)
  end subroutine hsmg_solve
end module hsmg
