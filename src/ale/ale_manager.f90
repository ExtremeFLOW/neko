!> ALE Manager: Handles Mesh Motion, Solver, and History
module ale_manager
  use num_types,          only : rp, dp
  use json_module,        only : json_file
  use json_utils,         only : json_get, json_get_or_default, json_extract_item
  use field,              only : field_t
  use coefs,              only : coef_t
  use ax_product,         only : ax_t, ax_helm_factory
  use krylov,             only : ksp_t, ksp_monitor_t, krylov_solver_factory
  use precon,             only : pc_t, precon_factory, precon_destroy
  use bc_list,            only : bc_list_t
  use zero_dirichlet,     only : zero_dirichlet_t
  use gather_scatter,     only : gs_t, GS_OP_ADD
  use mesh,               only : mesh_t
  use dofmap,             only : dofmap_t
  use jacobi,             only : jacobi_t
  use device_jacobi,      only : device_jacobi_t
  use sx_jacobi,          only : sx_jacobi_t
  use profiler,           only : profiler_start_region, profiler_end_region
  use math,               only : masked_copy_0
  use file,               only : file_t 
  use logger,             only : neko_log, LOG_SIZE
  use ale_motion,         only : ale_config_t, update_ale_mesh_velocity, &
                                 update_ale_mesh, update_ale_mass_history, &
                                 init_ale_pivot, pivot_state_t
  use utils,              only : neko_error 
  use mpi_f08,            only : MPI_WTIME, MPI_Barrier
  use comm
  use registry,           only : neko_registry
  use field_series,       only : field_series_t
  use time_state,         only : time_state_t

  implicit none
  private

  type, public :: ale_manager_t
     class(ax_t),        allocatable :: Ax
     class(ksp_t),       allocatable :: ksp
     class(pc_t),        allocatable :: pc
     type(ksp_monitor_t)             :: monitor(1)

     type(bc_list_t)            :: bc_list
     type(zero_dirichlet_t)     :: bc_moving 
     type(zero_dirichlet_t)     :: bc_fixed 
     
     type(ale_config_t) :: config  
     logical            :: has_moving_boundary = .false.
     logical            :: active = .false.
     
     real(kind=rp) :: abstol       = 1.0e-10_rp
     integer       :: ksp_max_iter = 10000
     character(len=:), allocatable :: ksp_solver
     character(len=:), allocatable :: precon_type

     !> Mesh velocity fields (Registered in neko_registry)
     type(field_t), pointer :: wm_x => null()
     type(field_t), pointer :: wm_y => null()
     type(field_t), pointer :: wm_z => null()

     !> History buffers for mesh velocity components (for AB/BDF schemes)
     type(field_series_t) :: wm_lag_x
     type(field_series_t) :: wm_lag_y
     type(field_series_t) :: wm_lag_z
     type(pivot_state_t) :: ale_pivot  
     !> Base shape field for mesh movement (Laplace solution)
     type(field_t) :: base_shape 

     !> Lagged B matrices for ALE (geometry history)
     real(kind=rp), pointer :: Blag(:,:,:,:) => null()
     real(kind=rp), pointer :: Blaglag(:,:,:,:) => null()

  contains
     procedure, pass(this) :: init  => ale_manager_init
     procedure, pass(this) :: free  => ale_manager_free
     procedure, pass(this) :: solve_base_displacement
     
     procedure, pass(this) :: advance_mesh
     procedure, pass(this) :: update_mesh_velocity
  end type ale_manager_t

contains

  !> Initialize ALE Manager
  !> Sets up solver, registers fields, solves for base shape, and inits history.
  subroutine ale_manager_init(this, coef, json)
    class(ale_manager_t), intent(inout) :: this
    type(coef_t),         intent(inout) :: coef
    type(json_file),      intent(inout) :: json

    type(json_file) :: bc_subdict
    type(time_state_t) :: t_init
    integer, allocatable :: zone_indices(:)
    character(len=:), allocatable :: bc_type
    character(len=:), allocatable :: tmp_str
    real(kind=rp), allocatable :: tmp_vec(:) 
    character(len=LOG_SIZE) :: log_buf
    integer :: i, j, n_bcs, n
    logical :: tmp_logical = .true.
    integer :: time_order

    ! Check Activation
    if (json%valid_path('case.fluid.ale.ale_active')) then
       call json%get('case.fluid.ale.ale_active', this%active)
    end if    
    
    if (.not. this%active) return

    call profiler_start_region('ALE_INIT')
    call neko_log%section("ALE Initialization")

    n = coef%dof%size()
    call json_get_or_default(json, 'case.numerics.time_order', time_order, 3)

    ! Register Mesh Velocity Fields
    call neko_registry%add_field(coef%dof, 'wm_x')
    call neko_registry%add_field(coef%dof, 'wm_y')
    call neko_registry%add_field(coef%dof, 'wm_z')
    this%wm_x => neko_registry%get_field('wm_x')
    this%wm_y => neko_registry%get_field('wm_y')
    this%wm_z => neko_registry%get_field('wm_z')

    ! Initialize Base Shape
    call this%base_shape%init(coef%dof, "base_shape")
    this%base_shape%x = 0.0_rp

    ! Linear Solver Configuration 
    if (.not. allocated(this%ksp_solver))  this%ksp_solver  = 'cg'
    if (.not. allocated(this%precon_type)) this%precon_type = 'jacobi'

    if (json%valid_path('case.fluid.ale.solver.type')) then
       call json%get('case.fluid.ale.solver.type', this%ksp_solver)
    end if
    if (json%valid_path('case.fluid.ale.solver.preconditioner.type')) then
       call json%get('case.fluid.ale.solver.preconditioner.type', this%precon_type)
    end if
    if (json%valid_path('case.fluid.ale.solver.absolute_tolerance')) then
       call json%get('case.fluid.ale.solver.absolute_tolerance', this%abstol)
    end if
    if (json%valid_path('case.fluid.ale.solver.max_iterations')) then
       call json%get('case.fluid.ale.solver.max_iterations', this%ksp_max_iter)
    end if
    if (json%valid_path('case.fluid.ale.solver.output_base_shape')) then
        call json%get('case.fluid.ale.solver.output_base_shape', tmp_logical)
        this%config%if_output_phi = tmp_logical
    end if

    ! Initialize Matrix and Boundary Conditions
    call ax_helm_factory(this%Ax, full_formulation = .false.)
    call this%bc_moving%init_from_components(coef)
    call this%bc_fixed%init_from_components(coef)

    if (json%valid_path('case.fluid.boundary_conditions')) then
       call json%info('case.fluid.boundary_conditions', n_children = n_bcs)
       do i = 1, n_bcs
          call json_extract_item(json, 'case.fluid.boundary_conditions', i, bc_subdict)
          call json_get(bc_subdict, 'type', bc_type)
          call json_get(bc_subdict, 'zone_indices', zone_indices)
          select case (trim(bc_type))
          case ('moving_boundary')
             do j = 1, size(zone_indices)
                call this%bc_moving%mark_zone(coef%msh%labeled_zones(zone_indices(j)))
             end do
             this%has_moving_boundary = .true.
          case default
             do j = 1, size(zone_indices)
                call this%bc_fixed%mark_zone(coef%msh%labeled_zones(zone_indices(j)))
             end do
          end select
       end do
    end if

    call this%bc_moving%finalize()
    call this%bc_fixed%finalize()
    call this%bc_list%init()
    call this%bc_list%append(this%bc_moving)
    call this%bc_list%append(this%bc_fixed)

    ! Stiffness Configuration
    if (json%valid_path('case.fluid.ale.solver.mesh_stiffness.type')) then
       call json%get('case.fluid.ale.solver.mesh_stiffness.type', tmp_str)
       this%config%stiffness_type = tmp_str
       if (.not. (trim(tmp_str) == 'built-in')) then
        call neko_error("ALE: stiffness_type must be 'built-in'")
       end if
    end if

    select case (trim(this%config%stiffness_type))
    case ('built-in')
       if (json%valid_path('case.fluid.ale.solver.mesh_stiffness.params.gain')) then
          call json%get('case.fluid.ale.solver.mesh_stiffness.params.gain', this%config%stiffness_gain)
       end if
       if (json%valid_path('case.fluid.ale.solver.mesh_stiffness.params.decay')) then
          call json%get('case.fluid.ale.solver.mesh_stiffness.params.decay', this%config%stiffness_decay)
       end if
    end select
    
    ! Motion Configuration
    if (json%valid_path('case.fluid.ale.motion.oscillation_amplitude')) then
        call json%get('case.fluid.ale.motion.oscillation_amplitude', tmp_vec)
        if (size(tmp_vec) == 3) then
           this%config%osc_amp = tmp_vec
        else
           call neko_error("ALE: oscillation_amplitude must be size 3")
        end if
    end if
    if (json%valid_path('case.fluid.ale.motion.oscillation_frequency')) then
        call json%get('case.fluid.ale.motion.oscillation_frequency', tmp_vec)
        if (size(tmp_vec) == 3) then
           this%config%osc_freq = tmp_vec
        else
           call neko_error("ALE: oscillation_frequency must be size 3")
        end if
    end if

    if (json%valid_path('case.fluid.ale.motion.rotation_center_type')) then
        call json%get('case.fluid.ale.motion.rotation_center_type', tmp_str)
        if (.not. (trim(tmp_str) == 'fixed' .or. &
        &trim(tmp_str) == 'relative_sin' .or. &
        &trim(tmp_str) == 'relative')) then
             call neko_error("ALE: rotation_center_type must be 'fixed' or &
             &'relative' or 'relative_sin'.")
        end if
        this%config%rotation_center_type = tmp_str
    end if

    if (json%valid_path('case.fluid.ale.motion.rotation_amplitude_deg')) then
        call json%get('case.fluid.ale.motion.rotation_amplitude_deg', tmp_vec)
        if (size(tmp_vec) == 3) then
           this%config%rot_amp_degree = tmp_vec
        else
           call neko_error("ALE: rotation_amplitude_deg must be size 3")
        end if
    end if
    if (json%valid_path('case.fluid.ale.motion.rotation_frequency')) then
        call json%get('case.fluid.ale.motion.rotation_frequency', tmp_vec)
        if (size(tmp_vec) == 3) then
           this%config%rot_freq = tmp_vec
        else
           call neko_error("ALE: rotation_frequency must be size 3")
        end if
    end if
    if (json%valid_path('case.fluid.ale.motion.rotation_center')) then
        call json%get('case.fluid.ale.motion.rotation_center', tmp_vec)
        if (size(tmp_vec) == 3) then
           this%config%rot_center = tmp_vec
        else
           call neko_error("ALE: rotation_center must be size 3")
        end if
    end if

    ! --- LOGGING ---
    call neko_log%message('Solver Type  : ('// trim(this%ksp_solver) // &
         ', ' // trim(this%precon_type) // ')')
    write(log_buf, '(A,ES13.6)') 'Abs tol      :', this%abstol
    call neko_log%message(log_buf)

    call neko_log%message('Mesh Stiffness    : ' // trim(this%config%stiffness_type))
    write(log_buf, '(A,F8.2)')   '  Gain       :', this%config%stiffness_gain
    call neko_log%message(log_buf)
    write(log_buf, '(A,F8.2)')   '  Decay      :', this%config%stiffness_decay
    call neko_log%message(log_buf)

    call neko_log%message('Motion       :')

    if (any(abs(this%config%osc_amp) > 0.0_rp)) then
       write(log_buf, '(A,3(F8.3,1X))') '  Osc Amp    :', this%config%osc_amp
       call neko_log%message(log_buf)
       write(log_buf, '(A,3(F8.3,1X))') '  Osc Freq   :', this%config%osc_freq
       call neko_log%message(log_buf)
    else
       call neko_log%message('  Oscillation: None')
    end if

    if (any(abs(this%config%rot_amp_degree) > 0.0_rp)) then
       call neko_log%message('  Pivot Type : ' // trim(this%config%rotation_center_type))
       if (trim(this%config%rotation_center_type) == 'fixed') then
          write(log_buf, '(A,3(F8.3,1X))') '  Fixed Pivot Point:', this%config%rot_center
          call neko_log%message(log_buf)
       else if (trim(this%config%rotation_center_type) == 'relative' .or. &
         &trim(this%config%rotation_center_type) == 'relative_sin') then
          write(log_buf, '(A,3(F8.3,1X))') '  Initial Pivot Point:', this%config%rot_center
          call neko_log%message(log_buf)
       end if
       write(log_buf, '(A,3(F8.3,1X))') '  Rot Amp(deg):', this%config%rot_amp_degree
       call neko_log%message(log_buf)
       write(log_buf, '(A,3(F8.3,1X))') '  Rot Freq    :', this%config%rot_freq
       call neko_log%message(log_buf)
    else
       call neko_log%message('  Rotation: None')
    end if

    ! Create KSP Solver
    call krylov_solver_factory(this%ksp, n, this%ksp_solver, &
         this%ksp_max_iter, this%abstol)
    call ale_precon_factory(this%pc, this%ksp, coef, coef%dof, &
         coef%gs_h, this%bc_list, this%precon_type)

    ! Calculate Base Shape (Laplace Solve)
    call this%solve_base_displacement(this%base_shape, coef)

    call init_ale_pivot(this%ale_pivot, this%config)
    ! Initialize Mesh Velocities at t=t_init
    t_init%t = 0.0_rp
    t_init%tstep = 0
    call update_ale_mesh_velocity(this%wm_x, this%wm_y, this%wm_z, &
                             this%base_shape, coef, &
                             t_init, this%config, this%ale_pivot, time_order)

    call this%wm_lag_x%init(this%wm_x, time_order)
    call this%wm_lag_y%init(this%wm_y, time_order)
    call this%wm_lag_z%init(this%wm_z, time_order)
    
    do i = 1, time_order
        this%wm_lag_x%lf(i)%x = 0.0_rp
        this%wm_lag_y%lf(i)%x = 0.0_rp
        this%wm_lag_z%lf(i)%x = 0.0_rp
    end do

    allocate(this%Blag(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(this%Blaglag(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    this%Blag = coef%B
    this%Blaglag = coef%B

    call profiler_end_region('ALE_INIT')
    call neko_log%end_section()
  end subroutine ale_manager_init

  !> Solves the Laplace equation to determine the base shape (phi)
  !> phi = 1 on moving boundaries, phi = 0 on fixed boundaries.
  subroutine solve_base_displacement(this, phi, coef)
    class(ale_manager_t), intent(inout) :: this
    type(field_t),        intent(inout) :: phi
    type(coef_t),         intent(inout) :: coef
    type(file_t)  :: phi_file
    type(field_t) :: rhs_field
    type(field_t) :: corr_field
    real(kind=dp) :: sample_start_time, sample_end_time
    real(kind=dp) :: sample_time
    character(len=LOG_SIZE) :: log_buf
    integer :: n, i, m, k, ierr
    
    if (.not. this%has_moving_boundary) return

    call profiler_start_region('ALE_SOLVE', 99)

    call MPI_Barrier(NEKO_COMM, ierr)
    sample_start_time = MPI_WTIME()

    call neko_log%message(" ")
    call neko_log%message("ALE: Starting base mesh motion solve...")
    n = coef%dof%size()

    call rhs_field%init(coef%dof)
    call corr_field%init(coef%dof)

    ! Compute stiffness
    call compute_stiffness_ale(coef, this%config)

    ! Lift boundary conditions
    phi%x = 0.0_rp
    m = this%bc_moving%msk(0)
    do i = 1, m
       k = this%bc_moving%msk(i)
       phi%x(k,1,1,1) = 1.0_rp
    end do
    call this%bc_fixed%apply_scalar(phi%x, n)

    ! Compute RHS
    call this%Ax%compute(rhs_field%x, phi%x, coef, coef%msh, coef%Xh)
    rhs_field%x = -rhs_field%x

    ! Mask RHS
    call this%bc_list%apply_scalar(rhs_field%x, n)
    call coef%gs_h%op(rhs_field, GS_OP_ADD)

    ! Solve
    corr_field%x = 0.0_rp
    call this%pc%update()

    this%monitor(1) = this%ksp%solve(this%Ax, corr_field, rhs_field%x, &
                                     n, coef, this%bc_list, coef%gs_h)

    ! Final Solution
    phi%x = phi%x + corr_field%x
    
    call MPI_Barrier(NEKO_COMM, ierr)
    sample_end_time = MPI_WTIME()
    sample_time = sample_end_time - sample_start_time
    write(log_buf, '(A, F10.4, A)') "  Laplace solve took ", sample_time, " (s)"
    call neko_log%message(log_buf)
    
    if (this%config%if_output_phi) then
       call neko_log%message("ALE: Writing base_shape to phi0.f00000")
       call phi_file%init('phi.fld')
       call phi_file%write(phi)
       call phi_file%free()
    end if

    call rhs_field%free()
    call corr_field%free()
    call profiler_end_region('ALE_SOLVE', 99)
  end subroutine solve_base_displacement

  !> Main routine to advance the mesh in time
  !> Updates history (B), moves mesh points, and updates mesh velocity history
  subroutine advance_mesh(this, coef, time, nadv)
     class(ale_manager_t), intent(inout) :: this
     type(coef_t), intent(inout) :: coef
     type(time_state_t), intent(in) :: time
     integer, intent(in) :: nadv

     if (.not. this%active) return
     if (.not. this%has_moving_boundary) return

     ! Update lagged B terms (geometry history)
     call update_ale_mass_history(coef%B, this%Blag, this%Blaglag, coef%dof%size())
     
     ! Update ALE mesh coordinates
     call update_ale_mesh(coef, this%wm_x, this%wm_y, this%wm_z, &
                          this%wm_lag_x, this%wm_lag_y, this%wm_lag_z, &
                          time, nadv, "ab")
     
     ! Update internal history of mesh velocity
     call this%wm_lag_x%update()
     call this%wm_lag_y%update()
     call this%wm_lag_z%update()
  end subroutine advance_mesh

  !> Updates the mesh velocity field based on current time and kinematics
  subroutine update_mesh_velocity(this, coef, time_s, nadv)
     class(ale_manager_t), intent(inout) :: this
     type(coef_t), intent(in) :: coef
     type(time_state_t), intent(in) :: time_s
     integer, intent(in) :: nadv

     if (.not. this%active) return
     if (.not. this%has_moving_boundary) return

     call update_ale_mesh_velocity(this%wm_x, this%wm_y, this%wm_z, &
                              this%base_shape, coef, time_s, this%config, &
                              this%ale_pivot, nadv)
                              
  end subroutine update_mesh_velocity

  subroutine ale_manager_free(this)
    class(ale_manager_t), intent(inout) :: this
    
    if (.not. this%active) return

    if (allocated(this%Ax)) deallocate(this%Ax)
    if (allocated(this%ksp)) then
       call this%ksp%free()
       deallocate(this%ksp)
    end if
    if (allocated(this%pc)) then
       call precon_destroy(this%pc)
       deallocate(this%pc)
    end if
    if (allocated(this%ksp_solver)) deallocate(this%ksp_solver)
    if (allocated(this%precon_type)) deallocate(this%precon_type)
    call this%bc_moving%free()
    call this%bc_fixed%free()
    call this%bc_list%free()

    ! wm_x, wm_y, wm_z are managed by registry, no need to deallocate here (right?)
    call this%base_shape%free()
    call this%wm_lag_x%free()
    call this%wm_lag_y%free()
    call this%wm_lag_z%free()
    if (associated(this%Blag)) deallocate(this%Blag)
    if (associated(this%Blaglag)) deallocate(this%Blaglag)

  end subroutine ale_manager_free

  !> Factory for ALE Preconditioner
  subroutine ale_precon_factory(pc, ksp, coef, dof, gs, bclst, pctype)
    class(pc_t), allocatable, target, intent(inout)    :: pc
    class(ksp_t), target, intent(inout)                :: ksp
    type(coef_t), target, intent(in)                   :: coef
    type(dofmap_t), target, intent(in)                 :: dof
    type(gs_t), target, intent(inout)                  :: gs
    type(bc_list_t), target, intent(inout)             :: bclst
    character(len=*), intent(in)                       :: pctype
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
  end subroutine ale_precon_factory

  !> Compute base mesh stiffness
  subroutine compute_stiffness_ale(coef, params)
    type(coef_t), intent(inout) :: coef
    type(ale_config_t), intent(in) :: params
    integer :: i, n
    real(kind=rp) :: x, y, rr, arg

    n = coef%dof%size()

    select case (trim(params%stiffness_type))
    case ('built-in')
       do concurrent (i = 1:n)
          x = coef%dof%x(i,1,1,1)
          y = coef%dof%y(i,1,1,1)
          rr = (x-0.0_rp)**2 + (y-0.0_rp)**2
          arg = -rr / (params%stiffness_decay**2)
          coef%h1(i,1,1,1) = 1.0_rp + params%stiffness_gain * exp(arg)
          coef%h2(i,1,1,1) = 0.0_rp
       end do
    case default
       call neko_error("ALE Manager: Unknown stiffness type")
    end select
    
    coef%ifh2 = .false.
  end subroutine compute_stiffness_ale

end module ale_manager