!> Implements `boundary_profile_t` to capture and output instantaneous boundary profiles.
module boundary_profile
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use registry, only : neko_registry
  use scratch_registry, only : neko_scratch_registry
  use field, only : field_t
  use case, only : case_t
  use json_utils, only : json_get_or_default, json_get_or_lookup
  use time_state, only : time_state_t
  use coefs, only : coef_t
  use dirichlet, only : dirichlet_t
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error
  use vector, only : vector_t
  use time_based_controller, only : time_based_controller_t
  use drag_torque, only : setup_normals
  use operators, only : strain_rate
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_rank, pe_size
  use mpi_f08, only : MPI_Comm_rank, MPI_Comm_size, MPI_Gather, MPI_Gatherv, &
       MPI_Exscan, MPI_INTEGER, MPI_SUM
  use math, only : masked_gather_copy_0
  use ale_manager, only : neko_ale
  use file, only : file_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use device_math, only : device_masked_gather_copy_0
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST

  implicit none
  private

  !> A simulation component for boundary profiles (pressure, wall shear)
  type, public, extends(simulation_component_t) :: boundary_profile_t
     type(field_t), pointer :: u => null()
     type(field_t), pointer :: v => null()
     type(field_t), pointer :: w => null()
     type(field_t), pointer :: p => null()
     type(field_t), pointer :: mu => null()
     type(coef_t), pointer :: coef => null()

     type(dirichlet_t) :: bc
     integer :: zone_id
     character(len=60) :: zone_name
     real(kind=rp) :: start_time = 0.0_rp
     logical :: ale_enabled = .false.
     logical :: format_short = .true.

     ! Output variable flags
     logical :: out_p = .true.
     logical :: out_tx = .true.
     logical :: out_ty = .true.
     logical :: out_tz = .true.

     integer, allocatable :: gll_id(:)
     type(vector_t) :: x_ref, y_ref, z_ref
     type(vector_t) :: r1, r2, r3, n1, n2, n3

     type(vector_t) :: p_work, mu_work
     type(vector_t) :: s11msk, s22msk, s33msk, s12msk, s13msk, s23msk

     character(len=256) :: data_filename

   contains
     procedure, pass(this) :: init => boundary_profile_init_from_json
     generic :: init_from_components => init_from_controllers, init_from_controllers_properties
     procedure, pass(this) :: init_from_controllers => boundary_profile_init_from_controllers
     procedure, pass(this) :: init_from_controllers_properties => boundary_profile_init_from_controllers_properties
     procedure, private, pass(this) :: init_common => boundary_profile_init_common
     procedure, pass(this) :: free => boundary_profile_free
     procedure, pass(this) :: compute_ => boundary_profile_compute
     procedure, pass(this) :: write_profile => boundary_profile_write
  end type boundary_profile_t

contains

  subroutine boundary_profile_init_from_json(this, json, case)
    class(boundary_profile_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case

    character(len=:), allocatable :: name, fluid_name, zone_name, var_str
    integer :: zone_id
    real(kind=rp) :: start_time
    logical :: format_long

    call this%free()
    call json_get_or_default(json, "name", name, "boundary_profile")
    call this%init_base(json, case)

    call json_get_or_default(json, 'fluid_name', fluid_name, 'fluid')
    call json_get_or_lookup(json, 'zone_id', zone_id)
    call json_get_or_default(json, 'zone_name', zone_name, 'profile_zone')
    call json_get_or_default(json, 'start_time', start_time, 0.0_rp)

    call json_get_or_default(json, 'format_long', format_long, .false.)
    if (format_long) then
       this%format_short = .false.
    else
       this%format_short = .true.
    end if

    call json_get_or_default(json, 'variables', var_str, 'all')
    this%out_p = (index(var_str, 'p') > 0 .or. index(var_str, 'all') > 0)
    this%out_tx = (index(var_str, 'tau_x') > 0 .or. index(var_str, 'all') > 0)
    this%out_ty = (index(var_str, 'tau_y') > 0 .or. index(var_str, 'all') > 0)
    this%out_tz = (index(var_str, 'tau_z') > 0 .or. index(var_str, 'all') > 0)

    call this%init_common(name, fluid_name, zone_id, zone_name, &
         start_time, case%fluid%c_Xh)

  end subroutine boundary_profile_init_from_json

  subroutine boundary_profile_init_from_controllers(this, name, case, order, &
        preprocess_controller, compute_controller, output_controller, &
        fluid_name, zone_id, zone_name, start_time, coef)
    class(boundary_profile_t), intent(inout) :: this
    character(len=*), intent(in) :: name, fluid_name, zone_name
    class(case_t), intent(inout), target :: case
    integer, intent(in) :: order, zone_id
    real(kind=rp), intent(in) :: start_time
    type(time_based_controller_t), intent(in) :: preprocess_controller, compute_controller, output_controller
    type(coef_t), intent(inout), target :: coef

    call this%free()
    call this%init_base_from_components(case, order, preprocess_controller, compute_controller, output_controller)
    call this%init_common(name, fluid_name, zone_id, zone_name, start_time, coef)
  end subroutine boundary_profile_init_from_controllers

  subroutine boundary_profile_init_from_controllers_properties(this, name, &
        case, order, preprocess_control, preprocess_value, compute_control, &
        compute_value, output_control, output_value, fluid_name, zone_id, &
        zone_name, start_time, coef)
    class(boundary_profile_t), intent(inout) :: this
    character(len=*), intent(in) :: name, preprocess_control, compute_control, output_control
    character(len=*), intent(in) :: fluid_name, zone_name
    class(case_t), intent(inout), target :: case
    integer, intent(in) :: order, zone_id
    real(kind=rp), intent(in) :: preprocess_value, compute_value, output_value, start_time
    type(coef_t), intent(inout), target :: coef

    call this%free()
    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, output_value)
    call this%init_common(name, fluid_name, zone_id, zone_name, start_time, coef)
  end subroutine boundary_profile_init_from_controllers_properties

  subroutine boundary_profile_init_common(this, name, fluid_name, zone_id, zone_name, &
       start_time, coef)
    class(boundary_profile_t), intent(inout) :: this
    character(len=*), intent(in) :: name, fluid_name, zone_name
    integer, intent(in) :: zone_id
    real(kind=rp), intent(in) :: start_time
    type(coef_t), intent(inout), target :: coef
    integer :: n_pts, i, init_rank, init_ierr, init_offset, total_n
    integer, allocatable :: recvcounts(:), displs(:)
    integer, allocatable :: local_id_buf(:), global_id_buf(:)
    real(kind=rp), allocatable :: local_mesh_buf(:,:), global_mesh_buf(:,:)
    character(len=256) :: filename_data, filename_mesh
    character(len=1000) :: log_buf
    character(len=512) :: header_str
    integer :: io_unit

    this%name = name
    this%zone_id = zone_id
    this%zone_name = zone_name
    this%start_time = start_time
    this%coef => coef
    this%ale_enabled = .false.
    if (associated(neko_ale)) then
       if (neko_ale%active) this%ale_enabled = .true.
    end if

    this%u => neko_registry%get_field_by_name("u")
    this%v => neko_registry%get_field_by_name("v")
    this%w => neko_registry%get_field_by_name("w")
    this%p => neko_registry%get_field_by_name("p")
    this%mu => neko_registry%get_field_by_name(trim(fluid_name) // '_mu_tot')

    call this%bc%init_base(this%coef)
    call this%bc%mark_zone(this%case%msh%labeled_zones(this%zone_id))
    call this%bc%finalize()

    n_pts = this%bc%msk(0)

    if (n_pts .gt. 0) then
       call this%r1%init(n_pts)
       call this%r2%init(n_pts)
       call this%r3%init(n_pts)
       call this%n1%init(n_pts)
       call this%n2%init(n_pts)
       call this%n3%init(n_pts)

       if (this%out_p) call this%p_work%init(n_pts)

       if (this%out_tx .or. this%out_ty .or. this%out_tz) then
          call this%mu_work%init(n_pts)
          call this%s11msk%init(n_pts)
          call this%s22msk%init(n_pts)
          call this%s33msk%init(n_pts)
          call this%s12msk%init(n_pts)
          call this%s13msk%init(n_pts)
          call this%s23msk%init(n_pts)
       end if

       call this%x_ref%init(n_pts)
       call this%y_ref%init(n_pts)
       call this%z_ref%init(n_pts)

       call setup_normals(this%coef, this%bc%msk, this%bc%facet, this%n1%x, this%n2%x, this%n3%x, n_pts)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_masked_gather_copy_0(this%r1%x_d, &
               this%coef%dof%x_d, this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_gather_copy_0(this%r2%x_d, &
               this%coef%dof%y_d, this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_gather_copy_0(this%r3%x_d, &
               this%coef%dof%z_d, this%bc%msk_d, this%u%size(), n_pts)

          call device_memcpy(this%r1%x, this%r1%x_d, n_pts, DEVICE_TO_HOST, .false.)
          call device_memcpy(this%r2%x, this%r2%x_d, n_pts, DEVICE_TO_HOST, .false.)
          call device_memcpy(this%r3%x, this%r3%x_d, n_pts, DEVICE_TO_HOST, .false.)

          call device_memcpy(this%n1%x, this%n1%x_d, n_pts, HOST_TO_DEVICE, .false.)
          call device_memcpy(this%n2%x, this%n2%x_d, n_pts, HOST_TO_DEVICE, .false.)
          call device_memcpy(this%n3%x, this%n3%x_d, n_pts, HOST_TO_DEVICE, .true.)
       else
          call masked_gather_copy_0(this%r1%x, this%coef%dof%x, &
             this%bc%msk, this%u%size(), n_pts)
          call masked_gather_copy_0(this%r2%x, this%coef%dof%y, &
             this%bc%msk, this%u%size(), n_pts)
          call masked_gather_copy_0(this%r3%x, this%coef%dof%z, &
             this%bc%msk, this%u%size(), n_pts)
       end if

       this%x_ref%x = this%r1%x
       this%y_ref%x = this%r2%x
       this%z_ref%x = this%r3%x

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(this%x_ref%x, this%x_ref%x_d, n_pts, HOST_TO_DEVICE, .false.)
          call device_memcpy(this%y_ref%x, this%y_ref%x_d, n_pts, HOST_TO_DEVICE, .false.)
          call device_memcpy(this%z_ref%x, this%z_ref%x_d, n_pts, HOST_TO_DEVICE, .true.)
       end if
    end if

    init_rank = pe_rank
    init_offset = 0
    call MPI_Exscan(n_pts, init_offset, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, init_ierr)
    if (init_rank == 0) init_offset = 0

    allocate(this%gll_id(max(0, n_pts)))
    do i = 1, n_pts
       this%gll_id(i) = init_offset + i
    end do

    if (pe_rank == 0) allocate(recvcounts(pe_size), displs(pe_size))
    call MPI_Gather(n_pts, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, 0, NEKO_COMM, init_ierr)

    init_offset = 0
    call MPI_Exscan(n_pts, init_offset, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, init_ierr)
    if (pe_rank == 0) init_offset = 0
    call MPI_Gather(init_offset, 1, MPI_INTEGER, displs, 1, MPI_INTEGER, 0, NEKO_COMM, init_ierr)

    allocate(local_id_buf(max(1, n_pts)))
    allocate(local_mesh_buf(6, max(1, n_pts)))
    do i = 1, n_pts
       local_id_buf(i) = this%gll_id(i)
       local_mesh_buf(1, i) = this%x_ref%x(i)
       local_mesh_buf(2, i) = this%y_ref%x(i)
       local_mesh_buf(3, i) = this%z_ref%x(i)
       local_mesh_buf(4, i) = -this%n1%x(i)
       local_mesh_buf(5, i) = -this%n2%x(i)
       local_mesh_buf(6, i) = -this%n3%x(i)
    end do

    if (pe_rank == 0) then
       total_n = sum(recvcounts)
       allocate(global_id_buf(max(1, total_n)))
       allocate(global_mesh_buf(6, max(1, total_n)))
    else
       allocate(global_id_buf(1), global_mesh_buf(1, 1))
    end if

    call MPI_Gatherv(local_id_buf, n_pts, MPI_INTEGER, &
         global_id_buf, recvcounts, displs, MPI_INTEGER, &
         0, NEKO_COMM, init_ierr)

    if (pe_rank == 0) then
       recvcounts = recvcounts * 6
       displs = displs * 6
    end if

    call MPI_Gatherv(local_mesh_buf, n_pts * 6, MPI_REAL_PRECISION, &
         global_mesh_buf, recvcounts, displs, MPI_REAL_PRECISION, &
         0, NEKO_COMM, init_ierr)

    if (pe_rank == 0 .and. total_n > 0) then
       write(filename_mesh, '("boundary_profile_",A,"_mesh.csv")') trim(this%name)
       open(newunit=io_unit, file=this%case%output_directory // trim(filename_mesh), &
            status='replace', action='write')

       write(io_unit, '(A)') "gll_id,x_ref,y_ref,z_ref,-nrm_x_ref,-nrm_y_ref,-nrm_z_ref"

       do i = 1, total_n
          write(io_unit, '(I0,A,6(G0.15,A))') &
             global_id_buf(i), ',', &
             global_mesh_buf(1, i), ',', global_mesh_buf(2, i), ',', global_mesh_buf(3, i), ',', &
             global_mesh_buf(4, i), ',', global_mesh_buf(5, i), ',', global_mesh_buf(6, i)
       end do
       close(io_unit)
    end if

    deallocate(local_id_buf, global_id_buf)
    deallocate(local_mesh_buf, global_mesh_buf)
    if (pe_rank == 0) deallocate(recvcounts, displs)

    write(this%data_filename, '("boundary_profile_",A,"_data.csv")') trim(this%name)

    if (pe_rank == 0) then
       open(newunit=io_unit, file=this%case%output_directory // trim(this%data_filename), &
            status='replace', action='write')

       ! Dynamically build the header
       header_str = "t,gll_id"
       if (this%ale_enabled) header_str = trim(header_str) // ",x,y,z,-nrm_x,-nrm_y,-nrm_z"
       if (this%out_p) header_str = trim(header_str) // ",p"
       if (this%out_tx) header_str = trim(header_str) // ",tau_x"
       if (this%out_ty) header_str = trim(header_str) // ",tau_y"
       if (this%out_tz) header_str = trim(header_str) // ",tau_z"

       write(io_unit, '(A)') trim(header_str)
       close(io_unit)
    end if

    call neko_log%section("Boundary Profile Operation")
    write(log_buf, '(A,A)') "Name: ", trim(this%name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,I0,A,A)') "Tracking Zone: ", this%zone_id, " ", trim(this%zone_name)
    call neko_log%message(log_buf)

    ! Log active variables
    header_str = "Output Variables: "
    if (this%out_p) header_str = trim(header_str) // " p"
    if (this%out_tx) header_str = trim(header_str) // " tau_x"
    if (this%out_ty) header_str = trim(header_str) // " tau_y"
    if (this%out_tz) header_str = trim(header_str) // " tau_z"
    call neko_log%message(trim(header_str))

    call neko_log%end_section()

  end subroutine boundary_profile_init_common

  subroutine boundary_profile_free(this)
    class(boundary_profile_t), intent(inout) :: this

    call this%bc%free()
    call this%r1%free(); call this%r2%free(); call this%r3%free()
    call this%n1%free(); call this%n2%free(); call this%n3%free()

    if (this%out_p) call this%p_work%free()
    if (this%out_tx .or. this%out_ty .or. this%out_tz) then
       call this%mu_work%free()
       call this%s11msk%free()
       call this%s22msk%free()
       call this%s33msk%free()
       call this%s12msk%free()
       call this%s13msk%free()
       call this%s23msk%free()
    end if

    call this%x_ref%free()
    call this%y_ref%free()
    call this%z_ref%free()
    if (allocated(this%gll_id)) deallocate(this%gll_id)

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)
    nullify(this%mu)
    nullify(this%coef)

    call this%free_base()
  end subroutine boundary_profile_free

  subroutine boundary_profile_compute(this, time)
    class(boundary_profile_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp) :: dt

    if (time%t < this%start_time) return

    if (this%compute_controller%check(time)) then
       dt = time%dt
       if (dt > 0.0_rp) then
          call this%write_profile(time)
       end if
       call this%compute_controller%register_execution()
    end if
  end subroutine boundary_profile_compute

  subroutine boundary_profile_write(this, time)
    class(boundary_profile_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    integer :: i, j, local_n, total_n, rank, num_procs, ierr, local_offset, num_cols
    integer :: temp_indices(6), idx
    integer, allocatable :: recvcounts(:), displs(:)
    integer, allocatable :: local_id_buf(:), global_id_buf(:)
    real(kind=rp), allocatable :: local_buffer(:,:), global_buffer(:,:)
    type(field_t), pointer :: s11, s22, s33, s12, s13, s23
    real(kind=rp) :: total_visc_x, total_visc_y, total_visc_z, dot_v_n, mu_val, p_val
    integer :: io_unit
    logical :: compute_tau
    character(len=256) :: fmt_str

    rank = pe_rank
    num_procs = pe_size
    local_n = this%bc%msk(0)
    compute_tau = this%out_tx .or. this%out_ty .or. this%out_tz

    num_cols = 0
    if (this%ale_enabled) num_cols = 6
    if (this%out_p) num_cols = num_cols + 1
    if (this%out_tx) num_cols = num_cols + 1
    if (this%out_ty) num_cols = num_cols + 1
    if (this%out_tz) num_cols = num_cols + 1

    if (compute_tau) then
       call neko_scratch_registry%request_field(s11, temp_indices(1), .false.)
       call neko_scratch_registry%request_field(s12, temp_indices(2), .false.)
       call neko_scratch_registry%request_field(s13, temp_indices(3), .false.)
       call neko_scratch_registry%request_field(s22, temp_indices(4), .false.)
       call neko_scratch_registry%request_field(s23, temp_indices(5), .false.)
       call neko_scratch_registry%request_field(s33, temp_indices(6), .false.)
       call strain_rate(s11%x, s22%x, s33%x, s12%x, s13%x, s23%x, this%u, this%v, this%w, this%coef)
    end if


    if (this%ale_enabled .and. local_n > 0) then
       ! get normals
       call setup_normals(this%coef, this%bc%msk, this%bc%facet, this%n1%x, this%n2%x, this%n3%x, local_n)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(this%n1%x, this%n1%x_d, local_n, HOST_TO_DEVICE, .false.)
          call device_memcpy(this%n2%x, this%n2%x_d, local_n, HOST_TO_DEVICE, .false.)
          call device_memcpy(this%n3%x, this%n3%x_d, local_n, HOST_TO_DEVICE, .false.)

          call device_masked_gather_copy_0(this%r1%x_d, this%coef%dof%x_d, &
               this%bc%msk_d, this%u%size(), local_n)
          call device_masked_gather_copy_0(this%r2%x_d, this%coef%dof%y_d, &
               this%bc%msk_d, this%u%size(), local_n)
          call device_masked_gather_copy_0(this%r3%x_d, this%coef%dof%z_d, &
               this%bc%msk_d, this%u%size(), local_n)


          call device_memcpy(this%r1%x, this%r1%x_d, local_n, HOST_TO_DEVICE, .false.)
          call device_memcpy(this%r2%x, this%r2%x_d, local_n, HOST_TO_DEVICE, .false.)
          call device_memcpy(this%r3%x, this%r3%x_d, local_n, HOST_TO_DEVICE, .false.)
       else
          call masked_gather_copy_0(this%r1%x, this%coef%dof%x, &
               this%bc%msk, this%u%size(), local_n)
          call masked_gather_copy_0(this%r2%x, this%coef%dof%y, &
               this%bc%msk, this%u%size(), local_n)
          call masked_gather_copy_0(this%r3%x, this%coef%dof%z, &
               this%bc%msk, this%u%size(), local_n)
       end if
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       if (this%out_p .and. local_n > 0) then
          call device_masked_gather_copy_0(this%p_work%x_d, this%p%x_d, this%bc%msk_d, this%p%size(), local_n)
          call device_memcpy(this%p_work%x, this%p_work%x_d, local_n, DEVICE_TO_HOST, .true.)
       end if

       if (compute_tau .and. local_n > 0) then
          call device_masked_gather_copy_0(this%mu_work%x_d, this%mu%x_d, this%bc%msk_d, this%mu%size(), local_n)
          call device_masked_gather_copy_0(this%s11msk%x_d, s11%x_d, this%bc%msk_d, s11%size(), local_n)
          call device_masked_gather_copy_0(this%s22msk%x_d, s22%x_d, this%bc%msk_d, s22%size(), local_n)
          call device_masked_gather_copy_0(this%s33msk%x_d, s33%x_d, this%bc%msk_d, s33%size(), local_n)
          call device_masked_gather_copy_0(this%s12msk%x_d, s12%x_d, this%bc%msk_d, s12%size(), local_n)
          call device_masked_gather_copy_0(this%s13msk%x_d, s13%x_d, this%bc%msk_d, s13%size(), local_n)
          call device_masked_gather_copy_0(this%s23msk%x_d, s23%x_d, this%bc%msk_d, s23%size(), local_n)

          call device_memcpy(this%mu_work%x, this%mu_work%x_d, local_n, DEVICE_TO_HOST, .false.)
          call device_memcpy(this%s11msk%x, this%s11msk%x_d, local_n, DEVICE_TO_HOST, .false.)
          call device_memcpy(this%s22msk%x, this%s22msk%x_d, local_n, DEVICE_TO_HOST, .false.)
          call device_memcpy(this%s33msk%x, this%s33msk%x_d, local_n, DEVICE_TO_HOST, .false.)
          call device_memcpy(this%s12msk%x, this%s12msk%x_d, local_n, DEVICE_TO_HOST, .false.)
          call device_memcpy(this%s13msk%x, this%s13msk%x_d, local_n, DEVICE_TO_HOST, .false.)
          call device_memcpy(this%s23msk%x, this%s23msk%x_d, local_n, DEVICE_TO_HOST, .true.)
       end if
    else
       if (this%out_p .and. local_n > 0) then
          call masked_gather_copy_0(this%p_work%x, this%p%x, this%bc%msk, this%p%size(), local_n)
       end if

       if (compute_tau .and. local_n > 0) then
          call masked_gather_copy_0(this%mu_work%x, this%mu%x, this%bc%msk, this%mu%size(), local_n)
          call masked_gather_copy_0(this%s11msk%x, s11%x, this%bc%msk, s11%size(), local_n)
          call masked_gather_copy_0(this%s22msk%x, s22%x, this%bc%msk, s22%size(), local_n)
          call masked_gather_copy_0(this%s33msk%x, s33%x, this%bc%msk, s33%size(), local_n)
          call masked_gather_copy_0(this%s12msk%x, s12%x, this%bc%msk, s12%size(), local_n)
          call masked_gather_copy_0(this%s13msk%x, s13%x, this%bc%msk, s13%size(), local_n)
          call masked_gather_copy_0(this%s23msk%x, s23%x, this%bc%msk, s23%size(), local_n)
       end if
    end if

    allocate(local_id_buf(max(1, local_n)))
    allocate(local_buffer(max(1, num_cols), max(1, local_n)))

    do i = 1, local_n
       idx = 1
       local_id_buf(i) = this%gll_id(i)

       if (this%ale_enabled) then
          local_buffer(idx, i) = this%r1%x(i)
          idx = idx + 1
          local_buffer(idx, i) = this%r2%x(i)
          idx = idx + 1
          local_buffer(idx, i) = this%r3%x(i)
          idx = idx + 1
          local_buffer(idx, i) = -this%n1%x(i)
          idx = idx + 1
          local_buffer(idx, i) = -this%n2%x(i)
          idx = idx + 1
          local_buffer(idx, i) = -this%n3%x(i)
          idx = idx + 1
       end if

       if (this%out_p) then
          local_buffer(idx, i) = this%p_work%x(i)
          idx = idx + 1
       end if

       if (compute_tau) then
          mu_val = this%mu_work%x(i)

          total_visc_x = mu_val * 2.0_rp * &
          (this%s11msk%x(i)*this%n1%x(i) + this%s12msk%x(i)*this%n2%x(i) + this%s13msk%x(i)*this%n3%x(i))
          total_visc_y = mu_val * 2.0_rp * &
          (this%s12msk%x(i)*this%n1%x(i) + this%s22msk%x(i)*this%n2%x(i) + this%s23msk%x(i)*this%n3%x(i))
          total_visc_z = mu_val * 2.0_rp * &
          (this%s13msk%x(i)*this%n1%x(i) + this%s23msk%x(i)*this%n2%x(i) + this%s33msk%x(i)*this%n3%x(i))

          dot_v_n = total_visc_x*this%n1%x(i) + total_visc_y*this%n2%x(i) + total_visc_z*this%n3%x(i)

          if (this%out_tx) then
             local_buffer(idx, i) = (total_visc_x - dot_v_n*this%n1%x(i)); idx = idx + 1
          end if
          if (this%out_ty) then
             local_buffer(idx, i) = (total_visc_y - dot_v_n*this%n2%x(i)); idx = idx + 1
          end if
          if (this%out_tz) then
             local_buffer(idx, i) = (total_visc_z - dot_v_n*this%n3%x(i)); idx = idx + 1
          end if
       end if
    end do

    if (compute_tau) call neko_scratch_registry%relinquish_field(temp_indices)

    if (rank == 0) allocate(recvcounts(num_procs), displs(num_procs))
    call MPI_Gather(local_n, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)

    local_offset = 0
    call MPI_Exscan(local_n, local_offset, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    if (rank == 0) local_offset = 0
    call MPI_Gather(local_offset, 1, MPI_INTEGER, displs, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)

    if (rank == 0) then
       total_n = sum(recvcounts)
       allocate(global_id_buf(max(1, total_n)))
       allocate(global_buffer(max(1, num_cols), max(1, total_n)))
    else
       allocate(global_id_buf(1), global_buffer(1, 1))
    end if

    call MPI_Gatherv(local_id_buf, local_n, MPI_INTEGER, &
         global_id_buf, recvcounts, displs, MPI_INTEGER, &
         0, NEKO_COMM, ierr)

    if (rank == 0) then
       recvcounts = recvcounts * max(1, num_cols)
       displs = displs * max(1, num_cols)
    end if

    if (num_cols > 0) then
       call MPI_Gatherv(local_buffer, local_n * num_cols, MPI_REAL_PRECISION, &
            global_buffer, recvcounts, displs, MPI_REAL_PRECISION, &
            0, NEKO_COMM, ierr)
    end if

    if (rank == 0 .and. total_n > 0) then
       fmt_str = "(G0.15,A,I0"
       if (num_cols > 0) then
          fmt_str = trim(fmt_str) // ",A"
          do j = 1, num_cols
             if (this%ale_enabled .and. j <= 3) then
                fmt_str = trim(fmt_str) // ",G0.15" ! Coords always long
             else if (this%ale_enabled .and. j <= 6) then
                if (this%format_short) then
                   fmt_str = trim(fmt_str) // ",G0.6"
                else
                   fmt_str = trim(fmt_str) // ",G0.15"
                end if
             else
                if (this%format_short) then
                   fmt_str = trim(fmt_str) // ",G0.6"
                else
                   fmt_str = trim(fmt_str) // ",G0.15"
                end if
             end if

             if (j < num_cols) fmt_str = trim(fmt_str) // ",A"
          end do
       end if
       fmt_str = trim(fmt_str) // ")"

       open(newunit=io_unit, file=this%case%output_directory // trim(this%data_filename), &
            position='append', action='write')

       do i = 1, total_n
          if (num_cols > 0) then
             write(io_unit, fmt_str) time%t, ',', global_id_buf(i), ',', &
                 (global_buffer(j, i), ',', j=1, num_cols-1), global_buffer(num_cols, i)
          else
             write(io_unit, fmt_str) time%t, ',', global_id_buf(i)
          end if
       end do
       close(io_unit)
    end if

    deallocate(local_id_buf, global_id_buf)
    deallocate(local_buffer, global_buffer)
    if (rank == 0) deallocate(recvcounts, displs)

  end subroutine boundary_profile_write
end module boundary_profile
