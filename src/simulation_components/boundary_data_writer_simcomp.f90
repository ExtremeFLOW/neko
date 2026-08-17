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
!> Implements the `boundary_data_writer_t` type.
module boundary_data_writer_simcomp
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use simulation_component, only : simulation_component_t
  use time_based_controller, only : time_based_controller_t
  use time_state, only : time_state_t
  use case, only : case_t
  use field, only : field_t
  use field_list, only : field_list_t
  use registry, only : neko_registry
  use coefs, only : coef_t
  use vector, only : vector_t
  use device, only : DEVICE_TO_HOST
  use math, only : copy
  use tensor, only : trsp
  use matrix, only : matrix_t
  use boundary_data, only : boundary_data_t
  use file, only : file_t
  use csv_file, only : csv_file_t
  use hdf5_file, only : hdf5_file_t
  use ale_manager, only : neko_ale
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error, neko_warning, NEKO_VARNAME_LEN, &
       NEKO_FNAME_LEN, filename_suffix, filename_suffix_pos
  use comm, only : NEKO_COMM, pe_rank, pe_size, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_Gather, MPI_Gatherv, MPI_Exscan, &
       MPI_INTEGER, MPI_SUM
  implicit none
  private

  !> A simulation component that writes registry fields sampled on labelled
  !! boundary zones.
  !! The component builds a mask over the requested labelled zones and, at
  !! every output, samples the requested quantities at the boundary GLL points.
  !! The result is written either as CSV or as HDF5.
  type, public, extends(simulation_component_t) :: boundary_data_writer_t
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> Boundary points of the chosen zones, with their geometry.
     type(boundary_data_t) :: bdata
     !> Labelled zones to include.
     integer, allocatable :: zone_indices(:)
     !> Fields sampled directly.
     type(field_list_t) :: fields
     !> Names of the output columns.
     character(len=NEKO_VARNAME_LEN), allocatable :: col_names(:)
     !> Reusable gathered values.
     type(vector_t) :: work
     !> Number of masked points on this rank.
     integer :: n_local = 0
     !> Total number of masked points.
     integer :: n_global = 0
     !> Number of data columns, excluding time.
     integer :: n_cols = 0
     !> Local sample buffer (point, column).
     real(kind=rp), allocatable :: local_buffer(:,:)
     !> Column-blocked transpose of `local_buffer`.
     real(kind=rp), allocatable :: local_buffer_t(:,:)
     !> Global sample buffer on rank zero.
     real(kind=rp), allocatable :: global_buffer(:,:)
     !> Receive counts for the gather onto rank zero.
     integer, allocatable :: recvcounts(:)
     !> Displacements for the gather onto rank zero.
     integer, allocatable :: displs(:)
     integer, allocatable :: recvcounts_c(:), displs_c(:)
     !> Output matrix.
     type(matrix_t) :: mat_out
     !> Output vector, used for the HDF5.
     type(vector_t) :: vec_out
     !> The output file.
     type(file_t) :: fout
     !> Time after which to start writing.
     real(kind=rp) :: start_time = 0.0_rp
     !> Whether the mesh geometry changes during the run.
     logical :: mesh_has_changed = .false.
     !> Whether the geometry is part of every sample.
     logical :: geometry_in_data = .false.
     !> Whether to output the unit normals.
     logical :: output_normals = .false.
     !> Whether to output the surface quadrature weights (area).
     logical :: output_area = .false.
     !> Append every sample into shared datasets (.true.), or write
     !! each sample into its own `Step_N` subgroup (.false.)
     logical :: append_out = .true.
   contains
     !> Constructors.
     procedure, pass(this) :: init => boundary_data_writer_init_from_json
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     procedure, pass(this) :: init_from_controllers => &
          boundary_data_writer_init_from_controllers
     procedure, pass(this) :: init_from_controllers_properties => &
          boundary_data_writer_init_from_controllers_properties
     procedure, private, pass(this) :: init_common => &
          boundary_data_writer_init_common
     !> Destructor.
     procedure, pass(this) :: free => boundary_data_writer_free
     !> Sample the boundary and write.
     procedure, pass(this) :: compute_ => boundary_data_writer_compute
     !> Gather a masked quantity into a column of the local buffer.
     procedure, private, pass(this) :: gather_column => &
          boundary_data_writer_gather_column
  end type boundary_data_writer_t

contains

  !> Constructor from json.
  !! @param json JSON object with the parameters.
  !! @param case The case object.
  subroutine boundary_data_writer_init_from_json(this, json, case)
    class(boundary_data_writer_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: name, output_filename
    character(len=NEKO_VARNAME_LEN), allocatable :: fields(:)
    integer, allocatable :: zone_indices(:)
    real(kind=rp) :: start_time
    logical :: output_normals
    logical :: output_area
    logical :: append_out
    logical :: user_set_compute

    call this%free()

    call json_get_or_default(json, "name", name, "boundary_data_writer")

    user_set_compute = json%valid_path("compute_control") .or. &
         json%valid_path("compute_value")

    call this%init_base(json, case)

    if (user_set_compute) then
       call neko_warning("boundary_data_writer: 'compute_control' and " // &
            "'compute_value' are ignored. Use 'output_control' and " // &
            "'output_value' to set how often the data is written.")
    end if
    call this%compute_controller%init(case%time%start_time, &
         case%time%end_time, "tsteps", 1.0_rp)

    call json_get(json, "zone_indices", zone_indices)
    call json_get(json, "output_filename", output_filename)
    call json_get_or_default(json, "output_normals", output_normals, .true.)
    call json_get_or_default(json, "output_area", output_area, .true.)
    call json_get_or_default(json, "append_output", append_out, .true.)
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)

    call json_get(json, "fields", fields)

    call this%init_common(name, case%fluid%c_Xh, zone_indices, fields, &
         output_filename, output_normals, output_area, append_out, &
         start_time)

  end subroutine boundary_data_writer_init_from_json

  !> Constructor from components, passing controllers.
  !! @param name The unique name of the simcomp.
  !! @param case The simulation case object.
  !! @param order The execution order priority of the simcomp.
  !! @param preprocess_controller The controller for running preprocessing.
  !! @param compute_controller The controller for running compute.
  !! @param output_controller The controller for producing output.
  !! @param coef The SEM coefficients.
  !! @param zone_indices Labelled zones to include.
  !! @param fields Names of the registry fields to sample.
  !! @param output_filename The output file, either `.csv` or `.h5`.
  !! @param output_normals Whether to output the unit normals.
  !! @param output_area Whether to output the surface quadrature weights.
  !! @param append_out Whether to append samples to shared datasets.
  !! @param start_time Time after which to start writing.
  subroutine boundary_data_writer_init_from_controllers(this, name, case, &
       order, preprocess_controller, compute_controller, output_controller, &
       coef, zone_indices, fields, output_filename, output_normals, &
       output_area, append_out, start_time)
    class(boundary_data_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer, intent(in) :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    type(coef_t), intent(inout), target :: coef
    integer, intent(in) :: zone_indices(:)
    character(len=*), intent(in) :: fields(:)
    character(len=*), intent(in) :: output_filename
    logical, intent(in) :: output_normals
    logical, intent(in) :: output_area
    logical, intent(in) :: append_out
    real(kind=rp), intent(in) :: start_time

    call this%free()

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)
    call this%init_common(name, coef, zone_indices, fields, &
         output_filename, output_normals, output_area, append_out, start_time)

  end subroutine boundary_data_writer_init_from_controllers

  !> Constructor from components, passing properties to the
  !! `time_based_controller` components in the base type.
  !! @param name The unique name of the simcomp.
  !! @param case The simulation case object.
  !! @param order The execution order priority of the simcomp.
  !! @param preprocess_control Control mode for preprocessing.
  !! @param preprocess_value Value parameter for preprocessing.
  !! @param compute_control Control mode for computing.
  !! @param compute_value Value parameter for computing.
  !! @param output_control Control mode for output.
  !! @param output_value Value parameter for output.
  !! @param coef The SEM coefficients.
  !! @param zone_indices Labelled zones to include.
  !! @param fields Names of the registry fields to sample.
  !! @param output_filename The output file, either `.csv` or `.h5`.
  !! @param output_normals Whether to output the unit normals.
  !! @param output_area Whether to output the surface quadrature weights.
  !! @param append_out Whether to append samples to shared datasets.
  !! @param start_time Time after which to start writing.
  subroutine boundary_data_writer_init_from_controllers_properties(this, &
       name, case, order, preprocess_control, preprocess_value, &
       compute_control, compute_value, output_control, output_value, coef, &
       zone_indices, fields, output_filename, output_normals, output_area, &
       append_out, start_time)
    class(boundary_data_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer, intent(in) :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    type(coef_t), intent(inout), target :: coef
    integer, intent(in) :: zone_indices(:)
    character(len=*), intent(in) :: fields(:)
    character(len=*), intent(in) :: output_filename
    logical, intent(in) :: output_normals
    logical, intent(in) :: output_area
    logical, intent(in) :: append_out
    real(kind=rp), intent(in) :: start_time

    call this%free()

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)
    call this%init_common(name, coef, zone_indices, fields, &
         output_filename, output_normals, output_area, append_out, start_time)

  end subroutine boundary_data_writer_init_from_controllers_properties

  !> Common part of all constructors.
  !! @param name The unique name of the simcomp.
  !! @param coef The SEM coefficients.
  !! @param zone_indices Labelled zones to include.
  !! @param fields Names of the registry fields to sample.
  !! @param output_filename The output file, either `.csv` or `.h5`.
  !! @param output_normals Whether to output the unit normals.
  !! @param output_area Whether to output the surface quadrature weights.
  !! @param append_out Whether to append samples to shared datasets.
  !! @param start_time Time after which to start writing.
  subroutine boundary_data_writer_init_common(this, name, coef, zone_indices, &
       fields, output_filename, output_normals, output_area, append_out, &
       start_time)
    class(boundary_data_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(coef_t), intent(inout), target :: coef
    integer, intent(in) :: zone_indices(:)
    character(len=*), intent(in) :: fields(:)
    character(len=*), intent(in) :: output_filename
    logical, intent(in) :: output_normals
    logical, intent(in) :: output_area
    logical, intent(in) :: append_out
    real(kind=rp), intent(in) :: start_time
    character(len=LOG_SIZE) :: log_buf
    character(len=80) :: suffix
    integer :: i, col, ierr, offset, n_geom
    logical :: ale_enabled

    this%name = name
    this%coef => coef
    this%start_time = start_time
    this%output_normals = output_normals
    this%output_area = output_area
    this%append_out = append_out

    ! Zones
    allocate(this%zone_indices(size(zone_indices)))
    this%zone_indices = zone_indices

    ! Fields
    if (size(fields) .eq. 0) then
       call neko_error("boundary_data_writer requires at least one entry " // &
            "in 'fields'")
    end if

    call this%fields%init(size(fields))
    do i = 1, size(fields)
       this%fields%items(i)%ptr => &
            neko_registry%get_field_by_name(trim(fields(i)))
    end do

    ale_enabled = .false.
    if (associated(neko_ale)) then
       if (neko_ale%active) ale_enabled = .true.
    end if

    ! Whether the mesh geometry varies during the run.
    this%mesh_has_changed = .false.
    if (ale_enabled) then
       this%mesh_has_changed = .true.
    end if

    ! When the mesh moves the geometry accompanies every sample, so that
    ! each output carries the coordinates it corresponds to.
    this%geometry_in_data = this%mesh_has_changed

    ! only_facets
    call this%bdata%init(this%coef, this%zone_indices)

    this%n_local = this%bdata%n_local
    this%n_global = this%bdata%n_global

    call this%work%init(this%n_local)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call this%bdata%x%copy_from(DEVICE_TO_HOST, .false.)
       call this%bdata%y%copy_from(DEVICE_TO_HOST, .false.)
       call this%bdata%z%copy_from(DEVICE_TO_HOST, .false.)
       call this%bdata%n_x%copy_from(DEVICE_TO_HOST, .false.)
       call this%bdata%n_y%copy_from(DEVICE_TO_HOST, .false.)
       call this%bdata%n_z%copy_from(DEVICE_TO_HOST, .false.)
       call this%bdata%area%copy_from(DEVICE_TO_HOST, .true.)
    end if

    ! Column layout
    n_geom = 0
    if (this%geometry_in_data) then
       n_geom = 3
       if (this%output_normals) n_geom = n_geom + 3
       if (this%output_area) n_geom = n_geom + 1
    end if

    this%n_cols = n_geom + this%fields%size()

    allocate(this%col_names(this%n_cols))
    col = 0
    if (this%geometry_in_data) then
       col = col + 1
       this%col_names(col) = "x"
       col = col + 1
       this%col_names(col) = "y"
       col = col + 1
       this%col_names(col) = "z"
       if (this%output_normals) then
          col = col + 1
          this%col_names(col) = "n_x"
          col = col + 1
          this%col_names(col) = "n_y"
          col = col + 1
          this%col_names(col) = "n_z"
       end if
       if (this%output_area) then
          col = col + 1
          this%col_names(col) = "area"
       end if
    end if
    do i = 1, this%fields%size()
       col = col + 1
       this%col_names(col) = trim(fields(i))
    end do

    allocate(this%local_buffer(this%n_local, this%n_cols))
    this%local_buffer = 0.0_rp

    ! For CSV.
    allocate(this%local_buffer_t(this%n_cols, this%n_local))
    this%local_buffer_t = 0.0_rp

    ! Offsets for the gather onto rank zero.
    allocate(this%recvcounts(pe_size))
    allocate(this%displs(pe_size))
    this%recvcounts = 0
    this%displs = 0

    call MPI_Gather(this%n_local, 1, MPI_INTEGER, this%recvcounts, 1, &
         MPI_INTEGER, 0, NEKO_COMM, ierr)

    offset = 0
    call MPI_Exscan(this%n_local, offset, 1, MPI_INTEGER, MPI_SUM, &
         NEKO_COMM, ierr)
    if (pe_rank .eq. 0) offset = 0

    call MPI_Gather(offset, 1, MPI_INTEGER, this%displs, 1, MPI_INTEGER, &
         0, NEKO_COMM, ierr)

    allocate(this%recvcounts_c(pe_size))
    allocate(this%displs_c(pe_size))
    this%recvcounts_c = this%recvcounts * this%n_cols
    this%displs_c = this%displs * this%n_cols

    if (pe_rank .eq. 0) then
       allocate(this%global_buffer(max(1, this%n_cols), max(1, this%n_global)))
    else
       allocate(this%global_buffer(1, 1))
    end if
    this%global_buffer = 0.0_rp

    ! Output file
    call filename_suffix(output_filename, suffix)
    if (trim(suffix) .ne. "csv" .and. trim(suffix) .ne. "h5" .and. &
         trim(suffix) .ne. "hdf5") then
       call neko_error("boundary_data_writer: output_filename must end in " // &
            "'.csv', '.h5' or '.hdf5'")
    end if

    call this%fout%init(this%case%output_directory // trim(output_filename))

    n_geom = 3
    if (this%output_normals) n_geom = n_geom + 3
    if (this%output_area) n_geom = n_geom + 1

    select type (ft => this%fout%file_type)
    type is (csv_file_t)
       call boundary_data_writer_setup_output_csv(this, output_filename, &
            n_geom)
    class is (hdf5_file_t)
       call boundary_data_writer_setup_output_hdf5(this, ft, n_geom)
    class default
       call neko_error("boundary_data_writer: expected csv_file_t or " // &
            "hdf5_file_t")
    end select

    ! Log
    call neko_log%section("Boundary data writer")
    write(log_buf, '(A,A)') "Name: ", trim(this%name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,*(I0,:,", "))') "Zone indices: ", this%zone_indices
    call neko_log%message(log_buf)
    write(log_buf, '(A,I0)') "Global number of masked points: ", this%n_global
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') "Output file: ", trim(output_filename)
    call neko_log%message(log_buf)
    call neko_log%message("Columns:")
    do i = 1, this%n_cols
       write(log_buf, '(A,A)') "  ", trim(this%col_names(i))
       call neko_log%message(log_buf)
    end do
    write(log_buf, '(A,L1)') "Moving mesh: ", this%mesh_has_changed
    call neko_log%message(log_buf)
    call neko_log%end_section()

  end subroutine boundary_data_writer_init_common

  !> Set up a CSV output.
  !! @param output_filename The base output filename.
  !! @param n_geom Number of geometry entries per point in the companion file.
  subroutine boundary_data_writer_setup_output_csv(this, output_filename, &
       n_geom)
    class(boundary_data_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: output_filename
    integer, intent(in) :: n_geom
    character(len=1024) :: header_line
    character(len=NEKO_FNAME_LEN) :: mesh_filename
    type(matrix_t) :: mat_geom
    type(file_t) :: fmesh
    integer :: i, suffix_pos

    header_line = "time"
    do i = 1, this%n_cols
       header_line = trim(header_line) // "," // trim(this%col_names(i))
    end do
    call this%fout%set_header(trim(header_line))

    if (pe_rank .eq. 0) then
       call this%mat_out%init(this%n_global, this%n_cols)
    else
       call this%mat_out%init(1, 1)
    end if

    ! The geometry goes into a companion file. For a moving mesh the
    ! geometry appears in every sample, but the initial
    ! geometry is additionally written once as `_initial_mesh`.
    call boundary_data_writer_gather_geometry(this, mat_geom, n_geom)

    if (pe_rank .eq. 0) then
       header_line = "x,y,z"
       if (this%output_normals) header_line = &
            trim(header_line) // ",n_x,n_y,n_z"
       if (this%output_area) header_line = trim(header_line) // ",area"

       mesh_filename = this%case%output_directory // &
            trim(output_filename)
       suffix_pos = filename_suffix_pos(mesh_filename)
       if (this%geometry_in_data) then
          mesh_filename = trim(mesh_filename(1:suffix_pos-1)) // &
               "_initial_mesh" // trim(mesh_filename(suffix_pos:))
       else
          mesh_filename = trim(mesh_filename(1:suffix_pos-1)) // "_mesh" // &
               trim(mesh_filename(suffix_pos:))
       end if

       call fmesh%init(trim(mesh_filename), &
            header = trim(header_line), overwrite = .true.)
       call fmesh%write(mat_geom)
       call fmesh%free()
    end if
    call mat_geom%free()

  end subroutine boundary_data_writer_setup_output_csv

  !> Set up an HDF5 output.
  !! @param ft The HDF5 file handle.
  !! @param n_geom Number of geometry entries per point.
  subroutine boundary_data_writer_setup_output_hdf5(this, ft, n_geom)
    class(boundary_data_writer_t), intent(inout) :: this
    class(hdf5_file_t), intent(inout) :: ft
    integer, intent(in) :: n_geom
    type(matrix_t) :: mat_geom
    integer :: out_int
    logical :: attr_exist

    call ft%set_overwrite(.not. this%append_out)

    call this%vec_out%init(max(0, this%n_local), "value")

    call ft%open("w")
    call ft%set_active_group("boundary_data")

    ! If the file already carries samples we are restarting into it, so the
    ! geometry must not be appended a second time. Similar to probes.
    call ft%read_attribute("NSteps", out_int, attr_exist)
    if (attr_exist) then
       this%output_controller%nexecutions = out_int
    else
       out_int = this%n_global
       call ft%write_attribute("NPoints", out_int)
       call boundary_data_writer_local_geometry(this, mat_geom, n_geom)
       if (this%geometry_in_data) then
          mat_geom%name = "initial_coordinates"
       else
          mat_geom%name = "coordinates"
       end if
       call ft%write_dataset(mat_geom)
       call mat_geom%free()
    end if
    call ft%close()

  end subroutine boundary_data_writer_setup_output_hdf5

  !> Destructor.
  subroutine boundary_data_writer_free(this)
    class(boundary_data_writer_t), intent(inout) :: this

    call this%bdata%free()

    call this%work%free()
    call this%vec_out%free()
    call this%mat_out%free()

    call this%fields%free()
    call this%fout%free()

    if (allocated(this%zone_indices)) deallocate(this%zone_indices)
    if (allocated(this%col_names)) deallocate(this%col_names)
    if (allocated(this%local_buffer)) deallocate(this%local_buffer)
    if (allocated(this%local_buffer_t)) deallocate(this%local_buffer_t)
    if (allocated(this%global_buffer)) deallocate(this%global_buffer)
    if (allocated(this%recvcounts)) deallocate(this%recvcounts)
    if (allocated(this%displs)) deallocate(this%displs)
    if (allocated(this%recvcounts_c)) deallocate(this%recvcounts_c)
    if (allocated(this%displs_c)) deallocate(this%displs_c)

    nullify(this%coef)

    call this%free_base()

  end subroutine boundary_data_writer_free


  !> Gather a field at the masked points into a column of the local buffer.
  !! @param f The source field.
  !! @param col The destination column.
  subroutine boundary_data_writer_gather_column(this, f, col)
    class(boundary_data_writer_t), intent(inout) :: this
    type(field_t), intent(in) :: f
    integer, intent(in) :: col

    if (this%n_local .le. 0) return

    call this%bdata%get(f, this%work)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call this%work%copy_from(DEVICE_TO_HOST, .true.)
    end if

    call copy(this%local_buffer(:, col), this%work%x, this%n_local)

  end subroutine boundary_data_writer_gather_column

  !> Sample the boundary and write.
  !! @param time The current time state.
  subroutine boundary_data_writer_compute(this, time)
    class(boundary_data_writer_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i, col, n_rows, out_int
    character(len=80) :: group_name
    real(kind=rp) :: time_
    type(vector_t) :: vec_time

    if (time%t .lt. this%start_time) then
       call this%output_controller%set_counter(time)
       return
    end if
    if (.not. this%output_controller%check(time)) return

    ! Re-gather the boundary geometry only when the mesh has moved.
    if (this%mesh_has_changed) then
       call this%bdata%update_geometry()

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call this%bdata%x%copy_from(DEVICE_TO_HOST, .false.)
          call this%bdata%y%copy_from(DEVICE_TO_HOST, .false.)
          call this%bdata%z%copy_from(DEVICE_TO_HOST, .false.)
          call this%bdata%n_x%copy_from(DEVICE_TO_HOST, .false.)
          call this%bdata%n_y%copy_from(DEVICE_TO_HOST, .false.)
          call this%bdata%n_z%copy_from(DEVICE_TO_HOST, .false.)
          call this%bdata%area%copy_from(DEVICE_TO_HOST, .true.)
       end if
    end if

    ! Fill the local buffer
    col = 0
    if (this%geometry_in_data) then
       call copy(this%local_buffer(:, col + 1), this%bdata%x%x, this%n_local)
       call copy(this%local_buffer(:, col + 2), this%bdata%y%x, this%n_local)
       call copy(this%local_buffer(:, col + 3), this%bdata%z%x, this%n_local)
       col = col + 3

       if (this%output_normals) then
          call copy(this%local_buffer(:, col + 1), this%bdata%n_x%x, &
               this%n_local)
          call copy(this%local_buffer(:, col + 2), this%bdata%n_y%x, &
               this%n_local)
          call copy(this%local_buffer(:, col + 3), this%bdata%n_z%x, &
               this%n_local)
          col = col + 3
       end if

       if (this%output_area) then
          call copy(this%local_buffer(:, col + 1), this%bdata%area%x, &
               this%n_local)
          col = col + 1
       end if
    end if

    do i = 1, this%fields%size()
       col = col + 1
       call this%gather_column(this%fields%items(i)%ptr, col)
    end do

    ! Write
    select type (ft => this%fout%file_type)
    type is (csv_file_t)
       call boundary_data_writer_gather_to_root(this)

       n_rows = this%n_global

       if (pe_rank .eq. 0) then
          do i = 1, n_rows
             this%mat_out%x(i, 1:this%n_cols) = &
                  this%global_buffer(1:this%n_cols, i)
          end do
          call this%fout%write(this%mat_out, time%t)
       end if

    class is (hdf5_file_t)
       out_int = this%output_controller%nexecutions + 1

       call ft%open("w")
       call ft%set_active_group("boundary_data")
       call ft%write_attribute("NSteps", out_int)

       if (.not. this%append_out) then
          write(group_name, '(A,I0)') "boundary_data/Step_", out_int
          call ft%set_active_group(trim(group_name))
       end if

       do i = 1, this%n_cols
          call copy(this%vec_out%x, this%local_buffer(:, i), &
               this%vec_out%size())
          this%vec_out%name = trim(this%col_names(i))
          call ft%write_dataset(this%vec_out)
       end do

       if (this%append_out) then
          ! Time is a growing dataset alongside the samples.
          if (pe_rank .eq. 0) then
             call vec_time%init(1, "time")
             vec_time%x(1) = time%t
          else
             call vec_time%init(0, "time")
          end if
          call ft%write_dataset(vec_time)
          call vec_time%free()
       else
          ! Time is an attribute on the Step_N subgroup.
          time_ = time%t
          call ft%write_attribute("time", time_)
       end if

       call ft%close()
    end select

    call this%output_controller%register_execution()

  end subroutine boundary_data_writer_compute

  !> Gather the local sample buffer onto rank zero.
  subroutine boundary_data_writer_gather_to_root(this)
    class(boundary_data_writer_t), intent(inout) :: this
    integer :: ierr

    call trsp(this%local_buffer_t, this%n_cols, this%local_buffer, this%n_local)

    call MPI_Gatherv(this%local_buffer_t, this%n_local * this%n_cols, &
         MPI_REAL_PRECISION, this%global_buffer, this%recvcounts_c, &
         this%displs_c, MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)

  end subroutine boundary_data_writer_gather_to_root

  !> Collect the reference geometry of the local points into a matrix.
  !! @param mat The matrix to fill.
  !! @param n_geom The number of geometry entries per point.
  subroutine boundary_data_writer_local_geometry(this, mat, n_geom)
    class(boundary_data_writer_t), intent(inout) :: this
    type(matrix_t), intent(inout) :: mat
    integer, intent(in) :: n_geom
    integer :: i, k

    call mat%init(n_geom, max(0, this%n_local), "coordinates")

    do i = 1, this%n_local
       mat%x(1, i) = this%bdata%x%x(i)
       mat%x(2, i) = this%bdata%y%x(i)
       mat%x(3, i) = this%bdata%z%x(i)
       k = 3
       if (this%output_normals) then
          mat%x(4, i) = this%bdata%n_x%x(i)
          mat%x(5, i) = this%bdata%n_y%x(i)
          mat%x(6, i) = this%bdata%n_z%x(i)
          k = 6
       end if
       if (this%output_area) mat%x(k + 1, i) = this%bdata%area%x(i)
    end do

  end subroutine boundary_data_writer_local_geometry

  !> Gather the reference geometry of all points onto rank zero.
  !! @param mat The matrix to fill.
  !! @param n_geom The number of geometry entries per point.
  subroutine boundary_data_writer_gather_geometry(this, mat, n_geom)
    class(boundary_data_writer_t), intent(inout) :: this
    type(matrix_t), intent(inout) :: mat
    integer, intent(in) :: n_geom
    real(kind=rp), allocatable :: sendbuf(:,:), recvbuf(:,:)
    integer, allocatable :: counts(:), disp(:)
    integer :: i, k, ierr

    allocate(sendbuf(n_geom, this%n_local))
    sendbuf = 0.0_rp

    do i = 1, this%n_local
       sendbuf(1, i) = this%bdata%x%x(i)
       sendbuf(2, i) = this%bdata%y%x(i)
       sendbuf(3, i) = this%bdata%z%x(i)
       k = 3
       if (this%output_normals) then
          sendbuf(4, i) = this%bdata%n_x%x(i)
          sendbuf(5, i) = this%bdata%n_y%x(i)
          sendbuf(6, i) = this%bdata%n_z%x(i)
          k = 6
       end if
       if (this%output_area) sendbuf(k + 1, i) = this%bdata%area%x(i)
    end do

    if (pe_rank .eq. 0) then
       allocate(recvbuf(n_geom, this%n_global))
    else
       allocate(recvbuf(n_geom, 1))
    end if
    recvbuf = 0.0_rp

    allocate(counts(pe_size), disp(pe_size))
    counts = this%recvcounts * n_geom
    disp = this%displs * n_geom

    call MPI_Gatherv(sendbuf, this%n_local * n_geom, MPI_REAL_PRECISION, &
         recvbuf, counts, disp, MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)

    if (pe_rank .eq. 0) then
       call mat%init(this%n_global, n_geom)
       do i = 1, this%n_global
          mat%x(i, 1:n_geom) = recvbuf(1:n_geom, i)
       end do
    else
       call mat%init(1, 1)
    end if

    deallocate(sendbuf, recvbuf, counts, disp)

  end subroutine boundary_data_writer_gather_geometry
end module boundary_data_writer_simcomp
