! Copyright (c) 2020-2023, The Neko Authors
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
!> Implements probes.
!! @note This modules uses functions from `gslib`, namely `findpts_setup`,
!! `findpts`, and `findpts_eval`. A full description of these subroutines can
!! be found at https://github.com/Nek5000/gslib/blob/master/src/findpts.c
module probes
  use num_types, only: rp
  use matrix, only: matrix_t
  use logger, only: neko_log, LOG_SIZE, NEKO_LOG_DEBUG
  use utils, only: neko_error, nonlinear_index
  use field_list, only: field_list_t
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use dofmap, only: dofmap_t
  use json_module, only : json_file, json_value, json_core
  use json_utils, only : json_get, json_extract_item, json_get_or_default
  use global_interpolation, only: global_interpolation_t
  use tensor, only: trsp
  use point_zone, only: point_zone_t
  use point_zone_registry, only: neko_point_zone_registry
  use comm
  use device
  use file, only : file_t, file_free
  use csv_file, only : csv_file_t
  use case, only : case_t
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public, extends(simulation_component_t) :: probes_t
     !> Number of output fields
     integer :: n_fields = 0
     type(global_interpolation_t) :: global_interp
     !> Global number of probes (needed for i/o)
     integer :: n_global_probes
     !> global offset for writing
     integer :: n_probes_offset
     !> x,y,z coordinates
     real(kind=rp), allocatable :: xyz(:,:)
     !> Interpolated values
     real(kind=rp), allocatable :: out_values(:,:)
     type(c_ptr), allocatable :: out_values_d(:)
     real(kind=rp), allocatable :: out_vals_trsp(:,:)
     !> Number of local elements per rank
     integer :: n_local_probes
     !> Fields to be probed
     type(field_list_t) :: sampled_fields
     character(len=20), allocatable :: which_fields(:)
     !> Allocated on rank 0
     integer, allocatable :: n_local_probes_tot_offset(:)
     integer, allocatable :: n_local_probes_tot(:)
     !>  For output on rank 0
     logical :: seq_io
     real(kind=rp), allocatable :: global_output_values(:,:)
     !> Output variables
     type(file_t) :: fout
     type(matrix_t) :: mat_out
   contains
     !> Initialize from json
     procedure, pass(this) :: init => probes_init_from_json
     ! Actual constructor
     procedure, pass(this) :: init_from_attributes => &
          probes_init_from_attributes
     !> Destructor
     procedure, pass(this) :: free => probes_free
     !> Setup offset for I/O when using sequential write/read from rank 0
     procedure, pass(this) :: setup_offset => probes_setup_offset
     !> Interpolate each probe from its `r,s,t` coordinates.
     procedure, pass(this) :: compute_ => probes_evaluate_and_write

     ! ----------------------------------------------------------------------- !
     ! Private methods

     !> Reader for file type points
     procedure, private, pass(this) :: read_file
     !> Reader for point type points
     procedure, private, pass(this) :: read_point
     !> Reader for line type points
     procedure, private, pass(this) :: read_line
     !> Reader for circle type points
     procedure, private, pass(this) :: read_circle
     !> Reader for point zone type points
     procedure, private, pass(this) :: read_point_zone

     !> Append a new list of points to the exsiting list.
     procedure, private, pass(this) :: add_points
  end type probes_t

contains

  !> Constructor from json.
  subroutine probes_init_from_json(this, json, case)
    class(probes_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: output_file
    character(len=:), allocatable :: input_file
    integer :: i, ierr

    ! JSON variables
    character(len=:), allocatable :: point_type
    type(json_value), pointer :: json_point_list
    type(json_file) :: json_point
    type(json_core) :: core
    integer :: idx, n_point_children

    ! Initialize the base class
    call this%free()
    call this%init_base(json, case)

    !> Read from case file
    call json%info('fields', n_children = this%n_fields)
    call json_get(json, 'fields', this%which_fields)
    call json_get(json, 'output_file', output_file)

    call this%sampled_fields%init(this%n_fields)
    do i = 1, this%n_fields

       call this%sampled_fields%assign(i, &
            & neko_field_registry%get_field(trim(this%which_fields(i))))
    end do

    ! Setup the required arrays and initialize variables.
    this%n_local_probes = 0
    this%n_global_probes = 0

    ! Read the legacy point specification from the points file.
    if (json%valid_path('points_file')) then

       ! Todo: We should add a deprecation notice here
       call json_get(json, 'points_file', input_file)

       ! This is distributed as to make it similar to parallel file
       ! formats later, Reads all into rank 0
       call read_probe_locations(this, this%xyz, this%n_local_probes, &
                                 this%n_global_probes, input_file)
    end if

    ! Go through the points list and construct the probe list
    call json%get('points', json_point_list)
    call json%info('points', n_children = n_point_children)

    do idx = 1, n_point_children
       call json_extract_item(core, json_point_list, idx, json_point)

       call json_get_or_default(json_point, 'type', point_type, 'none')
       select case (point_type)

         case ('file')
          call this%read_file(json_point)
         case ('points')
          call this%read_point(json_point)
         case ('line')
          call this%read_line(json_point)
         case ('plane')
          call neko_error('Plane probes not implemented yet.')
         case ('circle')
          call this%read_circle(json_point)
         case ('point_zone')
          call this%read_point_zone(json_point, case%fluid%dm_Xh)
         case ('none')
          call json_point%print()
          call neko_error('No point type specified.')
         case default
          call neko_error('Unknown region type ' // point_type)
       end select
    end do

    call mpi_allreduce(this%n_local_probes, this%n_global_probes, 1, &
                       MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    call probes_show(this)
    call this%init_from_attributes(case%fluid%dm_Xh, output_file)

  end subroutine probes_init_from_json

  ! ========================================================================== !
  ! Readers for different point types

  !> Read a list of points from a csv file.
  !! @note The points are expected to be in the form of a list of coordinates.
  !! @param[inout] this The probes object.
  !! @param[inout] json The json file object.
  subroutine read_file(this, json)
    class(probes_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: input_file
    real(kind=rp), dimension(:,:), allocatable :: point_list

    integer :: n_local, n_global

    if (pe_rank .ne. 0) return

    call json_get(json, 'file_name', input_file)

    call read_probe_locations(this, point_list, n_local, n_global, input_file)

    call this%add_points(point_list)
  end subroutine read_file

  !> Read a list of points from the json file.
  !! @note The points are expected to be in the form of a list of coordinates.
  !! @param[inout] this The probes object.
  !! @param[inout] json The json file object.
  subroutine read_point(this, json)
    class(probes_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    real(kind=rp), dimension(:,:), allocatable :: point_list
    real(kind=rp), dimension(:), allocatable :: rp_list_reader

    ! Ensure only rank 0 reads the coordinates.
    if (pe_rank .ne. 0) return
    call json_get(json, 'coordinates', rp_list_reader)

    if (mod(size(rp_list_reader), 3) /= 0) then
       call neko_error('Invalid number of coordinates.')
    end if

    ! Allocate list of points and reshape the input array
    allocate(point_list(3, size(rp_list_reader)/3))
    point_list = reshape(rp_list_reader, [3, size(rp_list_reader)/3])

    call this%add_points(point_list)
  end subroutine read_point

  !> Construct a list of points from a line.
  !! @param[inout] this The probes object.
  !! @param[inout] json The json file object.
  subroutine read_line(this, json)
    class(probes_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    real(kind=rp), dimension(:,:), allocatable :: point_list
    real(kind=rp), dimension(:), allocatable :: start, end
    real(kind=rp), dimension(3) :: direction
    real(kind=rp) :: t

    integer :: n_points, i

    ! Ensure only rank 0 reads the coordinates.
    if (pe_rank .ne. 0) return
    call json_get(json, "start", start)
    call json_get(json, "end", end)
    call json_get(json, "amount", n_points)

    ! If either start or end is not of length 3, error out
    if (size(start) /= 3 .or. size(end) /= 3) then
       call neko_error('Invalid start or end coordinates.')
    end if

    ! Calculate the number of points
    allocate(point_list(3, n_points))

    ! Calculate the direction vector
    direction = end - start
    do i = 1, n_points
       t = real(i - 1, kind = rp) / real(n_points - 1, kind = rp)
       point_list(:, i) = start + direction * t
    end do

    call this%add_points(point_list)
  end subroutine read_line

  !> Construct a list of points from a circle.
  !! @details The general structure of the circle is defined by a center point,
  !! a radius and a normal to the plane it lies on. The circle is then
  !! discretized into a number of points, based on the `amount` parameter.
  !! The points are added clockwise starting from a chosen axis, which is
  !! defined by the `axis` parameter.
  !! @note The axis must be one of the following: `x`, `y`, or `z`.
  !! @param[inout] this The probes object.
  !! @param[inout] json The json file object.
  subroutine read_circle(this, json)
    class(probes_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    real(kind=rp), dimension(:,:), allocatable :: point_list
    real(kind=rp), dimension(:), allocatable :: center, normal
    real(kind=rp) :: radius
    real(kind=rp) :: angle
    integer :: n_points, i
    character(len=:), allocatable :: axis

    real(kind=rp), dimension(3) :: zero_line, cross_line, temp
    real(kind=rp) :: pi

    ! Ensure only rank 0 reads the coordinates.
    if (pe_rank .ne. 0) return
    call json_get(json, "center", center)
    call json_get(json, "normal", normal)
    call json_get(json, "radius", radius)
    call json_get(json, "amount", n_points)
    call json_get(json, "axis", axis)

    ! If either center or normal is not of length 3, error out
    if (size(center) /= 3 .or. size(normal) /= 3) then
       call neko_error('Invalid center or normal coordinates.')
    end if
    if (axis /= 'x' .and. axis /= 'y' .and. axis /= 'z') then
       call neko_error('Invalid axis.')
    end if
    if (radius <= 0) then
       call neko_error('Invalid radius.')
    end if
    if (n_points <= 0) then
       call neko_error('Invalid number of points.')
    end if

    ! Normalize the normal vector
    normal = normal / norm2(normal)

    ! Set the zero line
    if (axis .eq. 'x') zero_line = [1.0, 0.0, 0.0]
    if (axis .eq. 'y') zero_line = [0.0, 1.0, 0.0]
    if (axis .eq. 'z') zero_line = [0.0, 0.0, 1.0]

    if (1.0_rp - dot_product(zero_line, normal) .le. 1e-6) then
       call neko_error('Invalid axis and normal.')
    end if

    zero_line = zero_line - dot_product(zero_line, normal) * normal
    zero_line = zero_line / norm2(zero_line)

    cross_line(1) = normal(2) * zero_line(3) - normal(3) * zero_line(2)
    cross_line(2) = normal(3) * zero_line(1) - normal(1) * zero_line(3)
    cross_line(3) = normal(1) * zero_line(2) - normal(2) * zero_line(1)

    ! Calculate the number of points
    allocate(point_list(3, n_points))

    pi = 4.0_rp * atan(1.0_rp)

    ! Calculate the points
    do i = 1, n_points
       angle = 2.0_rp * pi * real(i - 1, kind = rp) / real(n_points, kind = rp)
       temp = cos(angle) * zero_line + sin(angle) * cross_line

       point_list(:, i) = center + radius * temp
    end do

    call this%add_points(point_list)
  end subroutine read_circle

  !> Construct a list of points from a point zone.
  !! @details The GLL points are read from the point zone and added to the
  !! probe list.
  !! @param[inout] this The probes object.
  !! @param[inout] json The json file object.
  subroutine read_point_zone(this, json, dof)
    class(probes_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(dofmap_t), intent(in) :: dof

    real(kind=rp), dimension(:,:), allocatable :: point_list
    character(len=:), allocatable :: point_zone_name
    class(point_zone_t), pointer :: zone
    integer :: i, idx, lx, nlindex(4)
    real(kind=rp) :: x, y, z

    ! Ensure only rank 0 reads the coordinates.
    if (pe_rank .ne. 0) return

    call json_get(json, "name", point_zone_name)
    zone => neko_point_zone_registry%get_point_zone(point_zone_name)

    ! Allocate list of points and reshape the input array
    allocate(point_list(3, zone%size))

    lx = dof%Xh%lx
    do i = 1, zone%size
       idx = zone%mask(i)

       nlindex = nonlinear_index(idx, lx, lx, lx)
       x = dof%x(nlindex(1), nlindex(2), nlindex(3), nlindex(4))
       y = dof%y(nlindex(1), nlindex(2), nlindex(3), nlindex(4))
       z = dof%z(nlindex(1), nlindex(2), nlindex(3), nlindex(4))

       point_list(:, i) = [x, y, z]
    end do

    call this%add_points(point_list)
  end subroutine read_point_zone

  ! ========================================================================== !
  ! Supporting routines

  !> Append a new list of points to the exsiting list.
  !! @param[inout] this The probes object.
  !! @param[in] new_points The new points to be appended.
  subroutine add_points(this, new_points)
    class(probes_t), intent(inout) :: this
    real(kind=rp), dimension(:,:), intent(in) :: new_points

    real(kind=rp), dimension(:,:), allocatable :: temp
    integer :: n_old, n_new

    ! Get the current number of points
    n_old = this%n_local_probes
    n_new = size(new_points, 2)

    ! Move current points to a temporary array
    if (allocated(this%xyz)) then
       call move_alloc(this%xyz, temp)
    end if

    ! Allocate the new array and copy the full list of points
    allocate(this%xyz(3, n_old + n_new))
    if (allocated(temp)) then
       this%xyz(:, 1:n_old) = temp
    end if
    this%xyz(:, n_old+1:n_old+n_new) = new_points

    this%n_local_probes = this%n_local_probes + n_new
  end subroutine add_points

  ! ========================================================================== !
  ! General initialization routine

  !> Initialize without json things
  !! @param dof Dofmap to probe
  !! @output_file Name of output file, current must be CSV
  subroutine probes_init_from_attributes(this, dof, output_file)
    class(probes_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dof
    character(len=:), allocatable, intent(inout) :: output_file
    character(len=1024) :: header_line
    real(kind=rp), allocatable :: global_output_coords(:,:)
    integer :: i, ierr
    type(matrix_t) :: mat_coords

    !> Init interpolator
    call this%global_interp%init(dof)

    !> find probes and redistribute them
    call this%global_interp%find_points_and_redist(this%xyz, &
                                                   this%n_local_probes)

    !> Allocate output array
    allocate(this%out_values(this%n_local_probes, this%n_fields))
    allocate(this%out_values_d(this%n_fields))
    allocate(this%out_vals_trsp(this%n_fields, this%n_local_probes))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, this%n_fields
          this%out_values_d(i) = c_null_ptr
          call device_map(this%out_values(:,i), this%out_values_d(i), &
                          this%n_local_probes)
       end do
    end if

    !> Initialize the output file
    this%fout = file_t(trim(output_file))

    select type (ft => this%fout%file_type)
      type is (csv_file_t)

       this%seq_io = .true.

       ! Build the header
       write(header_line, '(I0,A,I0)') this%n_global_probes, ",", this%n_fields
       do i = 1, this%n_fields
          header_line = trim(header_line) // "," // trim(this%which_fields(i))
       end do
       call this%fout%set_header(header_line)

       !> Necessary for not-parallel csv format...
       !! offsets and n points per pe
       !! Needed at root for sequential csv i/o
       allocate(this%n_local_probes_tot(pe_size))
       allocate(this%n_local_probes_tot_offset(pe_size))
       call this%setup_offset()
       if (pe_rank .eq. 0) then
          allocate(global_output_coords(3, this%n_global_probes))
          call this%mat_out%init(this%n_global_probes, this%n_fields)
          allocate(this%global_output_values(this%n_fields, &
                                             this%n_global_probes))
          call mat_coords%init(this%n_global_probes,3)
       end if
       call MPI_Gatherv(this%xyz, 3*this%n_local_probes, &
                        MPI_DOUBLE_PRECISION, global_output_coords, &
                        3*this%n_local_probes_tot, &
                        3*this%n_local_probes_tot_offset, &
                        MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)
       if (pe_rank .eq. 0) then
          call trsp(mat_coords%x, this%n_global_probes, &
                    global_output_coords, 3)
          !! Write the data to the file
          call this%fout%write(mat_coords)
       end if
      class default
       call neko_error("Invalid data. Expected csv_file_t.")
    end select

  end subroutine probes_init_from_attributes

  !> Destructor
  subroutine probes_free(this)
    class(probes_t), intent(inout) :: this

    if (allocated(this%xyz)) then
       deallocate(this%xyz)
    end if

    if (allocated(this%out_values)) then
       deallocate(this%out_values)
    end if

    if (allocated(this%out_vals_trsp)) then
       deallocate(this%out_vals_trsp)
    end if

    call this%sampled_fields%free()

    if (allocated(this%n_local_probes_tot)) then
       deallocate(this%n_local_probes_tot)
    end if

    if (allocated(this%n_local_probes_tot_offset)) then
       deallocate(this%n_local_probes_tot_offset)
    end if

    if (allocated(this%global_output_values)) then
       deallocate(this%global_output_values)
    end if

    call this%global_interp%free()

  end subroutine probes_free

  !> Print current probe status, with number of probes and coordinates
  subroutine probes_show(this)
    class(probes_t), intent(in) :: this
    character(len=LOG_SIZE) :: log_buf ! For logging status
    integer :: i

    ! Probes summary
    call neko_log%section('Probes')
    write(log_buf, '(A,I6)') "Number of probes: ", this%n_global_probes
    call neko_log%message(log_buf)

    call neko_log%message("xyz-coordinates:", lvl = NEKO_LOG_DEBUG)
    do i = 1, this%n_local_probes
       write(log_buf, '("(",F10.6,",",F10.6,",",F10.6,")")') this%xyz(:,i)
       call neko_log%message(log_buf, lvl = NEKO_LOG_DEBUG)
    end do

    ! Field summary
    write(log_buf, '(A,I6)') "Number of fields: ", this%n_fields
    call neko_log%message(log_buf)
    do i = 1, this%n_fields
       write(log_buf, '(A,I6,A,A)') &
            "Field: ", i, " ", trim(this%which_fields(i))
       call neko_log%message(log_buf, lvl = NEKO_LOG_DEBUG)
    end do
    call neko_log%end_section()
    call neko_log%newline()

  end subroutine probes_show

  !> Show the status of processor/element owner and error code for each point
  subroutine probes_debug(this)
    class(probes_t) :: this

    character(len=LOG_SIZE) :: log_buf ! For logging status
    integer :: i

    do i = 1, this%n_local_probes
       write (log_buf, *) pe_rank, "/", this%global_interp%proc_owner(i), &
            "/" , this%global_interp%el_owner(i), &
            "/", this%global_interp%error_code(i)
       call neko_log%message(log_buf)
       write(log_buf, '(A5,"(",F10.6,",",F10.6,",",F10.6,")")') &
            "rst: ", this%global_interp%rst(:,i)
       call neko_log%message(log_buf)
    end do
  end subroutine probes_debug

  !> Setup offset for rank 0
  subroutine probes_setup_offset(this)
    class(probes_t) :: this
    integer :: ierr
    this%n_local_probes_tot = 0
    this%n_local_probes_tot_offset = 0
    this%n_probes_offset = 0
    call MPI_Gather(this%n_local_probes, 1, MPI_INTEGER, &
                    this%n_local_probes_tot, 1, MPI_INTEGER, &
                    0, NEKO_COMM, ierr)

    call MPI_Exscan(this%n_local_probes, this%n_probes_offset, 1, &
                    MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    call MPI_Gather(this%n_probes_offset, 1, MPI_INTEGER, &
                    this%n_local_probes_tot_offset, 1, MPI_INTEGER, &
                    0, NEKO_COMM, ierr)



  end subroutine probes_setup_offset

  !> Interpolate each probe from its `r,s,t` coordinates.
  !! @note The final interpolated field is only available on rank 0.
  !! @param t Current simulation time.
  !! @param tstep Current time step.
  subroutine probes_evaluate_and_write(this, t, tstep)
    class(probes_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: i, ierr

    !> Check controller to determine if we must write
    do i = 1, this%n_fields
       call this%global_interp%evaluate(this%out_values(:,i), &
                                        this%sampled_fields%items(i)%ptr%x)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, this%n_fields
          call device_memcpy(this%out_values(:,i), this%out_values_d(i), &
                             this%n_local_probes, DEVICE_TO_HOST, sync = .true.)
       end do
    end if

    if (this%output_controller%check(t, tstep)) then
       ! Gather all values to rank 0
       ! If io is only done at root
       if (this%seq_io) then
          call trsp(this%out_vals_trsp, this%n_fields, &
                    this%out_values, this%n_local_probes)
          call MPI_Gatherv(this%out_vals_trsp, &
                           this%n_fields*this%n_local_probes, &
                           MPI_DOUBLE_PRECISION, this%global_output_values, &
                           this%n_fields*this%n_local_probes_tot, &
                           this%n_fields*this%n_local_probes_tot_offset, &
                           MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)
          if (pe_rank .eq. 0) then
             call trsp(this%mat_out%x, this%n_global_probes, &
                       this%global_output_values, this%n_fields)
             call this%fout%write(this%mat_out, t)
          end if
       else
          call neko_error('probes sim comp, parallel io need implementation')
       end if

       !! Register the execution of the activity
       call this%output_controller%register_execution()
    end if

  end subroutine probes_evaluate_and_write

  !> Initialize the physical coordinates from a `csv` input file
  !! @param points_file A csv file containing probes.
  subroutine read_probe_locations(this, xyz, n_local_probes, n_global_probes, &
                                  points_file)
    class(probes_t), intent(inout) :: this
    character(len=:), allocatable :: points_file
    real(kind=rp), allocatable :: xyz(:,:)
    integer, intent(inout) :: n_local_probes, n_global_probes

    !> Supporting variables
    type(file_t) :: file_in

    file_in = file_t(trim(points_file))
    !> Reads on rank 0 and distributes the probes across the different ranks
    select type (ft => file_in%file_type)
      type is (csv_file_t)
       call read_xyz_from_csv(xyz, n_local_probes, n_global_probes, ft)
       this%seq_io = .true.
      class default
       call neko_error("Invalid data. Expected csv_file_t.")
    end select

    !> Close the file
    call file_free(file_in)

  end subroutine read_probe_locations

  !> Read and initialize the number of probes from a `csv` input file
  !! @param xyz xyz coordinates of the probes
  !! @param n_local_probes The number of probes local to this process
  !! @param n_global_probes The number of total probes on all processes
  !! @param f The csv file we read from
  subroutine read_xyz_from_csv(xyz, n_local_probes, n_global_probes, f)
    type(csv_file_t), intent(inout) :: f
    real(kind=rp), allocatable :: xyz(:,:)
    integer, intent(inout) :: n_local_probes, n_global_probes
    type(matrix_t) :: mat_in, mat_in2
    integer :: n_lines

    n_lines = f%count_lines()

    ! Update the number of probes
    n_global_probes = n_lines

    ! Initialize the temporal array
    if (pe_rank .eq. 0) then
       n_local_probes = n_global_probes
       allocate(xyz(3, n_local_probes))
       call mat_in%init(n_global_probes,3)
       call mat_in2%init(3, n_global_probes)
       call f%read(mat_in)
       call trsp(xyz, 3, mat_in%x, n_global_probes)
    else
       n_local_probes = 0
       allocate(xyz(3, n_local_probes))
    end if

  end subroutine read_xyz_from_csv
end module probes
