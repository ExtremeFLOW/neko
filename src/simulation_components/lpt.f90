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
!> Implements `lpt_t`. (Lagrangian Particle Tracking)
module lagrangian_particle_tracking
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use registry, only : neko_registry
  use field, only : field_t
  use case, only : case_t
  use json_utils, only : json_get, json_get_or_default, &
       json_get_subdict_or_empty
  use time_state, only : time_state_t
  use global_interpolation, only : global_interpolation_t
  use point, only : point_t
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error
  use file, only : file_t
  use matrix, only : matrix_t
  use tensor, only : trsp
  use comm, only : pe_rank, pe_size, NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_Exscan, MPI_Gather, MPI_Gatherv, &
       MPI_INTEGER, MPI_SUM
  use csv_file, only : csv_file_t
  implicit none
  private

  !> A simulation component for passive Lagrangian particle tracking.
  !! Particles are redistributed to the rank that owns their current location.
  type, public, extends(simulation_component_t) :: lpt_t
     !> X velocity field.
     type(field_t), pointer :: u => null()
     !> Y velocity field.
     type(field_t), pointer :: v => null()
     !> Z velocity field.
     type(field_t), pointer :: w => null()
     !> Point interpolation helper used to evaluate the carrier velocity.
     type(global_interpolation_t) :: global_interp
     !> Particle coordinates local to this MPI rank.
     real(kind=rp), allocatable :: xyz_particles(:,:)
     !> Optional CSV output file.
     type(file_t) :: output_file
     !> Whether trajectory output is enabled.
     logical :: output_enabled = .false.
     !> Whether to emit informational logs.
     logical :: log = .true.
     !> Time after which tracking should start.
     real(kind=rp) :: start_time = -huge(0.0_rp)
     !> Number of particles currently local to this rank.
     integer :: n_particles = 0
     !> Total number of particles across all ranks.
     integer :: n_global_particles = 0
   contains
     !> Construct the component from a case-file JSON object.
     procedure, pass(this) :: init => lpt_init_from_json
     !> Free the component.
     procedure, pass(this) :: free => lpt_free
     !> Advance particles and write output.
     procedure, pass(this) :: compute_ => lpt_compute
     !> Build initial particle coordinates from the case file.
     procedure, private, pass(this) :: read_particles_json
     !> Read particle coordinates from a CSV file.
     procedure, private, pass(this) :: read_particles_csv
     !> Redistribute particles to the owning MPI rank.
     procedure, private, pass(this) :: redistribute_particles
     !> Interpolate the carrier velocity at the local particles.
     procedure, private, pass(this) :: evaluate_velocity
     !> Write a trajectory snapshot.
     procedure, private, pass(this) :: write_output
     !> Emit a short setup summary.
     procedure, private, pass(this) :: log_status
  end type lpt_t

contains

  !> Construct from JSON.
  !! Supported particle input is either a flat `coordinates` array
  !! `[x1,y1,z1,x2,y2,z2,...]` or a three-column CSV `points_file`.
  subroutine lpt_init_from_json(this, json, case)
    class(lpt_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    type(json_file) :: interp_subdict
    character(len=:), allocatable :: name
    character(len=:), allocatable :: output_filename

    call this%free()

    call json_get_or_default(json, "name", name, "lpt")
    call json_get_or_default(json, "log", this%log, .true.)
    call json_get_or_default(json, "start_time", this%start_time, &
         -huge(0.0_rp))
    call this%init_base(json, case)

    this%name = name
    this%u => neko_registry%get_field_by_name("u")
    this%v => neko_registry%get_field_by_name("v")
    this%w => neko_registry%get_field_by_name("w")

    call this%read_particles_json(json)

    call json_get_subdict_or_empty(json, "interpolation", interp_subdict)
    call this%global_interp%init(case%fluid%dm_Xh, &
         params_subdict = interp_subdict)
    call this%redistribute_particles()

    call json_get_or_default(json, "output_filename", output_filename, &
         trim(this%name) // ".csv")
    call this%output_file%init(this%case%output_directory // &
         trim(output_filename), &
         header = "tstep,time,particle_id,x,y,z,u,v,w", overwrite = .true.)
    this%output_enabled = .true.

    call this%log_status()
  end subroutine lpt_init_from_json

  !> Read particle coordinates from JSON.
  subroutine read_particles_json(this, json)
    class(lpt_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    real(kind=rp), allocatable :: coords(:)

    if (pe_rank .eq. 0) then
       if (json%valid_path("coordinates")) then
          call json_get(json, "coordinates", coords)
          if (mod(size(coords), 3) .ne. 0) then
             call neko_error("lpt coordinates must contain 3 values per " // &
                  "particle")
          end if
          this%n_particles = size(coords) / 3
          allocate(this%xyz_particles(3, this%n_particles))
          this%xyz_particles = reshape(coords, [3, this%n_particles])
       else if (json%valid_path("points_file")) then
          call this%read_particles_csv(json)
       else
          call neko_error("lpt requires either coordinates or points_file")
       end if
    else
       this%n_particles = 0
       allocate(this%xyz_particles(3, 0))
    end if

    this%n_global_particles = this%n_particles
  end subroutine read_particles_json

  !> Read particle coordinates from a three-column CSV file.
  subroutine read_particles_csv(this, json)
    class(lpt_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: points_file
    type(file_t) :: file_in
    type(matrix_t) :: mat_in

    if (pe_rank .ne. 0) return

    call json_get(json, "points_file", points_file)
    call file_in%init(trim(points_file))

    select type (ft => file_in%file_type)
    type is (csv_file_t)
       this%n_particles = ft%count_lines()
       call mat_in%init(this%n_particles, 3)
       call ft%read(mat_in)
       allocate(this%xyz_particles(3, this%n_particles))
       this%xyz_particles = transpose(mat_in%x)
    class default
       call neko_error("lpt points_file must be a csv file")
    end select

    call mat_in%free()
    call file_in%free()
  end subroutine read_particles_csv

  !> Redistribute particles to the rank that owns their current location.
  subroutine redistribute_particles(this)
    class(lpt_t), intent(inout) :: this
    integer :: ierr

    call this%global_interp%find_points_and_redist(this%xyz_particles, &
         this%n_particles)
    call MPI_Allreduce(this%n_particles, this%n_global_particles, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
  end subroutine redistribute_particles

  !> Interpolate the carrier velocity at the local particles.
  subroutine evaluate_velocity(this, vel)
    class(lpt_t), intent(inout) :: this
    real(kind=rp), allocatable, intent(out) :: vel(:,:)
    logical :: do_interp_on_host

    allocate(vel(3, this%n_particles))
    if (this%n_particles .eq. 0) return

    do_interp_on_host = .false.
    call this%global_interp%evaluate(vel(1,:), this%u%x, do_interp_on_host)
    call this%global_interp%evaluate(vel(2,:), this%v%x, do_interp_on_host)
    call this%global_interp%evaluate(vel(3,:), this%w%x, do_interp_on_host)
  end subroutine evaluate_velocity

  !> Advance particles with an explicit Euler update and optionally write them.
  subroutine lpt_compute(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), allocatable :: vel(:,:)
    integer :: i
    type(point_t) :: x_p

    if (time%t .lt. this%start_time) return

    call this%redistribute_particles()
    call this%evaluate_velocity(vel)

    if (this%output_enabled) then
      if (this%output_controller%check(time)) then
         call this%write_output(time, vel)
         call this%output_controller%register_execution()
      end if
    end if

    do i = 1, this%n_particles
       this%xyz_particles(:,i) = this%xyz_particles(:,i) + time%dt * vel(:,i)
    end do

    if (allocated(vel)) deallocate(vel)
  end subroutine lpt_compute

  !> Write one trajectory snapshot to CSV by gathering local particle data to
  !! rank 0.
  subroutine write_output(this, time, vel)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: vel(:,:)
    type(matrix_t) :: block
    real(kind=rp), allocatable :: local_data(:,:)
    real(kind=rp), allocatable :: global_data(:,:)
    integer, allocatable :: n_local_particles_per_rank(:)
    integer, allocatable :: particle_offsets(:)
    integer, allocatable :: recvcounts(:)
    integer, allocatable :: displs(:)
    integer :: n_local
    integer :: particle_offset
    integer :: total_particles
    integer :: i
    integer :: ierr

    n_local = this%n_particles
    particle_offset = 0
    call MPI_Exscan(n_local, particle_offset, 1, MPI_INTEGER, MPI_SUM, &
         NEKO_COMM, ierr)
    if (pe_rank .eq. 0) particle_offset = 0

    allocate(local_data(9, n_local))
    do i = 1, n_local
       local_data(1,i) = real(time%tstep, rp)
       local_data(2,i) = time%t
       local_data(3,i) = real(particle_offset + i, rp)
       local_data(4,i) = this%xyz_particles(1,i)
       local_data(5,i) = this%xyz_particles(2,i)
       local_data(6,i) = this%xyz_particles(3,i)
       local_data(7,i) = vel(1,i)
       local_data(8,i) = vel(2,i)
       local_data(9,i) = vel(3,i)
    end do

    if (pe_rank .eq. 0) then
       allocate(n_local_particles_per_rank(pe_size))
    else
       allocate(n_local_particles_per_rank(0))
    end if
    call MPI_Gather(n_local, 1, MPI_INTEGER, n_local_particles_per_rank, 1, &
         MPI_INTEGER, 0, NEKO_COMM, ierr)

    if (pe_rank .eq. 0) then
       allocate(particle_offsets(pe_size))
       allocate(recvcounts(pe_size))
       allocate(displs(pe_size))
       particle_offsets = 0
       recvcounts = 0
       displs = 0

       total_particles = 0
       do i = 1, pe_size
          particle_offsets(i) = total_particles
          total_particles = total_particles + n_local_particles_per_rank(i)
          recvcounts(i) = 9 * n_local_particles_per_rank(i)
          displs(i) = 9 * particle_offsets(i)
       end do

       allocate(global_data(9, total_particles))
    else
       allocate(particle_offsets(0))
       allocate(recvcounts(0))
       allocate(displs(0))
       allocate(global_data(9, 0))
    end if

    call MPI_Gatherv(local_data, 9 * n_local, MPI_REAL_PRECISION, global_data, &
         recvcounts, displs, MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)

    if (pe_rank .eq. 0) then
       call block%init(total_particles, 9)
       call trsp(block%x, total_particles, global_data, 9)
       call this%output_file%write(block)
       call block%free()
    end if

    deallocate(global_data)
    deallocate(n_local_particles_per_rank)
    deallocate(particle_offsets)
    deallocate(recvcounts)
    deallocate(displs)
    deallocate(local_data)
  end subroutine write_output

  !> Free the component.
  subroutine lpt_free(this)
    class(lpt_t), intent(inout) :: this

    if (allocated(this%xyz_particles)) deallocate(this%xyz_particles)
    call this%global_interp%free()
    call this%output_file%free()

    this%u => null()
    this%v => null()
    this%w => null()
    this%output_enabled = .false.
    this%log = .true.
    this%start_time = -huge(0.0_rp)
    this%n_particles = 0
    this%n_global_particles = 0
    call this%free_base()
  end subroutine lpt_free

  !> Emit a setup summary.
  subroutine log_status(this)
    class(lpt_t), intent(in) :: this
    character(len=LOG_SIZE) :: log_buf

    if (.not. this%log) return

    call neko_log%section("Lagrangian particle tracking")
    write(log_buf, '(A,A)') "Name: ", trim(this%name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,I0)') "Global seeded particles: ", &
         this%n_global_particles
    call neko_log%message(log_buf)
    write(log_buf, '(A,I0)') "Local particles on rank 0 at init: ", &
         this%n_particles
    if (pe_rank .eq. 0) call neko_log%message(log_buf)
    call neko_log%message("Input supported from doc/pages/user-guide/" // &
         "case-file.md semantics for simcomp JSON; particle seeding here " // &
         "uses coordinates or points_file.")
    call neko_log%end_section()
  end subroutine log_status

end module lagrangian_particle_tracking
