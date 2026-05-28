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
  use comm, only : pe_rank, NEKO_COMM
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_SUM
  use csv_file, only : csv_file_t
  implicit none
  private

  !> A simulation component for passive Lagrangian particle tracking.
  !! Particles are seeded once, then advected with an explicit Euler update
  !! using a velocity field interpolated from registered flow fields.
  type, public, extends(simulation_component_t) :: lpt_t
     !> X velocity field.
     type(field_t), pointer :: u => null()
     !> Y velocity field.
     type(field_t), pointer :: v => null()
     !> Z velocity field.
     type(field_t), pointer :: w => null()
     !> Point interpolation helper used to evaluate the carrier velocity.
     type(global_interpolation_t) :: global_interp
     !> Names of the velocity fields to sample.
     character(len=80) :: field_names(3)
     !> Particle coordinates stored on rank 0.
     type(point_t), allocatable :: particles(:)
     !> Whether a particle is still active.
     logical, allocatable :: active(:)
     !> Optional CSV output file.
     type(file_t) :: output_file
     !> Whether trajectory output is enabled.
     logical :: output_enabled = .false.
     !> Whether to emit informational logs.
     logical :: log = .true.
     !> Time after which tracking should start.
     real(kind=rp) :: start_time = -huge(0.0_rp)
     !> Number of seeded particles on this rank.
     integer :: n_particles = 0
     !> Number of particles lost on this rank.
     integer :: n_lost = 0
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
     !> Convert particle coordinates to an xyz array.
     procedure, private, pass(this) :: particles_to_xyz
     !> Resolve active particles and interpolate the carrier velocity.
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
    character(len=80), allocatable :: field_names_json(:)

    call this%free()

    call json_get_or_default(json, "name", name, "lpt")
    call json_get_or_default(json, "log", this%log, .true.)
    call json_get_or_default(json, "start_time", this%start_time, &
         -huge(0.0_rp))
    call this%init_base(json, case)

    this%name = name
    this%field_names = [character(len=80) :: "u", "v", "w"]
    if (json%valid_path("field_names")) then
       call json_get(json, "field_names", field_names_json)
       if (size(field_names_json) .ne. 3) then
          call neko_error("lpt requires exactly three field_names")
       end if
       this%field_names = field_names_json
    end if

    this%u => neko_registry%get_field_by_name(trim(this%field_names(1)))
    this%v => neko_registry%get_field_by_name(trim(this%field_names(2)))
    this%w => neko_registry%get_field_by_name(trim(this%field_names(3)))

    call this%read_particles_json(json)

    call json_get_subdict_or_empty(json, "interpolation", interp_subdict)
    call this%global_interp%init(case%fluid%dm_Xh, &
         params_subdict = interp_subdict)

    block
      real(kind=rp), allocatable :: xyz(:,:)
      if (this%n_particles .gt. 0) then
         call this%particles_to_xyz(xyz)
         call this%global_interp%find_points(xyz, this%n_particles)
         deallocate(xyz)
      end if
    end block

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
    integer :: i

    if (pe_rank .eq. 0) then
       if (json%valid_path("coordinates")) then
          call json_get(json, "coordinates", coords)
          if (mod(size(coords), 3) .ne. 0) then
             call neko_error("lpt coordinates must contain 3 values per " // &
                  "particle")
          end if
          this%n_particles = size(coords) / 3
          allocate(this%particles(this%n_particles))
          do i = 1, this%n_particles
             this%particles(i)%x = real(coords(3*i-2:3*i), &
                  kind(this%particles(i)%x(1)))
          end do
       else if (json%valid_path("points_file")) then
          call this%read_particles_csv(json)
       else
          call neko_error("lpt requires either coordinates or points_file")
       end if

       allocate(this%active(this%n_particles))
       this%active = .true.
    else
       this%n_particles = 0
       allocate(this%particles(0))
       allocate(this%active(0))
    end if
  end subroutine read_particles_json

  !> Read particle coordinates from a three-column CSV file.
  subroutine read_particles_csv(this, json)
    class(lpt_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: points_file
    type(file_t) :: file_in
    type(matrix_t) :: mat_in
    integer :: i

    if (pe_rank .ne. 0) return

    call json_get(json, "points_file", points_file)
    call file_in%init(trim(points_file))

    select type (ft => file_in%file_type)
    type is (csv_file_t)
       this%n_particles = ft%count_lines()
       call mat_in%init(this%n_particles, 3)
       call ft%read(mat_in)
       allocate(this%particles(this%n_particles))
       do i = 1, this%n_particles
          this%particles(i)%x = real(mat_in%x(i,:), &
               kind(this%particles(i)%x(1)))
       end do
    class default
       call neko_error("lpt points_file must be a csv file")
    end select

    call mat_in%free()
    call file_in%free()
  end subroutine read_particles_csv

  !> Convert stored particles to an xyz array compatible with interpolation.
  subroutine particles_to_xyz(this, xyz)
    class(lpt_t), intent(in) :: this
    real(kind=rp), allocatable, intent(out) :: xyz(:,:)
    integer :: i

    allocate(xyz(3, this%n_particles))
    do i = 1, this%n_particles
      xyz(:,i) = real(this%particles(i)%x, rp)
    end do
  end subroutine particles_to_xyz

  !> Resolve active particles, interpolate their velocity and optionally
  !! deactivate particles that have left the mesh.
  subroutine evaluate_velocity(this, particle_ids, xyz, vel, n_active)
    class(lpt_t), intent(inout) :: this
    integer, allocatable, intent(out) :: particle_ids(:)
    real(kind=rp), allocatable, intent(out) :: xyz(:,:)
    real(kind=rp), allocatable, intent(out) :: vel(:,:)
    integer, intent(out) :: n_active
    real(kind=rp), allocatable :: u_vals(:), v_vals(:), w_vals(:)
    real(kind=rp), allocatable :: trial_xyz(:,:)
    integer, allocatable :: trial_ids(:)
    logical :: do_interp_on_host
    integer :: i, j, n_trial

    n_active = 0
    allocate(particle_ids(0))
    allocate(xyz(3, 0))
    allocate(vel(3, 0))

    if (pe_rank .ne. 0) return
    if (this%n_particles .eq. 0) return

    n_trial = count(this%active)
    if (n_trial .eq. 0) return

    allocate(trial_ids(n_trial))
    allocate(trial_xyz(3, n_trial))
    j = 0
    do i = 1, this%n_particles
       if (this%active(i)) then
          j = j + 1
          trial_ids(j) = i
          trial_xyz(:,j) = real(this%particles(i)%x, rp)
       end if
    end do

    call this%global_interp%find_points(trial_xyz, n_trial)

    n_active = 0
    do i = 1, n_trial
       if (this%global_interp%pe_owner(i) .ge. 0 .and. &
            this%global_interp%el_owner0(i) .ge. 0) then
          n_active = n_active + 1
       else
          this%active(trial_ids(i)) = .false.
          this%n_lost = this%n_lost + 1
       end if
    end do

    if (n_active .eq. 0) then
       deallocate(trial_ids)
       deallocate(trial_xyz)
       return
    end if

    deallocate(particle_ids)
    deallocate(xyz)
    deallocate(vel)
    allocate(particle_ids(n_active))
    allocate(xyz(3, n_active))
    allocate(vel(3, n_active))

    j = 0
    do i = 1, n_trial
       if (this%active(trial_ids(i))) then
          j = j + 1
          particle_ids(j) = trial_ids(i)
          xyz(:,j) = trial_xyz(:,i)
       end if
    end do

    call this%global_interp%find_points(xyz, n_active)
    allocate(u_vals(n_active), v_vals(n_active), w_vals(n_active))
    do_interp_on_host = .false.
    call this%global_interp%evaluate(u_vals, this%u%x, do_interp_on_host)
    call this%global_interp%evaluate(v_vals, this%v%x, do_interp_on_host)
    call this%global_interp%evaluate(w_vals, this%w%x, do_interp_on_host)

    vel(1,:) = u_vals
    vel(2,:) = v_vals
    vel(3,:) = w_vals

    deallocate(u_vals)
    deallocate(v_vals)
    deallocate(w_vals)
    deallocate(trial_ids)
    deallocate(trial_xyz)
  end subroutine evaluate_velocity

  !> Advance particles with an explicit Euler update and optionally write them.
  subroutine lpt_compute(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer, allocatable :: particle_ids(:)
    real(kind=rp), allocatable :: xyz(:,:), vel(:,:)
    integer :: n_active, i
    type(point_t) :: x_p, u_p

    if (time%t .lt. this%start_time) return
    write(*,*) "LPT active particles: ", count(this%active), " lost particles: ", &
         this%n_lost
    call this%evaluate_velocity(particle_ids, xyz, vel, n_active)

    if (this%output_enabled) then
       if (this%output_controller%check(time)) then
          call this%write_output(time, particle_ids, xyz, vel, n_active)
          call this%output_controller%register_execution()
       end if
    end if

    if (pe_rank .eq. 0) then
       do i = 1, n_active
          x_p%x = real(xyz(:,i), kind(x_p%x(1)))
          u_p%x = real(vel(:,i), kind(u_p%x(1)))
          x_p%x = x_p%x + time%dt * u_p%x
          this%particles(particle_ids(i))%x = x_p%x
       end do
    end if

    if (allocated(particle_ids)) deallocate(particle_ids)
    if (allocated(xyz)) deallocate(xyz)
    if (allocated(vel)) deallocate(vel)
  end subroutine lpt_compute

  !> Write one trajectory snapshot to CSV.
  subroutine write_output(this, time, particle_ids, xyz, vel, n_active)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: particle_ids(:)
    real(kind=rp), intent(in) :: xyz(:,:)
    real(kind=rp), intent(in) :: vel(:,:)
    integer, intent(in) :: n_active
    type(matrix_t) :: block
    integer :: i

    if (pe_rank .ne. 0) return
    if (n_active .eq. 0) return

    call block%init(n_active, 9)
    do i = 1, n_active
       block%x(i,1) = time%tstep
       block%x(i,2) = time%t
       block%x(i,3) = particle_ids(i)
       block%x(i,4) = xyz(1,i)
       block%x(i,5) = xyz(2,i)
       block%x(i,6) = xyz(3,i)
       block%x(i,7) = vel(1,i)
       block%x(i,8) = vel(2,i)
       block%x(i,9) = vel(3,i)
    end do
    call this%output_file%write(block)
    call block%free()
  end subroutine write_output

  !> Free the component.
  subroutine lpt_free(this)
    class(lpt_t), intent(inout) :: this

    if (allocated(this%particles)) deallocate(this%particles)
    if (allocated(this%active)) deallocate(this%active)
    call this%global_interp%free()
    call this%output_file%free()

    this%u => null()
    this%v => null()
    this%w => null()
    this%output_enabled = .false.
    this%log = .true.
    this%start_time = -huge(0.0_rp)
    this%n_particles = 0
    this%n_lost = 0
    call this%free_base()
  end subroutine lpt_free

  !> Emit a setup summary.
  subroutine log_status(this)
    class(lpt_t), intent(in) :: this
    character(len=LOG_SIZE) :: log_buf
    integer :: n_active_global, n_active_local, ierr

    if (.not. this%log) return

    n_active_local = count(this%active)
    call MPI_Allreduce(n_active_local, n_active_global, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)

    call neko_log%section("Lagrangian particle tracking")
    write(log_buf, '(A,A)') "Name: ", trim(this%name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,A,", ",A,", ",A)') "Fields: ", &
         trim(this%field_names(1)), trim(this%field_names(2)), &
         trim(this%field_names(3))
    call neko_log%message(log_buf)
    write(log_buf, '(A,I0)') "Seeded particles: ", n_active_global
    call neko_log%message(log_buf)
    call neko_log%message("Input supported from doc/pages/user-guide/" // &
         "case-file.md semantics for simcomp JSON; particle seeding here " // &
         "uses coordinates or points_file.")
    call neko_log%end_section()
  end subroutine log_status

end module lagrangian_particle_tracking
