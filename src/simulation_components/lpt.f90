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
  use glb_intrp_comm, only : glb_intrp_comm_t
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error
  use file, only : file_t
  use matrix, only : matrix_t
  use math, only : add2s2
  use tensor, only : trsp
  use lpt_periodic_bc, only : lpt_periodic_bc_t
  use comm, only : pe_rank, pe_size, NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_Gather, MPI_Gatherv, &
       MPI_INTEGER, MPI_SUM
  use csv_file, only : csv_file_t
  implicit none
  private

  integer, parameter :: LPT_VEL_HISTORY_LEN = 2

  type, private :: particles_t
     real(kind=rp), allocatable :: xyz(:,:)
     integer, allocatable :: ids(:)
     real(kind=rp), allocatable :: vel_lag(:,:,:)
     integer :: n = 0
     integer :: n_global = 0
   contains
     procedure, pass(this) :: init => particles_init
     procedure, pass(this) :: free => particles_free
  end type particles_t

  !> A simulation component for passive Lagrangian particle tracking.
  !! Particles are redistributed to the rank that owns their current location.
  type, public, extends(simulation_component_t) :: lpt_t
     type(field_t), pointer :: u => null()
     type(field_t), pointer :: v => null()
     type(field_t), pointer :: w => null()
     type(global_interpolation_t) :: global_interp
     type(lpt_periodic_bc_t) :: periodic_bc
     type(particles_t) :: particles
     type(file_t) :: output_file
     logical :: output_enabled = .false.
     logical :: log = .true.
     real(kind=rp) :: start_time = -huge(0.0_rp)
   contains
     procedure, pass(this) :: init => lpt_init_from_json
     procedure, pass(this) :: free => lpt_free
     procedure, pass(this) :: compute_ => lpt_compute
     procedure, private, pass(this) :: read_particles_json
     procedure, private, pass(this) :: read_particles_csv
     procedure, private, pass(this) :: redistribute_particles
     procedure, private, pass(this) :: evaluate_velocity
     procedure, private, pass(this) :: write_output
     procedure, private, pass(this) :: redistribute_particle_ids
     procedure, private, pass(this) :: redistribute_velocity_history
     procedure, private, pass(this) :: log_status
  end type lpt_t

contains

  subroutine particles_init(this, xyz)
    class(particles_t), intent(inout) :: this
    real(kind=rp), intent(in) :: xyz(:,:)
    integer :: i

    if (size(xyz, 1) .ne. 3) then
       call neko_error("particles coordinates must have size 3 in dim 1")
    end if

    call this%free()

    this%n = size(xyz, 2)
    this%n_global = this%n
    allocate(this%xyz(3, this%n))
    allocate(this%ids(this%n))
    allocate(this%vel_lag(3, LPT_VEL_HISTORY_LEN, this%n))
    this%xyz = xyz
    this%vel_lag = 0.0_rp
    do i = 1, this%n
       this%ids(i) = i
    end do
  end subroutine particles_init

  subroutine particles_free(this)
    class(particles_t), intent(inout) :: this

    if (allocated(this%xyz)) deallocate(this%xyz)
    if (allocated(this%ids)) deallocate(this%ids)
    if (allocated(this%vel_lag)) deallocate(this%vel_lag)
    this%n = 0
    this%n_global = 0
  end subroutine particles_free

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
    real(kind=rp), allocatable :: vel(:,:)

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
    call this%periodic_bc%init(case%fluid%msh, case%fluid%dm_Xh, &
         case%fluid%c_Xh)
    call this%redistribute_particles()

    call json_get_or_default(json, "output_filename", output_filename, &
         trim(this%name) // ".csv")
    call this%output_file%init(this%case%output_directory // &
         trim(output_filename), &
         header = "tstep,time,particle_id,x,y,z,u,v,w", overwrite = .true.)
    this%output_enabled = .true.
    if (case%time%t .ge. this%start_time) then
       call this%evaluate_velocity(vel)
       call this%write_output(case%time, vel)
       if (allocated(vel)) deallocate(vel)
    end if

    call this%log_status()
  end subroutine lpt_init_from_json

  !> Read particle coordinates from JSON.
  subroutine read_particles_json(this, json)
    class(lpt_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    real(kind=rp), allocatable :: coords(:)
    real(kind=rp), allocatable :: empty_xyz(:,:)

    if (pe_rank .eq. 0) then
       if (json%valid_path("coordinates")) then
          call json_get(json, "coordinates", coords)
          if (mod(size(coords), 3) .ne. 0) then
             call neko_error("lpt coordinates must contain 3 values per " // &
                  "particle")
          end if
          call this%particles%init(reshape(coords, [3, size(coords) / 3]))
       else if (json%valid_path("points_file")) then
          call this%read_particles_csv(json)
       else
          call neko_error("lpt requires either coordinates or points_file")
       end if
    else
       allocate(empty_xyz(3, 0))
       call this%particles%init(empty_xyz)
    end if
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
       call mat_in%init(ft%count_lines(), 3)
       call ft%read(mat_in)
       call this%particles%init(transpose(mat_in%x))
    class default
       call neko_error("lpt points_file must be a csv file")
    end select

    call mat_in%free()
    call file_in%free()
  end subroutine read_particles_csv

  !> Redistribute particles to the rank that owns their current location.
  subroutine redistribute_particles(this)
    class(lpt_t), intent(inout) :: this
    type(glb_intrp_comm_t) :: redist_comm
    integer :: n_particles_old
    integer, allocatable :: particle_ids_local(:)
    real(kind=rp), allocatable :: vel_particles_lag_local(:,:,:)
    integer :: ierr

    n_particles_old = this%particles%n
    call this%periodic_bc%wrap(this%particles%xyz, this%particles%n)
    call this%global_interp%find_points_and_redist(this%particles%xyz, &
         this%particles%n)
    call this%global_interp%init_redist_comm(redist_comm)
    call this%redistribute_particle_ids(redist_comm, n_particles_old, &
         particle_ids_local)
    call this%redistribute_velocity_history(redist_comm, n_particles_old, &
         vel_particles_lag_local)
    call redist_comm%free()
    call move_alloc(particle_ids_local, this%particles%ids)
    call move_alloc(vel_particles_lag_local, this%particles%vel_lag)
    call MPI_Allreduce(this%particles%n, this%particles%n_global, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
  end subroutine redistribute_particles

  !> Interpolate the carrier velocity at the local particles.
  subroutine evaluate_velocity(this, vel)
    class(lpt_t), intent(inout) :: this
    real(kind=rp), allocatable, intent(out) :: vel(:,:)
    logical :: do_interp_on_host

    allocate(vel(3, this%particles%n))
    if (this%particles%n .eq. 0) return

    do_interp_on_host = .false.
    call this%global_interp%evaluate(vel(1,:), this%u%x, do_interp_on_host)
    call this%global_interp%evaluate(vel(2,:), this%v%x, do_interp_on_host)
    call this%global_interp%evaluate(vel(3,:), this%w%x, do_interp_on_host)
  end subroutine evaluate_velocity

  !> Advance particles with the fluid scheme's explicit multistep coefficients
  !! and optionally write them.
  subroutine lpt_compute(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), allocatable :: vel(:,:)
    real(kind=rp) :: dtc
    integer :: i
    integer :: j
    integer :: nadv

    if (time%t .lt. this%start_time) return

    call this%redistribute_particles()
    call this%evaluate_velocity(vel)
    nadv = this%case%fluid%ext_bdf%nadv

    if (nadv .gt. 1 .and. all(this%particles%vel_lag .eq. 0.0_rp)) then
       do j = 1, LPT_VEL_HISTORY_LEN
          this%particles%vel_lag(:, j, :) = vel
       end do
    end if

    dtc = time%dt * this%case%fluid%ext_bdf%advection_coeffs%x(1)
    do i = 1, 3
       call add2s2(this%particles%xyz(i, :), vel(i, :), dtc, this%particles%n)
    end do

    do j = 2, nadv
       dtc = time%dt * this%case%fluid%ext_bdf%advection_coeffs%x(j)
       do i = 1, 3
          call add2s2(this%particles%xyz(i, :), this%particles%vel_lag(i, j - 1, :), &
               dtc, this%particles%n)
       end do
    end do

    do j = LPT_VEL_HISTORY_LEN, 2, -1
       this%particles%vel_lag(:, j, :) = this%particles%vel_lag(:, j - 1, :)
    end do
    this%particles%vel_lag(:, 1, :) = vel

    if (this%output_enabled) then
      if (this%output_controller%check(time)) then
         call this%write_output(time, vel)
         call this%output_controller%register_execution()
      end if
    end if

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
    integer, allocatable :: recvcounts(:)
    integer, allocatable :: displs(:)
    integer :: n_local
    integer :: total_particles
    integer :: i
    integer :: ierr

    n_local = this%particles%n

    allocate(local_data(9, n_local))
    do i = 1, n_local
       local_data(1,i) = real(time%tstep, rp)
       local_data(2,i) = time%t
       local_data(3,i) = real(this%particles%ids(i), rp)
       local_data(4,i) = this%particles%xyz(1,i)
       local_data(5,i) = this%particles%xyz(2,i)
       local_data(6,i) = this%particles%xyz(3,i)
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
       allocate(recvcounts(pe_size))
       allocate(displs(pe_size))
       recvcounts = 0
       displs = 0

       total_particles = 0
       do i = 1, pe_size
          total_particles = total_particles + n_local_particles_per_rank(i)
          recvcounts(i) = 9 * n_local_particles_per_rank(i)
          displs(i) = 9 * (total_particles - n_local_particles_per_rank(i))
       end do

       allocate(global_data(9, total_particles))
    else
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
    deallocate(recvcounts)
    deallocate(displs)
    deallocate(local_data)
  end subroutine write_output

  !> Free the component.
  subroutine lpt_free(this)
    class(lpt_t), intent(inout) :: this

    call this%particles%free()
    call this%global_interp%free()
    call this%periodic_bc%free()
    call this%output_file%free()

    this%u => null()
    this%v => null()
    this%w => null()
    this%output_enabled = .false.
    this%log = .true.
    this%start_time = -huge(0.0_rp)
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
         this%particles%n_global
    call neko_log%message(log_buf)
    if (this%periodic_bc%periodic_enabled) then
       write(log_buf, '(A,I0)') "Periodic wrap directions: ", &
            this%periodic_bc%n_periodic_dirs
       call neko_log%message(log_buf)
    end if
    if (this%periodic_bc%rotational_periodic_enabled) then
       write(log_buf, '(A,3(ES13.5,A),ES13.5)') &
            "Rotational periodic sector: theta_min=", &
            this%periodic_bc%rotational_theta_min, ", theta_max=", &
            this%periodic_bc%rotational_theta_max, ", theta_len=", &
            this%periodic_bc%rotational_theta_len, ""
       call neko_log%message(log_buf)
    end if
    write(log_buf, '(A,I0)') "Local particles on rank 0 at init: ", &
         this%particles%n
    if (pe_rank .eq. 0) call neko_log%message(log_buf)
    call neko_log%end_section()
  end subroutine log_status

  subroutine redistribute_particle_ids(this, redist_comm, n_particles_old, &
       particle_ids_local)
    class(lpt_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: redist_comm
    integer, intent(in) :: n_particles_old
    integer, allocatable, intent(out) :: particle_ids_local(:)
    real(kind=rp), allocatable :: sendbuf(:)
    real(kind=rp), allocatable :: recvbuf(:)
    integer :: n_local

    n_local = this%particles%n
    allocate(particle_ids_local(n_local))
    particle_ids_local = 0

    if (n_particles_old .eq. 0 .and. n_local .eq. 0) return

    allocate(sendbuf(n_particles_old))
    allocate(recvbuf(n_local))
    sendbuf = real(this%particles%ids, rp)
    recvbuf = 0.0_rp
    call redist_comm%sendrecv(sendbuf, recvbuf, n_particles_old, n_local)
    particle_ids_local = nint(recvbuf)
  end subroutine redistribute_particle_ids

  subroutine redistribute_velocity_history(this, redist_comm, n_particles_old, &
       vel_particles_lag_local)
    class(lpt_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: redist_comm
    integer, intent(in) :: n_particles_old
    real(kind=rp), allocatable, intent(out) :: vel_particles_lag_local(:,:,:)
    real(kind=rp), allocatable :: sendbuf(:)
    real(kind=rp), allocatable :: recvbuf(:)
    integer :: n_local
    integer :: i
    integer :: j

    n_local = this%particles%n
    allocate(vel_particles_lag_local(3, LPT_VEL_HISTORY_LEN, n_local))
    vel_particles_lag_local = 0.0_rp

    if (n_particles_old .eq. 0 .and. n_local .eq. 0) return

    allocate(sendbuf(n_particles_old))
    allocate(recvbuf(n_local))
    do j = 1, LPT_VEL_HISTORY_LEN
       do i = 1, 3
          sendbuf = this%particles%vel_lag(i, j, :)
          recvbuf = 0.0_rp
          call redist_comm%sendrecv(sendbuf, recvbuf, n_particles_old, n_local)
          vel_particles_lag_local(i, j, :) = recvbuf
       end do
    end do
  end subroutine redistribute_velocity_history

end module lagrangian_particle_tracking
