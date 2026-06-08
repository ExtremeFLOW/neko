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
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use json_utils, only : json_get, json_get_or_default, &
       json_get_subdict_or_empty
  use time_state, only : time_state_t
  use global_interpolation, only : global_interpolation_t
  use glb_intrp_comm, only : glb_intrp_comm_t
  use point, only : point_t
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error
  use file, only : file_t
  use matrix, only : matrix_t
  use math, only : add2s2
  use stack, only : stack_i4_t
  use tensor, only : trsp
  use mesh, only : mesh_t
  use comm, only : pe_rank, pe_size, NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_Allgather, MPI_Gather, &
       MPI_Gatherv, MPI_INTEGER, MPI_SUM, MPI_MIN, MPI_MAX
  use csv_file, only : csv_file_t
  implicit none
  private

  real(kind=rp), parameter :: LPT_PERIODIC_TOL = 1.0e-8_rp
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
     type(particles_t) :: particles
     type(file_t) :: output_file
     logical :: output_enabled = .false.
     logical :: log = .true.
     logical :: periodic_enabled = .false.
     logical :: rotational_periodic_enabled = .false.
     real(kind=rp) :: start_time = -huge(0.0_rp)
     integer :: n_periodic_dirs = 0
     real(kind=rp) :: periodic_shift(3, 3) = 0.0_rp
     real(kind=rp) :: periodic_dir(3, 3) = 0.0_rp
     real(kind=rp) :: periodic_min(3) = 0.0_rp
     real(kind=rp) :: periodic_max(3) = 0.0_rp
     real(kind=rp) :: periodic_len(3) = 0.0_rp
     real(kind=rp) :: rotational_theta_min = 0.0_rp
     real(kind=rp) :: rotational_theta_max = 0.0_rp
     real(kind=rp) :: rotational_theta_len = 0.0_rp
     integer :: n_periodic_maps = 0
     integer, allocatable :: periodic_map_npts(:)
     real(kind=rp), allocatable :: periodic_src_center(:,:)
     real(kind=rp), allocatable :: periodic_src_basis(:,:,:)
     real(kind=rp), allocatable :: periodic_src_bounds_min(:,:)
     real(kind=rp), allocatable :: periodic_src_bounds_max(:,:)
     real(kind=rp), allocatable :: periodic_tgt_center(:,:)
     real(kind=rp), allocatable :: periodic_tgt_basis(:,:,:)
   contains
     procedure, pass(this) :: init => lpt_init_from_json
     procedure, pass(this) :: free => lpt_free
     procedure, pass(this) :: compute_ => lpt_compute
     procedure, private, pass(this) :: read_particles_json
     procedure, private, pass(this) :: read_particles_csv
     procedure, private, pass(this) :: init_periodic_wrap
     procedure, private, pass(this) :: init_periodic_maps
     procedure, private, pass(this) :: wrap_particles_periodic
     procedure, private, pass(this) :: init_rotational_periodic_wrap
     procedure, private, pass(this) :: init_translational_periodic_wrap
     procedure, private, pass(this) :: wrap_particles_rotational
     procedure, private, pass(this) :: wrap_particles_translational
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
    call this%init_periodic_maps(case%fluid%msh)
    call this%init_periodic_wrap(case%fluid%msh, case%fluid%dm_Xh, &
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
    call this%wrap_particles_periodic()
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

    if (this%output_enabled) then
      if (this%output_controller%check(time)) then
         call this%write_output(time, vel)
         call this%output_controller%register_execution()
      end if
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
    if (allocated(this%periodic_map_npts)) deallocate(this%periodic_map_npts)
    if (allocated(this%periodic_src_center)) deallocate(this%periodic_src_center)
    if (allocated(this%periodic_src_basis)) deallocate(this%periodic_src_basis)
    if (allocated(this%periodic_src_bounds_min)) deallocate(this%periodic_src_bounds_min)
    if (allocated(this%periodic_src_bounds_max)) deallocate(this%periodic_src_bounds_max)
    if (allocated(this%periodic_tgt_center)) deallocate(this%periodic_tgt_center)
    if (allocated(this%periodic_tgt_basis)) deallocate(this%periodic_tgt_basis)
    call this%global_interp%free()
    call this%output_file%free()

    this%u => null()
    this%v => null()
    this%w => null()
    this%output_enabled = .false.
    this%log = .true.
    this%start_time = -huge(0.0_rp)
    this%n_periodic_maps = 0
    this%periodic_enabled = .false.
    this%rotational_periodic_enabled = .false.
    this%n_periodic_dirs = 0
    this%periodic_shift = 0.0_rp
    this%periodic_dir = 0.0_rp
    this%periodic_min = 0.0_rp
    this%periodic_max = 0.0_rp
    this%periodic_len = 0.0_rp
    this%rotational_theta_min = 0.0_rp
    this%rotational_theta_max = 0.0_rp
    this%rotational_theta_len = 0.0_rp
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
    if (this%periodic_enabled) then
       write(log_buf, '(A,I0)') "Periodic wrap directions: ", &
            this%n_periodic_dirs
       call neko_log%message(log_buf)
    end if
    if (this%rotational_periodic_enabled) then
       write(log_buf, '(A,3(ES13.5,A),ES13.5)') &
            "Rotational periodic sector: theta_min=", &
            this%rotational_theta_min, ", theta_max=", &
            this%rotational_theta_max, ", theta_len=", &
            this%rotational_theta_len, ""
       call neko_log%message(log_buf)
    end if
    if (this%n_periodic_maps .gt. 0) then
       write(log_buf, '(A,I0)') "Periodic facet transforms: ", &
            this%n_periodic_maps
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

  !> Initialize periodic wrapping from the mesh periodic facet pairs.
  !! For cyclic sectors, the angular extent is inferred directly from the mesh.
  subroutine init_periodic_wrap(this, msh, dm_Xh, coef)
    class(lpt_t), intent(inout) :: this
    type(mesh_t), intent(in), target :: msh
    type(dofmap_t), intent(in) :: dm_Xh
    type(coef_t), intent(in) :: coef
    this%periodic_enabled = .false.
    this%rotational_periodic_enabled = .false.
    this%n_periodic_dirs = 0
    this%periodic_shift = 0.0_rp
    this%periodic_dir = 0.0_rp
    this%periodic_min = 0.0_rp
    this%periodic_max = 0.0_rp
    this%periodic_len = 0.0_rp
    this%rotational_theta_min = 0.0_rp
    this%rotational_theta_max = 0.0_rp
    this%rotational_theta_len = 0.0_rp

    if (msh%periodic%size .eq. 0) return

    call this%init_rotational_periodic_wrap(msh, dm_Xh, coef)
    if (this%rotational_periodic_enabled) return

    call this%init_translational_periodic_wrap(msh, dm_Xh)
  end subroutine init_periodic_wrap

  !> Wrap particles back into the periodic domain before the ownership search.
  subroutine wrap_particles_periodic(this)
    class(lpt_t), intent(inout) :: this
    integer :: i
    integer :: k
    integer :: n_iters
    logical :: mapped

    if (this%particles%n .eq. 0) return

    if (this%rotational_periodic_enabled) then
       call this%wrap_particles_rotational()
       return
    end if

    if (this%n_periodic_maps .gt. 0) then
       do i = 1, this%particles%n
          do n_iters = 1, 3
             mapped = .false.
             do k = 1, this%n_periodic_maps
                if (lpt_apply_periodic_map_if_needed(this, i, k)) then
                   mapped = .true.
                   exit
                end if
             end do
             if (.not. mapped) exit
          end do
       end do
    end if

    call this%wrap_particles_translational()
  end subroutine wrap_particles_periodic

  subroutine init_rotational_periodic_wrap(this, msh, dm_Xh, coef)
    class(lpt_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    type(dofmap_t), intent(in) :: dm_Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp) :: theta_min
    real(kind=rp) :: theta_max
    real(kind=rp) :: pi
    integer :: ierr

    if (.not. coef%cyclic .or. msh%gdim .lt. 2) return

    pi = acos(-1.0_rp)
    if (msh%nelv .gt. 0) then
       theta_min = minval(modulo(atan2(dm_Xh%y, dm_Xh%x) + &
            2.0_rp * pi, 2.0_rp * pi))
       theta_max = maxval(modulo(atan2(dm_Xh%y, dm_Xh%x) + &
            2.0_rp * pi, 2.0_rp * pi))
    else
       theta_min = huge(0.0_rp)
       theta_max = -huge(0.0_rp)
    end if

    call MPI_Allreduce(theta_min, this%rotational_theta_min, 1, &
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
    call MPI_Allreduce(theta_max, this%rotational_theta_max, 1, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)

    this%rotational_theta_len = this%rotational_theta_max - &
         this%rotational_theta_min
    this%rotational_periodic_enabled = &
         this%rotational_theta_len .gt. LPT_PERIODIC_TOL
  end subroutine init_rotational_periodic_wrap

  subroutine init_translational_periodic_wrap(this, msh, dm_Xh)
    class(lpt_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    type(dofmap_t), intent(in) :: dm_Xh
    integer :: i
    integer :: idx
    integer :: ierr
    real(kind=rp) :: local_min(3)
    real(kind=rp) :: local_max(3)
    real(kind=rp) :: global_min(3)
    real(kind=rp) :: global_max(3)

    if (msh%nelv .gt. 0) then
       local_min(1) = minval(dm_Xh%x)
       local_max(1) = maxval(dm_Xh%x)
       local_min(2) = minval(dm_Xh%y)
       local_max(2) = maxval(dm_Xh%y)
       local_min(3) = minval(dm_Xh%z)
       local_max(3) = maxval(dm_Xh%z)
    else
       local_min = huge(0.0_rp)
       local_max = -huge(0.0_rp)
    end if

    call MPI_Allreduce(local_min, global_min, 3, MPI_REAL_PRECISION, MPI_MIN, &
         NEKO_COMM, ierr)
    call MPI_Allreduce(local_max, global_max, 3, MPI_REAL_PRECISION, MPI_MAX, &
         NEKO_COMM, ierr)

    do i = 1, msh%periodic%size
       idx = lpt_periodic_dir_from_facet(msh%gdim, msh%periodic%facet_el(i)%x(1))
       if (idx .eq. 0 .or. this%periodic_len(idx) .gt. 0.0_rp) cycle
       this%n_periodic_dirs = max(this%n_periodic_dirs, idx)
       this%periodic_dir(:, idx) = 0.0_rp
       this%periodic_dir(idx, idx) = 1.0_rp
       this%periodic_min(idx) = global_min(idx)
       this%periodic_max(idx) = global_max(idx)
       this%periodic_len(idx) = global_max(idx) - global_min(idx)
       this%periodic_shift(:, idx) = 0.0_rp
       this%periodic_shift(idx, idx) = this%periodic_len(idx)
    end do

    this%periodic_enabled = this%n_periodic_dirs .gt. 0
  end subroutine init_translational_periodic_wrap

  subroutine wrap_particles_rotational(this)
    class(lpt_t), intent(inout) :: this
    integer :: i
    real(kind=rp) :: radius
    real(kind=rp) :: theta
    real(kind=rp) :: pi

    pi = acos(-1.0_rp)
    do i = 1, this%particles%n
       radius = norm2(this%particles%xyz(1:2, i))
       theta = modulo(atan2(this%particles%xyz(2, i), &
            this%particles%xyz(1, i)) + 2.0_rp * pi, 2.0_rp * pi)

       do while (theta .lt. this%rotational_theta_min - LPT_PERIODIC_TOL)
          theta = theta + this%rotational_theta_len
       end do

       do while (theta .gt. this%rotational_theta_max + LPT_PERIODIC_TOL)
          theta = theta - this%rotational_theta_len
       end do

       this%particles%xyz(1, i) = radius * cos(theta)
       this%particles%xyz(2, i) = radius * sin(theta)
    end do
  end subroutine wrap_particles_rotational

  subroutine wrap_particles_translational(this)
    class(lpt_t), intent(inout) :: this
    integer :: i
    integer :: j
    real(kind=rp) :: coord

    if (.not. this%periodic_enabled) return

    do i = 1, this%particles%n
       do j = 1, this%n_periodic_dirs
          coord = dot_product(this%particles%xyz(:, i), this%periodic_dir(:, j))

          do while (coord .lt. this%periodic_min(j) - LPT_PERIODIC_TOL)
             this%particles%xyz(:, i) = this%particles%xyz(:, i) + &
                  this%periodic_shift(:, j)
             coord = coord + this%periodic_len(j)
          end do

          do while (coord .gt. this%periodic_max(j) + LPT_PERIODIC_TOL)
             this%particles%xyz(:, i) = this%particles%xyz(:, i) - &
                  this%periodic_shift(:, j)
             coord = coord - this%periodic_len(j)
          end do
       end do
    end do
  end subroutine wrap_particles_translational

  !> Map a periodic facet number to its Cartesian wrap direction.
  integer function lpt_periodic_dir_from_facet(gdim, facet) result(idx)
    integer, intent(in) :: gdim
    integer, intent(in) :: facet

    idx = 0
    select case (facet)
    case (1, 2)
       idx = 1
    case (3, 4)
       idx = 2
    case (5, 6)
       if (gdim .eq. 3) idx = 3
    end select
  end function lpt_periodic_dir_from_facet

  !> Build local periodic facet transforms, including cyclic mappings.
  subroutine init_periodic_maps(this, msh)
    class(lpt_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    integer :: n_local
    integer :: n_global
    integer :: i
    integer :: j
    integer :: idx
    integer :: ierr
    integer :: max_local
    integer, allocatable :: counts(:)
    integer, allocatable :: global_meta(:)
    integer, allocatable :: padded_meta(:)
    integer, allocatable :: gathered_meta(:)
    real(kind=rp), allocatable :: global_geom(:)
    real(kind=rp), allocatable :: padded_geom(:)
    real(kind=rp), allocatable :: gathered_geom(:)
    real(kind=rp) :: src_pts(3, 4)
    real(kind=rp) :: tgt_pts(3, 4)
    real(kind=rp) :: src_centroid(3)
    real(kind=rp) :: tgt_centroid(3)
    integer :: match_idx
    integer :: npts

    this%n_periodic_maps = 0
    if (allocated(this%periodic_map_npts)) deallocate(this%periodic_map_npts)
    if (allocated(this%periodic_src_center)) deallocate(this%periodic_src_center)
    if (allocated(this%periodic_src_basis)) deallocate(this%periodic_src_basis)
    if (allocated(this%periodic_src_bounds_min)) deallocate(this%periodic_src_bounds_min)
    if (allocated(this%periodic_src_bounds_max)) deallocate(this%periodic_src_bounds_max)
    if (allocated(this%periodic_tgt_center)) deallocate(this%periodic_tgt_center)
    if (allocated(this%periodic_tgt_basis)) deallocate(this%periodic_tgt_basis)

    n_local = msh%periodic%size
    if (n_local .eq. 0) return

    allocate(counts(pe_size))
    call MPI_Allgather(n_local, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, &
         NEKO_COMM, ierr)
    n_global = sum(counts)
    max_local = max(1, maxval(counts))

    allocate(padded_meta(7 * max_local))
    allocate(padded_geom(15 * max_local))
    padded_meta = 0
    padded_geom = 0.0_rp
    do i = 1, n_local
       npts = merge(4, 2, msh%gdim .eq. 3)
       padded_meta(7 * (i - 1) + 1) = msh%periodic%facet_el(i)%x(1)
       padded_meta(7 * (i - 1) + 2) = msh%periodic%facet_el(i)%x(2)
       padded_meta(7 * (i - 1) + 3) = npts
       padded_meta(7 * (i - 1) + 4:7 * (i - 1) + 7) = msh%periodic%org_ids(i)%x

       call lpt_get_facet_points(msh, msh%periodic%facet_el(i)%x(2), &
            msh%periodic%facet_el(i)%x(1), src_pts, src_centroid)
       padded_geom(15 * (i - 1) + 1:15 * (i - 1) + 12) = reshape(src_pts, [12])
       padded_geom(15 * (i - 1) + 13:15 * (i - 1) + 15) = src_centroid
    end do

    allocate(gathered_meta(7 * max_local * pe_size))
    allocate(gathered_geom(15 * max_local * pe_size))
    call MPI_Allgather(padded_meta, 7 * max_local, MPI_INTEGER, gathered_meta, &
         7 * max_local, MPI_INTEGER, NEKO_COMM, ierr)
    call MPI_Allgather(padded_geom, 15 * max_local, MPI_REAL_PRECISION, &
         gathered_geom, 15 * max_local, MPI_REAL_PRECISION, NEKO_COMM, ierr)

    allocate(global_meta(7 * n_global))
    allocate(global_geom(15 * n_global))
    idx = 0
    do i = 1, pe_size
       if (counts(i) .gt. 0) then
          global_meta(7 * idx + 1:7 * (idx + counts(i))) = &
               gathered_meta(7 * max_local * (i - 1) + 1: &
               7 * max_local * (i - 1) + 7 * counts(i))
          global_geom(15 * idx + 1:15 * (idx + counts(i))) = &
               gathered_geom(15 * max_local * (i - 1) + 1: &
               15 * max_local * (i - 1) + 15 * counts(i))
          idx = idx + counts(i)
       end if
    end do

    allocate(this%periodic_map_npts(n_local))
    allocate(this%periodic_src_center(3, n_local))
    allocate(this%periodic_src_basis(3, 3, n_local))
    allocate(this%periodic_src_bounds_min(2, n_local))
    allocate(this%periodic_src_bounds_max(2, n_local))
    allocate(this%periodic_tgt_center(3, n_local))
    allocate(this%periodic_tgt_basis(3, 3, n_local))

    do i = 1, n_local
       match_idx = 0
       do j = 1, n_global
          if (global_meta(7 * (j - 1) + 1) .eq. msh%periodic%p_facet_el(i)%x(1) &
               .and. global_meta(7 * (j - 1) + 2) .eq. &
               msh%periodic%p_facet_el(i)%x(2)) then
             match_idx = j
             exit
          end if
       end do
       if (match_idx .eq. 0) cycle

       npts = merge(4, 2, msh%gdim .eq. 3)
       call lpt_get_facet_points(msh, msh%periodic%facet_el(i)%x(2), &
            msh%periodic%facet_el(i)%x(1), src_pts, src_centroid)
       call lpt_reorder_facet_points( &
            global_meta(7 * (match_idx - 1) + 4:7 * (match_idx - 1) + 7), &
            reshape(global_geom(15 * (match_idx - 1) + 1:15 * (match_idx - 1) + 12), [3, 4]), &
            msh%periodic%p_ids(i)%x, npts, tgt_pts)
       tgt_centroid = global_geom(15 * (match_idx - 1) + 13:15 * (match_idx - 1) + 15)

       idx = this%n_periodic_maps + 1
       this%n_periodic_maps = idx
       this%periodic_map_npts(idx) = npts
       call lpt_build_facet_basis(src_pts, src_centroid, npts, &
            this%periodic_src_center(:, idx), this%periodic_src_basis(:, :, idx))
       call lpt_get_facet_bounds(src_pts, this%periodic_src_center(:, idx), &
            this%periodic_src_basis(:, :, idx), npts, &
            this%periodic_src_bounds_min(:, idx), &
            this%periodic_src_bounds_max(:, idx))
       call lpt_build_facet_basis(tgt_pts, tgt_centroid, npts, &
            this%periodic_tgt_center(:, idx), this%periodic_tgt_basis(:, :, idx))
    end do

    deallocate(global_meta)
    deallocate(global_geom)
    deallocate(padded_meta)
    deallocate(gathered_meta)
    deallocate(padded_geom)
    deallocate(gathered_geom)
    deallocate(counts)
  end subroutine init_periodic_maps

  logical function lpt_apply_periodic_map_if_needed(this, i_particle, i_map) result(mapped)
    class(lpt_t), intent(inout) :: this
    integer, intent(in) :: i_particle
    integer, intent(in) :: i_map
    real(kind=rp) :: rel(3)
    real(kind=rp) :: xi(3)

    mapped = .false.
    rel = this%particles%xyz(:, i_particle) - this%periodic_src_center(:, i_map)
    if (this%periodic_map_npts(i_map) .eq. 4) then
       xi(1) = dot_product(rel, this%periodic_src_basis(:, 1, i_map))
       xi(2) = dot_product(rel, this%periodic_src_basis(:, 2, i_map))
       xi(3) = dot_product(rel, this%periodic_src_basis(:, 3, i_map))
       if (xi(1) .ge. this%periodic_src_bounds_min(1, i_map) - &
            LPT_PERIODIC_TOL .and. &
            xi(1) .le. this%periodic_src_bounds_max(1, i_map) + &
            LPT_PERIODIC_TOL .and. &
            xi(2) .ge. this%periodic_src_bounds_min(2, i_map) - &
            LPT_PERIODIC_TOL .and. &
            xi(2) .le. this%periodic_src_bounds_max(2, i_map) + &
            LPT_PERIODIC_TOL .and. xi(3) .lt. -LPT_PERIODIC_TOL) then
          this%particles%xyz(:, i_particle) = this%periodic_tgt_center(:, i_map) + &
               xi(1) * this%periodic_tgt_basis(:, 1, i_map) + &
               xi(2) * this%periodic_tgt_basis(:, 2, i_map) + &
               xi(3) * this%periodic_tgt_basis(:, 3, i_map)
          mapped = .true.
       end if
    else
       xi(1) = dot_product(rel, this%periodic_src_basis(:, 1, i_map))
       xi(2) = dot_product(rel, this%periodic_src_basis(:, 3, i_map))
       if (xi(1) .ge. this%periodic_src_bounds_min(1, i_map) - &
            LPT_PERIODIC_TOL .and. &
            xi(1) .le. this%periodic_src_bounds_max(1, i_map) + &
            LPT_PERIODIC_TOL .and. xi(2) .lt. -LPT_PERIODIC_TOL) then
          this%particles%xyz(:, i_particle) = this%periodic_tgt_center(:, i_map) + &
               xi(1) * this%periodic_tgt_basis(:, 1, i_map) + &
               xi(2) * this%periodic_tgt_basis(:, 3, i_map)
          mapped = .true.
       end if
    end if
  end function lpt_apply_periodic_map_if_needed

  subroutine lpt_reorder_facet_points(ids_in, pts_in, ids_target, npts, pts_out)
    integer, intent(in) :: ids_in(4)
    real(kind=rp), intent(in) :: pts_in(3, 4)
    integer, intent(in) :: ids_target(4)
    integer, intent(in) :: npts
    real(kind=rp), intent(out) :: pts_out(3, 4)
    integer :: i
    integer :: j

    pts_out = 0.0_rp
    do i = 1, npts
       do j = 1, npts
          if (ids_in(j) .eq. ids_target(i)) then
             pts_out(:, i) = pts_in(:, j)
             exit
          end if
       end do
    end do
  end subroutine lpt_reorder_facet_points

  subroutine lpt_get_facet_points(msh, el, facet, pts, el_centroid)
    type(mesh_t), intent(in) :: msh
    integer, intent(in) :: el
    integer, intent(in) :: facet
    real(kind=rp), intent(out) :: pts(3, 4)
    real(kind=rp), intent(out) :: el_centroid(3)
    type(point_t) :: centroid
    integer, dimension(4, 6) :: face_nodes = reshape([ &
         1,5,7,3, &
         2,6,8,4, &
         1,2,6,5, &
         3,4,8,7, &
         1,2,4,3, &
         5,6,8,7], [4, 6])
    integer, dimension(2, 4) :: edge_nodes = reshape([ &
         1,3, &
         2,4, &
         1,2, &
         3,4], [2, 4])
    integer :: i

    pts = 0.0_rp
    centroid = msh%elements(el)%e%centroid()
    el_centroid = real(centroid%x, rp)
    if (msh%gdim .eq. 3) then
       do i = 1, 4
          pts(:, i) = real(msh%elements(el)%e%pts(face_nodes(i, facet))%p%x, rp)
       end do
    else
       do i = 1, 2
          pts(:, i) = real(msh%elements(el)%e%pts(edge_nodes(i, facet))%p%x, rp)
       end do
    end if
  end subroutine lpt_get_facet_points

  subroutine lpt_build_facet_basis(pts, el_centroid, npts, center, basis)
    real(kind=rp), intent(in) :: pts(3, 4)
    real(kind=rp), intent(in) :: el_centroid(3)
    integer, intent(in) :: npts
    real(kind=rp), intent(out) :: center(3)
    real(kind=rp), intent(out) :: basis(3, 3)
    real(kind=rp) :: v1(3)
    real(kind=rp) :: v2(3)
    real(kind=rp) :: n(3)

    center = 0.0_rp
    basis = 0.0_rp
    center = sum(pts(:, 1:npts), dim = 2) / real(npts, rp)

    v1 = pts(:, 2) - pts(:, 1)
    call lpt_normalize(v1)
    basis(:, 1) = v1

    if (npts .eq. 4) then
       v2 = pts(:, 4) - pts(:, 1)
       n = lpt_cross(v1, v2)
       call lpt_normalize(n)
       if (dot_product(el_centroid - center, n) .lt. 0.0_rp) n = -n
       basis(:, 3) = n
       basis(:, 2) = lpt_cross(n, v1)
       call lpt_normalize(basis(:, 2))
    else
       n = el_centroid - center
       call lpt_normalize(n)
       basis(:, 3) = n
    end if
  end subroutine lpt_build_facet_basis

  subroutine lpt_get_facet_bounds(pts, center, basis, npts, bounds_min, bounds_max)
    real(kind=rp), intent(in) :: pts(3, 4)
    real(kind=rp), intent(in) :: center(3)
    real(kind=rp), intent(in) :: basis(3, 3)
    integer, intent(in) :: npts
    real(kind=rp), intent(out) :: bounds_min(2)
    real(kind=rp), intent(out) :: bounds_max(2)
    real(kind=rp) :: rel(3)
    real(kind=rp) :: xi1
    real(kind=rp) :: xi2
    integer :: i

    bounds_min = huge(0.0_rp)
    bounds_max = -huge(0.0_rp)
    do i = 1, npts
       rel = pts(:, i) - center
       xi1 = dot_product(rel, basis(:, 1))
       bounds_min(1) = min(bounds_min(1), xi1)
       bounds_max(1) = max(bounds_max(1), xi1)
       if (npts .eq. 4) then
          xi2 = dot_product(rel, basis(:, 2))
          bounds_min(2) = min(bounds_min(2), xi2)
          bounds_max(2) = max(bounds_max(2), xi2)
       end if
    end do
    if (npts .eq. 2) then
       bounds_min(2) = 0.0_rp
       bounds_max(2) = 0.0_rp
    end if
  end subroutine lpt_get_facet_bounds

  subroutine lpt_push_comp_range(s, offset, n_comp)
    type(stack_i4_t), intent(inout) :: s
    integer, intent(in) :: offset
    integer, intent(in) :: n_comp
    integer :: i

    do i = 1, n_comp
       call s%push(offset + i)
    end do
  end subroutine lpt_push_comp_range

  function lpt_cross(a, b) result(c)
    real(kind=rp), intent(in) :: a(3)
    real(kind=rp), intent(in) :: b(3)
    real(kind=rp) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end function lpt_cross

  subroutine lpt_normalize(v)
    real(kind=rp), intent(inout) :: v(3)
    real(kind=rp) :: vnorm

    vnorm = norm2(v)
    if (vnorm .gt. LPT_PERIODIC_TOL) v = v / vnorm
  end subroutine lpt_normalize

end module lagrangian_particle_tracking
