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
  use neko_config, only : NEKO_BCKND_DEVICE
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use registry, only : neko_registry
  use field, only : field_t
  use case, only : case_t
  use dofmap, only : dofmap_t
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
  use device, only : device_map, device_memcpy, HOST_TO_DEVICE
  use mesh, only : mesh_t
  use comm, only : pe_rank, pe_size, NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_Allgather, MPI_Exscan, MPI_Gather, &
       MPI_Gatherv, MPI_INTEGER, MPI_SUM, MPI_MIN, MPI_MAX, MPI_Request, &
       MPI_Isend, MPI_Irecv, MPI_Waitall, MPI_STATUSES_IGNORE
  use csv_file, only : csv_file_t
  implicit none
  private

  real(kind=rp), parameter :: LPT_PERIODIC_TOL = 1.0e-8_rp
  integer, parameter :: LPT_VEL_HISTORY_LEN = 2

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
     !> Stable particle ids that travel with the particle during redistribution.
     integer, allocatable :: particle_ids(:)
     !> Lagged carrier velocities, kept with each particle for multistep time
     !! integration.
     real(kind=rp), allocatable :: vel_particles_lag(:,:,:)
     !> Optional CSV output file.
     type(file_t) :: output_file
     !> Whether trajectory output is enabled.
     logical :: output_enabled = .false.
     !> Whether to emit informational logs.
     logical :: log = .true.
     !> Whether periodic wrapping is active.
     logical :: periodic_enabled = .false.
     !> Time after which tracking should start.
     real(kind=rp) :: start_time = -huge(0.0_rp)
     !> Number of particles currently local to this rank.
     integer :: n_particles = 0
     !> Total number of particles across all ranks.
     integer :: n_global_particles = 0
     !> Number of unique periodic translation directions.
     integer :: n_periodic_dirs = 0
     !> Translation vectors from the low to the high periodic boundary.
     real(kind=rp) :: periodic_shift(3, 3) = 0.0_rp
     !> Unit vectors associated with the periodic translations.
     real(kind=rp) :: periodic_dir(3, 3) = 0.0_rp
     !> Low-side coordinate along each periodic direction.
     real(kind=rp) :: periodic_min(3) = 0.0_rp
     !> High-side coordinate along each periodic direction.
     real(kind=rp) :: periodic_max(3) = 0.0_rp
     !> Period length along each periodic direction.
     real(kind=rp) :: periodic_len(3) = 0.0_rp
     !> Number of local periodic facet transforms.
     integer :: n_periodic_maps = 0
     !> Owning element for each local periodic facet transform.
     integer, allocatable :: periodic_map_el(:)
     !> Local facet index for each periodic transform.
     integer, allocatable :: periodic_map_facet(:)
     !> Number of facet points for each transform (2 in 2D, 4 in 3D).
     integer, allocatable :: periodic_map_npts(:)
     !> Source facet centers.
     real(kind=rp), allocatable :: periodic_src_center(:,:)
     !> Source facet bases: tangent 1, tangent 2, inward normal.
     real(kind=rp), allocatable :: periodic_src_basis(:,:,:)
     !> Target facet centers.
     real(kind=rp), allocatable :: periodic_tgt_center(:,:)
     !> Target facet bases: tangent 1, tangent 2, inward normal.
     real(kind=rp), allocatable :: periodic_tgt_basis(:,:,:)
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
     !> Initialize periodic wrapping metadata from the mesh.
     procedure, private, pass(this) :: init_periodic_wrap
     !> Build local cyclic facet transforms from periodic facet pairings.
     procedure, private, pass(this) :: init_periodic_maps
     !> Wrap local particles back into the periodic domain.
     procedure, private, pass(this) :: wrap_particles_periodic
     !> Redistribute particles to the owning MPI rank.
     procedure, private, pass(this) :: redistribute_particles
     !> Interpolate the carrier velocity at the local particles.
     procedure, private, pass(this) :: evaluate_velocity
     !> Write a trajectory snapshot.
     procedure, private, pass(this) :: write_output
     !> Redistribute per-particle ids with the interpolation ownership.
     procedure, private, pass(this) :: redistribute_particle_ids
     !> Redistribute lagged per-particle data with the interpolation ownership.
     procedure, private, pass(this) :: redistribute_velocity_history
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
    call this%init_periodic_wrap(case%fluid%msh, case%fluid%dm_Xh)
    call this%init_periodic_maps(case%fluid%msh)
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
    integer :: i

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

    allocate(this%particle_ids(this%n_particles))
    do i = 1, this%n_particles
       this%particle_ids(i) = i
    end do

    allocate(this%vel_particles_lag(3, LPT_VEL_HISTORY_LEN, this%n_particles))
    this%vel_particles_lag = 0.0_rp

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
    integer, allocatable :: particle_ids_local(:)
    real(kind=rp), allocatable :: vel_particles_lag_local(:,:,:)
    integer :: ierr

    call this%wrap_particles_periodic()
    call this%global_interp%find_points(this%xyz_particles, this%n_particles)
    call this%redistribute_particle_ids(particle_ids_local)
    call this%redistribute_velocity_history(vel_particles_lag_local)

    this%n_particles = this%global_interp%n_points_local

    if (allocated(this%xyz_particles)) deallocate(this%xyz_particles)
    allocate(this%xyz_particles(3, this%n_particles))
    this%xyz_particles = this%global_interp%xyz_local

    call this%global_interp%free_points()
    this%global_interp%n_points = this%n_particles
    allocate(this%global_interp%pe_owner(this%n_particles))
    allocate(this%global_interp%el_owner0(this%n_particles))
    allocate(this%global_interp%xyz(3, this%n_particles))
    allocate(this%global_interp%rst(3, this%n_particles))
    this%global_interp%pe_owner = this%global_interp%pe_rank
    this%global_interp%el_owner0 = this%global_interp%el_owner0_local
    this%global_interp%xyz = this%global_interp%xyz_local
    this%global_interp%rst = this%global_interp%rst_local
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%global_interp%el_owner0, &
            this%global_interp%el_owner0_d, this%n_particles)
       call device_memcpy(this%global_interp%el_owner0, &
            this%global_interp%el_owner0_d, this%n_particles, HOST_TO_DEVICE, &
            sync = .true.)
    end if
    this%global_interp%all_points_local = .true.

    call move_alloc(particle_ids_local, this%particle_ids)
    call move_alloc(vel_particles_lag_local, this%vel_particles_lag)
    call MPI_Allreduce(this%n_particles, this%n_global_particles, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
  end subroutine redistribute_particles

  !> Redistribute stable particle ids to the owning rank, matching the
  !! ordering produced by `global_interpolation_t`.
  subroutine redistribute_particle_ids(this, particle_ids_local)
    class(lpt_t), intent(inout) :: this
    integer, allocatable, intent(out) :: particle_ids_local(:)
    integer, allocatable :: sendbuf(:)
    integer, allocatable :: recvbuf(:)
    integer, allocatable :: send_offsets(:)
    integer, allocatable :: recv_offsets(:)
    type(MPI_Request), allocatable :: requests(:)
    integer, pointer :: point_ids(:)
    integer :: rank
    integer :: i
    integer :: idx
    integer :: n_local
    integer :: nreq
    integer :: ierr

    n_local = this%global_interp%n_points_local
    allocate(particle_ids_local(n_local))
    particle_ids_local = 0

    if (this%n_particles .eq. 0 .and. n_local .eq. 0) return

    allocate(send_offsets(0:pe_size - 1))
    allocate(recv_offsets(0:pe_size - 1))
    send_offsets = 0
    recv_offsets = 0

    do rank = 0, pe_size - 1
       send_offsets(rank) = this%global_interp%n_points_offset_pe(rank)
       recv_offsets(rank) = this%global_interp%n_points_offset_pe_local(rank)
    end do

    allocate(sendbuf(this%n_particles))
    sendbuf = 0
    idx = 0
    do rank = 0, pe_size - 1
       if (this%global_interp%n_points_pe(rank) .le. 0) cycle
       point_ids => this%global_interp%points_at_pe(rank)%array()
       do i = 1, this%global_interp%n_points_pe(rank)
          idx = idx + 1
          sendbuf(idx) = this%particle_ids(point_ids(i))
       end do
    end do

    allocate(recvbuf(n_local))
    recvbuf = 0

    nreq = 0
    do rank = 0, pe_size - 1
       if (rank .ne. pe_rank .and. this%global_interp%n_points_pe_local(rank) .gt. 0) then
          nreq = nreq + 1
       end if
       if (rank .ne. pe_rank .and. this%global_interp%n_points_pe(rank) .gt. 0) then
          nreq = nreq + 1
       end if
    end do
    allocate(requests(max(1, nreq)))

    idx = 0
    do rank = 0, pe_size - 1
       if (rank .eq. pe_rank) cycle
       if (this%global_interp%n_points_pe_local(rank) .gt. 0) then
          idx = idx + 1
          call MPI_Irecv(recvbuf(recv_offsets(rank) + 1), &
               this%global_interp%n_points_pe_local(rank), MPI_INTEGER, rank, 0, &
               NEKO_COMM, requests(idx), ierr)
       end if
    end do

    do rank = 0, pe_size - 1
       if (this%global_interp%n_points_pe(rank) .le. 0) cycle
       if (rank .eq. pe_rank) then
          recvbuf(recv_offsets(rank) + 1:recv_offsets(rank) + &
               this%global_interp%n_points_pe(rank)) = &
               sendbuf(send_offsets(rank) + 1:send_offsets(rank) + &
               this%global_interp%n_points_pe(rank))
       else
          idx = idx + 1
          call MPI_Isend(sendbuf(send_offsets(rank) + 1), &
               this%global_interp%n_points_pe(rank), MPI_INTEGER, rank, 0, &
               NEKO_COMM, requests(idx), ierr)
       end if
    end do

    if (nreq .gt. 0) then
       call MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE, ierr)
    end if

    particle_ids_local = recvbuf

    deallocate(recvbuf)
    deallocate(sendbuf)
    deallocate(send_offsets)
    deallocate(recv_offsets)
    deallocate(requests)
  end subroutine redistribute_particle_ids

  !> Redistribute the lagged per-particle velocities to the owning rank,
  !! matching the ordering produced by `global_interpolation_t`.
  subroutine redistribute_velocity_history(this, vel_particles_lag_local)
    class(lpt_t), intent(inout) :: this
    real(kind=rp), allocatable, intent(out) :: vel_particles_lag_local(:,:,:)
    real(kind=rp), allocatable :: sendbuf(:)
    real(kind=rp), allocatable :: recvbuf(:)
    integer, allocatable :: send_offsets(:)
    integer, allocatable :: recv_offsets(:)
    type(MPI_Request), allocatable :: requests(:)
    integer, pointer :: point_ids(:)
    integer :: rank
    integer :: i
    integer :: j
    integer :: idx
    integer :: n_local
    integer :: nreq
    integer :: ierr

    n_local = this%global_interp%n_points_local
    allocate(vel_particles_lag_local(3, LPT_VEL_HISTORY_LEN, n_local))
    vel_particles_lag_local = 0.0_rp

    if (this%n_particles .eq. 0 .and. n_local .eq. 0) return

    allocate(send_offsets(0:pe_size - 1))
    allocate(recv_offsets(0:pe_size - 1))
    send_offsets = 0
    recv_offsets = 0

    do rank = 0, pe_size - 1
       send_offsets(rank) = 3 * LPT_VEL_HISTORY_LEN * &
            this%global_interp%n_points_offset_pe(rank)
       recv_offsets(rank) = 3 * LPT_VEL_HISTORY_LEN * &
            this%global_interp%n_points_offset_pe_local(rank)
    end do

    allocate(sendbuf(3 * LPT_VEL_HISTORY_LEN * this%n_particles))
    sendbuf = 0.0_rp
    idx = 0
    do rank = 0, pe_size - 1
       if (this%global_interp%n_points_pe(rank) .le. 0) cycle
       point_ids => this%global_interp%points_at_pe(rank)%array()
       do i = 1, this%global_interp%n_points_pe(rank)
          do j = 1, LPT_VEL_HISTORY_LEN
             sendbuf(idx + 1:idx + 3) = &
                  this%vel_particles_lag(:, j, point_ids(i))
             idx = idx + 3
          end do
       end do
    end do

    allocate(recvbuf(3 * LPT_VEL_HISTORY_LEN * n_local))
    recvbuf = 0.0_rp

    nreq = 0
    do rank = 0, pe_size - 1
       if (rank .ne. pe_rank .and. this%global_interp%n_points_pe_local(rank) .gt. 0) then
          nreq = nreq + 1
       end if
       if (rank .ne. pe_rank .and. this%global_interp%n_points_pe(rank) .gt. 0) then
          nreq = nreq + 1
       end if
    end do
    allocate(requests(max(1, nreq)))

    idx = 0
    do rank = 0, pe_size - 1
       if (rank .eq. pe_rank) cycle
       if (this%global_interp%n_points_pe_local(rank) .gt. 0) then
          idx = idx + 1
          call MPI_Irecv(recvbuf(recv_offsets(rank) + 1), &
               3 * LPT_VEL_HISTORY_LEN * this%global_interp%n_points_pe_local(rank), &
               MPI_REAL_PRECISION, rank, 0, NEKO_COMM, requests(idx), ierr)
       end if
    end do

    do rank = 0, pe_size - 1
       if (this%global_interp%n_points_pe(rank) .le. 0) cycle
       if (rank .eq. pe_rank) then
          recvbuf(recv_offsets(rank) + 1:recv_offsets(rank) + &
               3 * LPT_VEL_HISTORY_LEN * this%global_interp%n_points_pe(rank)) = &
               sendbuf(send_offsets(rank) + 1:send_offsets(rank) + &
               3 * LPT_VEL_HISTORY_LEN * this%global_interp%n_points_pe(rank))
       else
          idx = idx + 1
          call MPI_Isend(sendbuf(send_offsets(rank) + 1), &
               3 * LPT_VEL_HISTORY_LEN * this%global_interp%n_points_pe(rank), &
               MPI_REAL_PRECISION, rank, 0, NEKO_COMM, requests(idx), ierr)
       end if
    end do

    if (nreq .gt. 0) then
      call MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE, ierr)
    end if

    idx = 0
    do i = 1, n_local
       do j = 1, LPT_VEL_HISTORY_LEN
          vel_particles_lag_local(:, j, i) = recvbuf(idx + 1:idx + 3)
          idx = idx + 3
       end do
    end do

    deallocate(recvbuf)
    deallocate(sendbuf)
    deallocate(send_offsets)
    deallocate(recv_offsets)
    deallocate(requests)
  end subroutine redistribute_velocity_history

  !> Initialize periodic wrapping from the mesh periodic facet pairs.
  subroutine init_periodic_wrap(this, msh, dm_Xh)
    class(lpt_t), intent(inout) :: this
    type(mesh_t), intent(in), target :: msh
    type(dofmap_t), intent(in) :: dm_Xh
    integer :: i
    integer :: idx
    integer :: ierr
    real(kind=rp) :: local_min(3)
    real(kind=rp) :: local_max(3)
    real(kind=rp) :: global_min(3)
    real(kind=rp) :: global_max(3)

    this%periodic_enabled = .false.
    this%n_periodic_dirs = 0
    this%periodic_shift = 0.0_rp
    this%periodic_dir = 0.0_rp
    this%periodic_min = 0.0_rp
    this%periodic_max = 0.0_rp
    this%periodic_len = 0.0_rp

    if (msh%periodic%size .eq. 0) return

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
       if (idx .eq. 0) then
          cycle
       end if
       if (this%periodic_len(idx) .le. 0.0_rp) then
          this%n_periodic_dirs = max(this%n_periodic_dirs, idx)
          this%periodic_dir(:, idx) = 0.0_rp
          this%periodic_dir(idx, idx) = 1.0_rp
          this%periodic_min(idx) = global_min(idx)
          this%periodic_max(idx) = global_max(idx)
          this%periodic_len(idx) = global_max(idx) - global_min(idx)
          this%periodic_shift(:, idx) = 0.0_rp
          this%periodic_shift(idx, idx) = this%periodic_len(idx)
       end if
    end do

    this%periodic_enabled = this%n_periodic_dirs .gt. 0
  end subroutine init_periodic_wrap

  !> Wrap particles back into the periodic domain before the ownership search.
  subroutine wrap_particles_periodic(this)
    class(lpt_t), intent(inout) :: this
    integer :: i
    integer :: j
    integer :: k
    integer :: n_iters
    real(kind=rp) :: coord

    if (this%n_particles .eq. 0) return

    if (this%n_periodic_maps .gt. 0) then
       if (allocated(this%global_interp%el_owner0_local)) then
          do i = 1, this%n_particles
             do n_iters = 1, 3
                do k = 1, this%n_periodic_maps
                   if (this%global_interp%el_owner0_local(i) .eq. &
                        this%periodic_map_el(k)) then
                      call lpt_apply_periodic_map_if_needed(this, i, k)
                   end if
                end do
             end do
          end do
       end if
    end if

    if (.not. this%periodic_enabled) return

    do i = 1, this%n_particles
       do j = 1, this%n_periodic_dirs
          coord = dot_product(this%xyz_particles(:, i), this%periodic_dir(:, j))

          do while (coord .lt. this%periodic_min(j) - LPT_PERIODIC_TOL)
             this%xyz_particles(:, i) = this%xyz_particles(:, i) + &
                  this%periodic_shift(:, j)
             coord = coord + this%periodic_len(j)
          end do

          do while (coord .gt. this%periodic_max(j) + LPT_PERIODIC_TOL)
             this%xyz_particles(:, i) = this%xyz_particles(:, i) - &
                  this%periodic_shift(:, j)
             coord = coord - this%periodic_len(j)
          end do
       end do
    end do
  end subroutine wrap_particles_periodic

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

  !> Advance particles with the fluid scheme's explicit multistep coefficients
  !! and optionally write them.
  subroutine lpt_compute(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), allocatable :: vel(:,:)
    integer :: i
    integer :: j
    integer :: nadv

    if (time%t .lt. this%start_time) return

    call this%redistribute_particles()
    call this%evaluate_velocity(vel)
    nadv = this%case%fluid%ext_bdf%nadv

    if (nadv .gt. 1 .and. all(this%vel_particles_lag .eq. 0.0_rp)) then
       do j = 1, LPT_VEL_HISTORY_LEN
          this%vel_particles_lag(:, j, :) = vel
       end do
    end if

    if (this%output_enabled) then
      if (this%output_controller%check(time)) then
         call this%write_output(time, vel)
         call this%output_controller%register_execution()
      end if
    end if

    do i = 1, this%n_particles
       this%xyz_particles(:, i) = this%xyz_particles(:, i) + time%dt * &
            this%case%fluid%ext_bdf%advection_coeffs%x(1) * vel(:, i)
       do j = 2, nadv
          this%xyz_particles(:, i) = this%xyz_particles(:, i) + time%dt * &
               this%case%fluid%ext_bdf%advection_coeffs%x(j) * &
               this%vel_particles_lag(:, j - 1, i)
       end do
    end do

    do j = LPT_VEL_HISTORY_LEN, 2, -1
       this%vel_particles_lag(:, j, :) = this%vel_particles_lag(:, j - 1, :)
    end do
    this%vel_particles_lag(:, 1, :) = vel

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

    n_local = this%n_particles

    allocate(local_data(9, n_local))
    do i = 1, n_local
       local_data(1,i) = real(time%tstep, rp)
       local_data(2,i) = time%t
       local_data(3,i) = real(this%particle_ids(i), rp)
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

    if (allocated(this%xyz_particles)) deallocate(this%xyz_particles)
    if (allocated(this%particle_ids)) deallocate(this%particle_ids)
    if (allocated(this%vel_particles_lag)) deallocate(this%vel_particles_lag)
    if (allocated(this%periodic_map_el)) deallocate(this%periodic_map_el)
    if (allocated(this%periodic_map_facet)) deallocate(this%periodic_map_facet)
    if (allocated(this%periodic_map_npts)) deallocate(this%periodic_map_npts)
    if (allocated(this%periodic_src_center)) deallocate(this%periodic_src_center)
    if (allocated(this%periodic_src_basis)) deallocate(this%periodic_src_basis)
    if (allocated(this%periodic_tgt_center)) deallocate(this%periodic_tgt_center)
    if (allocated(this%periodic_tgt_basis)) deallocate(this%periodic_tgt_basis)
    call this%global_interp%free()
    call this%output_file%free()

    this%u => null()
    this%v => null()
    this%w => null()
    this%output_enabled = .false.
    this%log = .true.
    this%periodic_enabled = .false.
    this%start_time = -huge(0.0_rp)
    this%n_particles = 0
    this%n_global_particles = 0
    this%n_periodic_dirs = 0
    this%periodic_shift = 0.0_rp
    this%periodic_dir = 0.0_rp
    this%periodic_min = 0.0_rp
    this%periodic_max = 0.0_rp
    this%periodic_len = 0.0_rp
    this%n_periodic_maps = 0
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
    if (this%periodic_enabled) then
       write(log_buf, '(A,I0)') "Periodic wrap directions: ", &
            this%n_periodic_dirs
       call neko_log%message(log_buf)
    end if
    if (this%n_periodic_maps .gt. 0) then
       write(log_buf, '(A,I0)') "Periodic facet transforms: ", &
            this%n_periodic_maps
       call neko_log%message(log_buf)
    end if
    write(log_buf, '(A,I0)') "Local particles on rank 0 at init: ", &
         this%n_particles
    if (pe_rank .eq. 0) call neko_log%message(log_buf)
    call neko_log%message("Input supported from doc/pages/user-guide/" // &
         "case-file.md semantics for simcomp JSON; particle seeding here " // &
         "uses coordinates or points_file.")
    call neko_log%end_section()
  end subroutine log_status

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
    integer, allocatable :: displs(:)
    integer, allocatable :: local_meta(:)
    integer, allocatable :: global_meta(:)
    integer, allocatable :: padded_meta(:)
    integer, allocatable :: gathered_meta(:)
    real(kind=rp), allocatable :: local_geom(:)
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
    if (allocated(this%periodic_map_el)) deallocate(this%periodic_map_el)
    if (allocated(this%periodic_map_facet)) deallocate(this%periodic_map_facet)
    if (allocated(this%periodic_map_npts)) deallocate(this%periodic_map_npts)
    if (allocated(this%periodic_src_center)) deallocate(this%periodic_src_center)
    if (allocated(this%periodic_src_basis)) deallocate(this%periodic_src_basis)
    if (allocated(this%periodic_tgt_center)) deallocate(this%periodic_tgt_center)
    if (allocated(this%periodic_tgt_basis)) deallocate(this%periodic_tgt_basis)

    n_local = msh%periodic%size
    if (n_local .eq. 0) return

    allocate(counts(pe_size))
    allocate(displs(pe_size))
    call MPI_Allgather(n_local, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, &
         NEKO_COMM, ierr)
    displs(1) = 0
    do i = 2, pe_size
       displs(i) = displs(i - 1) + counts(i - 1)
    end do
    n_global = sum(counts)
    max_local = max(1, maxval(counts))

    allocate(local_meta(10 * n_local))
    allocate(local_geom(15 * n_local))
    do i = 1, n_local
       npts = merge(4, 2, msh%gdim .eq. 3)
       local_meta(10 * (i - 1) + 1) = msh%periodic%facet_el(i)%x(2)
       local_meta(10 * (i - 1) + 2) = msh%periodic%facet_el(i)%x(1)
       local_meta(10 * (i - 1) + 3) = npts
       local_meta(10 * (i - 1) + 4:10 * (i - 1) + 7) = msh%periodic%org_ids(i)%x
       local_meta(10 * (i - 1) + 8:10 * (i - 1) + 10) = 0

       call lpt_get_facet_points(msh, msh%periodic%facet_el(i)%x(2), &
            msh%periodic%facet_el(i)%x(1), src_pts, src_centroid)
       local_geom(15 * (i - 1) + 1:15 * (i - 1) + 12) = reshape(src_pts, [12])
       local_geom(15 * (i - 1) + 13:15 * (i - 1) + 15) = src_centroid
    end do

    allocate(padded_meta(10 * max_local))
    allocate(padded_geom(15 * max_local))
    padded_meta = 0
    padded_geom = 0.0_rp
    if (n_local .gt. 0) then
       padded_meta(1:10 * n_local) = local_meta
       padded_geom(1:15 * n_local) = local_geom
    end if

    allocate(gathered_meta(10 * max_local * pe_size))
    allocate(gathered_geom(15 * max_local * pe_size))
    call MPI_Allgather(padded_meta, 10 * max_local, MPI_INTEGER, gathered_meta, &
         10 * max_local, MPI_INTEGER, NEKO_COMM, ierr)
    call MPI_Allgather(padded_geom, 15 * max_local, MPI_REAL_PRECISION, &
         gathered_geom, 15 * max_local, MPI_REAL_PRECISION, NEKO_COMM, ierr)

    allocate(global_meta(10 * n_global))
    allocate(global_geom(15 * n_global))
    idx = 0
    do i = 1, pe_size
       if (counts(i) .gt. 0) then
          global_meta(10 * idx + 1:10 * (idx + counts(i))) = &
               gathered_meta(10 * max_local * (i - 1) + 1: &
               10 * max_local * (i - 1) + 10 * counts(i))
          global_geom(15 * idx + 1:15 * (idx + counts(i))) = &
               gathered_geom(15 * max_local * (i - 1) + 1: &
               15 * max_local * (i - 1) + 15 * counts(i))
          idx = idx + counts(i)
       end if
    end do

    allocate(this%periodic_map_el(n_local))
    allocate(this%periodic_map_facet(n_local))
    allocate(this%periodic_map_npts(n_local))
    allocate(this%periodic_src_center(3, n_local))
    allocate(this%periodic_src_basis(3, 3, n_local))
    allocate(this%periodic_tgt_center(3, n_local))
    allocate(this%periodic_tgt_basis(3, 3, n_local))

    do i = 1, n_local
       match_idx = 0
       do j = 1, n_global
          if (lpt_same_point_ids(global_meta(10 * (j - 1) + 4:10 * (j - 1) + 7), &
               msh%periodic%p_ids(i)%x, local_meta(10 * (i - 1) + 3))) then
             match_idx = j
             exit
          end if
       end do
       if (match_idx .eq. 0) cycle

       npts = local_meta(10 * (i - 1) + 3)
       src_pts = reshape(local_geom(15 * (i - 1) + 1:15 * (i - 1) + 12), [3, 4])
       src_centroid = local_geom(15 * (i - 1) + 13:15 * (i - 1) + 15)
       call lpt_reorder_facet_points( &
            global_meta(10 * (match_idx - 1) + 4:10 * (match_idx - 1) + 7), &
            reshape(global_geom(15 * (match_idx - 1) + 1:15 * (match_idx - 1) + 12), [3, 4]), &
            msh%periodic%p_ids(i)%x, npts, tgt_pts)
       tgt_centroid = global_geom(15 * (match_idx - 1) + 13:15 * (match_idx - 1) + 15)

       idx = this%n_periodic_maps + 1
       this%n_periodic_maps = idx
       this%periodic_map_el(idx) = msh%periodic%facet_el(i)%x(2)
       this%periodic_map_facet(idx) = msh%periodic%facet_el(i)%x(1)
       this%periodic_map_npts(idx) = npts
       call lpt_build_facet_basis(src_pts, src_centroid, npts, &
            this%periodic_src_center(:, idx), this%periodic_src_basis(:, :, idx))
       call lpt_build_facet_basis(tgt_pts, tgt_centroid, npts, &
            this%periodic_tgt_center(:, idx), this%periodic_tgt_basis(:, :, idx))
    end do

    deallocate(local_meta)
    deallocate(global_meta)
    deallocate(local_geom)
    deallocate(global_geom)
    deallocate(padded_meta)
    deallocate(gathered_meta)
    deallocate(padded_geom)
    deallocate(gathered_geom)
    deallocate(counts)
    deallocate(displs)
  end subroutine init_periodic_maps

  subroutine lpt_apply_periodic_map_if_needed(this, i_particle, i_map)
    class(lpt_t), intent(inout) :: this
    integer, intent(in) :: i_particle
    integer, intent(in) :: i_map
    real(kind=rp) :: rel(3)
    real(kind=rp) :: xi(3)

    rel = this%xyz_particles(:, i_particle) - this%periodic_src_center(:, i_map)
    if (this%periodic_map_npts(i_map) .eq. 4) then
       xi(1) = dot_product(rel, this%periodic_src_basis(:, 1, i_map))
       xi(2) = dot_product(rel, this%periodic_src_basis(:, 2, i_map))
       xi(3) = dot_product(rel, this%periodic_src_basis(:, 3, i_map))
       if (xi(3) .lt. -LPT_PERIODIC_TOL) then
          this%xyz_particles(:, i_particle) = this%periodic_tgt_center(:, i_map) + &
               xi(1) * this%periodic_tgt_basis(:, 1, i_map) + &
               xi(2) * this%periodic_tgt_basis(:, 2, i_map) + &
               xi(3) * this%periodic_tgt_basis(:, 3, i_map)
       end if
    else
       xi(1) = dot_product(rel, this%periodic_src_basis(:, 1, i_map))
       xi(2) = dot_product(rel, this%periodic_src_basis(:, 3, i_map))
       if (xi(2) .lt. -LPT_PERIODIC_TOL) then
          this%xyz_particles(:, i_particle) = this%periodic_tgt_center(:, i_map) + &
               xi(1) * this%periodic_tgt_basis(:, 1, i_map) + &
               xi(2) * this%periodic_tgt_basis(:, 3, i_map)
       end if
    end if
  end subroutine lpt_apply_periodic_map_if_needed

  logical function lpt_same_point_ids(ids_a, ids_b, npts) result(match)
    integer, intent(in) :: ids_a(4)
    integer, intent(in) :: ids_b(4)
    integer, intent(in) :: npts
    integer :: i
    integer :: j

    match = .true.
    do i = 1, npts
       if (.not. any(ids_a(i) .eq. ids_b(1:npts))) then
          match = .false.
          return
       end if
    end do
  end function lpt_same_point_ids

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
