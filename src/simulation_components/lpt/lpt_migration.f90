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
!> Particle redistribution support for LPT.
module lpt_migrate
  use num_types, only : rp
  use global_interpolation, only : global_interpolation_t
  use glb_intrp_comm, only : glb_intrp_comm_t
  use lpt_periodic_bc, only : lpt_periodic_bc_t
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_rank, pe_size
  use mpi_f08, only : MPI_Allreduce, MPI_Bcast, MPI_INTEGER, MPI_Scatterv, &
       MPI_SUM
  implicit none
  private

  integer, public, parameter :: LPT_MIGRATE_TO_OWNER = 1
  integer, public, parameter :: LPT_MIGRATE_NONE = 2

  type, public :: lpt_migrate_t
     integer :: lag_len = 0
     integer :: strategy = LPT_MIGRATE_TO_OWNER
   contains
     procedure, pass(this) :: init => lpt_migrate_init
     procedure, pass(this) :: free => lpt_migrate_free
     procedure, pass(this) :: initialize_particle_distribution
     procedure, pass(this) :: migrate_particles
     procedure, private, pass(this) :: distribute_particles_evenly
     procedure, private, pass(this) :: distribute_particle_ids
     procedure, private, pass(this) :: distribute_particle_field
     procedure, private, pass(this) :: distribute_particle_scalar
     procedure, private, pass(this) :: distribute_lags
     procedure, private, pass(this) :: localize_global_interpolation
     procedure, private, pass(this) :: migrate_particle_positions
     procedure, private, pass(this) :: migrate_particle_ids
     procedure, private, pass(this) :: migrate_particle_field
     procedure, private, pass(this) :: migrate_particle_scalar
     procedure, private, pass(this) :: migrate_lags
  end type lpt_migrate_t

contains

  subroutine lpt_migrate_init(this, lag_len, strategy)
    class(lpt_migrate_t), intent(inout) :: this
    integer, intent(in) :: lag_len
    integer, intent(in), optional :: strategy

    this%lag_len = lag_len
    this%strategy = LPT_MIGRATE_TO_OWNER
    if (present(strategy)) this%strategy = strategy
  end subroutine lpt_migrate_init

  subroutine lpt_migrate_free(this)
    class(lpt_migrate_t), intent(inout) :: this

    this%lag_len = 0
    this%strategy = LPT_MIGRATE_TO_OWNER
  end subroutine lpt_migrate_free

  subroutine build_even_particle_distribution(n_root, counts, offsets, n_local)
    integer, intent(in) :: n_root
    integer, allocatable, intent(out) :: counts(:)
    integer, allocatable, intent(out) :: offsets(:)
    integer, intent(out) :: n_local
    integer :: n_total
    integer :: rank
    integer :: ierr

    n_total = n_root
    call MPI_Bcast(n_total, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)

    allocate(counts(0:pe_size - 1))
    allocate(offsets(0:pe_size - 1))
    do rank = 0, pe_size - 1
       counts(rank) = n_total / pe_size
       if (rank .lt. mod(n_total, pe_size)) counts(rank) = counts(rank) + 1
    end do

    offsets(0) = 0
    do rank = 1, pe_size - 1
       offsets(rank) = offsets(rank - 1) + counts(rank - 1)
    end do
    n_local = counts(pe_rank)
  end subroutine build_even_particle_distribution

  subroutine initialize_particle_distribution(this, inertia, xyz, ids, &
       vel_lag, vel, acc, d, rho, acc_lag, n)
    class(lpt_migrate_t), intent(inout) :: this
    logical, intent(in) :: inertia
    real(kind=rp), allocatable, intent(inout) :: xyz(:, :)
    integer, allocatable, intent(inout) :: ids(:)
    real(kind=rp), allocatable, intent(inout) :: vel_lag(:, :, :)
    real(kind=rp), allocatable, intent(inout) :: vel(:, :)
    real(kind=rp), allocatable, intent(inout) :: acc(:, :)
    real(kind=rp), allocatable, intent(inout) :: d(:)
    real(kind=rp), allocatable, intent(inout) :: rho(:)
    real(kind=rp), allocatable, intent(inout) :: acc_lag(:, :, :)
    integer, intent(inout) :: n

    if (this%strategy .eq. LPT_MIGRATE_NONE) then
       call this%distribute_particles_evenly(inertia, xyz, ids, vel_lag, &
            vel, acc, d, rho, acc_lag, n)
    end if
  end subroutine initialize_particle_distribution

  !> Update particle ownership according to the selected migration strategy.
  subroutine migrate_particles(this, global_interp, periodic_bc, inertia, &
       xyz, ids, vel_lag, vel, acc, d, rho, acc_lag, n, n_global)
    class(lpt_migrate_t), intent(inout) :: this
    type(global_interpolation_t), intent(inout) :: global_interp
    type(lpt_periodic_bc_t), intent(inout) :: periodic_bc
    logical, intent(in) :: inertia
    real(kind=rp), allocatable, intent(inout) :: xyz(:, :)
    integer, allocatable, intent(inout) :: ids(:)
    real(kind=rp), allocatable, intent(inout) :: vel_lag(:, :, :)
    real(kind=rp), allocatable, intent(inout) :: vel(:, :)
    real(kind=rp), allocatable, intent(inout) :: acc(:, :)
    real(kind=rp), allocatable, intent(inout) :: d(:)
    real(kind=rp), allocatable, intent(inout) :: rho(:)
    real(kind=rp), allocatable, intent(inout) :: acc_lag(:, :, :)
    integer, intent(inout) :: n
    integer, intent(out) :: n_global
    type(glb_intrp_comm_t) :: migrate_comm
    integer :: n_particles_old
    integer :: n_particles_local
    integer, allocatable :: particle_ids_local(:)
    real(kind=rp), allocatable :: xyz_local(:, :)
    real(kind=rp), allocatable :: vel_local(:, :)
    real(kind=rp), allocatable :: acc_local(:, :)
    real(kind=rp), allocatable :: d_local(:)
    real(kind=rp), allocatable :: rho_local(:)
    real(kind=rp), allocatable :: vel_particles_lag_local(:, :, :)
    real(kind=rp), allocatable :: acc_particles_lag_local(:, :, :)
    integer :: ierr
    integer :: rank
    logical :: migration_needed

    call periodic_bc%wrap(xyz, n, vel, vel_lag, acc_lag)
    if (this%strategy .eq. LPT_MIGRATE_NONE) then
       call global_interp%find_points(xyz, n)
       call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, &
            ierr)
       return
    end if

    n_particles_old = n
    call global_interp%find_points(xyz, n)
    n_particles_local = global_interp%n_points_local

    migration_needed = .false.
    do rank = 0, global_interp%pe_size - 1
       if (rank .ne. pe_rank) then
          migration_needed = migration_needed .or. &
               global_interp%n_points_pe(rank) .gt. 0
          migration_needed = migration_needed .or. &
               global_interp%n_points_pe_local(rank) .gt. 0
       end if
    end do
    if (.not. migration_needed .and. n_particles_local .eq. n) then
       call this%localize_global_interpolation(global_interp, n)
       call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, &
            ierr)
       return
    end if

    call global_interp%init_redist_comm(migrate_comm)
    call this%migrate_particle_positions(migrate_comm, xyz, &
         n_particles_old, n_particles_local, xyz_local)
    call this%migrate_particle_ids(migrate_comm, ids, n_particles_old, &
         n_particles_local, particle_ids_local)
    call this%migrate_lags(migrate_comm, vel_lag, n_particles_old, &
         n_particles_local, &
         vel_particles_lag_local)
    if (inertia) then
       call this%migrate_particle_field(migrate_comm, vel, &
            n_particles_old, n_particles_local, vel_local)
       call this%migrate_particle_field(migrate_comm, acc, &
            n_particles_old, n_particles_local, acc_local)
       call this%migrate_particle_scalar(migrate_comm, d, &
            n_particles_old, n_particles_local, d_local)
       call this%migrate_particle_scalar(migrate_comm, rho, &
            n_particles_old, n_particles_local, rho_local)
       call this%migrate_lags(migrate_comm, acc_lag, n_particles_old, &
            n_particles_local, &
            acc_particles_lag_local)
    end if
    call migrate_comm%free()

    n = n_particles_local
    call move_alloc(xyz_local, xyz)
    call move_alloc(particle_ids_local, ids)
    call move_alloc(vel_particles_lag_local, vel_lag)
    if (inertia) then
      call move_alloc(vel_local, vel)
      call move_alloc(acc_local, acc)
      call move_alloc(d_local, d)
      call move_alloc(rho_local, rho)
      call move_alloc(acc_particles_lag_local, acc_lag)
    else
      if (allocated(vel)) deallocate(vel)
      if (allocated(acc)) deallocate(acc)
      if (allocated(d)) deallocate(d)
      if (allocated(rho)) deallocate(rho)
      if (allocated(acc_lag)) deallocate(acc_lag)
    end if

    call this%localize_global_interpolation(global_interp, n)
    call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
  end subroutine migrate_particles

  subroutine distribute_particles_evenly(this, inertia, xyz, ids, vel_lag, &
       vel, acc, d, rho, acc_lag, n)
    class(lpt_migrate_t), intent(inout) :: this
    logical, intent(in) :: inertia
    real(kind=rp), allocatable, intent(inout) :: xyz(:, :)
    integer, allocatable, intent(inout) :: ids(:)
    real(kind=rp), allocatable, intent(inout) :: vel_lag(:, :, :)
    real(kind=rp), allocatable, intent(inout) :: vel(:, :)
    real(kind=rp), allocatable, intent(inout) :: acc(:, :)
    real(kind=rp), allocatable, intent(inout) :: d(:)
    real(kind=rp), allocatable, intent(inout) :: rho(:)
    real(kind=rp), allocatable, intent(inout) :: acc_lag(:, :, :)
    integer, intent(inout) :: n
    integer, allocatable :: ids_local(:)
    real(kind=rp), allocatable :: xyz_local(:, :)
    real(kind=rp), allocatable :: vel_lag_local(:, :, :)
    real(kind=rp), allocatable :: vel_local(:, :)
    real(kind=rp), allocatable :: acc_local(:, :)
    real(kind=rp), allocatable :: d_local(:)
    real(kind=rp), allocatable :: rho_local(:)
    real(kind=rp), allocatable :: acc_lag_local(:, :, :)
    integer, allocatable :: counts(:)
    integer, allocatable :: offsets(:)
    integer :: n_local

    call build_even_particle_distribution(n, counts, offsets, n_local)

    call this%distribute_particle_field(xyz, counts, offsets, n_local, &
         xyz_local)
    call this%distribute_particle_ids(ids, counts, offsets, n_local, &
         ids_local)
    call this%distribute_lags(vel_lag, counts, offsets, n_local, &
         vel_lag_local)
    if (inertia) then
       call this%distribute_particle_field(vel, counts, offsets, n_local, &
            vel_local)
       call this%distribute_particle_field(acc, counts, offsets, n_local, &
            acc_local)
       call this%distribute_particle_scalar(d, counts, offsets, n_local, &
            d_local)
       call this%distribute_particle_scalar(rho, counts, offsets, n_local, &
            rho_local)
       call this%distribute_lags(acc_lag, counts, offsets, n_local, &
            acc_lag_local)
    end if

    n = n_local
    call move_alloc(xyz_local, xyz)
    call move_alloc(ids_local, ids)
    call move_alloc(vel_lag_local, vel_lag)
    if (inertia) then
       call move_alloc(vel_local, vel)
       call move_alloc(acc_local, acc)
       call move_alloc(d_local, d)
       call move_alloc(rho_local, rho)
       call move_alloc(acc_lag_local, acc_lag)
    else
       if (allocated(vel)) deallocate(vel)
       if (allocated(acc)) deallocate(acc)
       if (allocated(d)) deallocate(d)
       if (allocated(rho)) deallocate(rho)
       if (allocated(acc_lag)) deallocate(acc_lag)
    end if

    deallocate(counts)
    deallocate(offsets)
  end subroutine distribute_particles_evenly

  subroutine distribute_particle_ids(this, ids_old, counts, offsets, n_local, &
       ids_local)
    class(lpt_migrate_t), intent(inout) :: this
    integer, allocatable, intent(in) :: ids_old(:)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: offsets(0:)
    integer, intent(in) :: n_local
    integer, allocatable, intent(out) :: ids_local(:)
    integer, allocatable :: counts_i(:)
    integer, allocatable :: offsets_i(:)
    integer :: ierr

    allocate(ids_local(n_local))
    allocate(counts_i(0:pe_size - 1))
    allocate(offsets_i(0:pe_size - 1))
    counts_i = counts
    offsets_i = offsets
    call MPI_Scatterv(ids_old, counts_i, offsets_i, MPI_INTEGER, ids_local, &
         n_local, MPI_INTEGER, 0, NEKO_COMM, ierr)
    deallocate(counts_i)
    deallocate(offsets_i)
  end subroutine distribute_particle_ids

  subroutine distribute_particle_field(this, field_old, counts, offsets, &
       n_local, field_local)
    class(lpt_migrate_t), intent(inout) :: this
    real(kind=rp), allocatable, intent(in) :: field_old(:, :)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: offsets(0:)
    integer, intent(in) :: n_local
    real(kind=rp), allocatable, intent(out) :: field_local(:, :)
    integer, allocatable :: counts_r(:)
    integer, allocatable :: offsets_r(:)
    integer :: rank
    integer :: ierr

    allocate(field_local(3, n_local))
    allocate(counts_r(0:pe_size - 1))
    allocate(offsets_r(0:pe_size - 1))
    do rank = 0, pe_size - 1
       counts_r(rank) = 3 * counts(rank)
       offsets_r(rank) = 3 * offsets(rank)
    end do
    call MPI_Scatterv(field_old, counts_r, offsets_r, MPI_REAL_PRECISION, &
         field_local, 3 * n_local, MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)
    deallocate(counts_r)
    deallocate(offsets_r)
  end subroutine distribute_particle_field

  subroutine distribute_particle_scalar(this, scalar_old, counts, offsets, &
       n_local, scalar_local)
    class(lpt_migrate_t), intent(inout) :: this
    real(kind=rp), allocatable, intent(in) :: scalar_old(:)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: offsets(0:)
    integer, intent(in) :: n_local
    real(kind=rp), allocatable, intent(out) :: scalar_local(:)
    integer, allocatable :: counts_r(:)
    integer, allocatable :: offsets_r(:)
    integer :: ierr

    allocate(scalar_local(n_local))
    allocate(counts_r(0:pe_size - 1))
    allocate(offsets_r(0:pe_size - 1))
    counts_r = counts
    offsets_r = offsets
    call MPI_Scatterv(scalar_old, counts_r, offsets_r, MPI_REAL_PRECISION, &
         scalar_local, n_local, MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)
    deallocate(counts_r)
    deallocate(offsets_r)
  end subroutine distribute_particle_scalar

  subroutine distribute_lags(this, lag_old, counts, offsets, n_local, &
       lag_local)
    class(lpt_migrate_t), intent(inout) :: this
    real(kind=rp), allocatable, intent(in) :: lag_old(:, :, :)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: offsets(0:)
    integer, intent(in) :: n_local
    real(kind=rp), allocatable, intent(out) :: lag_local(:, :, :)
    integer, allocatable :: counts_r(:)
    integer, allocatable :: offsets_r(:)
    integer :: rank
    integer :: ierr

    allocate(lag_local(3, this%lag_len, n_local))
    allocate(counts_r(0:pe_size - 1))
    allocate(offsets_r(0:pe_size - 1))
    do rank = 0, pe_size - 1
       counts_r(rank) = 3 * this%lag_len * counts(rank)
       offsets_r(rank) = 3 * this%lag_len * offsets(rank)
    end do
    call MPI_Scatterv(lag_old, counts_r, offsets_r, MPI_REAL_PRECISION, &
         lag_local, 3 * this%lag_len * n_local, MPI_REAL_PRECISION, 0, &
         NEKO_COMM, ierr)
    deallocate(counts_r)
    deallocate(offsets_r)
  end subroutine distribute_lags

  subroutine localize_global_interpolation(this, global_interp, n_local)
    class(lpt_migrate_t), intent(inout) :: this
    type(global_interpolation_t), intent(inout) :: global_interp
    integer, intent(in) :: n_local

    global_interp%n_points = n_local
    global_interp%all_points_local = .true.
  end subroutine localize_global_interpolation

  subroutine migrate_particle_positions(this, migrate_comm, xyz_old, &
       n_particles_old, n_local, xyz_local)
    class(lpt_migrate_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: migrate_comm
    real(kind=rp), allocatable, intent(in) :: xyz_old(:, :)
    integer, intent(in) :: n_particles_old
    integer, intent(in) :: n_local
    real(kind=rp), allocatable, intent(out) :: xyz_local(:, :)

    call this%migrate_particle_field(migrate_comm, xyz_old, &
         n_particles_old, n_local, xyz_local)
  end subroutine migrate_particle_positions

  subroutine migrate_particle_ids(this, migrate_comm, ids_old, &
       n_particles_old, n_local, particle_ids_local)
    class(lpt_migrate_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: migrate_comm
    integer, allocatable, intent(in) :: ids_old(:)
    integer, intent(in) :: n_particles_old
    integer, intent(in) :: n_local
    integer, allocatable, intent(out) :: particle_ids_local(:)
    real(kind=rp), allocatable :: sendbuf(:)
    real(kind=rp), allocatable :: recvbuf(:)

    allocate(particle_ids_local(n_local))
    particle_ids_local = 0

    if (n_particles_old .eq. 0 .and. n_local .eq. 0) return

    allocate(sendbuf(n_particles_old))
    allocate(recvbuf(n_local))
    sendbuf = real(ids_old, rp)
    recvbuf = 0.0_rp
    call migrate_comm%sendrecv(sendbuf, recvbuf, n_particles_old, n_local)
    particle_ids_local = nint(recvbuf)
    deallocate(sendbuf)
    deallocate(recvbuf)
  end subroutine migrate_particle_ids

  subroutine migrate_lags(this, migrate_comm, lag_old, n_particles_old, &
       n_local, lag_local)
    class(lpt_migrate_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: migrate_comm
    real(kind=rp), allocatable, intent(in) :: lag_old(:, :, :)
    integer, intent(in) :: n_particles_old
    integer, intent(in) :: n_local
    real(kind=rp), allocatable, intent(out) :: lag_local(:, :, :)
    real(kind=rp), allocatable :: sendbuf(:)
    real(kind=rp), allocatable :: recvbuf(:)
    integer :: i
    integer :: j

    allocate(lag_local(3, this%lag_len, n_local))
    lag_local = 0.0_rp

    if (n_particles_old .eq. 0 .and. n_local .eq. 0) return

    allocate(sendbuf(n_particles_old))
    allocate(recvbuf(n_local))
    do j = 1, this%lag_len
       do i = 1, 3
          sendbuf = lag_old(i, j, :)
          recvbuf = 0.0_rp
          call migrate_comm%sendrecv(sendbuf, recvbuf, n_particles_old, n_local)
          lag_local(i, j, :) = recvbuf
       end do
    end do
    deallocate(sendbuf)
    deallocate(recvbuf)
  end subroutine migrate_lags

  subroutine migrate_particle_field(this, migrate_comm, field_old, &
       n_particles_old, n_local, field_local)
    class(lpt_migrate_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: migrate_comm
    real(kind=rp), allocatable, intent(in) :: field_old(:, :)
    integer, intent(in) :: n_particles_old
    integer, intent(in) :: n_local
    real(kind=rp), allocatable, intent(out) :: field_local(:, :)
    real(kind=rp), allocatable :: sendbuf(:)
    real(kind=rp), allocatable :: recvbuf(:)
    integer :: i

    allocate(field_local(3, n_local))
    field_local = 0.0_rp

    if (n_particles_old .eq. 0 .and. n_local .eq. 0) return

    allocate(sendbuf(n_particles_old))
    allocate(recvbuf(n_local))
    do i = 1, 3
       sendbuf = field_old(i, :)
       recvbuf = 0.0_rp
       call migrate_comm%sendrecv(sendbuf, recvbuf, n_particles_old, n_local)
       field_local(i, :) = recvbuf
    end do
    deallocate(sendbuf)
    deallocate(recvbuf)
  end subroutine migrate_particle_field

  subroutine migrate_particle_scalar(this, migrate_comm, scalar_old, &
       n_particles_old, n_local, scalar_local)
    class(lpt_migrate_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: migrate_comm
    real(kind=rp), allocatable, intent(in) :: scalar_old(:)
    integer, intent(in) :: n_particles_old
    integer, intent(in) :: n_local
    real(kind=rp), allocatable, intent(out) :: scalar_local(:)
    real(kind=rp), allocatable :: sendbuf(:)
    real(kind=rp), allocatable :: recvbuf(:)

    allocate(scalar_local(n_local))
    scalar_local = 0.0_rp

    if (n_particles_old .eq. 0 .and. n_local .eq. 0) return

    allocate(sendbuf(n_particles_old))
    allocate(recvbuf(n_local))
    sendbuf = scalar_old
    recvbuf = 0.0_rp
    call migrate_comm%sendrecv(sendbuf, recvbuf, n_particles_old, n_local)
    scalar_local = recvbuf
    deallocate(sendbuf)
    deallocate(recvbuf)
  end subroutine migrate_particle_scalar

end module lpt_migrate
