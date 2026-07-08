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
  use vector, only : vector_t
  use device, only : HOST_TO_DEVICE, DEVICE_TO_HOST
  use scratch_registry, only : neko_scratch_registry
  use particles, only : particles_t
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
     procedure, private, pass(this) :: distribute_particle_scalar
     procedure, private, pass(this) :: localize_global_interpolation
     procedure, private, pass(this) :: migrate_particle_ids
     procedure, private, pass(this) :: migrate_particle_scalar
  end type lpt_migrate_t

contains

  !> Initialise particle migration settings.
  !! @param lag_len Number of lagged RHS levels carried by particles.
  !! @param strategy Optional migration strategy identifier.
  subroutine lpt_migrate_init(this, lag_len, strategy)
    class(lpt_migrate_t), intent(inout) :: this
    integer, intent(in) :: lag_len
    integer, intent(in), optional :: strategy

    this%lag_len = lag_len
    this%strategy = LPT_MIGRATE_TO_OWNER
    if (present(strategy)) this%strategy = strategy
  end subroutine lpt_migrate_init

  !> Reset migration settings to defaults.
  subroutine lpt_migrate_free(this)
    class(lpt_migrate_t), intent(inout) :: this

    this%lag_len = 0
    this%strategy = LPT_MIGRATE_TO_OWNER
  end subroutine lpt_migrate_free

  !> Build an even rank distribution from the root-rank particle count.
  !! @param n_root Number of particles available on rank 0.
  !! @param counts Number of particles assigned to each rank.
  !! @param offsets Particle offsets into the root-rank arrays.
  !! @param n_local Number of particles assigned to this rank.
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

  !> Apply the initial particle distribution required by the strategy.
  !! @param inertia Whether particles carry inertial state.
  !! @param particles Particle storage to distribute.
  subroutine initialize_particle_distribution(this, inertia, particles)
    class(lpt_migrate_t), intent(inout) :: this
    logical, intent(in) :: inertia
    type(particles_t), intent(inout) :: particles

    if (this%strategy .eq. LPT_MIGRATE_NONE) then
       call this%distribute_particles_evenly(inertia, particles)
    end if
  end subroutine initialize_particle_distribution

  !> Update particle ownership according to the selected migration strategy.
  !! @param global_interp Interpolation object used to locate particle owners.
  !! @param periodic_bc Periodic wrapper applied before ownership lookup.
  !! @param inertia Whether particles carry inertial state.
  !! @param particles Particle storage to migrate.
  subroutine migrate_particles(this, global_interp, periodic_bc, inertia, &
       particles)
    class(lpt_migrate_t), intent(inout) :: this
    type(global_interpolation_t), intent(inout) :: global_interp
    type(lpt_periodic_bc_t), intent(inout) :: periodic_bc
    logical, intent(in) :: inertia
    type(particles_t), intent(inout) :: particles

    type(glb_intrp_comm_t) :: migrate_comm
    integer :: n_particles_old
    integer :: n_particles_local
    integer, allocatable :: particle_ids_local(:)
    type(vector_t), pointer :: x_local, y_local, z_local
    type(vector_t), pointer :: u_local, v_local, w_local
    type(vector_t), pointer :: acc_xlocal, acc_ylocal, acc_zlocal
    type(vector_t), pointer :: d_local, rho_local
    type(vector_t), pointer :: u_lag_local, v_lag_local, w_lag_local
    type(vector_t), pointer :: u_laglag_local, v_laglag_local, w_laglag_local
    type(vector_t), pointer :: acc_xlag_local, acc_ylag_local, acc_zlag_local
    type(vector_t), pointer :: acc_xlaglag_local, &
         acc_ylaglag_local, acc_zlaglag_local
    integer :: ierr
    integer :: rank
    logical :: migration_needed
    integer :: ind_basic(9)
    integer :: ind_inertia(14)

    associate( &
         x => particles%x, y => particles%y, z => particles%z, &
         u => particles%u, v => particles%v, w => particles%w, &
         acc_x => particles%acc_x, acc_y => particles%acc_y, &
         acc_z => particles%acc_z, d => particles%d, rho => particles%rho, &
         u_lag => particles%u_lag, v_lag => particles%v_lag, &
         w_lag => particles%w_lag, u_laglag => particles%u_laglag, &
         v_laglag => particles%v_laglag, w_laglag => particles%w_laglag, &
         acc_xlag => particles%acc_xlag, acc_ylag => particles%acc_ylag, &
         acc_zlag => particles%acc_zlag, acc_xlaglag => particles%acc_xlaglag, &
         acc_ylaglag => particles%acc_ylaglag, &
         acc_zlaglag => particles%acc_zlaglag, n => particles%n, &
         n_global => particles%n_global)

      if (inertia) then
         call periodic_bc%wrap(x, y, z, n, u, v, w, u_lag, &
              v_lag, w_lag, u_laglag, v_laglag, w_laglag, &
              acc_xlag, acc_ylag, acc_zlag, acc_xlaglag, &
              acc_ylaglag, acc_zlaglag)
      else
         call periodic_bc%wrap(x, y, z, n, u, v, w, u_lag, &
              v_lag, w_lag, u_laglag, v_laglag, w_laglag)
      end if
      if (this%strategy .eq. LPT_MIGRATE_NONE) then
         call global_interp%find_points(x%x, y%x, z%x, n)
         call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, &
              ierr)
         return
      end if

      n_particles_old = n
      call global_interp%find_points(x%x, y%x, z%x, n)
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

      allocate(particle_ids_local(n_particles_local))
      particle_ids_local = 0
      call neko_scratch_registry%request_vector(x_local, ind_basic(1), &
           n_particles_local, .false.)
      call neko_scratch_registry%request_vector(y_local, ind_basic(2), &
           n_particles_local, .false.)
      call neko_scratch_registry%request_vector(z_local, ind_basic(3), &
           n_particles_local, .false.)
      call neko_scratch_registry%request_vector(u_lag_local, ind_basic(4), &
           n_particles_local, .false.)
      call neko_scratch_registry%request_vector(v_lag_local, ind_basic(5), &
           n_particles_local, .false.)
      call neko_scratch_registry%request_vector(w_lag_local, ind_basic(6), &
           n_particles_local, .false.)
      call neko_scratch_registry%request_vector(u_laglag_local, ind_basic(7), &
           n_particles_local, .false.)
      call neko_scratch_registry%request_vector(v_laglag_local, ind_basic(8), &
           n_particles_local, .false.)
      call neko_scratch_registry%request_vector(w_laglag_local, ind_basic(9), &
           n_particles_local, .false.)
      if (inertia) then
         call neko_scratch_registry%request_vector(u_local, ind_inertia(1), &
              n_particles_local, .false.)
         call neko_scratch_registry%request_vector(v_local, ind_inertia(2), &
              n_particles_local, .false.)
         call neko_scratch_registry%request_vector(w_local, ind_inertia(3), &
              n_particles_local, .false.)
         call neko_scratch_registry%request_vector(acc_xlocal, ind_inertia(4), &
              n_particles_local, .false.)
         call neko_scratch_registry%request_vector(acc_ylocal, ind_inertia(5), &
              n_particles_local, .false.)
         call neko_scratch_registry%request_vector(acc_zlocal, ind_inertia(6), &
              n_particles_local, .false.)
         call neko_scratch_registry%request_vector(d_local, ind_inertia(7), &
              n_particles_local, .false.)
         call neko_scratch_registry%request_vector(rho_local, ind_inertia(8), &
              n_particles_local, .false.)
         call neko_scratch_registry%request_vector(acc_xlag_local, &
              ind_inertia(9), n_particles_local, .false.)
         call neko_scratch_registry%request_vector(acc_ylag_local, &
              ind_inertia(10), n_particles_local, .false.)
         call neko_scratch_registry%request_vector(acc_zlag_local, &
              ind_inertia(11), n_particles_local, .false.)
         call neko_scratch_registry%request_vector(acc_xlaglag_local, &
              ind_inertia(12), n_particles_local, .false.)
         call neko_scratch_registry%request_vector(acc_ylaglag_local, &
              ind_inertia(13), n_particles_local, .false.)
         call neko_scratch_registry%request_vector(acc_zlaglag_local, &
              ind_inertia(14), n_particles_local, .false.)
      end if

      call global_interp%init_redist_comm(migrate_comm)
      call this%migrate_particle_ids(migrate_comm, particles%ids, &
           n_particles_old, n_particles_local, particle_ids_local)
      call this%migrate_particle_scalar(migrate_comm, x, n_particles_old, &
           n_particles_local, x_local)
      call this%migrate_particle_scalar(migrate_comm, y, n_particles_old, &
           n_particles_local, y_local)
      call this%migrate_particle_scalar(migrate_comm, z, n_particles_old, &
           n_particles_local, z_local)
      call this%migrate_particle_scalar(migrate_comm, u_lag, n_particles_old, &
           n_particles_local, u_lag_local)
      call this%migrate_particle_scalar(migrate_comm, v_lag, n_particles_old, &
           n_particles_local, v_lag_local)
      call this%migrate_particle_scalar(migrate_comm, w_lag, n_particles_old, &
           n_particles_local, w_lag_local)
      call this%migrate_particle_scalar(migrate_comm, u_laglag, n_particles_old, &
           n_particles_local, u_laglag_local)
      call this%migrate_particle_scalar(migrate_comm, v_laglag, n_particles_old, &
           n_particles_local, v_laglag_local)
      call this%migrate_particle_scalar(migrate_comm, w_laglag, n_particles_old, &
           n_particles_local, w_laglag_local)
      if (inertia) then
         call this%migrate_particle_scalar(migrate_comm, u, n_particles_old, &
              n_particles_local, u_local)
         call this%migrate_particle_scalar(migrate_comm, v, n_particles_old, &
              n_particles_local, v_local)
         call this%migrate_particle_scalar(migrate_comm, w, n_particles_old, &
              n_particles_local, w_local)
         call this%migrate_particle_scalar(migrate_comm, acc_x, n_particles_old, &
              n_particles_local, acc_xlocal)
         call this%migrate_particle_scalar(migrate_comm, acc_y, n_particles_old, &
              n_particles_local, acc_ylocal)
         call this%migrate_particle_scalar(migrate_comm, acc_z, n_particles_old, &
              n_particles_local, acc_zlocal)
         call this%migrate_particle_scalar(migrate_comm, d, &
              n_particles_old, n_particles_local, d_local)
         call this%migrate_particle_scalar(migrate_comm, rho, &
              n_particles_old, n_particles_local, rho_local)
         call this%migrate_particle_scalar(migrate_comm, acc_xlag, &
              n_particles_old, n_particles_local, acc_xlag_local)
         call this%migrate_particle_scalar(migrate_comm, acc_ylag, &
              n_particles_old, n_particles_local, acc_ylag_local)
         call this%migrate_particle_scalar(migrate_comm, acc_zlag, &
              n_particles_old, n_particles_local, acc_zlag_local)
         call this%migrate_particle_scalar(migrate_comm, acc_xlaglag, &
              n_particles_old, n_particles_local, acc_xlaglag_local)
         call this%migrate_particle_scalar(migrate_comm, acc_ylaglag, &
              n_particles_old, n_particles_local, acc_ylaglag_local)
         call this%migrate_particle_scalar(migrate_comm, acc_zlaglag, &
              n_particles_old, n_particles_local, acc_zlaglag_local)
      end if
      call migrate_comm%free()

      n = n_particles_local
      x = x_local
      y = y_local
      z = z_local
      call move_alloc(particle_ids_local, particles%ids)
      u_lag = u_lag_local
      v_lag = v_lag_local
      w_lag = w_lag_local
      u_laglag = u_laglag_local
      v_laglag = v_laglag_local
      w_laglag = w_laglag_local
      if (inertia) then
         u = u_local
         v = v_local
         w = w_local
         acc_x = acc_xlocal
         acc_y = acc_ylocal
         acc_z = acc_zlocal
         d = d_local
         rho = rho_local
         acc_xlag = acc_xlag_local
         acc_ylag = acc_ylag_local
         acc_zlag = acc_zlag_local
         acc_xlaglag = acc_xlaglag_local
         acc_ylaglag = acc_ylaglag_local
         acc_zlaglag = acc_zlaglag_local
      end if

      if (allocated(particle_ids_local)) deallocate(particle_ids_local)
      call neko_scratch_registry%relinquish(ind_basic)
      if (inertia) call neko_scratch_registry%relinquish(ind_inertia)

      call this%localize_global_interpolation(global_interp, n)
      call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    end associate

  end subroutine migrate_particles

  !> Distribute particles from rank 0 evenly across all ranks.
  !! @param inertia Whether particles carry inertial state.
  !! @param particles Particle storage to redistribute.
  subroutine distribute_particles_evenly(this, inertia, particles)
    class(lpt_migrate_t), intent(inout) :: this
    logical, intent(in) :: inertia
    type(particles_t), intent(inout) :: particles

    integer, allocatable :: ids_local(:)
    type(vector_t), pointer :: x_local, y_local, z_local
    type(vector_t), pointer :: u_lag_local, v_lag_local, w_lag_local
    type(vector_t), pointer :: u_laglag_local, v_laglag_local, w_laglag_local
    type(vector_t), pointer :: u_local, v_local, w_local
    type(vector_t), pointer :: acc_xlocal, acc_ylocal, acc_zlocal
    type(vector_t), pointer :: d_local, rho_local
    type(vector_t), pointer :: acc_xlag_local, acc_ylag_local, acc_zlag_local
    type(vector_t), pointer :: acc_xlaglag_local, &
         acc_ylaglag_local, acc_zlaglag_local
    integer, allocatable :: counts(:)
    integer, allocatable :: offsets(:)
    integer :: n_local
    integer :: ind_basic(9), ind_inertia(14)

    associate( &
         x => particles%x, y => particles%y, z => particles%z, &
         u => particles%u, v => particles%v, w => particles%w, &
         acc_x => particles%acc_x, acc_y => particles%acc_y, &
         acc_z => particles%acc_z, d => particles%d, rho => particles%rho, &
         u_lag => particles%u_lag, v_lag => particles%v_lag, &
         w_lag => particles%w_lag, u_laglag => particles%u_laglag, &
         v_laglag => particles%v_laglag, w_laglag => particles%w_laglag, &
         acc_xlag => particles%acc_xlag, acc_ylag => particles%acc_ylag, &
         acc_zlag => particles%acc_zlag, acc_xlaglag => particles%acc_xlaglag, &
         acc_ylaglag => particles%acc_ylaglag, &
         acc_zlaglag => particles%acc_zlaglag, n => particles%n)

      call build_even_particle_distribution(n, counts, offsets, n_local)

      allocate(ids_local(n_local))
      call neko_scratch_registry%request_vector(x_local, ind_basic(1), &
           n_local, .false.)
      call neko_scratch_registry%request_vector(y_local, ind_basic(2), &
           n_local, .false.)
      call neko_scratch_registry%request_vector(z_local, ind_basic(3), &
           n_local, .false.)
      call neko_scratch_registry%request_vector(u_lag_local, ind_basic(4), &
           n_local, .false.)
      call neko_scratch_registry%request_vector(v_lag_local, ind_basic(5), &
           n_local, .false.)
      call neko_scratch_registry%request_vector(w_lag_local, ind_basic(6), &
           n_local, .false.)
      call neko_scratch_registry%request_vector(u_laglag_local, ind_basic(7), &
           n_local, .false.)
      call neko_scratch_registry%request_vector(v_laglag_local, ind_basic(8), &
           n_local, .false.)
      call neko_scratch_registry%request_vector(w_laglag_local, ind_basic(9), &
           n_local, .false.)
      if (inertia) then
         call neko_scratch_registry%request_vector(u_local, ind_inertia(1), &
              n_local, .false.)
         call neko_scratch_registry%request_vector(v_local, ind_inertia(2), &
              n_local, .false.)
         call neko_scratch_registry%request_vector(w_local, ind_inertia(3), &
              n_local, .false.)
         call neko_scratch_registry%request_vector(acc_xlocal, ind_inertia(4), &
              n_local, .false.)
         call neko_scratch_registry%request_vector(acc_ylocal, ind_inertia(5), &
              n_local, .false.)
         call neko_scratch_registry%request_vector(acc_zlocal, ind_inertia(6), &
              n_local, .false.)
         call neko_scratch_registry%request_vector(d_local, ind_inertia(7), &
              n_local, .false.)
         call neko_scratch_registry%request_vector(rho_local, ind_inertia(8), &
              n_local, .false.)
         call neko_scratch_registry%request_vector(acc_xlag_local, &
              ind_inertia(9), n_local, .false.)
         call neko_scratch_registry%request_vector(acc_ylag_local, &
              ind_inertia(10), n_local, .false.)
         call neko_scratch_registry%request_vector(acc_zlag_local, &
              ind_inertia(11), n_local, .false.)
         call neko_scratch_registry%request_vector(acc_xlaglag_local, &
              ind_inertia(12), n_local, .false.)
         call neko_scratch_registry%request_vector(acc_ylaglag_local, &
              ind_inertia(13), n_local, .false.)
         call neko_scratch_registry%request_vector(acc_zlaglag_local, &
              ind_inertia(14), n_local, .false.)
      end if

      call this%distribute_particle_ids(particles%ids, counts, offsets, n_local, &
           ids_local)
      call this%distribute_particle_scalar(x, counts, offsets, n_local, x_local)
      call this%distribute_particle_scalar(y, counts, offsets, n_local, y_local)
      call this%distribute_particle_scalar(z, counts, offsets, n_local, z_local)
      call this%distribute_particle_scalar(u_lag, counts, offsets, n_local, &
           u_lag_local)
      call this%distribute_particle_scalar(v_lag, counts, offsets, n_local, &
           v_lag_local)
      call this%distribute_particle_scalar(w_lag, counts, offsets, n_local, &
           w_lag_local)
      call this%distribute_particle_scalar(u_laglag, counts, offsets, n_local, &
           u_laglag_local)
      call this%distribute_particle_scalar(v_laglag, counts, offsets, n_local, &
           v_laglag_local)
      call this%distribute_particle_scalar(w_laglag, counts, offsets, n_local, &
           w_laglag_local)
      if (inertia) then
         call this%distribute_particle_scalar(u, counts, offsets, n_local, &
              u_local)
         call this%distribute_particle_scalar(v, counts, offsets, n_local, &
              v_local)
         call this%distribute_particle_scalar(w, counts, offsets, n_local, &
              w_local)
         call this%distribute_particle_scalar(acc_x, counts, offsets, n_local, &
              acc_xlocal)
         call this%distribute_particle_scalar(acc_y, counts, offsets, n_local, &
              acc_ylocal)
         call this%distribute_particle_scalar(acc_z, counts, offsets, n_local, &
              acc_zlocal)
         call this%distribute_particle_scalar(d, counts, offsets, n_local, &
              d_local)
         call this%distribute_particle_scalar(rho, counts, offsets, n_local, &
              rho_local)
         call this%distribute_particle_scalar(acc_xlag, counts, offsets, &
              n_local, acc_xlag_local)
         call this%distribute_particle_scalar(acc_ylag, counts, offsets, &
              n_local, acc_ylag_local)
         call this%distribute_particle_scalar(acc_zlag, counts, offsets, &
              n_local, acc_zlag_local)
         call this%distribute_particle_scalar(acc_xlaglag, counts, offsets, &
              n_local, acc_xlaglag_local)
         call this%distribute_particle_scalar(acc_ylaglag, counts, offsets, &
              n_local, acc_ylaglag_local)
         call this%distribute_particle_scalar(acc_zlaglag, counts, offsets, &
              n_local, acc_zlaglag_local)
      end if

      n = n_local
      x = x_local
      y = y_local
      z = z_local
      call move_alloc(ids_local, particles%ids)
      u_lag = u_lag_local
      v_lag = v_lag_local
      w_lag = w_lag_local
      u_laglag = u_laglag_local
      v_laglag = v_laglag_local
      w_laglag = w_laglag_local
      if (inertia) then
         u = u_local
         v = v_local
         w = w_local
         acc_x = acc_xlocal
         acc_y = acc_ylocal
         acc_z = acc_zlocal
         d = d_local
         rho = rho_local
         acc_xlag = acc_xlag_local
         acc_ylag = acc_ylag_local
         acc_zlag = acc_zlag_local
         acc_xlaglag = acc_xlaglag_local
         acc_ylaglag = acc_ylaglag_local
         acc_zlaglag = acc_zlaglag_local
      end if

      if (allocated(ids_local)) deallocate(ids_local)
      call neko_scratch_registry%relinquish(ind_basic)
      if (inertia) call neko_scratch_registry%relinquish(ind_inertia)

      deallocate(counts)
      deallocate(offsets)

    end associate

  end subroutine distribute_particles_evenly

  !> Scatter particle ids from rank 0 to the current rank.
  !! @param ids_old Root-rank particle ids.
  !! @param counts Number of ids assigned to each rank.
  !! @param offsets Offsets into `ids_old`.
  !! @param n_local Number of ids received by this rank.
  !! @param ids_local Local particle ids.
  subroutine distribute_particle_ids(this, ids_old, counts, offsets, n_local, &
       ids_local)
    class(lpt_migrate_t), intent(inout) :: this
    integer, allocatable, intent(in) :: ids_old(:)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: offsets(0:)
    integer, intent(in) :: n_local
    integer, allocatable, intent(inout) :: ids_local(:)
    integer, allocatable :: counts_i(:)
    integer, allocatable :: offsets_i(:)
    integer :: ierr

    allocate(counts_i(0:pe_size - 1))
    allocate(offsets_i(0:pe_size - 1))
    counts_i = counts
    offsets_i = offsets
    call MPI_Scatterv(ids_old, counts_i, offsets_i, MPI_INTEGER, ids_local, &
         n_local, MPI_INTEGER, 0, NEKO_COMM, ierr)
    deallocate(counts_i)
    deallocate(offsets_i)
  end subroutine distribute_particle_ids

  !> Scatter one particle scalar field from rank 0 to the current rank.
  !! @param scalar_old Root-rank scalar values.
  !! @param counts Number of values assigned to each rank.
  !! @param offsets Offsets into `scalar_old`.
  !! @param n_local Number of values received by this rank.
  !! @param scalar_local Local scalar values.
  subroutine distribute_particle_scalar(this, scalar_old, counts, offsets, &
       n_local, scalar_local)
    class(lpt_migrate_t), intent(inout) :: this
    type(vector_t), intent(inout) :: scalar_old
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: offsets(0:)
    integer, intent(in) :: n_local
    type(vector_t), intent(inout) :: scalar_local
    integer, allocatable :: counts_r(:)
    integer, allocatable :: offsets_r(:)
    integer :: ierr

    allocate(counts_r(0:pe_size - 1))
    allocate(offsets_r(0:pe_size - 1))
    counts_r = counts
    offsets_r = offsets
    call scalar_old%copy_from(DEVICE_TO_HOST, .true.)
    call MPI_Scatterv(scalar_old%x, counts_r, offsets_r, MPI_REAL_PRECISION, &
         scalar_local%x, n_local, MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)
    call scalar_local%copy_from(HOST_TO_DEVICE, .true.)
    deallocate(counts_r)
    deallocate(offsets_r)
  end subroutine distribute_particle_scalar

  !> Mark interpolation data as local after migration has completed.
  !! @param global_interp Interpolation object to update.
  !! @param n_local Number of local particles.
  subroutine localize_global_interpolation(this, global_interp, n_local)
    class(lpt_migrate_t), intent(inout) :: this
    type(global_interpolation_t), intent(inout) :: global_interp
    integer, intent(in) :: n_local

    global_interp%n_points = n_local
    global_interp%n_points_local = n_local
    call global_interp%temp_local%init(n_local)
    call global_interp%local_interp%init(global_interp%Xh, &
         global_interp%rst_local, n_local)
    global_interp%all_points_local = .true.
  end subroutine localize_global_interpolation

  !> Exchange particle ids according to global-interpolation ownership data.
  !! @param migrate_comm Redistribution communicator from interpolation.
  !! @param ids_old Particle ids before migration.
  !! @param n_particles_old Number of particles before migration.
  !! @param n_local Number of particles after migration.
  !! @param particle_ids_local Particle ids after migration.
  subroutine migrate_particle_ids(this, migrate_comm, ids_old, &
       n_particles_old, n_local, particle_ids_local)
    class(lpt_migrate_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: migrate_comm
    integer, allocatable, intent(in) :: ids_old(:)
    integer, intent(in) :: n_particles_old
    integer, intent(in) :: n_local
    integer, allocatable, intent(inout) :: particle_ids_local(:)
    real(kind=rp), allocatable :: sendbuf(:)
    real(kind=rp), allocatable :: recvbuf(:)
    integer :: n_sendbuf
    integer :: n_recvbuf

    if (n_particles_old .eq. 0 .and. n_local .eq. 0) return

    n_sendbuf = max(1, n_particles_old)
    n_recvbuf = max(1, n_local)
    allocate(sendbuf(n_sendbuf))
    allocate(recvbuf(n_recvbuf))
    sendbuf = 0.0_rp
    recvbuf = 0.0_rp
    if (n_particles_old .gt. 0) sendbuf = real(ids_old, rp)
    call migrate_comm%sendrecv(sendbuf, recvbuf, n_sendbuf, n_recvbuf)
    if (n_local .gt. 0) particle_ids_local = nint(recvbuf(1:n_local))
    if (allocated(sendbuf)) deallocate(sendbuf)
    if (allocated(recvbuf)) deallocate(recvbuf)
  end subroutine migrate_particle_ids

  !> Exchange one particle scalar according to interpolation ownership data.
  !! @param migrate_comm Redistribution communicator from interpolation.
  !! @param scalar_old Scalar values before migration.
  !! @param n_particles_old Number of particles before migration.
  !! @param n_local Number of particles after migration.
  !! @param scalar_local Scalar values after migration.
  subroutine migrate_particle_scalar(this, migrate_comm, scalar_old, &
       n_particles_old, n_local, scalar_local)
    class(lpt_migrate_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: migrate_comm
    type(vector_t), intent(inout) :: scalar_old
    integer, intent(in) :: n_particles_old
    integer, intent(in) :: n_local
    type(vector_t), intent(inout) :: scalar_local
    real(kind=rp), allocatable :: sendbuf(:)
    real(kind=rp), allocatable :: recvbuf(:)
    integer :: n_sendbuf
    integer :: n_recvbuf

    if (n_particles_old .eq. 0 .and. n_local .eq. 0) return

    n_sendbuf = max(1, n_particles_old)
    n_recvbuf = max(1, n_local)
    allocate(sendbuf(n_sendbuf))
    allocate(recvbuf(n_recvbuf))
    sendbuf = 0.0_rp
    recvbuf = 0.0_rp
    call scalar_old%copy_from(DEVICE_TO_HOST, .true.)
    if (n_particles_old .gt. 0) sendbuf = scalar_old%x
    call migrate_comm%sendrecv(sendbuf, recvbuf, n_sendbuf, n_recvbuf)
    call scalar_local%copy_from(HOST_TO_DEVICE, .true.)
    if (n_local .gt. 0) scalar_local = recvbuf(1:n_local)
    if (allocated(sendbuf)) deallocate(sendbuf)
    if (allocated(recvbuf)) deallocate(recvbuf)
  end subroutine migrate_particle_scalar

end module lpt_migrate
