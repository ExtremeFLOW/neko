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
     procedure, private, pass(this) :: distribute_particle_scalar
     procedure, private, pass(this) :: localize_global_interpolation
     procedure, private, pass(this) :: migrate_particle_ids
     procedure, private, pass(this) :: migrate_particle_scalar
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

  subroutine initialize_particle_distribution(this, inertia, x, y, z, ids, &
       u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, u, v, w, &
       acc_x, acc_y, acc_z, d, rho, acc_xlag, acc_ylag, acc_zlag, &
       acc_xlaglag, acc_ylaglag, acc_zlaglag, n)
    class(lpt_migrate_t), intent(inout) :: this
    integer, intent(inout) :: n
    logical, intent(in) :: inertia
    real(kind=rp), allocatable, intent(inout) :: x(:), y(:), z(:)
    integer, allocatable, intent(inout) :: ids(:)
    real(kind=rp), allocatable, intent(inout) :: u_lag(:), v_lag(:), w_lag(:)
    real(kind=rp), allocatable, intent(inout) :: u_laglag(:), v_laglag(:)
    real(kind=rp), allocatable, intent(inout) :: w_laglag(:)
    real(kind=rp), allocatable, intent(inout) :: u(:), v(:), w(:)
    real(kind=rp), allocatable, intent(inout) :: acc_x(:), acc_y(:), acc_z(:)
    real(kind=rp), allocatable, intent(inout) :: d(:)
    real(kind=rp), allocatable, intent(inout) :: rho(:)
    real(kind=rp), allocatable, intent(inout) :: acc_xlag(:), acc_ylag(:)
    real(kind=rp), allocatable, intent(inout) :: acc_zlag(:)
    real(kind=rp), allocatable, intent(inout) :: acc_xlaglag(:)
    real(kind=rp), allocatable, intent(inout) :: acc_ylaglag(:)
    real(kind=rp), allocatable, intent(inout) :: acc_zlaglag(:)

    if (this%strategy .eq. LPT_MIGRATE_NONE) then
       call this%distribute_particles_evenly(inertia, x, y, z, ids, &
            u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, &
            u, v, w, acc_x, acc_y, acc_z, d, rho, &
            acc_xlag, acc_ylag, acc_zlag, acc_xlaglag, acc_ylaglag, &
            acc_zlaglag, n)
    end if
  end subroutine initialize_particle_distribution

  !> Update particle ownership according to the selected migration strategy.
  subroutine migrate_particles(this, global_interp, periodic_bc, inertia, &
       x, y, z, ids, u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, &
       u, v, w, acc_x, acc_y, acc_z, d, rho, acc_xlag, acc_ylag, acc_zlag, &
       acc_xlaglag, acc_ylaglag, acc_zlaglag, n, n_global)
    class(lpt_migrate_t), intent(inout) :: this
    integer, intent(inout) :: n
    type(global_interpolation_t), intent(inout) :: global_interp
    type(lpt_periodic_bc_t), intent(inout) :: periodic_bc
    logical, intent(in) :: inertia
    real(kind=rp), allocatable, intent(inout) :: x(:), y(:), z(:)
    integer, allocatable, intent(inout) :: ids(:)
    real(kind=rp), allocatable, intent(inout) :: u_lag(:), v_lag(:), w_lag(:)
    real(kind=rp), allocatable, intent(inout) :: u_laglag(:), v_laglag(:)
    real(kind=rp), allocatable, intent(inout) :: w_laglag(:)
    real(kind=rp), allocatable, intent(inout) :: u(:), v(:), w(:)
    real(kind=rp), allocatable, intent(inout) :: acc_x(:), acc_y(:), acc_z(:)
    real(kind=rp), allocatable, intent(inout) :: d(:)
    real(kind=rp), allocatable, intent(inout) :: rho(:)
    real(kind=rp), allocatable, intent(inout) :: acc_xlag(:), acc_ylag(:)
    real(kind=rp), allocatable, intent(inout) :: acc_zlag(:)
    real(kind=rp), allocatable, intent(inout) :: acc_xlaglag(:)
    real(kind=rp), allocatable, intent(inout) :: acc_ylaglag(:)
    real(kind=rp), allocatable, intent(inout) :: acc_zlaglag(:)
    integer, intent(out) :: n_global
    type(glb_intrp_comm_t) :: migrate_comm
    integer :: n_particles_old
    integer :: n_particles_local
    integer, allocatable :: particle_ids_local(:)
    real(kind=rp), allocatable :: x_local(:)
    real(kind=rp), allocatable :: y_local(:)
    real(kind=rp), allocatable :: z_local(:)
    real(kind=rp), allocatable :: u_local(:)
    real(kind=rp), allocatable :: v_local(:)
    real(kind=rp), allocatable :: w_local(:)
    real(kind=rp), allocatable :: acc_xlocal(:)
    real(kind=rp), allocatable :: acc_ylocal(:)
    real(kind=rp), allocatable :: acc_zlocal(:)
    real(kind=rp), allocatable :: d_local(:)
    real(kind=rp), allocatable :: rho_local(:)
    real(kind=rp), allocatable :: u_lag_local(:), v_lag_local(:)
    real(kind=rp), allocatable :: w_lag_local(:), u_laglag_local(:)
    real(kind=rp), allocatable :: v_laglag_local(:), w_laglag_local(:)
    real(kind=rp), allocatable :: acc_xlag_local(:), acc_ylag_local(:)
    real(kind=rp), allocatable :: acc_zlag_local(:), acc_xlaglag_local(:)
    real(kind=rp), allocatable :: acc_ylaglag_local(:), acc_zlaglag_local(:)
    integer :: ierr
    integer :: rank
    logical :: migration_needed

    if (inertia) then
       call periodic_bc%wrap(x, y, z, n, u, v, w, u_lag, v_lag, w_lag, &
            u_laglag, v_laglag, w_laglag, acc_xlag, acc_ylag, acc_zlag, &
            acc_xlaglag, acc_ylaglag, acc_zlaglag)
    else
       call periodic_bc%wrap(x, y, z, n, u, v, w, u_lag, v_lag, w_lag, &
            u_laglag, v_laglag, w_laglag)
    end if
    if (this%strategy .eq. LPT_MIGRATE_NONE) then
       call global_interp%find_points(x, y, z, n)
       call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, &
            ierr)
       return
    end if

    n_particles_old = n
    call global_interp%find_points(x, y, z, n)
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
    call this%migrate_particle_scalar(migrate_comm, x, n_particles_old, &
         n_particles_local, x_local)
    call this%migrate_particle_scalar(migrate_comm, y, n_particles_old, &
         n_particles_local, y_local)
    call this%migrate_particle_scalar(migrate_comm, z, n_particles_old, &
         n_particles_local, z_local)
    call this%migrate_particle_ids(migrate_comm, ids, n_particles_old, &
         n_particles_local, particle_ids_local)
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
    call move_alloc(x_local, x)
    call move_alloc(y_local, y)
    call move_alloc(z_local, z)
    call move_alloc(particle_ids_local, ids)
    call move_alloc(u_lag_local, u_lag)
    call move_alloc(v_lag_local, v_lag)
    call move_alloc(w_lag_local, w_lag)
    call move_alloc(u_laglag_local, u_laglag)
    call move_alloc(v_laglag_local, v_laglag)
    call move_alloc(w_laglag_local, w_laglag)
    if (inertia) then
      call move_alloc(u_local, u)
      call move_alloc(v_local, v)
      call move_alloc(w_local, w)
      call move_alloc(acc_xlocal, acc_x)
      call move_alloc(acc_ylocal, acc_y)
      call move_alloc(acc_zlocal, acc_z)
      call move_alloc(d_local, d)
      call move_alloc(rho_local, rho)
      call move_alloc(acc_xlag_local, acc_xlag)
      call move_alloc(acc_ylag_local, acc_ylag)
      call move_alloc(acc_zlag_local, acc_zlag)
      call move_alloc(acc_xlaglag_local, acc_xlaglag)
      call move_alloc(acc_ylaglag_local, acc_ylaglag)
      call move_alloc(acc_zlaglag_local, acc_zlaglag)
    else
      if (allocated(u)) deallocate(u)
      if (allocated(v)) deallocate(v)
      if (allocated(w)) deallocate(w)
      if (allocated(acc_x)) deallocate(acc_x)
      if (allocated(acc_y)) deallocate(acc_y)
      if (allocated(acc_z)) deallocate(acc_z)
      if (allocated(d)) deallocate(d)
      if (allocated(rho)) deallocate(rho)
      if (allocated(acc_xlag)) deallocate(acc_xlag)
      if (allocated(acc_ylag)) deallocate(acc_ylag)
      if (allocated(acc_zlag)) deallocate(acc_zlag)
      if (allocated(acc_xlaglag)) deallocate(acc_xlaglag)
      if (allocated(acc_ylaglag)) deallocate(acc_ylaglag)
      if (allocated(acc_zlaglag)) deallocate(acc_zlaglag)
    end if

    call this%localize_global_interpolation(global_interp, n)
    call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
  end subroutine migrate_particles

  subroutine distribute_particles_evenly(this, inertia, x, y, z, ids, &
       u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, u, v, w, &
       acc_x, acc_y, acc_z, d, rho, acc_xlag, acc_ylag, acc_zlag, &
       acc_xlaglag, acc_ylaglag, acc_zlaglag, n)
    class(lpt_migrate_t), intent(inout) :: this
    logical, intent(in) :: inertia
    real(kind=rp), allocatable, intent(inout) :: x(:), y(:), z(:)
    integer, allocatable, intent(inout) :: ids(:)
    real(kind=rp), allocatable, intent(inout) :: u_lag(:), v_lag(:), w_lag(:)
    real(kind=rp), allocatable, intent(inout) :: u_laglag(:), v_laglag(:)
    real(kind=rp), allocatable, intent(inout) :: w_laglag(:)
    real(kind=rp), allocatable, intent(inout) :: u(:), v(:), w(:)
    real(kind=rp), allocatable, intent(inout) :: acc_x(:), acc_y(:), acc_z(:)
    real(kind=rp), allocatable, intent(inout) :: d(:)
    real(kind=rp), allocatable, intent(inout) :: rho(:)
    real(kind=rp), allocatable, intent(inout) :: acc_xlag(:), acc_ylag(:)
    real(kind=rp), allocatable, intent(inout) :: acc_zlag(:)
    real(kind=rp), allocatable, intent(inout) :: acc_xlaglag(:)
    real(kind=rp), allocatable, intent(inout) :: acc_ylaglag(:)
    real(kind=rp), allocatable, intent(inout) :: acc_zlaglag(:)
    integer, intent(inout) :: n
    integer, allocatable :: ids_local(:)
    real(kind=rp), allocatable :: x_local(:)
    real(kind=rp), allocatable :: y_local(:)
    real(kind=rp), allocatable :: z_local(:)
    real(kind=rp), allocatable :: u_lag_local(:), v_lag_local(:)
    real(kind=rp), allocatable :: w_lag_local(:), u_laglag_local(:)
    real(kind=rp), allocatable :: v_laglag_local(:), w_laglag_local(:)
    real(kind=rp), allocatable :: u_local(:)
    real(kind=rp), allocatable :: v_local(:)
    real(kind=rp), allocatable :: w_local(:)
    real(kind=rp), allocatable :: acc_xlocal(:)
    real(kind=rp), allocatable :: acc_ylocal(:)
    real(kind=rp), allocatable :: acc_zlocal(:)
    real(kind=rp), allocatable :: d_local(:)
    real(kind=rp), allocatable :: rho_local(:)
    real(kind=rp), allocatable :: acc_xlag_local(:), acc_ylag_local(:)
    real(kind=rp), allocatable :: acc_zlag_local(:), acc_xlaglag_local(:)
    real(kind=rp), allocatable :: acc_ylaglag_local(:), acc_zlaglag_local(:)
    integer, allocatable :: counts(:)
    integer, allocatable :: offsets(:)
    integer :: n_local

    call build_even_particle_distribution(n, counts, offsets, n_local)

    call this%distribute_particle_scalar(x, counts, offsets, n_local, x_local)
    call this%distribute_particle_scalar(y, counts, offsets, n_local, y_local)
    call this%distribute_particle_scalar(z, counts, offsets, n_local, z_local)
    call this%distribute_particle_ids(ids, counts, offsets, n_local, &
         ids_local)
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
       call this%distribute_particle_scalar(acc_xlag, counts, offsets, n_local, &
            acc_xlag_local)
       call this%distribute_particle_scalar(acc_ylag, counts, offsets, n_local, &
            acc_ylag_local)
       call this%distribute_particle_scalar(acc_zlag, counts, offsets, n_local, &
            acc_zlag_local)
       call this%distribute_particle_scalar(acc_xlaglag, counts, offsets, &
            n_local, acc_xlaglag_local)
       call this%distribute_particle_scalar(acc_ylaglag, counts, offsets, &
            n_local, acc_ylaglag_local)
       call this%distribute_particle_scalar(acc_zlaglag, counts, offsets, &
            n_local, acc_zlaglag_local)
    end if

    n = n_local
    call move_alloc(x_local, x)
    call move_alloc(y_local, y)
    call move_alloc(z_local, z)
    call move_alloc(ids_local, ids)
    call move_alloc(u_lag_local, u_lag)
    call move_alloc(v_lag_local, v_lag)
    call move_alloc(w_lag_local, w_lag)
    call move_alloc(u_laglag_local, u_laglag)
    call move_alloc(v_laglag_local, v_laglag)
    call move_alloc(w_laglag_local, w_laglag)
    if (inertia) then
       call move_alloc(u_local, u)
       call move_alloc(v_local, v)
       call move_alloc(w_local, w)
       call move_alloc(acc_xlocal, acc_x)
       call move_alloc(acc_ylocal, acc_y)
       call move_alloc(acc_zlocal, acc_z)
       call move_alloc(d_local, d)
       call move_alloc(rho_local, rho)
       call move_alloc(acc_xlag_local, acc_xlag)
       call move_alloc(acc_ylag_local, acc_ylag)
       call move_alloc(acc_zlag_local, acc_zlag)
       call move_alloc(acc_xlaglag_local, acc_xlaglag)
       call move_alloc(acc_ylaglag_local, acc_ylaglag)
       call move_alloc(acc_zlaglag_local, acc_zlaglag)
    else
       if (allocated(u)) deallocate(u)
       if (allocated(v)) deallocate(v)
       if (allocated(w)) deallocate(w)
       if (allocated(acc_x)) deallocate(acc_x)
       if (allocated(acc_y)) deallocate(acc_y)
       if (allocated(acc_z)) deallocate(acc_z)
       if (allocated(d)) deallocate(d)
       if (allocated(rho)) deallocate(rho)
       if (allocated(acc_xlag)) deallocate(acc_xlag)
       if (allocated(acc_ylag)) deallocate(acc_ylag)
       if (allocated(acc_zlag)) deallocate(acc_zlag)
       if (allocated(acc_xlaglag)) deallocate(acc_xlaglag)
       if (allocated(acc_ylaglag)) deallocate(acc_ylaglag)
       if (allocated(acc_zlaglag)) deallocate(acc_zlaglag)
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

  subroutine localize_global_interpolation(this, global_interp, n_local)
    class(lpt_migrate_t), intent(inout) :: this
    type(global_interpolation_t), intent(inout) :: global_interp
    integer, intent(in) :: n_local

    global_interp%n_points = n_local
    global_interp%all_points_local = .true.
  end subroutine localize_global_interpolation

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
