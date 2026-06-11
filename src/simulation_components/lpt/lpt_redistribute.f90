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
module lpt_redistribute
  use num_types, only : rp
  use global_interpolation, only : global_interpolation_t
  use glb_intrp_comm, only : glb_intrp_comm_t
  use lpt_periodic_bc, only : lpt_periodic_bc_t
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_SUM
  implicit none
  private

  type, public :: lpt_redistribute_t
     integer :: lag_len = 0
   contains
     procedure, pass(this) :: init => lpt_redistribute_init
     procedure, pass(this) :: free => lpt_redistribute_free
     procedure, pass(this) :: redistribute_particles
     procedure, private, pass(this) :: redistribute_particle_ids
     procedure, private, pass(this) :: redistribute_particle_field
     procedure, private, pass(this) :: redistribute_particle_scalar
     procedure, private, pass(this) :: redistribute_lags
  end type lpt_redistribute_t

contains

  subroutine lpt_redistribute_init(this, lag_len)
    class(lpt_redistribute_t), intent(inout) :: this
    integer, intent(in) :: lag_len

    this%lag_len = lag_len
  end subroutine lpt_redistribute_init

  subroutine lpt_redistribute_free(this)
    class(lpt_redistribute_t), intent(inout) :: this

    this%lag_len = 0
  end subroutine lpt_redistribute_free

  !> Redistribute particles to the rank that owns their current location.
  subroutine redistribute_particles(this, global_interp, periodic_bc, inertia, &
       xyz, ids, vel_lag, vel, acc, d, rho, acc_lag, n, n_global)
    class(lpt_redistribute_t), intent(inout) :: this
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
    type(glb_intrp_comm_t) :: redist_comm
    integer :: n_particles_old
    integer, allocatable :: particle_ids_local(:)
    real(kind=rp), allocatable :: vel_local(:, :)
    real(kind=rp), allocatable :: acc_local(:, :)
    real(kind=rp), allocatable :: d_local(:)
    real(kind=rp), allocatable :: rho_local(:)
    real(kind=rp), allocatable :: vel_particles_lag_local(:, :, :)
    real(kind=rp), allocatable :: acc_particles_lag_local(:, :, :)
    integer :: ierr

    n_particles_old = n
    call periodic_bc%wrap(xyz, n, vel, vel_lag, acc_lag)
    call global_interp%find_points_and_redist(xyz, n)
    call global_interp%init_redist_comm(redist_comm)
    call this%redistribute_particle_ids(redist_comm, ids, n_particles_old, n, &
         particle_ids_local)
    call this%redistribute_lags(redist_comm, vel_lag, n_particles_old, n, &
         vel_particles_lag_local)
    if (inertia) then
       call this%redistribute_particle_field(redist_comm, vel, &
            n_particles_old, n, vel_local)
       call this%redistribute_particle_field(redist_comm, acc, &
            n_particles_old, n, acc_local)
       call this%redistribute_particle_scalar(redist_comm, d, &
            n_particles_old, n, d_local)
       call this%redistribute_particle_scalar(redist_comm, rho, &
            n_particles_old, n, rho_local)
       call this%redistribute_lags(redist_comm, acc_lag, n_particles_old, n, &
            acc_particles_lag_local)
    end if
    call redist_comm%free()

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

    call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
  end subroutine redistribute_particles

  subroutine redistribute_particle_ids(this, redist_comm, ids_old, &
       n_particles_old, n_local, particle_ids_local)
    class(lpt_redistribute_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: redist_comm
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
    call redist_comm%sendrecv(sendbuf, recvbuf, n_particles_old, n_local)
    particle_ids_local = nint(recvbuf)
    deallocate(sendbuf)
    deallocate(recvbuf)
  end subroutine redistribute_particle_ids

  subroutine redistribute_lags(this, redist_comm, lag_old, n_particles_old, &
       n_local, lag_local)
    class(lpt_redistribute_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: redist_comm
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
          call redist_comm%sendrecv(sendbuf, recvbuf, n_particles_old, n_local)
          lag_local(i, j, :) = recvbuf
       end do
    end do
    deallocate(sendbuf)
    deallocate(recvbuf)
  end subroutine redistribute_lags

  subroutine redistribute_particle_field(this, redist_comm, field_old, &
       n_particles_old, n_local, field_local)
    class(lpt_redistribute_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: redist_comm
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
       call redist_comm%sendrecv(sendbuf, recvbuf, n_particles_old, n_local)
       field_local(i, :) = recvbuf
    end do
    deallocate(sendbuf)
    deallocate(recvbuf)
  end subroutine redistribute_particle_field

  subroutine redistribute_particle_scalar(this, redist_comm, scalar_old, &
       n_particles_old, n_local, scalar_local)
    class(lpt_redistribute_t), intent(inout) :: this
    type(glb_intrp_comm_t), intent(inout) :: redist_comm
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
    call redist_comm%sendrecv(sendbuf, recvbuf, n_particles_old, n_local)
    scalar_local = recvbuf
    deallocate(sendbuf)
    deallocate(recvbuf)
  end subroutine redistribute_particle_scalar

end module lpt_redistribute
