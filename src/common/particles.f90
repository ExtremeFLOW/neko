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
!> Defines a collection of Lagrangian particles.
module particles
  use num_types, only : rp
  use vector, only : vector_t
  use utils, only : neko_error
  implicit none
  private

  !> Particle positions, velocities, properties, and time-history data.
  type, public :: particles_t
     type(vector_t) :: x, y, z
     type(vector_t) :: u, v, w
     integer, allocatable :: ids(:)
     type(vector_t) :: u_lag, v_lag, w_lag
     type(vector_t) :: u_laglag, v_laglag, w_laglag
     integer :: time_order, lag_len
     integer :: n = 0
     integer :: n_global = 0
     logical :: inertia = .false.
     type(vector_t) :: acc_x, acc_y, acc_z
     type(vector_t) :: d
     type(vector_t) :: rho
     type(vector_t) :: acc_xlag, acc_ylag, acc_zlag
     type(vector_t) :: acc_xlaglag, acc_ylaglag, acc_zlaglag
   contains
     procedure, pass(this) :: init => particles_init
     procedure, pass(this) :: free => particles_free
     procedure, pass(this) :: device_sync => particles_device_sync
  end type particles_t

contains

  subroutine particles_init(this, x, y, z, time_order, u, v, w, diameter, &
       density)
    class(particles_t), intent(inout) :: this
    real(kind=rp), intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: time_order
    real(kind=rp), intent(in), optional :: u(:), v(:), w(:)
    real(kind=rp), intent(in), optional :: diameter(:)
    real(kind=rp), intent(in), optional :: density(:)

    integer :: i

    if (size(y) .ne. size(x) .or. size(z) .ne. size(x)) then
       call neko_error("particle coordinate arrays must have the same size")
    end if

    call this%free()

    this%n = size(x)
    this%n_global = this%n
    this%time_order = time_order
    this%lag_len = time_order - 1

    if (present(diameter) .neqv. present(density)) then
       call neko_error("particle diameter and density must both be provided")
    end if
    this%inertia = present(diameter)
    call this%x%init(this%n)
    call this%y%init(this%n)
    call this%z%init(this%n)
    call this%u%init(this%n)
    call this%v%init(this%n)
    call this%w%init(this%n)
    call this%acc_x%init(this%n)
    call this%acc_y%init(this%n)
    call this%acc_z%init(this%n)
    allocate(this%ids(this%n))
    call this%u_lag%init(this%n)
    call this%v_lag%init(this%n)
    call this%w_lag%init(this%n)
    call this%u_laglag%init(this%n)
    call this%v_laglag%init(this%n)
    call this%w_laglag%init(this%n)
    call this%acc_xlag%init(this%n)
    call this%acc_ylag%init(this%n)
    call this%acc_zlag%init(this%n)
    call this%acc_xlaglag%init(this%n)
    call this%acc_ylaglag%init(this%n)
    call this%acc_zlaglag%init(this%n)
    call this%d%init(this%n)
    call this%rho%init(this%n)
    this%x = x
    this%y = y
    this%z = z
    if (present(u) .and. present(v) .and. present(w)) then
       this%u = u
       this%v = v
       this%w = w
    else
       this%u = 0.0_rp
       this%v = 0.0_rp
       this%w = 0.0_rp
    end if
    if (present(diameter)) then
       this%d = diameter
    else
       this%d = 0.0_rp
    end if
    if (present(density)) then
       this%rho = density
    else
       this%rho = 0.0_rp
    end if
    do i = 1, this%n
       this%ids(i) = i
    end do
  end subroutine particles_init

  subroutine particles_free(this)
    class(particles_t), intent(inout) :: this

    call this%x%free()
    call this%y%free()
    call this%z%free()
    if (allocated(this%ids)) deallocate(this%ids)
    call this%u%free()
    call this%v%free()
    call this%w%free()
    call this%acc_x%free()
    call this%acc_y%free()
    call this%acc_z%free()
    call this%u_lag%free()
    call this%v_lag%free()
    call this%w_lag%free()
    call this%u_laglag%free()
    call this%v_laglag%free()
    call this%w_laglag%free()
    call this%acc_xlag%free()
    call this%acc_ylag%free()
    call this%acc_zlag%free()
    call this%acc_xlaglag%free()
    call this%acc_ylaglag%free()
    call this%acc_zlaglag%free()
    call this%d%free()
    call this%rho%free()
    this%n = 0
    this%n_global = 0
  end subroutine particles_free

  !> Synchronise the host and device data for particles.
  subroutine particles_device_sync(this, memdir)
    class(particles_t), intent(inout) :: this
    integer, intent(in) :: memdir

    call this%x%copy_from(memdir, .false.)
    call this%y%copy_from(memdir, .false.)
    call this%z%copy_from(memdir, .false.)
    call this%u%copy_from(memdir, .false.)
    call this%v%copy_from(memdir, .false.)
    call this%w%copy_from(memdir, .true.)
    call this%acc_x%copy_from(memdir, .false.)
    call this%acc_y%copy_from(memdir, .false.)
    call this%acc_z%copy_from(memdir, .true.)
    call this%u_lag%copy_from(memdir, .false.)
    call this%v_lag%copy_from(memdir, .false.)
    call this%w_lag%copy_from(memdir, .true.)
    call this%u_laglag%copy_from(memdir, .false.)
    call this%v_laglag%copy_from(memdir, .false.)
    call this%w_laglag%copy_from(memdir, .true.)
    call this%acc_xlag%copy_from(memdir, .false.)
    call this%acc_ylag%copy_from(memdir, .false.)
    call this%acc_zlag%copy_from(memdir, .true.)
    call this%acc_xlaglag%copy_from(memdir, .false.)
    call this%acc_ylaglag%copy_from(memdir, .false.)
    call this%acc_zlaglag%copy_from(memdir, .true.)
    call this%d%copy_from(memdir, .false.)
    call this%rho%copy_from(memdir, .true.)

  end subroutine particles_device_sync

end module particles
