! Copyright (c) 2021-2026, The Neko Authors
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
!> Defines format-independent checkpoint registration and restart state.
module checkpoint
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp, dp
  use field_series, only : field_series_t
  use checkpoint_payload, only : checkpoint_payload_t, &
       checkpoint_payload_ptr_t, checkpoint_array_t
  use space, only : space_t, operator(.ne.)
  use device, only : device_memcpy, DEVICE_TO_HOST, HOST_TO_DEVICE, &
       device_sync, glb_cmd_queue
  use field, only : field_t
  use utils, only : neko_error
  use mesh, only : mesh_t
  use global_interpolation, only : GLOB_INTERP_TOL
  use time_state, only : time_state_t
  use, intrinsic :: iso_c_binding, only : c_associated
  implicit none
  private

  !> Collection of live simulation data registered for checkpointing.
  type, public :: chkp_t
     !> Format-independent checkpoint payloads.
     type(checkpoint_payload_ptr_t), allocatable :: payloads(:)

     !> Restart time, valid after loading a checkpoint.
     real(kind=dp) :: t = 0d0
     !> Mesh on which an interpolated legacy checkpoint was written.
     type(mesh_t) :: previous_mesh
     !> Function space of a loaded legacy checkpoint.
     type(space_t) :: previous_Xh
     !> Point-search tolerance for legacy mesh-to-mesh restart.
     real(kind=dp) :: mesh2mesh_tol = GLOB_INTERP_TOL
   contains
     !> Initialize an empty checkpoint.
     procedure, pass(this) :: init => chkp_init
     !> Synchronise registered checkpoint data from device to host.
     procedure, pass(this) :: sync_host => chkp_sync_host
     !> Synchronise registered checkpoint data from host to device.
     procedure, pass(this) :: sync_device => chkp_sync_device
     !> Register the fluid fields.
     procedure, pass(this) :: add_fluid => chkp_add_fluid
     !> Register the lagged velocity fields.
     procedure, pass(this) :: add_lag => chkp_add_lag
     !> Add or return a checkpoint payload.
     procedure, pass(this) :: add_payload => chkp_add_payload
     !> Return a checkpoint payload by name.
     procedure, pass(this) :: get_payload => chkp_get_payload
     !> Return the number of registered checkpoint payloads.
     procedure, pass(this) :: payload_count => chkp_payload_count
     !> Return the number of registered scalar payloads.
     procedure, pass(this) :: scalar_payload_count => chkp_scalar_payload_count
     !> Register the simulation time history as a replicated payload.
     procedure, pass(this) :: add_time_state => chkp_add_time_state
     !> Return pointers to the simulation time history.
     procedure, pass(this) :: get_time_history => chkp_get_time_history
     !> Return the restart time from a loaded checkpoint.
     procedure, pass(this) :: restart_time => chkp_restart_time
     !> Restore a simulation time state from the checkpoint.
     procedure, pass(this) :: set_time_state => chkp_set_time_state
     !> Release all checkpoint payloads and restart metadata.
     procedure, pass(this) :: free => chkp_free
  end type chkp_t

contains

  !> Initialize an empty checkpoint.
  subroutine chkp_init(this)
    class(chkp_t), intent(inout) :: this

    ! Make sure the object is clean
    call this%free()

  end subroutine chkp_init

  !> Release all checkpoint payloads and restart metadata.
  subroutine chkp_free(this)
    class(chkp_t), intent(inout) :: this
    integer :: i

    this%t = 0d0

    if (allocated(this%payloads)) then
       do i = 1, size(this%payloads)
          if (associated(this%payloads(i)%ptr)) then
             call this%payloads(i)%ptr%free()
             deallocate(this%payloads(i)%ptr)
          end if
       end do
       deallocate(this%payloads)
    end if

    call this%previous_mesh%free()
    call this%previous_Xh%free()

  end subroutine chkp_free

  !> Synchronise registered checkpoint data from device to host.
  subroutine chkp_sync_host(this)
    class(chkp_t), intent(inout) :: this
    integer :: i, j

    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, this%payload_count()
          do j = 1, this%payloads(i)%ptr%field_count()
             call this%payloads(i)%ptr%fields(j)%ptr%copy_from( &
                  DEVICE_TO_HOST, sync = .false.)
          end do
          do j = 1, this%payloads(i)%ptr%series_count()
             block
               integer :: k
               do k = 1, this%payloads(i)%ptr%series(j)%ptr%size()
                  call this%payloads(i)%ptr%series(j)%ptr%lf(k)%copy_from( &
                       DEVICE_TO_HOST, sync = .false.)
               end do
             end block
          end do
          do j = 1, this%payloads(i)%ptr%array_count()
             if (c_associated( &
                  this%payloads(i)%ptr%arrays(j)%ptr%x_d)) then
                call device_memcpy( &
                     this%payloads(i)%ptr%arrays(j)%ptr%x, &
                     this%payloads(i)%ptr%arrays(j)%ptr%x_d, &
                     size(this%payloads(i)%ptr%arrays(j)%ptr%x), &
                     DEVICE_TO_HOST, .false.)
             end if
          end do
          do j = 1, this%payloads(i)%ptr%mesh_array_count()
             if (c_associated( &
                  this%payloads(i)%ptr%mesh_arrays(j)%ptr%x_d)) then
                call device_memcpy( &
                     this%payloads(i)%ptr%mesh_arrays(j)%ptr%x, &
                     this%payloads(i)%ptr%mesh_arrays(j)%ptr%x_d, &
                     size(this%payloads(i)%ptr%mesh_arrays(j)%ptr%x), &
                     DEVICE_TO_HOST, .false.)
             end if
          end do
       end do

       call device_sync(glb_cmd_queue)
    end if

  end subroutine chkp_sync_host

  !> Synchronise registered checkpoint data from host to device.
  subroutine chkp_sync_device(this)
    class(chkp_t), intent(inout) :: this
    integer :: i, j

    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, this%payload_count()
          do j = 1, this%payloads(i)%ptr%field_count()
             call this%payloads(i)%ptr%fields(j)%ptr%copy_from( &
                  HOST_TO_DEVICE, sync = .false.)
          end do
          do j = 1, this%payloads(i)%ptr%series_count()
             block
               integer :: k
               do k = 1, this%payloads(i)%ptr%series(j)%ptr%size()
                  call this%payloads(i)%ptr%series(j)%ptr%lf(k)%copy_from( &
                       HOST_TO_DEVICE, sync = .false.)
               end do
             end block
          end do
          do j = 1, this%payloads(i)%ptr%array_count()
             if (c_associated( &
                  this%payloads(i)%ptr%arrays(j)%ptr%x_d)) then
                call device_memcpy( &
                     this%payloads(i)%ptr%arrays(j)%ptr%x, &
                     this%payloads(i)%ptr%arrays(j)%ptr%x_d, &
                     size(this%payloads(i)%ptr%arrays(j)%ptr%x), &
                     HOST_TO_DEVICE, .false.)
             end if
          end do
          do j = 1, this%payloads(i)%ptr%mesh_array_count()
             if (c_associated( &
                  this%payloads(i)%ptr%mesh_arrays(j)%ptr%x_d)) then
                call device_memcpy( &
                     this%payloads(i)%ptr%mesh_arrays(j)%ptr%x, &
                     this%payloads(i)%ptr%mesh_arrays(j)%ptr%x_d, &
                     size(this%payloads(i)%ptr%mesh_arrays(j)%ptr%x), &
                     HOST_TO_DEVICE, .false.)
             end if
          end do
       end do

    end if

  end subroutine chkp_sync_device

  !> Register the fluid fields.
  !! @param u X component of the fluid velocity.
  !! @param v Y component of the fluid velocity.
  !! @param w Z component of the fluid velocity.
  !! @param p Fluid pressure.
  subroutine chkp_add_fluid(this, u, v, w, p)
    class(chkp_t), intent(inout) :: this
    type(field_t), target :: u
    type(field_t), target :: v
    type(field_t), target :: w
    type(field_t), target :: p
    type(checkpoint_payload_t), pointer :: payload

    ! Check that all velocity components are defined on the same
    ! function space
    if ( u%Xh .ne. v%Xh .or. &
         u%Xh .ne. w%Xh ) then
       call neko_error('Different function spaces for each velocity component')
    end if

    ! Check that both velocity and pressure is defined on the same mesh
    if ( u%msh%nelv .ne. p%msh%nelv ) then
       call neko_error('Velocity and pressure defined on different meshes')
    end if

    payload => this%add_payload("fluid")
    call payload%add_field(u)
    call payload%add_field(v)
    call payload%add_field(w)
    call payload%add_field(p)

  end subroutine chkp_add_fluid

  !> Register the lagged velocity fields.
  !! @param ulag Lagged X velocity fields.
  !! @param vlag Lagged Y velocity fields.
  !! @param wlag Lagged Z velocity fields.
  subroutine chkp_add_lag(this, ulag, vlag, wlag)
    class(chkp_t), intent(inout) :: this
    type(field_series_t), target :: ulag
    type(field_series_t), target :: vlag
    type(field_series_t), target :: wlag
    type(checkpoint_payload_t), pointer :: payload

    payload => this%add_payload("fluid")
    call payload%add_series(ulag)
    call payload%add_series(vlag)
    call payload%add_series(wlag)

  end subroutine chkp_add_lag

  !> Return the number of registered scalar payloads.
  !! @return Number of registered scalar payloads.
  pure function chkp_scalar_payload_count(this) result(n)
    class(chkp_t), intent(in) :: this
    integer :: i, n

    n = 0
    do i = 1, this%payload_count()
       if (index(trim(this%payloads(i)%ptr%name), "scalars/") .eq. 1) then
          n = n + 1
       end if
    end do

  end function chkp_scalar_payload_count

  !> Add or return a checkpoint payload.
  !! @param name Payload path to add or find.
  !! @return Pointer to the registered payload.
  function chkp_add_payload(this, name) result(payload)
    class(chkp_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(checkpoint_payload_t), pointer :: payload
    type(checkpoint_payload_ptr_t), allocatable :: tmp(:)
    integer :: i, n

    do i = 1, this%payload_count()
       if (trim(this%payloads(i)%ptr%name) .eq. trim(name)) then
          payload => this%payloads(i)%ptr
          return
       end if
    end do

    n = this%payload_count()
    allocate(tmp(n + 1))
    if (n .gt. 0) tmp(1:n) = this%payloads
    allocate(tmp(n + 1)%ptr)
    call tmp(n + 1)%ptr%init(name)
    call move_alloc(tmp, this%payloads)
    payload => this%payloads(n + 1)%ptr

  end function chkp_add_payload

  !> Return a checkpoint payload by name.
  !! @param name Payload path to find.
  !! @return Pointer to the registered payload.
  function chkp_get_payload(this, name) result(payload)
    class(chkp_t), intent(in) :: this
    character(len=*), intent(in) :: name
    type(checkpoint_payload_t), pointer :: payload
    integer :: i

    nullify(payload)
    do i = 1, this%payload_count()
       if (trim(this%payloads(i)%ptr%name) .eq. trim(name)) then
          payload => this%payloads(i)%ptr
          return
       end if
    end do

    call neko_error("Checkpoint payload '" // trim(name) // "' not found")

  end function chkp_get_payload

  !> Return the number of registered checkpoint payloads.
  !! @return Number of registered payloads.
  pure function chkp_payload_count(this) result(n)
    class(chkp_t), intent(in) :: this
    integer :: n

    if (allocated(this%payloads)) then
       n = size(this%payloads)
    else
       n = 0
    end if

  end function chkp_payload_count

  !> Register the simulation time history as a replicated payload.
  !! @param time_state Simulation time history to register.
  subroutine chkp_add_time_state(this, time_state)
    class(chkp_t), intent(inout) :: this
    type(time_state_t), target, intent(inout) :: time_state
    type(checkpoint_payload_t), pointer :: payload

    payload => this%add_payload("time")
    call payload%add_array("tlag", time_state%tlag, replicated = .true.)
    call payload%add_array("dtlag", time_state%dtlag, replicated = .true.)

  end subroutine chkp_add_time_state

  !> Return pointers to the simulation time history.
  !! @param tlag Pointer to the registered previous times.
  !! @param dtlag Pointer to the registered previous time-step sizes.
  subroutine chkp_get_time_history(this, tlag, dtlag)
    class(chkp_t), intent(in) :: this
    real(kind=rp), pointer, intent(out) :: tlag(:), dtlag(:)
    type(checkpoint_payload_t), pointer :: payload
    type(checkpoint_array_t), pointer :: array

    payload => this%get_payload("time")
    array => payload%find_array("tlag")
    tlag => array%x
    array => payload%find_array("dtlag")
    dtlag => array%x

  end subroutine chkp_get_time_history

  !> Return the restart time from a loaded checkpoint.
  !! @return Restart time stored in the checkpoint.
  pure function chkp_restart_time(this) result(rtime)
    class(chkp_t), intent(in) :: this
    real(kind=dp) :: rtime

    rtime = this%t
  end function chkp_restart_time

  !> Restore a simulation time state from the checkpoint.
  !! @param time_state Simulation time state to restore.
  subroutine chkp_set_time_state(this, time_state)
    class(chkp_t), intent(in) :: this
    type(time_state_t), intent(inout) :: time_state
    real(kind=rp), pointer :: tlag(:), dtlag(:)

    call this%get_time_history(tlag, dtlag)
    time_state%t = this%t
    time_state%dtlag = dtlag
    time_state%tlag = tlag
  end subroutine chkp_set_time_state

end module checkpoint
