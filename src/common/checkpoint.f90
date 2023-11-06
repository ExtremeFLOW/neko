! Copyright (c) 2021, The Neko Authors
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
!> Defines a checkpoint
module checkpoint
  use neko_config
  use num_types
  use field_series
  use space
  use device
  use field
  use space
  use utils
  use mesh, only: mesh_t
  implicit none
  private

  type, public :: chkp_t
     type(field_t), pointer :: u => null()
     type(field_t), pointer :: v => null()
     type(field_t), pointer :: w => null()
     type(field_t), pointer :: p => null()

     !
     ! Optional payload
     !
     type(field_series_t), pointer :: ulag => null()
     type(field_series_t), pointer :: vlag => null()
     type(field_series_t), pointer :: wlag => null()

     type(field_t), pointer :: s => null()

     real(kind=dp) :: t         !< Restart time (valid after load)
     type(mesh_t) :: previous_mesh
     real(kind=rp) :: mesh2mesh_tol = 1d-6

   contains
     procedure, pass(this) :: init => chkp_init
     procedure, pass(this) :: sync_host => chkp_sync_host
     procedure, pass(this) :: sync_device => chkp_sync_device
     procedure, pass(this) :: add_lag => chkp_add_lag
     procedure, pass(this) :: add_scalar => chkp_add_scalar
     procedure, pass(this) :: restart_time => chkp_restart_time
     final :: chkp_free
  end type chkp_t

contains

  !> Initialize checkpoint structure with mandatory data
  subroutine chkp_init(this, u, v, w, p)
    class(chkp_t), intent(inout) :: this
    type(field_t), intent(in), target :: u
    type(field_t), intent(in), target :: v
    type(field_t), intent(in), target :: w
    type(field_t), intent(in), target :: p

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
    
    this%u => u
    this%v => v
    this%w => w
    this%p => p

    this%t = 0d0
    
  end subroutine chkp_init

  !> Reset checkpoint
  subroutine chkp_free(this)
    type(chkp_t), intent(inout) :: this

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)

    nullify(this%ulag)
    nullify(this%vlag)
    nullify(this%wlag)
    
  end subroutine chkp_free

  !> Synchronize checkpoint with device
  subroutine chkp_sync_host(this)
    class(chkp_t), intent(inout) :: this

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(u=>this%u, v=>this%v, w=>this%w, &
            ulag=>this%ulag, vlag=>this%vlag, wlag=>this%wlag)

         if (associated(this%u) .and. associated(this%v) .and. &
              associated(this%w)) then
            call device_memcpy(u%x, u%x_d, u%dof%size(), DEVICE_TO_HOST)
            call device_memcpy(v%x, v%x_d, v%dof%size(), DEVICE_TO_HOST)
            call device_memcpy(w%x, w%x_d, w%dof%size(), DEVICE_TO_HOST)
         end if
         
         if (associated(this%ulag) .and. associated(this%vlag) .and. &
              associated(this%wlag)) then
            call device_memcpy(ulag%lf(1)%x, ulag%lf(1)%x_d, &
                               u%dof%size(), DEVICE_TO_HOST)
            call device_memcpy(ulag%lf(2)%x, ulag%lf(2)%x_d, &
                               u%dof%size(), DEVICE_TO_HOST)
  
            call device_memcpy(vlag%lf(1)%x, vlag%lf(1)%x_d, &
                               v%dof%size(), DEVICE_TO_HOST)
            call device_memcpy(vlag%lf(2)%x, vlag%lf(2)%x_d, &
                               v%dof%size(), DEVICE_TO_HOST)
    
            call device_memcpy(wlag%lf(1)%x, wlag%lf(1)%x_d, &
                               w%dof%size(), DEVICE_TO_HOST)
            call device_memcpy(wlag%lf(2)%x, wlag%lf(2)%x_d, &
                               w%dof%size(), DEVICE_TO_HOST)
         end if
         if (associated(this%s)) then
            call device_memcpy(this%s%x, this%s%x_d, &
                               this%s%dof%size(), DEVICE_TO_HOST)
         end if
       end associate
    end if
         
  end subroutine chkp_sync_host

  !> Synchronize device with checkpoint
  subroutine chkp_sync_device(this)
    class(chkp_t), intent(inout) :: this

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(u=>this%u, v=>this%v, w=>this%w, &
            ulag=>this%ulag, vlag=>this%vlag, wlag=>this%wlag)

         if (associated(this%u) .and. associated(this%v) .and. &
              associated(this%w)) then
            call device_memcpy(u%x, u%x_d, u%dof%size(), HOST_TO_DEVICE)
            call device_memcpy(v%x, v%x_d, v%dof%size(), HOST_TO_DEVICE)
            call device_memcpy(w%x, w%x_d, w%dof%size(), HOST_TO_DEVICE)
         end if
         
         if (associated(this%ulag) .and. associated(this%vlag) .and. &
              associated(this%wlag)) then
            call device_memcpy(ulag%lf(1)%x, ulag%lf(1)%x_d, &
                               u%dof%size(), HOST_TO_DEVICE)
            call device_memcpy(ulag%lf(2)%x, ulag%lf(2)%x_d, &
                               u%dof%size(), HOST_TO_DEVICE)
  
            call device_memcpy(vlag%lf(1)%x, vlag%lf(1)%x_d, &
                               v%dof%size(), HOST_TO_DEVICE)
            call device_memcpy(vlag%lf(2)%x, vlag%lf(2)%x_d, &
                               v%dof%size(), HOST_TO_DEVICE)
    
            call device_memcpy(wlag%lf(1)%x, wlag%lf(1)%x_d, &
                               w%dof%size(), HOST_TO_DEVICE)
            call device_memcpy(wlag%lf(2)%x, wlag%lf(2)%x_d, &
                               w%dof%size(), HOST_TO_DEVICE)
         end if
         if (associated(this%s)) then
            call device_memcpy(this%s%x, this%s%x_d, &
                               this%s%dof%size(), HOST_TO_DEVICE)
         end if
       end associate
    end if
         
  end subroutine chkp_sync_device

  !> Add lagged velocity terms
  subroutine chkp_add_lag(this, ulag, vlag, wlag)
    class(chkp_t), intent(inout) :: this    
    type(field_series_t), target :: ulag
    type(field_series_t), target :: vlag
    type(field_series_t), target :: wlag

    this%ulag => ulag
    this%vlag => vlag
    this%wlag => wlag
    
  end subroutine chkp_add_lag

  !> Add scalars
  subroutine chkp_add_scalar(this, s)
    class(chkp_t), intent(inout) :: this    
    type(field_t), target :: s

    this%s => s
    
  end subroutine chkp_add_scalar


  !> Return restart time from a loaded checkpoint
  pure function chkp_restart_time(this) result(rtime)
    class(chkp_t), intent(in) :: this
    real(kind=dp) :: rtime

    rtime = this%t
  end function chkp_restart_time
  
end module checkpoint
