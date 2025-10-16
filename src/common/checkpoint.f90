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
  use num_types, only : rp, dp
  use field_series, only : field_series_t
  use field_series_list, only : field_series_list_t
  use space, only : space_t, operator(.ne.)
  use device, only : device_memcpy, DEVICE_TO_HOST, HOST_TO_DEVICE, &
       device_sync, glb_cmd_queue
  use field, only : field_t, field_ptr_t
  use utils, only : neko_error, filename_suffix_pos
  use mesh, only: mesh_t
  use math, only : NEKO_EPS
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

     real(kind=rp), pointer :: tlag(:) => null()
     real(kind=rp), pointer :: dtlag(:) => null()

     !> for pnpn
     type(field_t), pointer :: abx1 => null()
     type(field_t), pointer :: abx2 => null()
     type(field_t), pointer :: aby1 => null()
     type(field_t), pointer :: aby2 => null()
     type(field_t), pointer :: abz1 => null()
     type(field_t), pointer :: abz2 => null()

     type(field_t), pointer :: s => null()
     type(field_series_t), pointer :: slag => null()
     type(field_t), pointer :: abs1 => null()
     type(field_t), pointer :: abs2 => null()

     type(field_series_list_t) :: scalar_lags

     !> Multi-scalar ABX fields
     type(field_ptr_t), allocatable :: scalar_abx1(:) !< ABX1 fields for each scalar
     type(field_ptr_t), allocatable :: scalar_abx2(:) !< ABX2 fields for each scalar

     real(kind=dp) :: t !< Restart time (valid after load)
     type(mesh_t) :: previous_mesh
     type(space_t) :: previous_Xh
     real(kind=rp) :: mesh2mesh_tol = NEKO_EPS*1e3_rp

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

    ! Scalar cleanup
    nullify(this%s)
    nullify(this%slag)
    nullify(this%abs1)
    nullify(this%abs2)

    ! Free scalar lag list if it was initialized
    if (allocated(this%scalar_lags%items)) then
       call this%scalar_lags%free()
    end if

    ! Free multi-scalar ABX field arrays
    if (allocated(this%scalar_abx1)) then
       deallocate(this%scalar_abx1)
    end if
    if (allocated(this%scalar_abx2)) then
       deallocate(this%scalar_abx2)
    end if

  end subroutine chkp_free

  !> Synchronize checkpoint with device
  subroutine chkp_sync_host(this)
    class(chkp_t), intent(inout) :: this
    integer :: i, j

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(u=>this%u, v=>this%v, w=>this%w, &
            ulag=>this%ulag, vlag=>this%vlag, wlag=>this%wlag, &
            p=>this%p)

         if (associated(this%u) .and. associated(this%v) .and. &
              associated(this%w) .and. associated(this%p)) then
            call device_memcpy(u%x, u%x_d, u%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(v%x, v%x_d, v%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(w%x, w%x_d, w%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(p%x, p%x_d, p%dof%size(), DEVICE_TO_HOST, sync=.false.)
         end if

         if (associated(this%ulag) .and. associated(this%vlag) .and. &
              associated(this%wlag)) then
            call device_memcpy(ulag%lf(1)%x, ulag%lf(1)%x_d, &
                 u%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(ulag%lf(2)%x, ulag%lf(2)%x_d, &
                 u%dof%size(), DEVICE_TO_HOST, sync=.false.)

            call device_memcpy(vlag%lf(1)%x, vlag%lf(1)%x_d, &
                 v%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(vlag%lf(2)%x, vlag%lf(2)%x_d, &
                 v%dof%size(), DEVICE_TO_HOST, sync=.false.)

            call device_memcpy(wlag%lf(1)%x, wlag%lf(1)%x_d, &
                 w%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(wlag%lf(2)%x, wlag%lf(2)%x_d, &
                 w%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(this%abx1%x, this%abx1%x_d, &
                 w%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(this%abx2%x, this%abx2%x_d, &
                 w%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(this%aby1%x, this%aby1%x_d, &
                 w%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(this%aby2%x, this%aby2%x_d, &
                 w%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(this%abz1%x, this%abz1%x_d, &
                 w%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(this%abz2%x, this%abz2%x_d, &
                 w%dof%size(), DEVICE_TO_HOST, sync=.false.)
         end if
         if (associated(this%s)) then
            call device_memcpy(this%s%x, this%s%x_d, &
                 this%s%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(this%slag%lf(1)%x, this%slag%lf(1)%x_d, &
                 this%s%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(this%slag%lf(2)%x, this%slag%lf(2)%x_d, &
                 this%s%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(this%abs1%x, this%abs1%x_d, &
                 w%dof%size(), DEVICE_TO_HOST, sync=.false.)
            call device_memcpy(this%abs2%x, this%abs2%x_d, &
                 w%dof%size(), DEVICE_TO_HOST, sync=.false.)
         end if

         ! Multi-scalar lag field synchronization
         if (allocated(this%scalar_lags%items) .and. this%scalar_lags%size() > 0) then
            do i = 1, this%scalar_lags%size()
               block
                 type(field_series_t), pointer :: slag
                 integer :: slag_size, dof_size
                 slag => this%scalar_lags%get(i)
                 slag_size = slag%size()
                 dof_size = slag%f%dof%size()
                 do j = 1, slag_size
                    call device_memcpy(slag%lf(j)%x, slag%lf(j)%x_d, &
                         dof_size, DEVICE_TO_HOST, sync=.false.)
                 end do
               end block
            end do
         end if

         ! Multi-scalar ABX field synchronization
         if (allocated(this%scalar_abx1) .and. allocated(this%scalar_abx2)) then
            do i = 1, size(this%scalar_abx1)
               call device_memcpy(this%scalar_abx1(i)%ptr%x, this%scalar_abx1(i)%ptr%x_d, &
                    this%scalar_abx1(i)%ptr%dof%size(), DEVICE_TO_HOST, sync=.false.)
               call device_memcpy(this%scalar_abx2(i)%ptr%x, this%scalar_abx2(i)%ptr%x_d, &
                    this%scalar_abx2(i)%ptr%dof%size(), DEVICE_TO_HOST, sync=.false.)
            end do
         end if
       end associate
       call device_sync(glb_cmd_queue)
    end if

  end subroutine chkp_sync_host

  !> Synchronize device with checkpoint
  subroutine chkp_sync_device(this)
    class(chkp_t), intent(inout) :: this
    integer :: i, j

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(u=>this%u, v=>this%v, w=>this%w, &
            ulag=>this%ulag, vlag=>this%vlag, wlag=>this%wlag,&
            p=>this%p)

         if (associated(this%u) .and. associated(this%v) .and. &
              associated(this%w)) then
            call device_memcpy(u%x, u%x_d, u%dof%size(), &
                 HOST_TO_DEVICE, sync=.false.)
            call device_memcpy(v%x, v%x_d, v%dof%size(), &
                 HOST_TO_DEVICE, sync=.false.)
            call device_memcpy(w%x, w%x_d, w%dof%size(), &
                 HOST_TO_DEVICE, sync=.false.)
            call device_memcpy(p%x, p%x_d, p%dof%size(), &
                 HOST_TO_DEVICE, sync=.false.)
         end if

         if (associated(this%ulag) .and. associated(this%vlag) .and. &
              associated(this%wlag)) then
            call device_memcpy(ulag%lf(1)%x, ulag%lf(1)%x_d, u%dof%size(), &
                 HOST_TO_DEVICE, sync=.false.)
            call device_memcpy(ulag%lf(2)%x, ulag%lf(2)%x_d, u%dof%size(), &
                 HOST_TO_DEVICE, sync=.false.)

            call device_memcpy(vlag%lf(1)%x, vlag%lf(1)%x_d, v%dof%size(), &
                 HOST_TO_DEVICE, sync=.false.)
            call device_memcpy(vlag%lf(2)%x, vlag%lf(2)%x_d, v%dof%size(), &
                 HOST_TO_DEVICE, sync=.false.)

            call device_memcpy(wlag%lf(1)%x, wlag%lf(1)%x_d, w%dof%size(), &
                 HOST_TO_DEVICE, sync=.false.)
            call device_memcpy(wlag%lf(2)%x, wlag%lf(2)%x_d, w%dof%size(), &
                 HOST_TO_DEVICE, sync=.false.)
         end if
         if (associated(this%s)) then
            call device_memcpy(this%s%x, this%s%x_d, this%s%dof%size(), &
                 HOST_TO_DEVICE, sync=.false.)

            call device_memcpy(this%slag%lf(1)%x, this%slag%lf(1)%x_d, &
                 this%s%dof%size(), HOST_TO_DEVICE, sync=.false.)
            call device_memcpy(this%slag%lf(2)%x, this%slag%lf(2)%x_d, &
                 this%s%dof%size(), HOST_TO_DEVICE, sync=.false.)
            call device_memcpy(this%abs1%x, this%abs1%x_d, &
                 w%dof%size(), HOST_TO_DEVICE, sync=.false.)
            call device_memcpy(this%abs2%x, this%abs2%x_d, &
                 w%dof%size(), HOST_TO_DEVICE, sync=.false.)
         end if

         ! Multi-scalar lag field synchronization
         if (allocated(this%scalar_lags%items) .and. this%scalar_lags%size() > 0) then
            do i = 1, this%scalar_lags%size()
               block
                 type(field_series_t), pointer :: slag
                 integer :: slag_size, dof_size
                 slag => this%scalar_lags%get(i)
                 slag_size = slag%size()
                 dof_size = slag%f%dof%size()
                 do j = 1, slag_size
                    call device_memcpy(slag%lf(j)%x, slag%lf(j)%x_d, &
                         dof_size, HOST_TO_DEVICE, sync=.false.)
                 end do
               end block
            end do
         end if

         ! Multi-scalar ABX field synchronization
         if (allocated(this%scalar_abx1) .and. allocated(this%scalar_abx2)) then
            do i = 1, size(this%scalar_abx1)
               call device_memcpy(this%scalar_abx1(i)%ptr%x, this%scalar_abx1(i)%ptr%x_d, &
                    this%scalar_abx1(i)%ptr%dof%size(), HOST_TO_DEVICE, sync=.false.)
               call device_memcpy(this%scalar_abx2(i)%ptr%x, this%scalar_abx2(i)%ptr%x_d, &
                    this%scalar_abx2(i)%ptr%dof%size(), HOST_TO_DEVICE, sync=.false.)
            end do
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



  !> Add a scalar to checkpointing
  subroutine chkp_add_scalar(this, s, slag, abs1, abs2)
    class(chkp_t), intent(inout) :: this
    type(field_t), target, intent(in) :: s
    type(field_series_t), target, intent(in) :: slag
    type(field_t), target, intent(in), optional :: abs1, abs2

    this%s => s
    this%slag => slag

    if (present(abs1)) this%abs1 => abs1
    if (present(abs2)) this%abs2 => abs2

  end subroutine chkp_add_scalar


  !> Return restart time from a loaded checkpoint
  pure function chkp_restart_time(this) result(rtime)
    class(chkp_t), intent(in) :: this
    real(kind=dp) :: rtime

    rtime = this%t
  end function chkp_restart_time

end module checkpoint
