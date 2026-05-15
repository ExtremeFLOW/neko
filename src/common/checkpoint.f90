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
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp, dp
  use field_series, only : field_series_t
  use field_series_list, only : field_series_list_t
  use space, only : space_t, operator(.ne.)
  use device, only : device_memcpy, DEVICE_TO_HOST, HOST_TO_DEVICE, &
       device_sync, glb_cmd_queue
  use field, only : field_t, field_ptr_t
  use utils, only : neko_error
  use mesh, only : mesh_t
  use math, only : NEKO_EPS
  use global_interpolation, only : GLOB_INTERP_TOL
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

     real(kind=dp) :: t = 0d0 !< Restart time (valid after load)
     type(mesh_t) :: previous_mesh
     type(space_t) :: previous_Xh
     real(kind=dp) :: mesh2mesh_tol = GLOB_INTERP_TOL

     !> ALE fields
     type(field_t), pointer :: wm_x => null()
     type(field_t), pointer :: wm_y => null()
     type(field_t), pointer :: wm_z => null()
     type(field_series_t), pointer :: wm_x_lag => null()
     type(field_series_t), pointer :: wm_y_lag => null()
     type(field_series_t), pointer :: wm_z_lag => null()
     real(kind=rp), pointer :: msh_x(:,:,:,:) => null()
     real(kind=rp), pointer :: msh_y(:,:,:,:) => null()
     real(kind=rp), pointer :: msh_z(:,:,:,:) => null()
     real(kind=rp), pointer :: pivot_pos(:) => null()
     real(kind=rp), pointer :: pivot_vel_lag(:,:) => null()
     real(kind=rp), pointer :: Blag(:,:,:,:) => null()
     real(kind=rp), pointer :: Blaglag(:,:,:,:) => null()
     real(kind=rp), pointer :: basis_pos(:) => null()
     real(kind=rp), pointer :: basis_vel_lag(:,:) => null()
   contains
     procedure, pass(this) :: init => chkp_init
     procedure, pass(this) :: sync_host => chkp_sync_host
     procedure, pass(this) :: sync_device => chkp_sync_device
     procedure, pass(this) :: add_fluid => chkp_add_fluid
     procedure, pass(this) :: add_lag => chkp_add_lag
     procedure, pass(this) :: add_scalar => chkp_add_scalar
     procedure, pass(this) :: add_ale => chkp_add_ale
     procedure, pass(this) :: restart_time => chkp_restart_time
     procedure, pass(this) :: free => chkp_free
  end type chkp_t

contains

  !> Initialize checkpoint structure with mandatory data
  subroutine chkp_init(this)
    class(chkp_t), intent(inout) :: this

    ! Make sure the object is clean
    call this%free()

  end subroutine chkp_init

  !> Reset checkpoint
  subroutine chkp_free(this)
    class(chkp_t), intent(inout) :: this

    this%t = 0d0

    if (associated(this%u)) nullify(this%u)
    if (associated(this%v)) nullify(this%v)
    if (associated(this%w)) nullify(this%w)
    if (associated(this%p)) nullify(this%p)

    if (associated(this%ulag)) nullify(this%ulag)
    if (associated(this%vlag)) nullify(this%vlag)
    if (associated(this%wlag)) nullify(this%wlag)

    if (associated(this%tlag)) nullify(this%tlag)
    if (associated(this%dtlag)) nullify(this%dtlag)

    ! ALE cleanup
    if (associated(this%wm_x)) nullify(this%wm_x)
    if (associated(this%wm_y)) nullify(this%wm_y)
    if (associated(this%wm_z)) nullify(this%wm_z)
    if (associated(this%wm_x_lag)) nullify(this%wm_x_lag)
    if (associated(this%wm_y_lag)) nullify(this%wm_y_lag)
    if (associated(this%wm_z_lag)) nullify(this%wm_z_lag)
    if (associated(this%basis_vel_lag)) nullify(this%basis_vel_lag)
    if (associated(this%msh_x)) nullify(this%msh_x)
    if (associated(this%msh_y)) nullify(this%msh_y)
    if (associated(this%msh_z)) nullify(this%msh_z)
    if (associated(this%pivot_pos)) nullify(this%pivot_pos)
    if (associated(this%pivot_vel_lag)) nullify(this%pivot_vel_lag)
    if (associated(this%Blag)) nullify(this%Blag)
    if (associated(this%Blaglag)) nullify(this%Blaglag)
    if (associated(this%basis_pos)) nullify(this%basis_pos)

    if (associated(this%abx1)) nullify(this%abx1)
    if (associated(this%abx2)) nullify(this%abx2)
    if (associated(this%aby1)) nullify(this%aby1)
    if (associated(this%aby2)) nullify(this%aby2)
    if (associated(this%abz1)) nullify(this%abz1)
    if (associated(this%abz2)) nullify(this%abz2)

    if (associated(this%s)) nullify(this%s)
    if (associated(this%slag)) nullify(this%slag)
    if (associated(this%abs1)) nullify(this%abs1)
    if (associated(this%abs2)) nullify(this%abs2)

    call this%scalar_lags%free()

    if (allocated(this%scalar_abx1)) then
       deallocate(this%scalar_abx1)
    end if

    if (allocated(this%scalar_abx2)) then
       deallocate(this%scalar_abx2)
    end if

    call this%previous_mesh%free()
    call this%previous_Xh%free()

  end subroutine chkp_free

  !> Synchronize checkpoint with device
  subroutine chkp_sync_host(this)
    class(chkp_t), intent(inout) :: this
    integer :: i, j

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(u => this%u, v => this%v, w => this%w, &
            ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
            p => this%p)

         if (associated(this%u) .and. associated(this%v) .and. &
              associated(this%w) .and. associated(this%p)) then
            call u%copy_from(DEVICE_TO_HOST, sync = .false.)
            call v%copy_from(DEVICE_TO_HOST, sync = .false.)
            call w%copy_from(DEVICE_TO_HOST, sync = .false.)
            call p%copy_from(DEVICE_TO_HOST, sync = .false.)
         end if

         if (associated(this%ulag) .and. associated(this%vlag) .and. &
              associated(this%wlag)) then
            call ulag%lf(1)%copy_from(DEVICE_TO_HOST, sync = .false.)
            call ulag%lf(2)%copy_from(DEVICE_TO_HOST, sync = .false.)

            call vlag%lf(1)%copy_from(DEVICE_TO_HOST, sync = .false.)
            call vlag%lf(2)%copy_from(DEVICE_TO_HOST, sync = .false.)

            call wlag%lf(1)%copy_from(DEVICE_TO_HOST, sync = .false.)
            call wlag%lf(2)%copy_from(DEVICE_TO_HOST, sync = .false.)
            call this%abx1%copy_from(DEVICE_TO_HOST, sync = .false.)
            call this%abx2%copy_from(DEVICE_TO_HOST, sync = .false.)
            call this%aby1%copy_from(DEVICE_TO_HOST, sync = .false.)
            call this%aby2%copy_from(DEVICE_TO_HOST, sync = .false.)
            call this%abz1%copy_from(DEVICE_TO_HOST, sync = .false.)
            call this%abz2%copy_from(DEVICE_TO_HOST, sync = .false.)
         end if

         if (associated(this%s)) then
            call this%s%copy_from(DEVICE_TO_HOST, sync = .false.)
            call this%slag%lf(1)%copy_from(DEVICE_TO_HOST, sync = .false.)
            call this%slag%lf(2)%copy_from(DEVICE_TO_HOST, sync = .false.)
            call this%abs1%copy_from(DEVICE_TO_HOST, sync = .false.)
            call this%abs2%copy_from(DEVICE_TO_HOST, sync = .false.)
         end if

         ! ALE field synchronization
         if (associated(this%wm_x) .and. associated(this%wm_y) .and. &
              associated(this%wm_z)) then
            call this%wm_x%copy_from(DEVICE_TO_HOST, sync = .false.)
            call this%wm_y%copy_from(DEVICE_TO_HOST, sync = .false.)
            call this%wm_z%copy_from(DEVICE_TO_HOST, sync = .false.)

            if (associated(this%wm_x_lag) .and. associated(this%wm_y_lag) &
                 .and. associated(this%wm_z_lag)) then

               call this%wm_x_lag%lf(1)%copy_from(DEVICE_TO_HOST, &
                    sync = .false.)
               call this%wm_x_lag%lf(2)%copy_from(DEVICE_TO_HOST, &
                    sync = .false.)

               call this%wm_y_lag%lf(1)%copy_from(DEVICE_TO_HOST, &
                    sync = .false.)
               call this%wm_y_lag%lf(2)%copy_from(DEVICE_TO_HOST, &
                    sync = .false.)

               call this%wm_z_lag%lf(1)%copy_from(DEVICE_TO_HOST, &
                    sync = .false.)
               call this%wm_z_lag%lf(2)%copy_from(DEVICE_TO_HOST, &
                    sync = .false.)
            end if
         end if

         ! Multi-scalar lag field synchronization
         do i = 1, this%scalar_lags%size()
            block
              type(field_series_t), pointer :: slag
              integer :: slag_size
              slag => this%scalar_lags%get(i)
              slag_size = slag%size()
              do j = 1, slag_size
                 call slag%lf(j)%copy_from(DEVICE_TO_HOST, sync = .false.)
              end do
            end block
         end do

         ! Multi-scalar ABX field synchronization
         if (allocated(this%scalar_abx1) .and. allocated(this%scalar_abx2)) then
            do i = 1, size(this%scalar_abx1)
               call this%scalar_abx1(i)%ptr%copy_from(DEVICE_TO_HOST, &
                    sync = .false.)
               call this%scalar_abx2(i)%ptr%copy_from(DEVICE_TO_HOST, &
                    sync = .false.)
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
       associate(u => this%u, v => this%v, w => this%w, &
            ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
            p => this%p)

         if (associated(this%u) .and. associated(this%v) .and. &
              associated(this%w)) then
            call u%copy_from(HOST_TO_DEVICE, sync = .false.)
            call v%copy_from(HOST_TO_DEVICE, sync = .false.)
            call w%copy_from(HOST_TO_DEVICE, sync = .false.)
            call p%copy_from(HOST_TO_DEVICE, sync = .false.)
         end if

         if (associated(this%ulag) .and. associated(this%vlag) .and. &
              associated(this%wlag)) then
            call ulag%lf(1)%copy_from(HOST_TO_DEVICE, sync = .false.)
            call ulag%lf(2)%copy_from(HOST_TO_DEVICE, sync = .false.)

            call vlag%lf(1)%copy_from(HOST_TO_DEVICE, sync = .false.)
            call vlag%lf(2)%copy_from(HOST_TO_DEVICE, sync = .false.)

            call wlag%lf(1)%copy_from(HOST_TO_DEVICE, sync = .false.)
            call wlag%lf(2)%copy_from(HOST_TO_DEVICE, sync = .false.)
         end if

         if (associated(this%s)) then
            call this%s%copy_from(HOST_TO_DEVICE, sync = .false.)

            call this%slag%lf(1)%copy_from(HOST_TO_DEVICE, sync = .false.)
            call this%slag%lf(2)%copy_from(HOST_TO_DEVICE, sync = .false.)
            call this%abs1%copy_from(HOST_TO_DEVICE, sync = .false.)
            call this%abs2%copy_from(HOST_TO_DEVICE, sync = .false.)
         end if

         ! ALE field synchronization
         if (associated(this%wm_x) .and. associated(this%wm_y) .and. &
              associated(this%wm_z)) then
            call this%wm_x%copy_from(HOST_TO_DEVICE, sync = .false.)
            call this%wm_y%copy_from(HOST_TO_DEVICE, sync = .false.)
            call this%wm_z%copy_from(HOST_TO_DEVICE, sync = .false.)
            if (associated(this%wm_x_lag) .and. associated(this%wm_y_lag) .and. &
                associated(this%wm_z_lag)) then
               call this%wm_x_lag%lf(1)%copy_from(HOST_TO_DEVICE, &
                    sync = .false.)
               call this%wm_x_lag%lf(2)%copy_from(HOST_TO_DEVICE, &
                    sync = .false.)

               call this%wm_y_lag%lf(1)%copy_from(HOST_TO_DEVICE, &
                    sync = .false.)
               call this%wm_y_lag%lf(2)%copy_from(HOST_TO_DEVICE, &
                    sync = .false.)

               call this%wm_z_lag%lf(1)%copy_from(HOST_TO_DEVICE, &
                    sync = .false.)
               call this%wm_z_lag%lf(2)%copy_from(HOST_TO_DEVICE, &
                    sync = .false.)
            end if
         end if

         ! Multi-scalar lag field synchronization
         if (allocated(this%scalar_lags%items) .and. &
              this%scalar_lags%size() > 0) then
            do i = 1, this%scalar_lags%size()
               block
                 type(field_series_t), pointer :: slag
                 integer :: slag_size, dof_size
                 slag => this%scalar_lags%get(i)
                 slag_size = slag%size()
                 dof_size = slag%f%dof%size()
                 do j = 1, slag_size
                    call slag%lf(j)%copy_from(HOST_TO_DEVICE, sync = .false.)
                 end do
               end block
            end do
         end if

         ! Multi-scalar ABX field synchronization
         if (allocated(this%scalar_abx1) .and. allocated(this%scalar_abx2)) then
            do i = 1, size(this%scalar_abx1)
               call this%scalar_abx1(i)%ptr%copy_from(HOST_TO_DEVICE, &
                    sync = .false.)
               call this%scalar_abx2(i)%ptr%copy_from(HOST_TO_DEVICE, &
                    sync = .false.)
            end do
         end if
       end associate
    end if

  end subroutine chkp_sync_device

  !> Add a fluid to the checkpoint
  subroutine chkp_add_fluid(this, u, v, w, p)
    class(chkp_t), intent(inout) :: this
    type(field_t), target :: u
    type(field_t), target :: v
    type(field_t), target :: w
    type(field_t), target :: p

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

  end subroutine chkp_add_fluid

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

  !> Add mesh velocity and other required variables to checkpointing for ALE
  subroutine chkp_add_ale(this, x, y, z, Blag, Blaglag, wm_x, wm_y, wm_z, &
       wm_x_lag, wm_y_lag, wm_z_lag, pivot_pos, pivot_vel_lag, basis_pos, &
       basis_vel_lag)
    class(chkp_t), intent(inout) :: this
    type(field_t), target, intent(in) :: wm_x, wm_y, wm_z
    real(kind=rp), intent(in), pointer :: pivot_pos(:), pivot_vel_lag(:,:)
    type(field_series_t), target, intent(in) :: wm_x_lag, wm_y_lag, wm_z_lag
    real(kind=rp), intent(in), pointer :: x(:,:,:,:), y(:,:,:,:), z(:,:,:,:)
    real(kind=rp), pointer, intent(in) :: Blag(:,:,:,:), Blaglag(:,:,:,:)
    real(kind=rp), intent(in), pointer :: basis_pos(:)
    real(kind=rp), intent(in), pointer :: basis_vel_lag(:,:)

    this%msh_x => x
    this%msh_y => y
    this%msh_z => z
    this%wm_x => wm_x
    this%wm_y => wm_y
    this%wm_z => wm_z
    this%wm_x_lag => wm_x_lag
    this%wm_y_lag => wm_y_lag
    this%wm_z_lag => wm_z_lag
    this%Blag => Blag
    this%Blaglag => Blaglag
    this%pivot_pos => pivot_pos
    this%pivot_vel_lag => pivot_vel_lag
    this%basis_pos => basis_pos
    this%basis_vel_lag => basis_vel_lag
  end subroutine chkp_add_ale

  !> Return restart time from a loaded checkpoint
  pure function chkp_restart_time(this) result(rtime)
    class(chkp_t), intent(in) :: this
    real(kind=dp) :: rtime

    rtime = this%t
  end function chkp_restart_time

end module checkpoint
