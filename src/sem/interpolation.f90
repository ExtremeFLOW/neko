! Copyright (c) 2021-2023, The Neko Authors
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
!> Routines to interpolate between different spaces
module interpolation
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use device, only : device_map, device_unmap, device_memcpy, HOST_TO_DEVICE
  use fast3d, only : setup_intp
  use field, only : field_t
  use tensor, only : tnsr3d
  use tensor_cpu, only : tnsr3d_cpu
  use tensor_device, only : tnsr3d_device
  use space, only : space_t, operator(.eq.), GL, GLL
  use utils, only : neko_error
  use logger, only : neko_log
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Interpolation between two \ref space::space_t.
  !! @details
  !! This type implements functionality to interpolate between a pair of spaces.
  !! Simply put, given some data of form (lx1, lx1, lx1, nelem) we can map it to
  !! (lx2, lx2, lx2, nelem), corresponding to a different polynomial order in
  !! each element.
  type, public :: interpolator_t
     !> First space.
     type(space_t), pointer :: Xh
     !> Second space.
     type(space_t), pointer :: Yh
     !> Interpolation weights from Xh to Yh.
     real(kind=rp), allocatable :: Xh_to_Yh(:,:), Xh_to_YhT(:,:)
     !> Interpolation weights from Yh to Xh.
     real(kind=rp), allocatable :: Yh_to_Xh(:,:), Yh_to_XhT(:,:)
     !> Device pointer for Xh_to_Yh.
     type(c_ptr) :: Xh_to_Yh_d = C_NULL_PTR
     !> Device pointer for Xh_to_YhT.
     type(c_ptr) :: Xh_to_YhT_d = C_NULL_PTR
     !> Device pointer for Yh_to_Xh.
     type(c_ptr) :: Yh_to_Xh_d = C_NULL_PTR
     !> Device pointer for Yh_to_XhT.
     type(c_ptr) :: Yh_to_XhT_d = C_NULL_PTR
   contains

     !> Constructor.
     procedure, pass(this) :: init => interpolator_init
     !> Destructor.
     procedure, pass(this) :: free => interpolator_free

     !> Generic interface to map an array to one of Xh or Yh.
     generic :: map => map_old, map_device, map_field
     !> Interpolate an array to one of Xh or Yh on the host.
     procedure, pass(this) :: map_host => interpolator_map_host
     !> Interpolate an array to one of Xh or Yh.
     procedure, pass(this) :: map_device => interpolator_map_device
     !> Interpolate an array to one of Xh or Yh.
     procedure, pass(this) :: map_field => interpolator_map_field

     !> Interpolate an array to one of Xh or Yh.
     !! @deprecated This is an older version of the interpolation routine, still
     !! relying on the implicit device mapping.
     procedure, pass(this) :: map_old => interpolator_map_old

  end type interpolator_t

contains

  !> Constructor to initialize with two different spaces.
  !> @param Xh The first space.
  !> @param Xh The second space.
  subroutine interpolator_init(this, Xh, Yh)
    class(interpolator_t), intent(inout), target :: this
    type(space_t), intent(inout), target :: Xh
    type(space_t), intent(inout), target :: Yh
    integer :: deg_derivate

    call this%free()

    allocate(this%Xh_to_Yh(Yh%lx,Xh%lx))
    allocate(this%Xh_to_YhT(Xh%lx,Yh%lx))
    allocate(this%Yh_to_Xh(Xh%lx,Yh%lx))
    allocate(this%Yh_to_XhT(Yh%lx,Xh%lx))

    if (Xh%t .eq. GLL .and. Yh%t .eq. GLL) then
    else if ((Xh%t .eq. GL .and. Yh%t .eq. GLL) .or. &
         (Yh%t .eq. GL .and. Xh%t .eq. GLL)) then
    else
       call neko_error('Unsupported interpolation')
    end if
    deg_derivate = 0
    call setup_intp(this%Xh_to_Yh, this%Xh_to_YhT, &
         Yh%zg, Xh%zg, Yh%lx, Xh%lx, deg_derivate)
    call setup_intp(this%Yh_to_Xh, this%Yh_to_XhT, &
         Xh%zg, Yh%zg, Xh%lx, Yh%lx, deg_derivate)

    this%Xh => Xh
    this%Yh => Yh
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%Xh_to_Yh, this%Xh_to_Yh_d, Yh%lx*Xh%lx)
       call device_map(this%Xh_to_YhT, this%Xh_to_YhT_d, Yh%lx*Xh%lx)
       call device_map(this%Yh_to_Xh, this%Yh_to_Xh_d, Yh%lx*Xh%lx)
       call device_map(this%Yh_to_XhT, this%Yh_to_XhT_d, Yh%lx*Xh%lx)

       call device_memcpy(this%Xh_to_Yh, this%Xh_to_Yh_d, Yh%lx*Xh%lx, &
            HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(this%Xh_to_YhT, this%Xh_to_YhT_d, Yh%lx*Xh%lx, &
            HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(this%Yh_to_Xh, this%Yh_to_Xh_d, Yh%lx*Xh%lx, &
            HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(this%Yh_to_XhT, this%Yh_to_XhT_d, Yh%lx*Xh%lx, &
            HOST_TO_DEVICE, sync=.true.)
    end if

  end subroutine interpolator_init

  subroutine interpolator_free(this)
    class(interpolator_t), intent(inout) :: this

    if (c_associated(this%Yh_to_Xh_d)) then
       call device_unmap(this%Yh_to_Xh, this%Yh_to_Xh_d)
    end if
    if (c_associated(this%Yh_to_XhT_d)) then
       call device_unmap(this%Yh_to_XhT, this%Yh_to_XhT_d)
    end if
    if (c_associated(this%Xh_to_Yh_d)) then
       call device_unmap(this%Xh_to_Yh, this%Xh_to_Yh_d)
    end if
    if (c_associated(this%Xh_to_YhT_d)) then
       call device_unmap(this%Xh_to_YhT, this%Xh_to_YhT_d)
    end if

    if (allocated(this%Xh_to_Yh)) then
       deallocate(this%Xh_to_Yh)
    end if
    if (allocated(this%Xh_to_YhT)) then
       deallocate(this%Xh_to_YhT)
    end if
    if (allocated(this%Yh_to_Xh)) then
       deallocate(this%Yh_to_Xh)
    end if
    if (allocated(this%Yh_to_XhT)) then
       deallocate(this%Yh_to_XhT)
    end if

    nullify(this%Xh)
    nullify(this%Yh)

  end subroutine interpolator_free

  !> Interpolates an array to one of Xh or Yh.
  !! @param x Original array.
  !! @param y Interpolated array.
  !! @param nel Number of elements in the mesh.
  !! @param to_space The space to interpolate to, must be either Xh or Yh.
  subroutine interpolator_map_old(this, y, x, nel, to_space)
    class(interpolator_t), intent(in) :: this
    integer, intent(in) :: nel
    type(space_t), intent(in) :: to_space
    real(kind=rp), intent(in) :: x(this%Xh%lx, this%Xh%ly, this%Xh%lz, nel)
    real(kind=rp), intent(inout) :: y(this%Yh%lx, this%Yh%ly, this%Yh%lz, nel)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_log%deprecated('Interpolator: map implicit device', '2.0.0', &
            'Please call the generic map or direct interfaces instead.')
    end if

    if (to_space .eq. this%Yh) then
       call tnsr3d(y, this%Yh%lx, x, &
            this%Xh%lx,this%Yh_to_XhT, &
            this%Yh_to_Xh, this%Yh_to_Xh, nel)
    else if (to_space .eq. this%Xh) then
       call tnsr3d(y, this%Xh%lx, x, &
            this%Yh%lx,this%Yh_to_Xh, &
            this%Yh_to_XhT, this%Yh_to_XhT, nel)
    else
       call neko_error('Invalid interpolation')
    end if
  end subroutine interpolator_map_old

  !> Interpolates an array to one of Xh or Yh on host.
  !! @param x Original array.
  !! @param y Interpolated array.
  !! @param nel Number of elements in the mesh.
  !! @param to_space The space to interpolate to, must be either Xh or Yh.
  !! @note Not optimized for performance, should only be used during init.
  subroutine interpolator_map_host(this, y, x, nel, to_space)
    class(interpolator_t), intent(in) :: this
    integer, intent(in) :: nel
    type(space_t), intent(in) :: to_space
    real(kind=rp), intent(in) :: x(this%Xh%lx, this%Xh%ly, this%Xh%lz, nel)
    real(kind=rp), intent(inout) :: y(this%Yh%lx, this%Yh%ly, this%Yh%lz, nel)
    if (to_space .eq. this%Yh) then
       call tnsr3d_cpu(y, this%Yh%lx, x, this%Xh%lx, &
            this%Yh_to_XhT, this%Yh_to_Xh, this%Yh_to_Xh, nel)
    else if (to_space .eq. this%Xh) then
       call tnsr3d_cpu(y, this%Xh%lx, x, this%Yh%lx, &
            this%Yh_to_Xh, this%Yh_to_XhT, this%Yh_to_XhT, nel)
    else
       call neko_error('Invalid interpolation')
    end if
  end subroutine interpolator_map_host

  !> Interpolates a field to one of Xh or Yh.
  !! @param x Original field.
  !! @param y Interpolated field.
  !! @param nel Number of elements in the mesh.
  !! @param to_space The space to interpolate to, must be either Xh or Yh.
  subroutine interpolator_map_field(this, y, x, nel, to_space)
    class(interpolator_t), intent(in) :: this
    integer, intent(in) :: nel
    type(space_t), intent(in) :: to_space
    type(field_t), intent(in) :: x
    type(field_t), intent(inout) :: y
    if (to_space .eq. this%Yh) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call tnsr3d_device(y%x_d, this%Yh%lx, x%x_d, this%Xh%lx, &
               this%Yh_to_XhT_d, this%Yh_to_Xh_d, this%Yh_to_Xh_d, nel)
       else
          call tnsr3d_cpu(y%x, this%Yh%lx, x%x, this%Xh%lx, &
               this%Yh_to_XhT, this%Yh_to_Xh, this%Yh_to_Xh, nel)
       end if
    else if (to_space .eq. this%Xh) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call tnsr3d_device(y%x_d, this%Xh%lx, x%x_d, this%Yh%lx, &
               this%Yh_to_Xh_d, this%Yh_to_XhT_d, this%Yh_to_XhT_d, nel)
       else
          call tnsr3d_cpu(y%x, this%Xh%lx, x%x, this%Yh%lx, &
               this%Yh_to_Xh, this%Yh_to_XhT, this%Yh_to_XhT, nel)
       end if
    else
       call neko_error('Invalid interpolation')
    end if
  end subroutine interpolator_map_field

  !> Interpolates an array to one of Xh or Yh on device.
  !! @param x Original array.
  !! @param y Interpolated array.
  !! @param nel Number of elements in the mesh.
  !! @param to_space The space to interpolate to, must be either Xh or Yh.
  subroutine interpolator_map_device(this, y_d, x_d, nel, to_space)
    class(interpolator_t), intent(in) :: this
    integer, intent(in) :: nel
    type(space_t), intent(in) :: to_space
    type(c_ptr), intent(in) :: x_d
    type(c_ptr), intent(inout) :: y_d
    if (to_space .eq. this%Yh) then
       call tnsr3d_device(y_d, this%Yh%lx, x_d, this%Xh%lx, &
            this%Yh_to_XhT_d, this%Yh_to_Xh_d, this%Yh_to_Xh_d, nel)
    else if (to_space .eq. this%Xh) then
       call tnsr3d_device(y_d, this%Xh%lx, x_d, this%Yh%lx, &
            this%Yh_to_Xh_d, this%Yh_to_XhT_d, this%Yh_to_XhT_d, nel)
    else
       call neko_error('Invalid interpolation')
    end if
  end subroutine interpolator_map_device

end module interpolation
