! Copyright (c) 2021-2022, The Neko Authors
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
  use speclib
  use device
  use utils
  use math
  use fast3d
  use tensor
  use tensor_cpu
  use space
  use device
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
     !> First space
     type(space_t), pointer :: Xh
     !> Second space
     type(space_t), pointer :: Yh
     !> Interpolation weights from Xh to Yh
     real(kind=rp), allocatable :: Xh_to_Yh(:,:), Xh_to_YhT(:,:)
     !> Interpolation weights from Yh to Xh
     real(kind=rp), allocatable :: Yh_to_Xh(:,:), Yh_to_XhT(:,:)
     !> Device pointer for Xh_to_Yh
     type(c_ptr) :: Xh_Yh_d = C_NULL_PTR
     !> Device pointer for Xh_to_YhT
     type(c_ptr) :: Xh_YhT_d = C_NULL_PTR
     !> Device pointer for Yh_to_Xh
     type(c_ptr) :: Yh_Xh_d = C_NULL_PTR
     !> Device pointer for Yh_to_XhT
     type(c_ptr) :: Yh_XhT_d = C_NULL_PTR

   contains
     !> Constructor
     procedure, pass(this) :: init => interpolator_init
     !> Destructor
     procedure, pass(this) :: free => interpolator_free
     !> Interpolate an array to one of Xh or Yh.
     procedure, pass(this) :: map => interpolator_map
     !> Interpolate an array to one of Xh or Yh on the host.
     procedure, pass(this) :: map_host => interpolator_map_host
  end type interpolator_t
  
contains
  
  !> Constructor 
  !> @param Xh The first space
  !> @param Xh The second space
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
       call device_map(this%Xh_to_Yh, this%Xh_Yh_d, Yh%lx*Xh%lx)
       call device_map(this%Xh_to_YhT, this%Xh_YhT_d, Yh%lx*Xh%lx)
       call device_map(this%Yh_to_Xh, this%Yh_Xh_d, Yh%lx*Xh%lx)
       call device_map(this%Yh_to_XhT, this%Yh_XhT_d, Yh%lx*Xh%lx)
       call device_memcpy(this%Xh_to_Yh, this%Xh_Yh_d, Yh%lx*Xh%lx, HOST_TO_DEVICE)
       call device_memcpy(this%Xh_to_YhT, this%Xh_YhT_d, Yh%lx*Xh%lx, HOST_TO_DEVICE)
       call device_memcpy(this%Yh_to_Xh, this%Yh_Xh_d, Yh%lx*Xh%lx, HOST_TO_DEVICE)
       call device_memcpy(this%Yh_to_XhT, this%Yh_XhT_d, Yh%lx*Xh%lx, HOST_TO_DEVICE)
    end if

  end subroutine interpolator_init

  subroutine interpolator_free(this)
    class(interpolator_t), intent(inout) :: this

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
    if (c_associated(this%Yh_Xh_d)) then
       call device_free(this%Yh_Xh_d)
    end if
    if (c_associated(this%Yh_XhT_d)) then
       call device_free(this%Yh_XhT_d)
    end if
    if (c_associated(this%Xh_Yh_d)) then
       call device_free(this%Xh_Yh_d)
    end if
    if (c_associated(this%Xh_YhT_d)) then
       call device_free(this%Xh_YhT_d)
    end if

  end subroutine interpolator_free

  !> Interpolates an array to one of Xh or Yh.
  !! @param x Original array.
  !! @param y Interpolated array.
  !! @param nel Number of elements in the mesh.
  !! @param to_space The space to interpolate to, must be either Xh or Yh.
  subroutine interpolator_map(this, y, x, nel, to_space)
    class(interpolator_t), intent(inout) :: this
    integer :: nel
    type(space_t) :: to_space
    real(kind=rp), intent(inout) :: x(1,nel)
    real(kind=rp), intent(inout) :: y(1,nel)
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
  end subroutine interpolator_map

  !> Interpolates an array to one of Xh or Yh on host.
  !! @param x Original array.
  !! @param y Interpolated array.
  !! @param nel Number of elements in the mesh.
  !! @param to_space The space to interpolate to, must be either Xh or Yh.
  subroutine interpolator_map_host(this, y, x, nel, to_space)
    class(interpolator_t), intent(inout) :: this
    integer :: nel
    type(space_t) :: to_space
    real(kind=rp), intent(inout) :: x(1,nel)
    real(kind=rp), intent(inout) :: y(1,nel)
    if (to_space .eq. this%Yh) then
       call tnsr3d_cpu(y, this%Yh%lx, x, &
                   this%Xh%lx,this%Yh_to_XhT, &
                   this%Yh_to_Xh, this%Yh_to_Xh, nel)
    else if (to_space .eq. this%Xh) then
       call tnsr3d_cpu(y, this%Xh%lx, x, &
                   this%Yh%lx,this%Yh_to_Xh, &
                   this%Yh_to_XhT, this%Yh_to_XhT, nel)
    else
       call neko_error('Invalid interpolation')
    end if
  end subroutine interpolator_map_host


end module interpolation
