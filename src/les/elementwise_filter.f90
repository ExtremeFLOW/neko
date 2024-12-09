! Copyright (c) 2024, The Neko Authors
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
!
!> Implements `explicit_filter_t`.
module elementwise_filter
  use num_types, only : rp
  use math, only : rzero, rone
  use field, only : field_t
  use utils, only : neko_error
  use neko_config, only : NEKO_BCKND_DEVICE
  use elementwise_filter_cpu
  use tensor, only : tnsr3d
  use device, only : device_map, device_free, c_ptr, &
                    C_NULL_PTR, device_memcpy, HOST_TO_DEVICE
  use device_math, only : device_cfill
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Implements the explicit filter for SEM.
  type, public :: elementwise_filter_t
     !> filter type:
     !> possible options: "Boyd", "nonBoyd"
     character(len=64) :: filter_type
     !> dimension
     integer :: nx
     !> filtered wavenumber
     integer :: nt
     !> matrix for 1d elementwise filtering
     real(kind=rp), allocatable :: fh(:,:), fht(:,:)
     type(c_ptr) :: fh_d = C_NULL_PTR
     type(c_ptr) :: fht_d = C_NULL_PTR
     !> transfer function
     real(kind=rp), allocatable :: trnsfr(:)
   contains
     !> Constructor.
     procedure, pass(this) :: init => elementwise_filter_init
     !> Destructor.
     procedure, pass(this) :: free => elementwise_filter_free
     !> Set up 1D filter inside an element.
     procedure, pass(this) :: build_1d
     !> Filter a 3D field
     procedure, pass(this) :: filter_3d => elementwise_field_filter_3d
  end type elementwise_filter_t

contains
  !> Constructor.
  !! @param nx number of points in an elements in one direction.
  !! @param filter_type possible options: "Boyd", "nonBoyd"
  subroutine elementwise_filter_init(this, nx, filter_type)
    class(elementwise_filter_t), intent(inout) :: this
    character(len=*) :: filter_type
    integer :: nx
    
    this%nx = nx
    this%nt = nx ! initialize as if nothing is filtered yet 
    this%filter_type = filter_type

    allocate(this%fh(nx, nx))
    allocate(this%fht(nx, nx))
    allocate(this%trnsfr(nx))

    call rzero(this%fh, nx*nx)
    call rzero(this%fht, nx*nx)
    call rone(this%trnsfr, nx) ! initialize as if nothing is filtered yet

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%fh, this%fh_d, this%nx * this%nx)
       call device_map(this%fht, this%fht_d, this%nx * this%nx)
       call device_cfill(this%fh_d, 0.0_rp, this%nx * this%nx)
       call device_cfill(this%fht_d, 0.0_rp, this%nx * this%nx)
    end if
    
  end subroutine elementwise_filter_init

  !> Destructor.
  subroutine elementwise_filter_free(this)
    class(elementwise_filter_t), intent(inout) :: this

    if (allocated(this%fh)) then
       deallocate(this%fh)
    end if

    if (allocated(this%fht)) then
       deallocate(this%fht)
    end if

    if (allocated(this%trnsfr)) then
       deallocate(this%trnsfr)
    end if

    if (c_associated(this%fh_d)) then
       call device_free(this%fh_d)
    end if

    if (c_associated(this%fht_d)) then
       call device_free(this%fht_d)
    end if

    this%filter_type = ""
    this%nx = 0
    this%nt = 0

  end subroutine elementwise_filter_free

  !> Build the 1d filter for an element.
  subroutine build_1d(this)
    class(elementwise_filter_t), intent(inout) :: this

    call build_1d_cpu(this%fh, this%fht, this%trnsfr, &
                               this%nx, this%filter_type)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%fh, this%fh_d, &
                          this%nx * this%nx, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%fht, this%fht_d, &
                          this%nx * this%nx, HOST_TO_DEVICE, sync = .false.)
    end if

  end subroutine build_1d

  !> Filter a 3D field.
  subroutine elementwise_field_filter_3d(this, v, u, nelv)
    class(elementwise_filter_t), intent(inout) :: this
    integer, intent(inout) :: nelv
    real(kind=rp), intent(inout), dimension(this%nx, this%nx, this%nx, nelv) :: u,v

    ! v = fh x fh x fh x u
    call tnsr3d(v, this%nx, u, this%nx, this%fh, this%fht, this%fht, nelv)

  end subroutine elementwise_field_filter_3d

end module elementwise_filter
