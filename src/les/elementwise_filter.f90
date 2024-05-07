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
  use math
  use field, only : field_t
  use utils, only : neko_error
  use neko_config, only : NEKO_BCKND_DEVICE
  use elementwise_filter_cpu
  use tensor, only : tnsr3d
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

    this%filter_type = ""
    this%nx = 0
    this%nt = 0

  end subroutine elementwise_filter_free

  !> Build the 1d filter for an element.
  subroutine build_1d(this)
    class(elementwise_filter_t), intent(inout) :: this

    if (NEKO_BCKND_DEVICE .eq. 1) then
        call neko_error("build_1d not implemented on accelarators.")
    else
        call build_1d_cpu(this%fh, this%fht, this%trnsfr, &
                               this%nx, this%filter_type)
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