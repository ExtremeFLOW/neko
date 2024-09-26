! Copyright (c) 2024, Gregor Weiss (HLRS)
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
!> Generic buffer that is extended with buffers of varying rank
module buffer_4d_npar
  use num_types
  use vector
  use adios2
  use buffer
  implicit none

  integer, private :: nthpar
  integer, private :: npar

  type, extends(buffer_t) :: buffer_4d_npar_t
     integer(kind=8), dimension(5) :: shape_dims, start_dims, count_dims
     real(kind=dp), private, allocatable :: data_dp(:,:,:,:,:)
     real(kind=sp), private, allocatable :: data_sp(:,:,:,:,:)
   contains
     procedure :: init => buffer_4d_npar_init
     procedure :: fill => buffer_4d_npar_fill
     procedure :: define => buffer_4d_npar_define
     procedure :: inquire => buffer_4d_npar_inquire
     procedure :: write => buffer_4d_npar_write
     procedure :: read => buffer_4d_npar_read
     procedure :: copy => buffer_4d_npar_copy
  end type buffer_4d_npar_t

contains

  subroutine buffer_4d_npar_init(this, precision, gdim, glb_nelv, offset_el, nelv, lx, ly, lz)
    class(buffer_4d_npar_t), intent(inout) :: this
    logical, intent(in) :: precision
    integer, intent(in) :: gdim, glb_nelv, offset_el, nelv, lx, ly, lz
    integer :: lxyz

    nthpar = 0
    npar = gdim
    lxyz = lx*ly*lz

    call buffer_set_precision(this, precision)

    if (this%dp_precision) then
       if (allocated(this%data_dp)) then
          deallocate(this%data_dp)
       end if
       allocate(this%data_dp(nelv, lx, ly, lz, npar))
    else
       if (allocated(this%data_sp)) then
          deallocate(this%data_sp)
       end if
       allocate(this%data_sp(nelv, lx, ly, lz, npar))
    end if

    this%shape_dims = [int(glb_nelv, i8), int(lx, i8), int(ly, i8), int(lz, i8), int(npar, i8)]
    this%start_dims = [int(offset_el, i8), int(0, i8), int(0, i8), int(0, i8), int(0, i8)]
    this%count_dims = [int(nelv, i8), int(lx, i8), int(ly, i8), int(lz, i8), int(npar, i8)]

  end subroutine buffer_4d_npar_init

  subroutine buffer_4d_npar_fill(this, x, n)
    class(buffer_4d_npar_t), intent(inout) :: this
    integer, intent(inout) :: n
    real(kind=rp), intent(inout) :: x(n)
    integer :: i, j, k, l, nelv, lx, ly, lz, index

    nthpar = nthpar + 1
    if (nthpar .le. npar) then

    nelv = this%count_dims(1)
    lx = this%count_dims(2)
    ly = this%count_dims(3)
    lz = this%count_dims(4)

    if (this%dp_precision) then
       do i = 1, nelv
          do j = 1, lz
             do k = 1, ly
                do l = 1, lx
                   index = (l-1) + lx*(k-1) + lx*ly*(j-1) + lx*ly*lz*(i-1) + 1
                   this%data_dp(i,l,k,j,nthpar) = real(x(index),dp)
                end do
             end do
          end do
       end do
    else
       do i = 1, nelv
          do j = 1, lz
             do k = 1, ly
                do l = 1, lx
                   index = (l-1) + lx*(k-1) + lx*ly*(j-1) + lx*ly*lz*(i-1) + 1
                   this%data_sp(i,l,k,j,nthpar) = real(x(index),sp)
                end do
             end do
          end do
       end do
    end if

    end if

  end subroutine buffer_4d_npar_fill

  subroutine buffer_4d_npar_define(this, variable, io, variable_name, ierr)
    class(buffer_4d_npar_t), intent(inout) :: this
    type(adios2_variable), intent(inout) :: variable
    type(adios2_io), intent(inout) :: io
    character(len=*), intent(in) :: variable_name
    integer, intent(inout) :: ierr
    integer :: adios2_type

    if (this%dp_precision) then
       adios2_type = adios2_type_dp
    else
       adios2_type = adios2_type_real
    end if

    call adios2_inquire_variable(variable, io, trim(variable_name), ierr)
    if (.not.variable%valid) then
       !> @todo could the shape and slice be fixed?
       call adios2_define_variable(variable, io, variable_name, adios2_type, &
            size(this%shape_dims), this%shape_dims, this%start_dims, &
            this%count_dims, .false., ierr)
    else
       call adios2_set_selection(variable, size(this%start_dims), &
            this%start_dims, this%count_dims, ierr)
    end if

  end subroutine buffer_4d_npar_define

  subroutine buffer_4d_npar_inquire(this, variable, io, variable_name, ierr)
    class(buffer_4d_npar_t), intent(inout) :: this
    type(adios2_variable), intent(inout) :: variable
    type(adios2_io), intent(inout) :: io
    character(len=*), intent(in) :: variable_name
    integer, intent(inout) :: ierr

    call adios2_inquire_variable(variable, io, trim(variable_name), ierr)
    if (variable%valid) then
       call adios2_set_selection(variable, size(this%start_dims), &
            this%start_dims, this%count_dims, ierr)
    end if

  end subroutine buffer_4d_npar_inquire

  subroutine buffer_4d_npar_write(this, engine, variable, ierr)
    class(buffer_4d_npar_t), intent(inout) :: this
    type(adios2_engine), intent(in) :: engine
    type(adios2_variable), intent(in) :: variable
    integer, intent(inout) :: ierr

    if (this%dp_precision) then
       call adios2_put(engine, variable, this%data_dp, adios2_mode_sync, ierr)
    else
       call adios2_put(engine, variable, this%data_sp, adios2_mode_sync, ierr)
    end if

  end subroutine buffer_4d_npar_write
  
  subroutine buffer_4d_npar_read(this, engine, variable, ierr)
    class(buffer_4d_npar_t), intent(inout) :: this
    type(adios2_engine), intent(in) :: engine
    type(adios2_variable), intent(in) :: variable
    integer, intent(inout) :: ierr

    if (this%dp_precision) then
       call adios2_get(engine, variable, this%data_dp, adios2_mode_sync, ierr)
    else
       call adios2_get(engine, variable, this%data_sp, adios2_mode_sync, ierr)
    end if

  end subroutine buffer_4d_npar_read

  subroutine buffer_4d_npar_copy(this, x)
    class(buffer_4d_npar_t), intent(inout) :: this
    type(vector_t), intent(inout) :: x
    integer :: i, j, k, l, nelv, lx, ly, lz, index

    nthpar = nthpar + 1
    if (nthpar .le. npar) then

    nelv = this%count_dims(1)
    lx = this%count_dims(2)
    ly = this%count_dims(3)
    lz = this%count_dims(4)

    if (this%dp_precision) then
       do i = 1, nelv
          do j = 1, lz
             do k = 1, ly
                do l = 1, lx
                   index = (l-1) + lx*(k-1) + lx*ly*(j-1) + lx*ly*lz*(i-1) + 1
                   x%x(index) = this%data_dp(i,l,k,j,nthpar)
                end do
             end do
          end do
       end do
    else
       do i = 1, nelv
          do j = 1, lz
             do k = 1, ly
                do l = 1, lx
                   index = (l-1) + lx*(k-1) + lx*ly*(j-1) + lx*ly*lz*(i-1) + 1
                   x%x(index) = this%data_sp(i,l,k,j,nthpar)
                end do
             end do
          end do
       end do
    end if

    end if

  end subroutine buffer_4d_npar_copy

end module buffer_4d_npar
