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
module buffer_1d
  use num_types
  use vector
#ifdef HAVE_ADIOS2_FORTRAN
  use adios2
#endif
  use buffer
  implicit none

  type, extends(buffer_t) :: buffer_1d_t
     integer(kind=8), dimension(1) :: shape_dims, start_dims, count_dims
     real(kind=dp), private, allocatable :: data_dp(:)
     real(kind=sp), private, allocatable :: data_sp(:)
   contains
     procedure :: init => buffer_1d_init
     procedure :: fill => buffer_1d_fill
#ifdef HAVE_ADIOS2_FORTRAN
     procedure :: define => buffer_1d_define
     procedure :: inquire => buffer_1d_inquire
     procedure :: write => buffer_1d_write
     procedure :: read => buffer_1d_read
#endif
     procedure :: copy => buffer_1d_copy
  end type buffer_1d_t

contains

  subroutine buffer_1d_init(this, precision, gdim, glb_nelv, offset_el, nelv, lx, ly, lz)
    class(buffer_1d_t), intent(inout) :: this
    logical, intent(in) :: precision
    integer, intent(in) :: gdim, glb_nelv, offset_el, nelv, lx, ly, lz
    integer :: lxyz

    lxyz = lx*ly*lz

    call buffer_set_precision(this, precision)

    if (this%dp_precision) then
       if (allocated(this%data_dp)) then
          deallocate(this%data_dp)
       end if
       allocate(this%data_dp(gdim*nelv*lxyz))
    else
       if (allocated(this%data_sp)) then
          deallocate(this%data_sp)
       end if
       allocate(this%data_sp(gdim*nelv*lxyz))
    end if

    this%shape_dims = (int(glb_nelv, i8) * int(lxyz, i8))
    this%start_dims = (int(offset_el, i8) * int(lxyz, i8))
    this%count_dims = (int(nelv, i8) * int(lxyz, i8))

  end subroutine buffer_1d_init

  subroutine buffer_1d_fill(this, x, n)
    class(buffer_1d_t), intent(inout) :: this
    integer, intent(inout) :: n
    real(kind=rp), intent(inout) :: x(n)
    integer :: i

    if (this%dp_precision) then
       do i = 1, n
          this%data_dp(i) = real(x(i),dp)
       end do
    else
       do i = 1, n
          this%data_sp(i) = real(x(i),sp)
       end do
    end if

  end subroutine buffer_1d_fill

#ifdef HAVE_ADIOS2_FORTRAN

  subroutine buffer_1d_define(this, variable, io, variable_name, ierr)
    class(buffer_1d_t), intent(inout) :: this
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

  end subroutine buffer_1d_define

  subroutine buffer_1d_inquire(this, variable, io, variable_name, ierr)
    class(buffer_1d_t), intent(inout) :: this
    type(adios2_variable), intent(inout) :: variable
    type(adios2_io), intent(inout) :: io
    character(len=*), intent(in) :: variable_name
    integer, intent(inout) :: ierr

    call adios2_inquire_variable(variable, io, trim(variable_name), ierr)
    if (variable%valid) then
       call adios2_set_selection(variable, size(this%start_dims), &
            this%start_dims, this%count_dims, ierr)
    end if

  end subroutine buffer_1d_inquire

  subroutine buffer_1d_write(this, engine, variable, ierr)
    class(buffer_1d_t), intent(inout) :: this
    type(adios2_engine), intent(in) :: engine
    type(adios2_variable), intent(in) :: variable
    integer, intent(inout) :: ierr

    if (this%dp_precision) then
       call adios2_put(engine, variable, this%data_dp, adios2_mode_sync, ierr)
    else
       call adios2_put(engine, variable, this%data_sp, adios2_mode_sync, ierr)
    end if

  end subroutine buffer_1d_write

  subroutine buffer_1d_read(this, engine, variable, ierr)
    class(buffer_1d_t), intent(inout) :: this
    type(adios2_engine), intent(in) :: engine
    type(adios2_variable), intent(in) :: variable
    integer, intent(inout) :: ierr

    if (this%dp_precision) then
       call adios2_get(engine, variable, this%data_dp, adios2_mode_sync, ierr)
    else
       call adios2_get(engine, variable, this%data_sp, adios2_mode_sync, ierr)
    end if

  end subroutine buffer_1d_read

#endif

  subroutine buffer_1d_copy(this, x)
    class(buffer_1d_t), intent(inout) :: this
    type(vector_t), intent(inout) :: x
    integer :: n, i

    n = x%n

    if (this%dp_precision) then
       do i = 1, n
          x%x(i) = this%data_dp(i)
       end do
    else
       do i = 1, n
          x%x(i) = this%data_sp(i)
       end do
    end if

  end subroutine buffer_1d_copy

end module buffer_1d
