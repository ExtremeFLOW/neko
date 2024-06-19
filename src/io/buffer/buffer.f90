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
module buffer
  use num_types
  use vector
  use adios2
  implicit none

  type, abstract :: buffer_t
     logical :: dp_precision = .false. !< Precision of output data
   contains
     procedure :: init => buffer_init
     procedure :: fill => buffer_fill
     procedure :: define => buffer_define
     procedure :: inquire => buffer_inquire
     procedure :: write => buffer_write
     procedure :: read => buffer_read
     procedure :: copy => buffer_copy
     procedure :: set_precision => buffer_set_precision
  end type buffer_t

contains

  subroutine buffer_init(this, precision, gdim, glb_nelv, offset_el, nelv, lx, ly, lz)
    class(buffer_t), intent(inout) :: this
    logical, intent(in) :: precision
    integer, intent(in) :: gdim, glb_nelv, offset_el, nelv, lx, ly, lz
  end subroutine buffer_init

  subroutine buffer_fill(this, x, n)
    class(buffer_t), intent(inout) :: this
    integer, intent(inout) :: n
    real(kind=rp), intent(inout) :: x(n)
  end subroutine buffer_fill

  subroutine buffer_define(this, variable, io, variable_name, ierr)
    class(buffer_t), intent(inout) :: this
    type(adios2_variable), intent(inout) :: variable
    type(adios2_io), intent(inout) :: io
    character(len=*), intent(in) :: variable_name
    integer, intent(inout) :: ierr
  end subroutine buffer_define

  subroutine buffer_inquire(this, variable, io, variable_name, ierr)
    class(buffer_t), intent(inout) :: this
    type(adios2_variable), intent(inout) :: variable
    type(adios2_io), intent(inout) :: io
    character(len=*), intent(in) :: variable_name
    integer, intent(inout) :: ierr
  end subroutine buffer_inquire

  subroutine buffer_write(this, engine, variable, ierr)
    class(buffer_t), intent(inout) :: this
    type(adios2_engine), intent(in) :: engine
    type(adios2_variable), intent(in) :: variable
    integer, intent(inout) :: ierr
  end subroutine buffer_write

  subroutine buffer_read(this, engine, variable, ierr)
    class(buffer_t), intent(inout) :: this
    type(adios2_engine), intent(in) :: engine
    type(adios2_variable), intent(in) :: variable
    integer, intent(inout) :: ierr
  end subroutine buffer_read

  subroutine buffer_copy(this, x)
    class(buffer_t), intent(inout) :: this
    type(vector_t), intent(inout) :: x
  end subroutine buffer_copy

  subroutine buffer_set_precision(this, precision)
    class(buffer_t) :: this
    logical, intent(in) :: precision

    !< Switches between double and single precision
    this%dp_precision = precision

  end subroutine buffer_set_precision

end module buffer
