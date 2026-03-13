! Copyright (c) 2020-2023, The Neko Authors
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
!> Defines an output
module output_m
  use num_types_m, only : rp
  use file_m, only : file_t
  implicit none
  private

  !> Abstract type defining an output type
  type, public, abstract :: output_t
     type(file_t) :: file_
   contains
     procedure, pass(this) :: init_base => output_init_base
     procedure, pass(this) :: free_base => output_free_base
     procedure, pass(this) :: set_counter => output_set_counter
     procedure, pass(this) :: set_start_counter => output_set_start_counter
     procedure(output_sample), pass(this), deferred :: sample
     procedure(output_free), pass(this), deferred :: free
  end type output_t

  !> Wrapper around an `output_t` pointer.
  type, public :: output_ptr_t
     class(output_t), pointer :: ptr => null()
  end type output_ptr_t

  !> Abstract interface for sampling an output type at time @a t
  abstract interface
     subroutine output_sample(this, t)
       import :: output_t
       import rp
       class(output_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t
     end subroutine output_sample

     subroutine output_free(this)
       import :: output_t
       class(output_t), intent(inout) :: this
     end subroutine output_free
  end interface

contains

  !> Output constructor.
  !! @param fname Name of the output file.
  !! @param precision Output precision (sp or dp).
  subroutine output_init_base(this, fname, precision, layout, overwrite)
    class(output_t), intent(inout) :: this
    character(len=*), intent(inout) :: fname
    integer, intent(in), optional :: precision
    integer, intent(in), optional :: layout
    logical, intent(in), optional :: overwrite

    call this%file_%init(fname, precision = precision, layout = layout, &
         overwrite = overwrite)

  end subroutine output_init_base

  !> Update the output's file counter
  subroutine output_set_counter(this, n)
    class(output_t), intent(inout) :: this
    integer, intent(in) :: n
    call this%file_%set_counter(n)
  end subroutine output_set_counter

  !> Update the start of output's file counter
  subroutine output_set_start_counter(this, n)
    class(output_t), intent(inout) :: this
    integer, intent(in) :: n
    call this%file_%set_start_counter(n)
  end subroutine output_set_start_counter

  !> Free the output
  subroutine output_free_base(this)
    class(output_t), intent(inout) :: this
    call this%file_%free()
  end subroutine output_free_base

end module output_m
