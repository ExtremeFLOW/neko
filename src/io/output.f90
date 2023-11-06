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
module output
  use num_types, only : rp
  use file, only : file_t
  implicit none

  !> Abstract type defining an output type
  type, public, abstract :: output_t
     type(file_t) :: file_
   contains
     procedure, pass(this) :: init => output_init
     procedure, pass(this) :: set_counter => output_set_counter
     procedure, pass(this) :: set_start_counter => output_set_start_counter
     procedure(output_sample), pass(this), deferred :: sample
  end type output_t

  !> Abstract interface for sampling an output type at time @a t
  abstract interface
     subroutine output_sample(this, t)
       import :: output_t
       import rp
       class(output_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t
     end subroutine output_sample
  end interface

contains

  !> Output constructor
  subroutine output_init(this, fname)
    class(output_t), intent(inout) :: this
    character(len=*), intent(inout) :: fname

    this%file_ = file_t(fname)
    
  end subroutine output_init

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
  
 
end module output
