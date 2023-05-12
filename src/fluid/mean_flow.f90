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
!> Defines a mean flow field
module mean_flow
  use mean_field, only : mean_field_t
  use field, only : field_t
  implicit none
  private

  type, public :: mean_flow_t
     type(mean_field_t) :: u
     type(mean_field_t) :: v
     type(mean_field_t) :: w
     type(mean_field_t) :: p
   contains
     procedure, pass(this) :: init => mean_flow_init
     procedure, pass(this) :: free => mean_flow_free
     procedure, pass(this) :: reset => mean_flow_reset
  end type mean_flow_t

contains

  !> Initialize a mean flow field
  subroutine mean_flow_init(this, u, v, w, p)
    class(mean_flow_t), intent(inout) :: this
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p

    call this%free()

    call this%u%init(u)
    call this%v%init(v)
    call this%w%init(w)
    call this%p%init(p)
    
  end subroutine mean_flow_init


  !> Deallocates a mean flow field
  subroutine mean_flow_free(this)
    class(mean_flow_t), intent(inout) :: this

    call this%u%free()
    call this%v%free()
    call this%w%free()
    call this%p%free()

  end subroutine mean_flow_free
 
  !> Resets a mean flow field
  subroutine mean_flow_reset(this)
    class(mean_flow_t), intent(inout) :: this

    call this%u%reset()
    call this%v%reset()
    call this%w%reset()
    call this%p%reset()

  end subroutine mean_flow_reset
   
end module mean_flow
