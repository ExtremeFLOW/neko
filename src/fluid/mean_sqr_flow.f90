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
!> Defines a mean squared flow field
module mean_sqr_flow
  use mean_sqr_field
  implicit none

  type :: mean_sqr_flow_t
     type(mean_sqr_field_t) :: uu
     type(mean_sqr_field_t) :: vv
     type(mean_sqr_field_t) :: ww
     type(mean_sqr_field_t) :: pp
   contains
     procedure, pass(this) :: init => mean_sqr_flow_init
     procedure, pass(this) :: free => mean_sqr_flow_free
  end type mean_sqr_flow_t
  
contains
  
  !> Initialize a mean squared flow field
  subroutine mean_sqr_flow_init(this, u, v, w, p)
    class(mean_sqr_flow_t), intent(inout) :: this
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p

    call this%free()

    call this%uu%init(u)
    call this%vv%init(v)
    call this%ww%init(w)
    call this%pp%init(p)
    
  end subroutine mean_sqr_flow_init


  !> Deallocates a mean squared flow field
  subroutine mean_sqr_flow_free(this)
    class(mean_sqr_flow_t), intent(inout) :: this

    call this%uu%free()
    call this%vv%free()
    call this%ww%free()
    call this%pp%free()

  end subroutine mean_sqr_flow_free
  
end module mean_sqr_flow
