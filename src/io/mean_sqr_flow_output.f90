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
!> Defines an output for a mean squared flow field
module mean_sqr_flow_output
  use mean_sqr_flow
  use num_types
  use output
  implicit none
  private

  type, public, extends(output_t) :: mean_sqr_flow_output_t
     type(mean_sqr_flow_t), pointer :: msqrf
     real(kind=rp) :: T_begin
   contains
     procedure, pass(this) :: sample => mean_sqr_flow_output_sample
  end type mean_sqr_flow_output_t

  interface mean_sqr_flow_output_t
     module procedure mean_sqr_flow_output_init
  end interface mean_sqr_flow_output_t

contains
  
  function mean_sqr_flow_output_init(msqrf, T_begin, name, path) result(this)
    type(mean_sqr_flow_t), intent(in), target ::msqrf
    real(kind=rp), intent(in) :: T_begin
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    type(mean_sqr_flow_output_t) :: this
    character(len=1024) :: fname
    
    if (present(name) .and. present(path)) then
       fname = trim(path) // trim(name) // '.fld'
    else if (present(name)) then
       fname = trim(name) // '.fld'
    else if (present(path)) then
       fname = trim(path) // 'mean_sqr_field.fld'
    else
       fname = 'mean_sqr_field.fld'
    end if

    call output_init(this, fname)
    this%msqrf => msqrf
    this%T_begin = T_begin
  end function mean_sqr_flow_output_init

  !> Sample a mean squared flow field at time @a t
  subroutine mean_sqr_flow_output_sample(this, t)
    class(mean_sqr_flow_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    if (t .ge. this%T_begin) then
       call this%file_%write(this%msqrf, t)
    end if

  end subroutine mean_sqr_flow_output_sample
  
end module mean_sqr_flow_output


