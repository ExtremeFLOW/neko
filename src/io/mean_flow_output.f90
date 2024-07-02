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
!> Defines an output for a mean flow field
module mean_flow_output
  use mean_flow, only : mean_flow_t
  use num_types, only : rp
  use device
  use output, only : output_t
  use map_1d, only : map_1d_t
  use matrix, only : matrix_t
  implicit none
  private

  type, public, extends(output_t) :: mean_flow_output_t
     type(mean_flow_t), pointer :: mf
     real(kind=rp) :: T_begin
     type(map_1d_t) :: map_1d
   contains
     procedure, pass(this) :: sample => mean_flow_output_sample
  end type mean_flow_output_t

  interface mean_flow_output_t
     module procedure mean_flow_output_init
  end interface mean_flow_output_t

contains

  function mean_flow_output_init(mf, T_begin, hom_dir, name, path) result(this)
    type(mean_flow_t), intent(in), target ::mf
    real(kind=rp), intent(in) :: T_begin
    character(len=*), intent(inout) :: hom_dir
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    type(mean_flow_output_t) :: this
    character(len=1024) :: fname

    if (trim(hom_dir) .eq. 'none') then
       if (present(name) .and. present(path)) then
          fname = trim(path) // trim(name) // '.fld'
       else if (present(name)) then
          fname = trim(name) // '.fld'
       else if (present(path)) then
          fname = trim(path) // 'mean_field.fld'
       else
          fname = 'mean_field.fld'
       end if
    else 
       if (present(name) .and. present(path)) then
          fname = trim(path) // trim(name) // '.csv'
       else if (present(name)) then
          fname = trim(name) // '.csv'
       else if (present(path)) then
          fname = trim(path) // 'mean_field.csv'
       else
          fname = 'mean_field.csv'
       end if
       call this%map_1d%init_char(mf%coef, hom_dir,1e-7_rp)
    end if


    call this%init_base(fname)
    this%mf => mf
    this%T_begin = T_begin
  end function mean_flow_output_init

  !> Sample a mean flow field at time @a t
  subroutine mean_flow_output_sample(this, t)
    class(mean_flow_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    type(matrix_t) :: avg_output_1d

    if (t .ge. this%T_begin) then
       call device_memcpy(this%mf%p%mf%x, this%mf%p%mf%x_d, this%mf%p%mf%dof%size(), &
                          DEVICE_TO_HOST, sync=.false.)
       call device_memcpy(this%mf%u%mf%x, this%mf%u%mf%x_d, this%mf%p%mf%dof%size(), &
                          DEVICE_TO_HOST, sync=.false.)
       call device_memcpy(this%mf%v%mf%x, this%mf%v%mf%x_d, this%mf%p%mf%dof%size(), &
                          DEVICE_TO_HOST, sync=.false.)
       call device_memcpy(this%mf%w%mf%x, this%mf%w%mf%x_d, this%mf%p%mf%dof%size(), &
                          DEVICE_TO_HOST, sync=.true.)
       if (allocated(this%map_1d%pt_lvl)) then
            call this%map_1d%average_planes(avg_output_1d, this%mf%list)
            call this%file_%write(avg_output_1d, t)
       else
          call this%file_%write(this%mf, t)
       end if
       call this%mf%reset()
    end if

  end subroutine mean_flow_output_sample

end module mean_flow_output


