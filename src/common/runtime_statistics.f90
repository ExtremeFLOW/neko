! Copyright (c) 2023, The Neko Authors
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
!> Runtime statistics
module runtime_stats
  use logger
  use num_types
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file
  implicit none
  private

  integer, public, parameter :: RT_STATS_FLUID=1, RT_STATS_SCALAR=2, &
                                RT_STATS_POST=3, RT_STATS_TOTAL=4
  
  character(len=13) :: RT_STATS_ID(4) = ['Fluid        ', &
                                         'Scalar       ', &
                                         'Post         ', &
                                         'Total        ']
  !> Runtime statistics backend
  type, public :: runtime_stats_t
     logical, private :: enabled_
     real(kind=dp), private, allocatable :: elapsed_time_(:)
     integer, private, allocatable :: nsamples_(:)
     integer, private :: start_step_
     integer, private :: total_steps_
     real(kind=dp), allocatable :: data(:,:)
   contains
     procedure, pass(this) :: init => runtime_stats_init
     procedure, pass(this) :: record => runtime_stats_record
     procedure, pass(this) :: report => runtime_stats_report
  end type runtime_stats_t

contains

  !> Initialize runtime statistics
  subroutine runtime_stats_init(this, params)
    class(runtime_stats_t), intent(inout) :: this
    type(json_file), target, intent(inout) :: params
    logical :: logical_val
    real(kind=rp) :: dt, end_time

    call json_get_or_default(params, 'case.runtime_statistics.enabled', &
                             this%enabled_, .false.)
    call json_get_or_default(params, 'case.runtime_statistics.start_step', &
                             this%start_step_, 1)

    
    if (this%enabled_) then
       call json_get(params, 'case.timestep', dt)
       call json_get(params, 'case.end_time', end_time)

       this%total_steps_ =  ceiling(end_time / dt)

       allocate(this%data(this%total_steps_, size(RT_STATS_ID)))
       allocate(this%elapsed_time_(size(RT_STATS_ID)))
       allocate(this%nsamples_(size(RT_STATS_ID)))

       this%data = 0d0
       this%elapsed_time_ = 0.0
       this%nsamples_ = 0
       
       
    end if

  end subroutine runtime_stats_init

  !> Record the elapsed time for a specific region at a given time-step
  subroutine runtime_stats_record(this, region_id, elapsed_time, time, tstep)
    class(runtime_stats_t), intent(inout) :: this
    integer, intent(in) :: region_id
    real(kind=dp), intent(in) :: elapsed_time
    real(kind=rp), intent(in) :: time
    integer, intent(in) :: tstep

    if (this%enabled_ .and. tstep .ge. this%start_step_) then

       this%elapsed_time_(region_id) = &
            this%elapsed_time_(region_id) + elapsed_time
       this%nsamples_(region_id) = this%nsamples_(region_id) + 1
       this%data(tstep, region_id) = elapsed_time
       
    end if

  end subroutine runtime_stats_record

  !> Report collected runtime statistics
  subroutine runtime_stats_report(this)
    class(runtime_stats_t), intent(inout) :: this
    character(len=LOG_SIZE) :: log_buf
    real(kind=dp) :: avg, std, sem
    integer :: i

    if (this%enabled_) then
       call neko_log%section('Runtime statistics')
       call neko_log%newline()
       write(log_buf, '(A,A,1x,A,A)') '               ','   Total time  ',&
                                     '     Avg. time ','    Range  +/-'
       call neko_log%message(log_buf)
       write(log_buf, '(A)') &
            '------------------------------------------------------------'
       call neko_log%message(log_buf)
       
       do i = 1, size(RT_STATS_ID)
          if (this%nsamples_(i) .gt. 0) then
             avg = this%elapsed_time_(i) / this%nsamples_(i)
             std = sum((this%data(this%start_step_: this%start_step_ &
                  + (this%nsamples_(i) - 1), i) - avg)**2) / this%nsamples_(i)
             sem = std / sqrt(real(this%nsamples_(i), dp))
             write(log_buf, '(A,E15.7,1x,1x,E15.7, E15.7)')  RT_STATS_ID(i), &
                  this%elapsed_time_(i), avg, 2.5758_rp * sem
             call neko_log%message(log_buf)
          end if
       end do

       call neko_log%newline()
       call neko_log%end_section()
    end if
  end subroutine runtime_stats_report
  
end module runtime_stats
