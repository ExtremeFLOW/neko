! Copyright (c) 2024, The Neko Authors
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
  use logger, only : neko_log, LOG_SIZE
  use stack, only : stack_r8_t, stack_i4r8t2_t
  use tuple, only : tuple_i4r8_t
  use num_types, only : dp
  use comm
  use mpi_f08
  implicit none
  private

  integer :: RT_STATS_MAX_REGIONS = 50
  
  type :: runtime_stats_t
     !> Name of measured region
     character(len=19), allocatable :: rt_stats_id(:)
     !> Elapsed time for each measured region
     type(stack_r8_t), allocatable :: elapsed_time_(:)
     !> Stack to hold current active region timestamps
     type(stack_i4r8t2_t) :: region_timestamp_
    contains
     procedure, pass(this) :: init => runtime_stats_init
     procedure, pass(this) :: free => runtime_stats_free
     procedure, pass(this) :: start_region => runtime_stats_start_region
     procedure, pass(this) :: end_region => runtime_stats_end_region
     procedure, pass(this) :: report => runtime_stats_report
  end type runtime_stats_t

  type(runtime_stats_t), public :: neko_rt_stats

contains

  !> Initialise runtime statistics
  subroutine runtime_stats_init(this)
    class(runtime_stats_t), intent(inout) :: this
    integer :: i
    
    call this%free()

    allocate(this%rt_stats_id(RT_STATS_MAX_REGIONS))

    this%rt_stats_id = ''

    allocate(this%elapsed_time_(RT_STATS_MAX_REGIONS))
    do i = 1, RT_STATS_MAX_REGIONS
       call this%elapsed_time_(i)%init()
    end do

    call this%region_timestamp_%init(100)
    
  end subroutine runtime_stats_init

  !> Destroy runtime statistics 
  subroutine runtime_stats_free(this)
    class(runtime_stats_t), intent(inout) :: this
    integer :: i
    
    if (allocated(this%rt_stats_id)) then
       deallocate(this%rt_stats_id)
    end if

    if (allocated(this%elapsed_time_)) then
       do i = 1, size(this%elapsed_time_)
          call this%elapsed_time_(i)%free()
       end do
       deallocate(this%elapsed_time_)
    end if

    call this%region_timestamp_%free()
    
  end subroutine runtime_stats_free

  !> Start measuring time for the region
  !! named @a name with id @a region_id
  subroutine runtime_stats_start_region(this, name, region_id)
    class(runtime_stats_t), intent(inout) :: this
    character(len=*) :: name
    integer, optional, intent(in) :: region_id
    type(tuple_i4r8_t) :: region_data

    if (present(region_id)) then
       if (region_id .le. RT_STATS_MAX_REGIONS) then
          if (len_trim(this%rt_stats_id(region_id)) .eq. 0) then
             this%rt_stats_id(region_id) = trim(name)
          end if
          region_data%y = MPI_Wtime()
          region_data%x = region_id
          call this%region_timestamp_%push(region_data)
       end if
    else
       region_data%x = 0
       region_data%y = MPI_Wtime()
       call this%region_timestamp_%push(region_data)
    end if
    
  end subroutine runtime_stats_start_region

  !> Compute elapsed time for the current region
  subroutine runtime_stats_end_region(this)
    class(runtime_stats_t), intent(inout) :: this
    real(kind=dp) :: end_time, elapsed_time
    type(tuple_i4r8_t) :: region_data

    end_time = MPI_Wtime()
    region_data = this%region_timestamp_%pop()    
    
    if (region_data%x .gt. 0) then
       elapsed_time = end_time - region_data%y
       call this%elapsed_time_(region_data%x)%push(elapsed_time)
    end if

  end subroutine runtime_stats_end_region

  !> Report runtime statistics for all recorded regions
  subroutine runtime_stats_report(this)
    class(runtime_stats_t), intent(inout) :: this    
    character(len=LOG_SIZE) :: log_buf
    real(kind=dp) :: avg, std, sem, total
    integer :: i, nsamples

    call neko_log%section('Runtime statistics')
    call neko_log%newline()
    write(log_buf, '(A,A,1x,A,1x,A)') '                  ',&
         '      Total time  ','     Avg. time ','     Range  +/-'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') &
         '--------------------------------------------------------------------'
    call neko_log%message(log_buf)

    do i = 1, size(this%elapsed_time_)
       if (len_trim(this%rt_stats_id(i)) .gt. 0) then
          nsamples = this%elapsed_time_(i)%size()
          if (nsamples .gt. 0) then
             select type (region_sample => this%elapsed_time_(i)%data)
             type is (double precision)
                total = sum(region_sample(1:nsamples))
                call MPI_Allreduce(MPI_IN_PLACE, total, 1, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, NEKO_COMM)
                total = total / pe_size
                avg = total / nsamples
                std = (total - avg)**2 / nsamples
                sem = std /sqrt(real(nsamples, dp))
             end select
             write(log_buf, '(A, E15.7,1x,1x,E15.7,1x,1x,E15.7)')  &
                  this%rt_stats_id(i), total, avg, 2.5758_dp * sem
             call neko_log%message(log_buf)
          end if
       end if
    end do
    
    call neko_log%newline()
    call neko_log%end_section()

  end subroutine runtime_stats_report
  
end module runtime_stats
