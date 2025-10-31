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
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_QUIET
  use stack, only : stack_r8_t, stack_i4r8t2_t
  use tuple, only : tuple_i4r8_t
  use num_types, only : dp
  use json_utils, only : json_get_or_default
  use json_module, only : json_file
  use file, only : file_t
  use matrix, only : matrix_t
  use utils, only : neko_error
  use comm, only : pe_rank, pe_size, NEKO_COMM
  use mpi_f08, only : MPI_Wtime, MPI_Allreduce, MPI_IN_PLACE, &
       MPI_DOUBLE_PRECISION, MPI_SUM
  implicit none
  private

  integer, parameter :: RT_STATS_MAX_REGIONS = 50
  integer, parameter :: RT_STATS_RESERVED_REGIONS = 25
  integer, parameter :: RT_STATS_MAX_NAME_LEN = 25

  type :: runtime_stats_t
     private
     !> Name of measured region
     character(len=RT_STATS_MAX_NAME_LEN), allocatable :: rt_stats_id(:)
     !> Elapsed time for each measured region
     type(stack_r8_t), allocatable :: elapsed_time(:)
     !> Stack to hold current active region timestamps
     type(stack_i4r8t2_t) :: region_timestamp
     logical :: enabled = .false.
     logical :: output_profile = .false.
   contains
     procedure, public, pass(this) :: init => runtime_stats_init
     procedure, public, pass(this) :: free => runtime_stats_free
     procedure, public, pass(this) :: start_region => runtime_stats_start_region
     procedure, public, pass(this) :: end_region => runtime_stats_end_region
     procedure, public, pass(this) :: report => runtime_stats_report

     procedure, pass(this) :: find_region_id => runtime_stats_find_region_id
  end type runtime_stats_t

  type(runtime_stats_t), public :: neko_rt_stats

contains

  !> Initialise runtime statistics
  subroutine runtime_stats_init(this, params)
    class(runtime_stats_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    integer :: i

    call this%free()

    call json_get_or_default(params, 'case.runtime_statistics.enabled', &
         this%enabled, .false.)
    call json_get_or_default(params, &
         'case.runtime_statistics.output_profile', &
         this%output_profile, .false.)

    if (this%enabled) then

       allocate(this%rt_stats_id(RT_STATS_MAX_REGIONS))

       this%rt_stats_id = ''

       allocate(this%elapsed_time(RT_STATS_MAX_REGIONS))
       do i = 1, RT_STATS_MAX_REGIONS
          call this%elapsed_time(i)%init()
       end do

       call this%region_timestamp%init(100)

    end if

  end subroutine runtime_stats_init

  !> Destroy runtime statistics
  subroutine runtime_stats_free(this)
    class(runtime_stats_t), intent(inout) :: this
    integer :: i

    if (allocated(this%rt_stats_id)) then
       deallocate(this%rt_stats_id)
    end if

    if (allocated(this%elapsed_time)) then
       do i = 1, size(this%elapsed_time)
          call this%elapsed_time(i)%free()
       end do
       deallocate(this%elapsed_time)
    end if

    call this%region_timestamp%free()

  end subroutine runtime_stats_free

  !> Start measuring time for the region
  !! named @a name with id @a region_id
  subroutine runtime_stats_start_region(this, name, region_id)
    class(runtime_stats_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, optional, intent(in) :: region_id
    type(tuple_i4r8_t) :: region_data
    integer :: id

    if (.not. this%enabled) return

    if (present(region_id)) then
       id = region_id
    else
       call this%find_region_id(name, id)
    end if

    if (id .gt. 0 .and. id .le. RT_STATS_MAX_REGIONS) then
       if (len_trim(this%rt_stats_id(id)) .eq. 0) then
          this%rt_stats_id(id) = trim(name)
       else if (trim(this%rt_stats_id(id)) .ne. trim(name)) then
          call neko_error('Profile region renamed')
       end if
       region_data%x = id
       region_data%y = MPI_Wtime()
       call this%region_timestamp%push(region_data)
    else
       call neko_error('Invalid profiling region id')
    end if

  end subroutine runtime_stats_start_region

  !> Compute elapsed time for the current region
  !! @param name Optional name of the region to close.
  !! @param region_id Optional id of the region to close.
  subroutine runtime_stats_end_region(this, name, region_id)
    class(runtime_stats_t), intent(inout) :: this
    character(len=*), optional, intent(in) :: name
    integer, optional, intent(in) :: region_id
    real(kind=dp) :: end_time, elapsed_time
    type(tuple_i4r8_t) :: region_data
    character(len=1024) :: error_msg
    integer :: id

    if (.not. this%enabled) return

    end_time = MPI_Wtime()
    region_data = this%region_timestamp%pop()

    if (region_data%x .le. 0) then
       call neko_error('Invalid profiling region closed')
    end if

    ! If we are given a name, check it matches the id
    if (present(name)) then
       if (present(region_id)) then
          id = region_id
       else
          call this%find_region_id(name, id)
       end if

       if (trim(this%rt_stats_id(id)) .ne. trim(name)) then
          write(error_msg, '(A,I0,A,A,A)') 'Invalid profiler region closed (', &
               id, ', expected: ', trim(this%rt_stats_id(id)), ')'
          call neko_error(trim(error_msg))

       else if (region_data%x .ne. id) then

          write(error_msg, '(A,A,A,A,A)') 'Invalid profiler region closed (', &
               trim(this%rt_stats_id(region_data%x)), ', expected: ', &
               trim(this%rt_stats_id(id)), ')'
          call neko_error(trim(error_msg))
       end if
    end if

    elapsed_time = end_time - region_data%y
    call this%elapsed_time(region_data%x)%push(elapsed_time)

  end subroutine runtime_stats_end_region

  !> Report runtime statistics for all recorded regions
  subroutine runtime_stats_report(this)
    class(runtime_stats_t), intent(inout) :: this
    character(len=LOG_SIZE) :: log_buf, fmt
    character(len=1250) :: hdr
    real(kind=dp) :: avg, std, sem, total
    integer :: i, nsamples, ncols, nrows, col_idx
    type(matrix_t) :: profile_data

    if (.not. this%enabled) return

    call neko_log%section('Runtime statistics', NEKO_LOG_QUIET)
    call neko_log%newline(NEKO_LOG_QUIET)
    write(fmt, '(A,I0,A)') '(', RT_STATS_MAX_NAME_LEN, 'x,1x,A15,2x,A15,2x,A15)'
    write(log_buf, fmt) 'Total time', 'Avg. time', 'Range +/-'
    call neko_log%message(log_buf, NEKO_LOG_QUIET)
    write(log_buf, '(A)') repeat('-', RT_STATS_MAX_NAME_LEN + 50)
    call neko_log%message(log_buf, NEKO_LOG_QUIET)

    ncols = 0
    nrows = 0
    hdr = ''
    do i = 1, RT_STATS_MAX_REGIONS
       if (len_trim(this%rt_stats_id(i)) .gt. 0) then
          nsamples = this%elapsed_time(i)%size()
          ncols = ncols + 1
          hdr = trim(hdr) // trim(this%rt_stats_id(i)) // ', '
          nrows = max(nrows, nsamples)
          if (nsamples .gt. 0) then
             select type (region_sample => this%elapsed_time(i)%data)
             type is (double precision)
                total = sum(region_sample(1:nsamples))
                call MPI_Allreduce(MPI_IN_PLACE, total, 1, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, NEKO_COMM)
                total = total / pe_size
                avg = total / nsamples
                std = (total - avg)**2 / nsamples
                sem = std / sqrt(real(nsamples, dp))
             end select
             write(fmt, '(A,I0,A)') '(A', RT_STATS_MAX_NAME_LEN, &
                  ',1x,E15.7,2x,E15.7,2x,E15.7)'
             write(log_buf, fmt) this%rt_stats_id(i), total, avg, &
                  2.5758_dp * sem
             call neko_log%message(log_buf, NEKO_LOG_QUIET)
          end if
       end if
    end do

    call neko_log%newline(NEKO_LOG_QUIET)

    if (this%output_profile) then
       col_idx = 0
       call profile_data%init(nrows, ncols)
       do i = 1, size(this%elapsed_time)
          if (len_trim(this%rt_stats_id(i)) .gt. 0) then
             nsamples = this%elapsed_time(i)%size()
             col_idx = col_idx + 1
             if (nsamples .gt. 0) then
                select type (region_sample => this%elapsed_time(i)%data)
                type is (double precision)
                   profile_data%x(1:nsamples,col_idx) = &
                        region_sample(1:nsamples)
                   call MPI_Allreduce(MPI_IN_PLACE, &
                        profile_data%x(1:nsamples,col_idx), nsamples, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, NEKO_COMM)
                   profile_data%x(1:nsamples, col_idx) = &
                        profile_data%x(1:nsamples, col_idx) / pe_size
                end select
             end if
          end if
       end do

       if (pe_rank .eq. 0) then
          block
            type(file_t) :: profile_file
            call profile_file%init('profile.csv')
            call profile_file%set_header(hdr)
            call profile_file%write(profile_data)
          end block
       end if
    end if
    call neko_log%end_section()

    call profile_data%free()

  end subroutine runtime_stats_report

  !> Find or allocate a region id for the named region @a name
  subroutine runtime_stats_find_region_id(this, name, region_id)
    class(runtime_stats_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(out) :: region_id
    integer :: i

    region_id = -1

    ! Look for the region name first
    do i = RT_STATS_RESERVED_REGIONS + 1, RT_STATS_MAX_REGIONS
       if (trim(this%rt_stats_id(i)) .eq. trim(name)) then
          region_id = i
          exit
       end if
    end do

    ! If found, return
    if (region_id .ne. -1) return

    ! Otherwise, look for an empty slot
    do i = RT_STATS_RESERVED_REGIONS + 1, RT_STATS_MAX_REGIONS
       if (len_trim(this%rt_stats_id(i)) .eq. 0) then
          region_id = i
          exit
       end if
    end do

    if (region_id .eq. -1) then
       call neko_error('Not enough profiling regions available')
    end if

  end subroutine runtime_stats_find_region_id

end module runtime_stats
