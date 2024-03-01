! Copyright (c) 2022, The Neko Authors
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
!> Built in profiling tool for neko
module neko_profiler
  use utils, only: neko_error
  use stack, only: stack_i4_t
  use comm, only: NEKO_COMM, pe_rank, pe_size
  use mpi_f08, only: mpi_wtime, mpi_ibcast, mpi_waitall
  use mpi_f08, only: MPI_INTEGER, MPI_DOUBLE, MPI_CHARACTER, &
    MPI_MAX, MPI_STATUSES_IGNORE, MPI_REQUEST
  implicit none

  private
  public :: init_profiler, free_profiler, start_region, end_region

  type :: profiler_region_t
     !> Region name
     character(len=80) :: name
     !> Region ID
     integer :: id
     !> List of execution times
     real(kind=8) :: time
     !> Number of executions
     integer :: n = 0

     !> Latest start time
     real(kind=8) :: start
     !> Logical value indicating if the region is currently active
     logical :: active = .false.
  end type profiler_region_t

  !> List of regions
  type(profiler_region_t), dimension(:), allocatable :: region_list
  !> Number of regions
  integer :: n_regions = 0
  !> Stack of active regions
  type(stack_i4_t) :: active_regions

contains

  ! ========================================================================== !
  ! Definition of the CPU based profiler.

  !> Initialize the profiler
  subroutine init_profiler()

    if (allocated(region_list)) then
       call neko_error("Error in profiler: double initialization")
    end if

    allocate(region_list(0))
    call active_regions%init(10)

    call start_region("Total Time")

  end subroutine init_profiler

  !> Finalize the profiler
  !! @note This will deallocate the region list and the stack of active regions
  subroutine free_profiler(filename)
    character(len=*), intent(in), optional :: filename
    integer :: i

    ! Stop the rest of the regions
    do while (.not. active_regions%is_empty())
       call end_region()
    end do
    call active_regions%free()

    ! Save the profiler log
    call print_regions(filename)

    deallocate(region_list)
    n_regions = 0

  end subroutine free_profiler

  ! ========================================================================== !
  ! Region manipulation

  !> Start a region
  !! @param name Name of the region
  !! @param id Optional ID of the region
  subroutine start_region(name, id)
    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in), optional :: id
    integer :: i_reg, i_stack, j, id_
    type(profiler_region_t) :: tmp_region

    ! If the profiler is not initialized, return
    if (.not. allocated(region_list)) return

    if (present(id)) then
       id_ = id
    else
       id_ = 0
    end if

    do i_reg = 1, n_regions
       if (region_list(i_reg)%name == trim(name) &
           .and. region_list(i_reg)%id == id_) then

          if (region_list(i_reg)%active) then
             call neko_error("Error in profiler: region is already active: " &
                             // name)
          end if

          region_list(i_reg)%active = .true.
          region_list(i_reg)%start = mpi_wtime()

          i_stack = i_reg
          call active_regions%push(i_stack)
          return
       end if
    end do

    tmp_region%name = trim(name)
    tmp_region%id = id_
    tmp_region%active = .true.
    tmp_region%n = 0
    tmp_region%start = mpi_wtime()
    tmp_region%time = 0.0d0

    region_list = [region_list, tmp_region]
    n_regions = size(region_list)

    call active_regions%push(n_regions)

  end subroutine start_region

  !> End a region
  !! @param name Name of the region
  !! @param id Optional ID of the region
  subroutine end_region(name, id)
    character(len=*), intent(in), optional :: name
    integer, intent(in), optional :: id
    integer :: i, id_, i_stack
    real(kind=8) :: t
    type(stack_i4_t) :: tmp_stack

    ! If the profiler is not initialized, return
    if (.not. allocated(region_list)) return

    if (present(name)) then

       id_ = 0
       if (present(id)) id_ = id

       do i = 1, n_regions
          if (region_list(i)%name == trim(name) &
              .and. region_list(i)%id == id_) then
             exit
          end if
       end do

       i_stack = active_regions%pop()
       do while (i_stack /= i)
          call tmp_stack%push(i_stack)
          i_stack = active_regions%pop()
       end do

       do while (.not. tmp_stack%is_empty())
          i_stack = tmp_stack%pop()
          call active_regions%push(i_stack)
       end do

    else

       if (active_regions%is_empty()) then
          call neko_error("Error in profiler: no active region to end")
       end if

       i = active_regions%pop()
    end if


    if (.not. region_list(i)%active) then
       call neko_error("Error in profiler: region is not active: " // name)
    end if

    t = mpi_wtime() - region_list(i)%start

    region_list(i)%time = region_list(i)%time + t
    region_list(i)%n = region_list(i)%n + 1
    region_list(i)%active = .false.

  end subroutine end_region

  !> Print the regions to a file, profiler.log
  subroutine print_regions(filename)
    character(len=*), intent(in), optional :: filename
    integer :: i, n, file_id, ierr
    real(kind=8) :: time
    character(len=80) :: name
    character(len=80) :: time_str
    character(len=80) :: n_str

    ! If the profiler is not initialized, return
    if (.not. allocated(region_list)) return

    do i = 1, n_regions
       if (region_list(i)%active) then
          call neko_error("Error in profiler: region is still active:" &
                          // region_list(i)%name)
       end if
    end do

    if (pe_size .gt. 1) then
       call synchronize_regions()
    end if
    if (pe_rank .ne. 0) return

    if (present(filename)) then
       open(newunit=file_id, file=filename, status="new", action="write")
    else
       open(newunit=file_id, file="profiler.log", status="new", action="write")
    end if

    write(file_id, *) 'Region, CPU Time [s], Evaluations'

    do i = 1, n_regions
       name = trim(region_list(i)%name)
       write (time_str, '(f14.7)') region_list(i)%time
       write (n_str, '(i14)') region_list(i)%n

       write(file_id, '(a,a1,a,a1,a)') &
         trim(adjustl(name)), ",", &
         trim(adjustl(time_str)), ",", &
         trim(adjustl(n_str))
    end do

    close(file_id)
  end subroutine print_regions

  !> Synchronize the regions across all ranks
  !! This will overwrite the local region list with the global region list
  !! and is called before printing the regions.
  subroutine synchronize_regions()
    integer, parameter :: name_len = 80

    ! Buffers for the region information
    integer :: buffer_size
    character(len=name_len), dimension(:,:), allocatable :: buffer_names
    real(kind=8), dimension(:,:), allocatable :: buffer_times
    integer, dimension(:,:), allocatable :: buffer_n, buffer_id

    ! Requests for the broadcasts
    TYPE(MPI_Request), dimension(:), allocatable :: req_names, req_id
    TYPE(MPI_Request), dimension(:), allocatable :: req_times, req_n

    ! Global list of regions, will overwrite the local list
    type(profiler_region_t), dimension(:), allocatable :: global_region_list
    type(profiler_region_t) :: temp_region

    ! Local variables
    character(len=name_len) :: name
    real(kind=8) :: time
    integer :: id, n, i_pe, i_reg, i, ierr

    ! Ensure the buffers are large enough to hold all regions
    call mpi_allreduce(n_regions, buffer_size, 1, MPI_INTEGER, MPI_MAX, &
                       NEKO_COMM, ierr)

    ! Allocate the buffers and requests
    allocate(buffer_names(buffer_size, pe_size))
    allocate(buffer_times(buffer_size, pe_size))
    allocate(buffer_id(buffer_size, pe_size))
    allocate(buffer_n(buffer_size, pe_size))

    allocate(req_names(buffer_size))
    allocate(req_id(buffer_size))
    allocate(req_times(buffer_size))
    allocate(req_n(buffer_size))

    ! Fill the buffers with the local region information
    do i_reg = 1, buffer_size

       if (i_reg > n_regions) then
          buffer_names(i_reg, pe_rank + 1) = ""
          buffer_id(i_reg, pe_rank + 1) = 0
          buffer_times(i_reg, pe_rank + 1) = 0.0
          buffer_n(i_reg, pe_rank + 1) = 0
       else
          buffer_names(i_reg, pe_rank + 1) = region_list(i_reg)%name
          buffer_id(i_reg, pe_rank + 1) = region_list(i_reg)%id

          buffer_times(i_reg, pe_rank + 1) = region_list(i_reg)%time
          buffer_n(i_reg, pe_rank + 1) = region_list(i_reg)%n
       end if

    end do

    ! Broadcast the information to everyone
    do i_pe = 0, pe_size - 1
       do i_reg = 1, buffer_size
          call mpi_ibcast(buffer_names(i_reg, i_pe + 1), name_len, MPI_CHARACTER, i_pe, &
                          NEKO_COMM, req_names(i_reg), ierr)

          call mpi_ibcast(buffer_id(i_reg, i_pe + 1), 1, MPI_INTEGER, i_pe, &
                          NEKO_COMM, req_id(i_reg), ierr)

          call mpi_ibcast(buffer_times(i_reg, i_pe + 1), 1, MPI_DOUBLE, i_pe, &
                          NEKO_COMM, req_times(i_reg), ierr)

          call mpi_ibcast(buffer_n(i_reg, i_pe + 1), 1, MPI_INTEGER, i_pe, &
                          NEKO_COMM, req_n(i_reg), ierr)
       end do

       ! Wait for the broadcasts to finish
       call mpi_waitall(buffer_size, req_names, MPI_STATUSES_IGNORE, ierr)
       call mpi_waitall(buffer_size, req_id, MPI_STATUSES_IGNORE, ierr)
       call mpi_waitall(buffer_size, req_times, MPI_STATUSES_IGNORE, ierr)
       call mpi_waitall(buffer_size, req_n, MPI_STATUSES_IGNORE, ierr)
    end do

    ! Build a new list of regions. Every entry must be unique and the same on all
    ! ranks. Two regions are considered the same if they have the same name and
    ! the same ID. If a region is not present on a rank, it is added with a time
    ! of 0.0 and a count of 0.

    ! First, build the global list from the local list
    allocate(global_region_list(0))

    ! Loop through the global lists and add any newly encountered regions
    ! to the global list
    loop_region: do i_reg = 1, buffer_size
       loop_rank: do i_pe = 1, pe_size
          name = buffer_names(i_reg, i_pe)
          id = buffer_id(i_reg, i_pe)
          time = buffer_times(i_reg, i_pe)
          n = buffer_n(i_reg, i_pe)

          if (trim(name) == "") cycle loop_rank

          ! Determine if the region is already in the global list
          do i = 1, size(global_region_list)
             if (global_region_list(i)%name == name .and. &
                 global_region_list(i)%id == id) then

                ! Save the time and count
                global_region_list(i)%time = global_region_list(i)%time + time
                global_region_list(i)%n = global_region_list(i)%n + n

                cycle loop_rank
             end if
          end do

          ! If the region is not in the global list, add it
          temp_region%name = name
          temp_region%id = id
          temp_region%time = time
          temp_region%n = n

          global_region_list = [global_region_list, temp_region]
       end do loop_rank
    end do loop_region

    ! Save the global list
    region_list = global_region_list

  end subroutine synchronize_regions


end module neko_profiler
