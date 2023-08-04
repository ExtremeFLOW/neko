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
!> Defines a container for all statistics
module stats
  use num_types
  use stats_quant
  use logger
  use comm
  implicit none

  !> Pointer to an arbitrary quantitiy
  type, private :: quantp_t
     class(stats_quant_t), pointer :: quantp
  end type quantp_t

  !> Statistics backend
  type :: stats_t
     type(quantp_t), allocatable :: quant_list(:)
     integer :: n
     integer :: size
     real(kind=rp) :: T_begin !< Start time for averaging
     real(kind=rp) :: t_diff !< Time since last sample
     integer :: samp_interval !< Number of time steps between samples
   contains
     procedure, pass(this) :: init => stats_init
     procedure, pass(this) :: free => stats_free
     procedure, pass(this) :: add => stats_add
     procedure, pass(this) :: eval => stats_eval
  end type stats_t

contains

  !> Initialize statistics, computed after @a T_begin
  subroutine stats_init(this, T_begin, samp_interval, size)
    class(stats_t), intent(inout) :: this
    real(kind=rp), intent(in) :: T_begin
    integer, intent(in) :: samp_interval
    integer, intent(inout), optional ::size
    integer :: n, i
    
    call this%free()

    if (present(size)) then
       n = size
    else
       n = 1
    end if

    allocate(this%quant_list(n))

    do i = 1, n
       this%quant_list(i)%quantp => null()
    end do

    this%n = 0
    this%size = n
    this%T_begin = T_begin
    this%samp_interval = samp_interval
    this%t_diff = 0.0
    
  end subroutine stats_init

  !> Deallocate
  subroutine stats_free(this)
    class(stats_t), intent(inout) :: this

    if (allocated(this%quant_list)) then
       deallocate(this%quant_list)
    end if

    this%n = 0
    this%size = 0    
  end subroutine stats_free

  !> Add a statistic quantitiy @a quant to the backend
  subroutine stats_add(this, quant)
    class(stats_t), intent(inout) :: this
    class(stats_quant_t), intent(inout), target :: quant
    type(quantp_t), allocatable :: tmp(:)

    if (this%n .ge. this%size) then
       allocate(tmp(this%size * 2))
       tmp(1:this%size) = this%quant_list
       call move_alloc(tmp, this%quant_list)
       this%size = this%size * 2
    end if

    this%n = this%n + 1
    this%quant_list(this%n)%quantp => quant
  end subroutine stats_add

  !> Evaluated all statistical quantities
  subroutine stats_eval(this, t, dt, tstep)
    class(stats_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: tstep
    integer :: i, ierr
    character(len=LOG_SIZE) :: log_buf    
    real(kind=rp) :: sample_start_time, sample_end_time
    real(kind=dp) :: sample_time

    if (t .ge. this%T_begin .and. this%n .gt. 0) then
       this%t_diff = this%t_diff + dt
       ! There is technically an issue here for the last sample if we
       ! reset the stats If the reset is not on a multiple of
       ! samp_interval the weight of the last sample is wrong.
       if (mod(tstep,this%samp_interval) .eq. 0) then
          call neko_log%section('Statistics')
          call MPI_Barrier(NEKO_COMM, ierr)
          sample_start_time = MPI_WTIME()
          do i = 1, this%n
             call this%quant_list(i)%quantp%update(this%t_diff)
          end do
          this%t_diff = 0.0_rp
          call MPI_Barrier(NEKO_COMM, ierr)
          sample_end_time = MPI_WTIME()
          sample_time = sample_end_time - sample_start_time
          write(log_buf,'(A17,1x,F10.6,A,F9.6)') 'Sampling at time:', t, &
          ' Sampling time (s): ', sample_time
          call neko_log%message(log_buf)
          call neko_log%end_section()
       end if
    end if
  end subroutine stats_eval

end module stats
