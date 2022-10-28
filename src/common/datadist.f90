! Copyright (c) 2019-2021, The Neko Authors
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
!> Defines practical data distributions
module datadist
  use mpi_f08
  use logger
  implicit none
  private

  type dist_t
     type(MPI_Comm) :: comm !< Communicator on which the dist. is defined
     integer :: pe_rank     !< Pe's rank in the given distribution
     integer :: pe_size     !< Size of communicator in the given dist.
     integer :: weight !< Weight of this rank 
                       !! (How large piece of the cake should it get?)
     integer :: start_weight !< Weight of all ranks before this
     integer :: total_weight !< Weight of all processes
     integer :: M !< Total size
     integer :: Ip !< Number of elements on this rank
     integer :: start !< Starting idx for data on this rank
  end type dist_t

  !> Load-balanced linear distribution \f$ M = PL + R \f$
  type, extends(dist_t) :: linear_dist_t
   contains 
     procedure :: num_local => linear_dist_Ip
     procedure :: num_global => linear_dist_M
     procedure :: start_idx => linear_dist_start
     procedure :: end_idx => linear_dist_end
  end type linear_dist_t

  interface linear_dist_t
     module procedure linear_dist_init
  end interface linear_dist_t

  public :: linear_dist_t

contains
  
  function linear_dist_init(n, rank, size, comm, weight) result(this)
    integer, intent(in) :: n    !< Total size
    integer :: rank             !< PE's rank to define the dist. over
    integer :: size             !< Size of comm where the dist. is def. on
    type(MPI_Comm) :: comm      !< comm. to define the dist. over
    type(linear_dist_t), target :: this
    integer, optional :: weight
    character(len=40) :: str_weight !< input string of long enough length
    integer :: L, R, ierr, len, status
    character(len=LOG_SIZE) :: log_buf

    this%M = n
    this%comm = comm
    this%pe_rank = rank
    this%pe_size = size
    if (present(weight)) then
       this%weight = weight
    else
       call get_environment_variable ("rank_weight", str_weight, len, status)
       if (status .eq. 0) then
          read(str_weight, *) this%weight
          print *, this%weight
       else
          this%weight = 1
       end if
    end if

    
    this%start_weight = 0
    call MPI_Exscan(this%weight, this%start_weight, 1, &
         MPI_INTEGER, MPI_SUM, comm, ierr)
    this%total_weight = this%weight
    call MPI_Allreduce(MPI_IN_PLACE, this%total_weight, &
                       1, MPI_INTEGER,MPI_SUM, comm, ierr)

    call neko_log%message("Mesh distribution")
    if (this%total_weight .eq. this%pe_size) then 
       L = floor(dble(this%M) / dble(this%pe_size))
       R = modulo(this%M, this%pe_size)
       this%Ip = floor((dble(this%M) + dble(this%pe_size) - &
            dble(this%pe_rank) - 1d0) / dble(this%pe_size))    
       this%start = this%pe_rank * L + min(this%pe_rank, R)
       write(log_buf,1) size
1      format('Mesh split uniformly across', i8, ' ranks')
       call neko_log%message(log_buf)
    else
       this%Ip = this%weight*floor((dble(this%M))/dble(this%total_weight))&
               + min(max(modulo(this%M, this%total_weight)-this%start_weight,0),this%weight)
       this%start = 0
       call MPI_Exscan(this%Ip,this%start,1, &
            MPI_INTEGER, MPI_SUM, comm, ierr)
       if (size -1 .eq. rank) then
          this%Ip = n - this%start
       end if
       write(log_buf,2) size
2      format('Mesh split non-uniformly across ', i8, ' ranks')
       call neko_log%message(log_buf)
       write (*,'(A,i6,A,i6,A,i9,A,i9)') 'Rank: ', rank, ' N el: ', this%Ip,&
             ' Start el: ', this%start, ' Cumulative el: ', this%start+this%Ip
    end if
  end function linear_dist_init

  pure function linear_dist_Ip(this) result(n)
    class(linear_dist_t), intent(in) :: this
    integer :: n
    n = this%Ip
  end function linear_dist_Ip

  pure function linear_dist_M(this) result(M)
    class(linear_dist_t), intent(in) :: this
    integer :: M
    M = this%M
  end function linear_dist_M

  pure function linear_dist_start(this) result(start)
    class(linear_dist_t), intent(in) :: this
    integer :: start
    start = this%start
  end function linear_dist_start

  function linear_dist_end(this) result(end)
    class(linear_dist_t), intent(inout) :: this
    integer :: end
    end = linear_dist_start(this) + (this%Ip - 1)
  end function linear_dist_end
end module datadist
