! Copyright (c) 2020-2022, The Neko Authors
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
!> File format for .dat files, used for any read/write operations involving
!> floating point data.
!! @details This module defines an interface for read/write operations of user 
!! specified floating-point data
module dat_file
  use probe
  use dat_outputs
  use logger
  use generic_file
  use comm
  use mpi_types
  use utils
  implicit none

  type, public, extends(generic_file_t) :: dat_file_t
   contains
     procedure :: write => dat_file_write
     procedure :: read => dat_file_read
  end type dat_file_t

contains

  !> Writes data to an output file
  subroutine dat_file_write(this, data, t)
    class (dat_file_t), intent(inout) :: this 
    class(*), target, intent(in) :: data  
    real(kind=rp), intent(in), optional :: t

    real(kind=rp) :: time
    type(vector_output_t), pointer :: vec => null()
    type(matrix_output_t), pointer :: mat => null()

    ! Select which type of output we have
    select type(data)
    type is (vector_output_t)

       if (.not. associated(data%data)) then
          call neko_error("Vector is not associated! Use &
          vector_output_t%init() to associate your array &
          with the vector_output_t object")
       end if
       vec => data

    type is (matrix_output_t)
      
       if (.not. associated(data%data)) then
          call neko_error("Matrix is not associated! Use &
          matrix_output_t%init() to associate your array &
          with the vector_output_t object")
       end if
       mat => data

    class default
       call neko_error("Invalid data. Expected vector_output_t or &
       matrix_output_t")
    end select    

    if (pe_rank .eq. 0) then

       ! Check if we have the time argument
       if ( present(t) ) then
          time = real(t, dp)
       else
          time = 0d0
       end if

       if (associated(vec)) then
          call dat_file_write_vector(this, vec, time)
       else if (associated(mat)) then
          call dat_file_write_matrix(this, mat, time)
       end if

    end if

  end subroutine dat_file_write

  !> Writes a @ vector_output_t object to an output file
  subroutine dat_file_write_vector(d, data, t)
    class(dat_file_t), intent(inout) :: d
    type(vector_output_t), intent(in) :: data
    real(kind=rp), intent(in), optional :: t

    character(len=1024) :: fname
    integer :: suffix_pos, file_unit, i, ierr

    suffix_pos = filename_suffix_pos(d%fname)
    fname = trim(d%fname(1:suffix_pos-1)) // "-out.dat"

    open(file=trim(fname), position="append", iostat=ierr, newunit=file_unit)

    if (present(t)) write (file_unit, '(E15.7)', advance="no") t

    do i=1, data%n - 1
       write (file_unit, '(E15.7)', advance="no") data%data(i)
    end do

    write (file_unit, '(E15.7)') data%data(data%n)

    close(file_unit)

  end subroutine dat_file_write_vector

  !> Writes a @ matrix_output_t object to an output file
  subroutine dat_file_write_matrix(d, data, t)
    class(dat_file_t), intent(inout) :: d
    type(matrix_output_t), intent(in) :: data
    real(kind=rp), intent(in), optional :: t

    character(len=1024) :: fname
    integer :: suffix_pos, file_unit, i,j, ierr

    ! The output file's name will be <input_filename>-out.dat
    suffix_pos = filename_suffix_pos(d%fname)
    fname = trim(d%fname(1:suffix_pos-1)) // "-out.dat"

    open(file=trim(fname), position="append", iostat=ierr, newunit=file_unit)

    do i = 1, data%ni
       if (present(t)) write (file_unit, '(E15.7)', advance="no") t

       do j = 1, data%nj - 1
          write (file_unit, '(E15.7)', advance="no") data%data(i,j)
       end do

       write (file_unit, '(E15.7)') data%data(i,data%nj)
    end do

    close(file_unit)

  end subroutine dat_file_write_matrix

  !> Reads data from an input file.
  !> @note Only possible to read probe formatted text files!
  subroutine dat_file_read(this, data)
    class(dat_file_t) :: this               !< dat_file_t object
    class(*), target, intent(inout) :: data !< data to read

    type(probe_t), pointer :: pb => null() 
    integer :: i, ierr, file_unit, npts
    character(len=LOG_SIZE) :: log_buf
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    real(kind=rp), dimension(:,:), allocatable :: temp_send_coords

    select type(data) 
    type is (probe_t) 
       pb => data
    class default
       call neko_error("Invalid data type for dat_file (expected: probe_t)")
    end select

    npts = -1

    ! == First, get npts from the file.
    if (pe_rank .eq. 0) then

      call neko_log%message("Reading dat file " // trim(this%fname))  
      open(file=trim(this%fname), status='old', iostat=ierr, newunit=file_unit)
      read(file_unit, *) npts
      
    end if ! end pe_rank
    
    ! Send npts to all other ranks
    call MPI_Bcast(npts, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
    
    ! Initialize the probe object with an array of npts for coordinates
    call probe_init(pb, npts)

    ! Allocate a temporary array that will store the points
    ! very ugly
    allocate(temp_send_coords(3,npts))

    ! == Now, read the coordinates
    if (pe_rank .eq. 0) then
      
      do i=1, npts
        read (file_unit, *) pb%xyz_coords(i)%x(1), pb%xyz_coords(i)%x(2), &
        pb%xyz_coords(i)%x(3)
      end do

      close(unit=file_unit)
      
      !very ugly
      temp_send_coords(1,:) = pb%xyz_coords(:)%x(1)
      temp_send_coords(2,:) = pb%xyz_coords(:)%x(2)
      temp_send_coords(3,:) = pb%xyz_coords(:)%x(3)
    end if

    ! ouch
    call MPI_Bcast(temp_send_coords, pb%npts, MPI_POINT, 0, NEKO_COMM, ierr)

    if (pe_rank .ne. 0) then
      !very ugly
      pb%xyz_coords(:)%x(1) = temp_send_coords(1,:)
      pb%xyz_coords(:)%x(2) = temp_send_coords(2,:)
      pb%xyz_coords(:)%x(3) = temp_send_coords(3,:)
    end if

  end subroutine dat_file_read


end module dat_file