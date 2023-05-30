! Copyright (c) 2020-2023, The Neko Authors
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
!> File format for .csv files, used for any read/write operations involving
!> floating point data.
!! @details This module defines an interface for read/write operations of user 
!! specified floating-point data
module csv_file
  use probe, only: probe_t, probe_init
  use vector, only: vector_t
  use matrix, only: matrix_t
  use logger
  use generic_file, only: generic_file_t
  use comm
  use mpi_types, only: MPI_POINT
  use utils, only: neko_error
  use num_types
  implicit none

  type, public, extends(generic_file_t) :: csv_file_t
   contains
     procedure :: write => csv_file_write
     procedure :: read => csv_file_read
  end type csv_file_t

contains

  !> Writes data to an output file
  !! @param this csv file to write in
  !! @param data Data to write, can be vector_t or matrix_t
  !! @param t Time
  subroutine csv_file_write(this, data, t)
    class (csv_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t

    real(kind=rp) :: time
    type(vector_t), pointer :: vec => null()
    type(matrix_t), pointer :: mat => null()

    select type(data)
    type is (vector_t)
       if (.not. allocated(data%x)) then
          call neko_error("Vector is not allocated! Use &
               vector_output_t%init() to associate your array &
               with the vector_output_t object")
       end if
       vec => data

    type is (matrix_t)
       if (.not. allocated(data%x)) then
          call neko_error("Matrix is not allocated! Use &
               matrix_output_t%init() to associate your array &
               with the vector_output_t object")
       end if
       mat => data

    class default
       call neko_error("Invalid data. Expected vector_t or &
            matrix_t")
    end select

    ! Write is performed on rank 0
    if (pe_rank .eq. 0) then

       if (associated(vec)) then
          call csv_file_write_vector(this, vec, t)
       else if (associated(mat)) then
          call csv_file_write_matrix(this, mat, t)
       end if

    end if

  end subroutine csv_file_write

  !> Writes a @ vector_t object to an output file
  !! @param f csv file in which to write
  !! @param data Vector to write
  !! @param Time
  subroutine csv_file_write_vector(f, data, t)
    class(csv_file_t), intent(inout) :: f
    type(vector_t), intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    character(len=1024) :: fname
    integer :: suffix_pos, file_unit, i, ierr

    open(file=trim(f%fname), position="append", iostat=ierr, newunit=file_unit)

    if (present(t)) write (file_unit, '(g0,",")', advance="no") t

    write (file_unit, '(*(g0,","))', advance="no") data%x(1:data%n-1)
    write(file_unit,'(g0)') data%x(data%n)

    close(file_unit)

  end subroutine csv_file_write_vector

  !> Writes a @ matrix_t object to an output file
  !! @param f csv file in which to write
  !! @param data Matrix to write
  !! @param Time
  subroutine csv_file_write_matrix(d, data, t)
    class(csv_file_t), intent(inout) :: d
    type(matrix_t), intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    integer :: file_unit, i,j, ierr

    open(file=trim(d%fname), position="append", iostat=ierr, newunit=file_unit)

    do i = 1, data%nrows
       if (present(t)) write (file_unit, '(g0,",")', advance="no") t
       write (file_unit, '(*(g0,","))', advance="no") data%x(i,1:data%ncols-1)
       write (file_unit, '(g0)') data%x(i,data%ncols)
    end do

    close(file_unit)

  end subroutine csv_file_write_matrix

  !> Reads data from an input file.
  !! @param this csv file in which to read
  !! @param data Data to read
  subroutine csv_file_read(this, data)
    class(csv_file_t) :: this
    class(*), target, intent(inout) :: data
    type(probe_t), pointer :: pb => null()
    type(vector_t), pointer :: vec => null()
    type(matrix_t), pointer :: mat => null()

    select type(data)
    type is (probe_t)
       pb => data
    type is (vector_t)
       vec => data
    type is (matrix_t)
       mat => data
    class default
       call neko_error("Invalid data type for csv_file (expected: probe_t)")
    end select

    if (associated(pb)) then
       call csv_file_read_probe(this, pb)
    else if (associated(vec)) then
       call csv_file_read_vector(this, vec)
    else if (associated(mat)) then
       call csv_file_read_matrix(this, mat)
    end if

  end subroutine csv_file_read

  !> Read a probe from a dat file
  !! @param d csv file in which to read
  !! @param pb Probe object in which to store the probes
  subroutine csv_file_read_probe(d, pb)
    type(csv_file_t), intent(in) :: d !< dat file object
    type(probe_t), pointer, intent(inout) :: pb !< probe
    real(kind=rp), dimension(:,:), allocatable :: temp_send_coords
    integer :: i, ierr, file_unit, npts

    !
    ! First, read npts from the file.
    !
    if (pe_rank .eq. 0) then

       call neko_log%message("Reading dat file " // trim(d%fname))
       open(file=trim(d%fname), status='old', iostat=ierr, newunit=file_unit)
       read(file_unit, *) npts

    end if

    ! Send npts to all other ranks
    call MPI_Bcast(npts, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)

    ! Initialize the probe object with an array of npts for coordinates
    call probe_init(pb, npts)

    ! Allocate a tempoqrary array that will store the points
    allocate(temp_send_coords(3,npts))

    !
    ! Now, read the coordinates from the file
    !
    if (pe_rank .eq. 0) then

       do i=1, npts
          read (file_unit, *) pb%xyz_coords(i)%x(1), pb%xyz_coords(i)%x(2), &
               pb%xyz_coords(i)%x(3)
       end do

       close(unit=file_unit)

       ! Store probes in temporary array
       temp_send_coords(1,:) = pb%xyz_coords(:)%x(1)
       temp_send_coords(2,:) = pb%xyz_coords(:)%x(2)
       temp_send_coords(3,:) = pb%xyz_coords(:)%x(3)
    end if

    !
    ! Send the probe coordintes to all other ranks
    ! I'm doing this because sending pb%xyz_coords
    ! throws an error from MPI_Bcast
    !
    call MPI_Bcast(temp_send_coords, pb%npts, MPI_POINT, 0, NEKO_COMM, ierr)

    !
    ! And for all other ranks load the probe coordinates
    !
    if (pe_rank .ne. 0) then
       pb%xyz_coords(:)%x(1) = temp_send_coords(1,:)
       pb%xyz_coords(:)%x(2) = temp_send_coords(2,:)
       pb%xyz_coords(:)%x(3) = temp_send_coords(3,:)
    end if

    deallocate(temp_send_coords)

  end subroutine csv_file_read_probe

  !> Read a vector from a csv file
  !! @param d csv file from which to read
  !! @param vec Vector object in which to store the file contents
  !> @TODO: Broadcast to other processes
  subroutine csv_file_read_vector(d, vec)
    type(csv_file_t), intent(in) :: d
    type(vector_t), intent(inout) :: vec
    integer :: ierr, file_unit

    if ( pe_rank .eq. 0 ) then
       open(file=trim(d%fname), status='old', iostat=ierr, newunit=file_unit)
       read (file_unit,*) vec%x
       close(unit=file_unit)
    end if

  end subroutine csv_file_read_vector

  !> Read a vector from a csv file
  !! @param d csv file from which to read
  !! @param vec Matrix object in which to store the file contents
  !! @TODO: Broadcast to other processes
  subroutine csv_file_read_matrix(d, mat)
    type(csv_file_t), intent(in) :: d
    type(matrix_t), intent(inout) :: mat
    integer :: ierr, file_unit, i

    if ( pe_rank .eq. 0 ) then
       open(file=trim(d%fname), status='old', iostat=ierr, newunit=file_unit)
       do i=1, mat%nrows
          read (file_unit,*) mat%x(i,:)
       end do
       close(unit=file_unit)
    end if

  end subroutine csv_file_read_matrix


end module csv_file
