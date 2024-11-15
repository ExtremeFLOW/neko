! Copyright (c) 2020-2024, The Neko Authors
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
!! floating point data.
module csv_file
  use vector, only: vector_t
  use matrix, only: matrix_t
  use generic_file, only: generic_file_t
  use utils, only: neko_error
  use num_types, only: rp
  use logger, only: neko_log, log_size
  use comm
  implicit none

  type, public, extends(generic_file_t) :: csv_file_t
     character(len=1024) :: header = ""     !< Contains header of file.
     logical :: header_is_written = .false. !< Has header already been written?
   contains
     !> Writes data to an output file.
     procedure :: write => csv_file_write
     !> Reads data from an input file.
     procedure :: read => csv_file_read
     !> Sets the header for a csv file.
     procedure :: set_header => csv_file_set_header
     !> Count the number of lines in a file
     procedure :: count_lines => csv_file_count_lines
  end type csv_file_t

contains

  !> Writes data to an output file.
  !! @param this csv file to write in.
  !! @param data Data to write, can be vector_t or matrix_t.
  !! @param t Time.
  subroutine csv_file_write(this, data, t)
    class (csv_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t

    type(vector_t), pointer :: vec
    type(matrix_t), pointer :: mat

    nullify(vec)
    nullify(mat)

    select type (data)
    type is (vector_t)
       if (.not. allocated(data%x)) then
          call neko_error("Vector is not allocated. Use &
&vector%init() to associate your array &
&with a vector_t object")
       end if
       vec => data

    type is (matrix_t)
       if (.not. allocated(data%x)) then
          call neko_error("Matrix is not allocated. Use &
&matrix%init() to associate your array &
&with a matrix_t object")
       end if
       mat => data

    class default
       call neko_error("Invalid data. Expected vector_t or &
&matrix_t")
    end select

    ! Write is performed on rank 0
    if (pe_rank .eq. 0) then

       call neko_log%message("Writing to " // trim(this%fname))
       if (associated(vec)) then
          call csv_file_write_vector(this, vec, t)
       else if (associated(mat)) then
          call csv_file_write_matrix(this, mat, t)
       end if

    end if

  end subroutine csv_file_write

  !> Writes a `vector_t` object to an output file, in a row format.
  !! If the parameter `t` is present, it will be appended at the start
  !! of the row, followed by the contents of the vector.
  !! @param f csv file in which to write.
  !! @param data Vector to write.
  !! @param Time.
  !! @note Equivalent to writing a matrix with `nrows = 1` and `ncols = N`.
  !! To write a vector in a column format, use a `matrix_t` with `ncols = 1`.
  subroutine csv_file_write_vector(f, data, t)
    class(csv_file_t), intent(inout) :: f
    type(vector_t), intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    integer :: file_unit, ierr

    open(file = trim(f%fname), position = "append", iostat = ierr, &
         newunit = file_unit)
    if (ierr .ne. 0) call neko_error("Error while opening " // trim(f%fname))

    ! write header if not empty and if not already written
    if (f%header .ne. "" .and. .not. f%header_is_written) then
       write (file_unit, '(A)') trim(f%header)
       f%header_is_written = .true.
    end if

    ! Add time at the beginning if specified
    if (present(t)) write (file_unit, '(g0,",")', advance = "no") t

    write (file_unit, '(*(g0,","))', advance = "no") data%x(1:data%n-1)
    write (file_unit,'(g0)') data%x(data%n)

    close(file_unit)

  end subroutine csv_file_write_vector

  !> Writes a `matrix_t` object to an output file. If the parameter `t`
  !! is present, it will be appended at the start of **each row**.
  !! @param f csv file in which to write.
  !! @param data Matrix to write.
  !! @param Time.
  subroutine csv_file_write_matrix(f, data, t)
    class(csv_file_t), intent(inout) :: f
    type(matrix_t), intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    integer :: file_unit, i, ierr

    open(file = trim(f%fname), position = "append", iostat = ierr, &
         newunit = file_unit)
    if (ierr .ne. 0) call neko_error("Error while opening " // trim(f%fname))

    ! write header if not empty and if not already written
    if (f%header .ne. "" .and. .not. f%header_is_written) then
       write (file_unit, '(A)') trim(f%header)
       f%header_is_written = .true.
    end if

    do i = 1, data%nrows
       if (present(t)) write (file_unit, '(g0,",")', advance = "no") t
       write (file_unit, '(*(g0,","))', advance = "no") &
            data%x(i, 1:data%ncols-1)
       write (file_unit, '(g0)') data%x(i, data%ncols)
    end do

    close(file_unit)

  end subroutine csv_file_write_matrix

  !> Reads data from an input file.
  !! @param this csv file in which to read.
  !! @param data `matrix_t` or `vector_t`, will contain the read data.
  subroutine csv_file_read(this, data)
    class(csv_file_t) :: this
    class(*), target, intent(inout) :: data
    type(vector_t), pointer :: vec
    type(matrix_t), pointer :: mat

    call this%check_exists()

    nullify(vec)
    nullify(mat)

    select type (data)
    type is (vector_t)
       vec => data
       if (.not. allocated(data%x)) then
          call neko_error("Vector is not allocated. Use &
&vector%init() to associate your array &
&with a vector_t object")
       end if

    type is (matrix_t)
       mat => data
       if (.not. allocated(data%x)) then
          call neko_error("Matrix is not allocated. Use &
&matrix%init() to associate your array &
&with a matrix_t object")
       end if


    class default
       call neko_error("Invalid data type for csv_file (expected: vector_t, &
&matrix_t)")
    end select

    if (pe_rank .eq. 0) then

       call neko_log%newline()
       call neko_log%message("Reading csv file " // trim(this%fname))
       if (associated(vec)) then
          call csv_file_read_vector(this, vec)
       else if (associated(mat)) then
          call csv_file_read_matrix(this, mat)
       end if

    end if

  end subroutine csv_file_read

  !> Read a vector (i.e. data on a single row) from a csv file
  !! @param d csv file from which to read.
  !! @param vec Vector object in which to store the file contents.
  !! @note - Equivalent to reading a matrix with `nrows = 1` and `ncols = N`.
  !! To read a vector in a column format, use a `matrix_t` with `ncols = 1`.
  !! @note - If the number of lines in the file is larger than 1,
  !! it will be assumed that a one-line header is present.
  subroutine csv_file_read_vector(f, vec)
    type(csv_file_t), intent(inout) :: f
    type(vector_t), intent(inout) :: vec
    integer :: ierr, file_unit, n_lines
    character(len=80) :: tmp

    n_lines = f%count_lines()

    open(file = trim(f%fname), status = 'old', newunit = file_unit, &
         iostat = ierr)
    if (ierr .ne. 0) call neko_error("Error while opening " // trim(f%fname))

    ! If there is more than 1 line, assume that means there is a header
    if (n_lines .gt. 1) then
       read (file_unit, '(A)') tmp
       f%header = trim(tmp)
    end if

    read (file_unit,*) vec%x
    close(unit = file_unit)


  end subroutine csv_file_read_vector

  !> Read a matrix from a csv file.
  !! @param d csv file from which to read.
  !! @param vec Matrix object in which to store the file contents.
  !! @note If the number of lines in the file is larger than the number
  !! of rows of `mat`, it will be assumed that a one-line header is present.
  subroutine csv_file_read_matrix(f, mat)
    type(csv_file_t), intent(inout) :: f
    type(matrix_t), intent(inout) :: mat
    integer :: ierr, file_unit, i, n_lines
    character(len=80) :: tmp

    n_lines = f%count_lines()

    open(file = trim(f%fname), status = 'old', newunit = file_unit, &
         iostat = ierr)
    if (ierr .ne. 0) call neko_error("Error while opening " // trim(f%fname))

    ! If the number of lines is larger than the number of rows in the
    ! matrix, assume that means there is a header
    if (n_lines .gt. mat%nrows) then
       read (file_unit, '(A)') tmp
       f%header = trim(tmp)
    end if

    do i = 1, mat%nrows
       read (file_unit,*) mat%x(i,:)
    end do
    close(unit = file_unit)

  end subroutine csv_file_read_matrix

  !> Sets the header for a csv file. For example: `hd = "u,v,w,p"`.
  !! @param hd Header.
  !! @note The header will be written "as is", meaning there will be no
  !! checks performed on the header separators, number of columns, etc.
  subroutine csv_file_set_header(this, hd)
    class(csv_file_t), intent(inout) :: this
    character(len=*), intent(in) :: hd

    this%header = trim(hd)

  end subroutine csv_file_set_header

  !> Count the number of lines in a file by going through it entirely
  !! until the end is reached.
  function csv_file_count_lines(this) result(n)
    class(csv_file_t), intent(in) :: this

    integer :: n
    integer :: ierr, file_unit

    call this%check_exists()

    open(file = trim(this%fname), status = 'old', newunit = file_unit, &
         iostat = ierr)
    if (ierr .ne. 0) call neko_error("Error while opening " // trim(this%fname))
    rewind(file_unit)

    n = 0

    ! Keep reading (ierr = 0) until we reach the end (ierr != 0)
    do
       read (file_unit, *, iostat = ierr)
       if (ierr .ne. 0) exit
       n = n + 1
    end do
    rewind(file_unit)
    close(unit = file_unit)

  end function csv_file_count_lines


end module csv_file
