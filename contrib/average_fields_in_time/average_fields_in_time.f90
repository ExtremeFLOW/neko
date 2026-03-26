!> Program to sum up averaged fields computed for statistics and mean field
!! Martin Karp 27/01-23
!! Victor Baconnet 26/03/26
program average_fields_in_time
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: in_fname, out_fname, inputchar
  character(len=80) :: extension
  real(kind=rp) :: start_time
  integer :: argc, i

  argc = command_argument_count()

  if ((argc .lt. 3) .or. (argc .gt. 3)) call usage()

  call neko_init

  call get_command_argument(1, inputchar)
  read(inputchar, *) in_fname
  call get_command_argument(2, inputchar)
  read(inputchar, *) start_time
  call get_command_argument(3, inputchar)
  read(inputchar, *) out_fname

  call filename_suffix(in_fname, extension)
  select case (trim(extension))
  case ("csv")
     call avg_flds_in_time_csv(in_fname, out_fname, start_time)
  case ("fld")
     call avg_flds_in_time_fld(in_fname, out_fname, start_time)
  case default
     block
       character(len=LOG_SIZE) :: log_buf
       write (log_buf,*) extension, "file format is not supported!"
       call neko_log%message(log_buf)
     end block

     call usage()
  end select

  call neko_finalize

contains

  subroutine usage()
    if (pe_rank .eq. 0) then
       write(*,*) 'Usage: ./average_fields_in_time field_series_name.fld start_time output_name.fld'
       write(*,*) 'Example command: ./average_fields_in_time mean_field104.fld 103.2 mean_field_avg.fld'
       write(*,*) 'Computes the average field over the fld files described in mean_field104.nek5000'
       write(*,*) 'The start time is the time at which the first file startsto collect stats'
       write(*,*) 'The files need to be aranged chronological order.'
       write(*,*) 'The average field is then stored in a fld series, i.e. output_name.nek5000 and output_name.f00000'
    end if
    stop
  end subroutine usage

  subroutine avg_flds_in_time_fld(fld_fname, output_fname, start_time)
    character(len=NEKO_FNAME_LEN), intent(in) :: fld_fname, output_fname
    real(kind=rp), intent(in) :: start_time

    type(file_t) :: fld_file, output_file
    type(fld_file_data_t) :: fld_data, fld_data_avg
    character(len=LOG_SIZE) :: log_buf
    integer :: i

    call fld_file%init(trim(fld_fname))

    call fld_data%init()
    call fld_data_avg%init()

    call fld_file%read(fld_data_avg)

    write (log_buf, '(A, g0)') "dt: ", fld_data_avg%time - start_time
    call neko_log%message(log_buf)

    call fld_data_avg%scale(fld_data_avg%time-start_time)

    do i = 1, fld_data_avg%meta_nsamples-1
       call fld_file%read(fld_data)
       call fld_data%scale(fld_data%time-fld_data_avg%time)
       call fld_data_avg%add(fld_data)

       write (log_buf, '(A, g0)') "dt: ", fld_data%time - fld_data_avg%time
       call neko_log%message(log_buf)

       fld_data_avg%time = fld_data%time
    end do
    call fld_data_avg%scale(1.0_rp/(fld_data_avg%time-start_time))

    call output_file%init(trim(output_fname))

    call neko_log%message('Writing file: ' // trim(output_fname))
    call output_file%write(fld_data_avg, fld_data_avg%time)
    call neko_log%message('Done')

  end subroutine avg_flds_in_time_fld

  !> Do the time averaging on a csv file
  subroutine avg_flds_in_time_csv(in_fname, out_fname, start_time)
    character(len=NEKO_FNAME_LEN), intent(in) :: in_fname, out_fname
    real(kind=rp), intent(in) :: start_time

    type(file_t) :: out_file
    type(matrix_t) :: data, avg_data
    character(len=LOG_SIZE) :: log_buf
    integer :: sample_size, n_samples

    ! This will initialize and populate the matrix data
    call populate_matrix(trim(in_fname), data)

    ! Compute the size of a statistics sample
    sample_size = determine_size_of_csv_sample(data)
    n_samples = data%get_nrows() / sample_size
    write (log_buf, '(A,I5)') "Size of each sample: ", sample_size
    call neko_log%message(log_buf)
    write (log_buf, '(A,I5)') "Number of samples  : ", n_samples
    call neko_log%message(log_buf)

    ! Initialize the matrix for the averaged data
    call avg_data%init(sample_size, data%get_ncols())

    block
      integer :: i
      real(kind=rp) :: dt, ta, tb

      tb = data%x(1,1) ! First time stamp
      ta = start_time ! User defined
      dt = tb - ta ! First "sampling length"
      write (log_buf, '(A,I0,A,I0,A,g0)') "Length of sample ", 1, "/", &
           n_samples, ": ", dt
      call neko_log%message(log_buf)

      ! First sample, scaled by the starting time
      avg_data%x(:,1:) = dt*data%x(1:sample_size, :)

      do i = 1, n_samples-1

         ta = data%x(sample_size*(i-1) + 1, 1) ! Previous time stamp
         tb = data%x(sample_size* i + 1, 1) ! Current time stamp
         dt = tb - ta ! Sampling time for the current sample

         write (log_buf, '(A,I0,A,I0,A,g0)') "Length of sample ", i+1, "/", &
              n_samples, ": ", dt
         call neko_log%message(log_buf)

         avg_data%x = avg_data%x + dt * data%x(sample_size*i + 1:sample_size * (i+1), :)
      end do

      ! Final rescaling with the total length of signal
      avg_data%x = avg_data%x / (data%x(data%get_nrows(), 1) - start_time)

      ! Set the final time
      avg_data%x(:,1) = data%x(data%get_nrows(), 1)
    end block

    ! Write the final averaged fields
    call out_file%init(trim(out_fname))
    call out_file%set_overwrite(.true.)
    call out_file%write(avg_data)
    call file_free(out_file)

  end subroutine avg_flds_in_time_csv

  !> Finds the size of a statistics sample based on when the time stamp changes
  function determine_size_of_csv_sample(m) result(n)
    type(matrix_t), intent(in) :: m
    integer :: n

    n = 1
    ! Increment n while the time stamps (in the first column) have the same value
    do while ( abscmp(m%x(n,1), m%x(n+1,1)) )
       n = n + 1
       if (n+1 > m%get_nrows()) &
            call neko_error("Out of bounds while searching for sampling size!")
    end do

  end function determine_size_of_csv_sample

  !> Read a csv file, determine the # of rows and columns, allocate the matrix
  !! accordingly and populate the matrix with the data from the csv file.
  subroutine populate_matrix(file_name, m)
    character(len=*), intent(in) :: file_name
    type(matrix_t), intent(inout) :: m

    type(file_t) :: f
    integer :: unit, i, ierr, num_columns, num_lines, ios
    character(len=1000) :: line
    character(len=1) :: delimiter
    character(len=LOG_SIZE) :: log_buf
    delimiter = ','

    if (pe_rank .eq. 0) then

       !
       ! Count lines and columns for initialization
       !
       open(unit=unit, file=trim(file_name), status='old', action='read', iostat=ios)
       if (ios /= 0) then
          call neko_error("Error opening file " // trim(file_name))
       end if

       num_columns = 1
       num_lines = 0

       ! Read the file line by line
       do
          read(unit, '(A)', iostat=ios) line
          if (ios /= 0) exit

          ! If it's the first line, count the columns
          if (num_lines .eq. 0) then

             ! Count the number of delimiters in the line
             do i = 1, len_trim(line)
                if (line(i:i) == delimiter) then
                   num_columns = num_columns + 1
                end if
             end do

          end if ! if num_columns .eq. 1

          num_lines = num_lines + 1
       end do

       ! Close the file
       close(unit)

    end if

    write (log_buf, '(A,A,A,I5,I5)') "Reading file ", trim(file_name), " of size ", num_lines, num_columns
    call neko_log%message(log_buf)

    !
    ! Reading of the inflow profile
    !
    call m%init(num_lines, num_columns)

    ! Read csv file and broadcast to all ranks since csv is only read in serial
    call f%init(trim(file_name))
    call f%read(m)
    call file_free(f)

  end subroutine populate_matrix


end program average_fields_in_time
