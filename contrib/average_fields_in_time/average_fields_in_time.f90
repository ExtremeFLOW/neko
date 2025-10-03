!> Program to sum up averaged fields computed for statistics and mean field
!! Martin Karp 27/01-23
program average_fields_in_time
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, mesh_fname, fld_fname, output_fname
  type(file_t) :: fld_file, part_file, output_file
  real(kind=rp) :: start_time
  type(fld_file_data_t) :: fld_data, fld_data_avg
  integer :: argc, i

  character(len=LOG_SIZE) :: log_buf

  argc = command_argument_count()

  if ((argc .lt. 3) .or. (argc .gt. 3)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./average_fields_in_time field_series_name.fld start_time output_name.fld'
        write(*,*) 'Example command: ./average_fields_in_time mean_field104.fld 103.2 mean_field_avg.fld'
        write(*,*) 'Computes the average field over the fld files described in mean_field104.nek5000'
        write(*,*) 'The start time is the time at which the first file startsto collect stats'
        write(*,*) 'The files need to be aranged chronological order.'
        write(*,*) 'The average field is then stored in a fld series, i.e. output_name.nek5000 and output_name.f00000'
     end if
     stop
  end if

  call neko_init

  call get_command_argument(1, inputchar)
  read(inputchar, *) fld_fname
  call fld_file%init(trim(fld_fname))
  call get_command_argument(2, inputchar)
  read(inputchar, *) start_time
  call get_command_argument(3, inputchar)
  read(inputchar, *) output_fname

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

  call neko_finalize

end program average_fields_in_time
