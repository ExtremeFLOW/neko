!> Program to sum up averaged fields computed for statistics and mean field
!! Martin Karp 27/01-23
program average_fields
  use neko
  implicit none
  
  character(len=NEKO_FNAME_LEN) :: inputchar, mesh_fname, fld_fname, output_fname
  type(file_t) :: fld_file, part_file, output_file
  real(kind=rp) :: start_time
  type(fld_file_data_t) :: fld_data, fld_data_avg
  integer :: argc, i
  
  argc = command_argument_count()

  if ((argc .lt. 3) .or. (argc .gt. 3)) then
     if (pe_rank .eq. 0) then
     write(*,*) 'Usage: ./average_fields field_series_name.fld start_time output_name.fld' 
     write(*,*) 'Example command: ./average_fields mean_field104.fld 103.2 mean_field_avg.fld'
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
  fld_file = file_t(trim(fld_fname))
  call get_command_argument(2, inputchar) 
  read(inputchar, *) start_time
  call get_command_argument(3, inputchar) 
  read(inputchar, *) output_fname

  call fld_data%init()
  call fld_data_avg%init()

  call fld_file%read(fld_data_avg)
  call fld_data_avg%scale(fld_data_avg%time-start_time)

  do i = 1, fld_data_avg%meta_nsamples-1
     call fld_file%read(fld_data)
     if (pe_rank .eq. 0) write(*,*) 'Reading file:', i+1
     call fld_data%scale(fld_data%time-fld_data_avg%time)
     call fld_data_avg%add(fld_data)

     if (pe_rank .eq. 0) write(*,*) 'dt', fld_data%time - fld_data_avg%time
     fld_data_avg%time = fld_data%time
  end do
  call fld_data_avg%scale(1.0_rp/(fld_data_avg%time-start_time))

  output_file = file_t(trim(output_fname))
  
  call output_file%write(fld_data_avg, fld_data_avg%time)
  
  call neko_finalize

end program average_fields
