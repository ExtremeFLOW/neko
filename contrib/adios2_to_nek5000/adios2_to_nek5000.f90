!> Program to convert from adios2 to nek5000 file format.
program adios2_to_nek5000
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, field_fname, output_fname
  type(file_t) :: field_file, output_file
  type(fld_file_data_t) :: field_data
  integer :: argc, file_precision, i, layout
  logical :: dp_precision

  argc = command_argument_count()

  if ((argc .lt. 4) .or. (argc .gt. 4)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./adios2_to_nek5000 field0.adios2 field0.nek5000 precision'
        write(*,*) 'Example command: ./map_to_equidistant_1d field0.adios2 field0.nek5000 .true.'
        write(*,*) 'Convert adios2 to nek5000 file format'
     end if
     stop
  end if

  call neko_init

  call get_command_argument(1, inputchar)
  read(inputchar, fmt='(A)') field_fname
  call get_command_argument(2, inputchar)
  read(inputchar, fmt='(A)') output_fname
  call get_command_argument(3, inputchar)
  read(inputchar, *) dp_precision
  call get_command_argument(4, inputchar)
  read(inputchar, *) layout

  if (dp_precision) then
     file_precision = dp
  else
     file_precision = sp
  end if

  field_file = file_t(trim(field_fname),precision=file_precision)

  call field_data%init()
  if (pe_rank .eq. 0) write(*,*) 'Reading file:', 1
  call field_file%read(field_data)

  ! output at t=0
  output_file = file_t(trim(output_fname),precision=file_precision, layout=layout)
  call output_file%write(field_data, field_data%time)

  ! output for t>0
  do i = 1, field_data%meta_nsamples-1
     if (pe_rank .eq. 0) write(*,*) 'Reading file:', i+1
     call field_file%read(field_data)
     ! output for t>0
     call output_file%write(field_data, field_data%time)
  end do

  if (pe_rank .eq. 0) write(*,*) 'Done'

  call neko_finalize

end program adios2_to_nek5000
