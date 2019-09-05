program rea2nbin
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: fname, output
  type(mesh_t) :: msh
  type(file_t) :: rea_file, nmsh_file
  
  if (command_argument_count() .lt. 2) then
     write(*,*) 'Usage: ./rea2nbin <reafile> <neko mesh>'
     stop
  end if
  
  call neko_init 
  
  call get_command_argument(1, fname) 
  call get_command_argument(2, output)
  
  rea_file = file_t(fname)
  
  call rea_file%read(msh)
  
  nmsh_file = file_t(output)
  
  call nmsh_file%write(msh)
  
  call mesh_free(msh)
  
  call neko_finalize

end program rea2nbin
