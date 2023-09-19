program rea2nbin
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: fname, output_
  type(mesh_t) :: msh
  type(file_t) :: rea_file, nmsh_file
  integer :: argc
  character(len=80) :: suffix
  
  argc = command_argument_count()

  if ((argc .lt. 1) .or. (argc .gt. 2)) then
     write(*,*) 'Usage: ./rea2nbin <reafile> <neko mesh>'
     stop
  end if
  
  call neko_init 
  
  call get_command_argument(1, fname) 

  if (argc .eq. 2) then
     call get_command_argument(2, output_)
     call filename_suffix(output_, suffix)

     if (suffix .ne. "nmsh") then
        call neko_error("Invalid output format")
     end if
  else
     call filename_chsuffix(fname, output_, 'nmsh')
  end if
     
  
  rea_file = file_t(fname)
    
  msh%lgenc = .false.
  
  call rea_file%read(msh)
  
  nmsh_file = file_t(output_)
  
  call nmsh_file%write(msh)
  
  call msh%free()
  
  call neko_finalize

end program rea2nbin
