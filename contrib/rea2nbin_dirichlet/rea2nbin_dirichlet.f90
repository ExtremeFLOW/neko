program rea2nbin_dirichlet
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: fname_in, fname_out
  type(mesh_t) :: msh
  type(file_t) :: nmsh_file, rea_dirichlet_file
  integer :: argc
  character(len=80) :: suffix
  logical :: remove_at_end = .false.

  argc = command_argument_count()

  if ((argc .lt. 1) .or. (argc .gt. 2)) then
     write(*,*) 'Usage: rea2nbin_dirichlet <reafile> <neko mesh>'
     stop
  end if
  
  call neko_init 
  
  !
  ! Parse first argument: rea file
  ! If the user gives a .rea file, we work
  ! from a temporary copy with extension .readirichlet
  !
  call get_command_argument(1, fname_in) 
  call filename_suffix(fname_in, suffix)

  if (suffix .eq. "rea") then
     write (*,*) ".rea format specified, creating temporary copy .readirichlet"
     call execute_command_line("cp "//trim(fname_in)//" temp.readirichlet")
     fname_in = "temp.readirichlet"
     remove_at_end = .true.
  else if (suffix .ne. "readirichlet") then
     call neko_error("Invalid input format")
  end if

  ! 
  ! Parse second argument: nmsh file
  ! If not specified, default to dirichlet.nmsh
  !
  if (argc .eq. 2) then
     call get_command_argument(2, fname_out)
     call filename_suffix(fname_out, suffix)

     if (suffix .ne. "nmsh") then
        call neko_error("Invalid output format")
     end if
  else
     fname_out = "dirichlet.nmsh"
  end if
 
  ! 
  ! Read nmsh file with labels
  !
  rea_dirichlet_file = file_t(fname_in)
  msh%lgenc = .false.  
  call rea_dirichlet_file%read(msh)
  
  !
  ! Write nmsh file with dirichlet
  !
  nmsh_file = file_t(fname_out)
  call nmsh_file%write(msh)
  call msh%free()
  
  call neko_finalize

  if (remove_at_end) call execute_command_line("rm temp.readirichlet")

end program rea2nbin_dirichlet
