program rea2nbin
  use num_types
  use mesh
  use file
  implicit none

  integer :: i, j, ierr
  integer :: ndim, nparam, nskip, nlogic
  integer :: nelgs, nelgv
  character(len=80) :: fname, output, opt
  type(mesh_t) :: msh
  type(file_t) :: rea_file
  
  if (command_argument_count() .lt. 2) then
     write(*,*) 'Usage: ./rea2nbin <reafile> <neko mesh>'
     stop
  end if
  
  call get_command_argument(1, fname)
  call get_command_argument(2, output)

  rea_file = file_t(fname)

  call rea_file%read(msh)

  call mesh_free(msh)

  close(9)

end program rea2nbin
