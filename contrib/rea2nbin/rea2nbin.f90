program rea2nbin
  use num_types
  use mesh
  implicit none

  integer :: i, j, ierr
  integer :: ndim, nparam, nskip, nlogic
  integer :: nelgs, nelgv
  character(len=80) :: fname, output, opt
  type(mesh_t) :: msh
  
  if (iargc() .lt. 2) then
     write(*,*) 'Usage: ./rea2nbin <reafile> <neko mesh>'
     stop
  end if
  
  call getarg(1, fname)
  call getarg(2, output)

  open(unit=9,file=trim(fname), status='old', iostat=ierr)
  write(*, '(A,A)') ' Reading ', fname
  
  read(9, *)
  read(9, *)
  read(9, *) ndim
  read(9, *) nparam
  
  ! Skip parameters
  do i = 1, nparam
     read(9, *)
  end do

  ! Skip passive scalars
  read(9, *) nskip
  do i = 1, nskip
     read(9, *)
  end do

  ! Skip logic switches
  read(9, *) nlogic
  do i = 1, nlogic
     read(9, *)
  end do

  ! Read mesh info
  read(9, *)
  read(9, *)
  read(9, *) nelgs,ndim, nelgv
  write(*,*) nelgs, ndim, nelgv

  call mesh_init_coordinates(msh, ndim, nelgv)

  do i = 1, nelgv
     read(9, *)
     if (ndim .eq. 2) then
        read(9, *) (msh%xc(j, i),j=1,4)
        read(9, *) (msh%yc(j, i),j=1,4)
     else if (ndim .eq. 3) then
        read(9, *) (msh%xc(j, i),j=1,4)
        read(9, *) (msh%yc(j, i),j=1,4)
        read(9, *) (msh%zc(j, i),j=1,4)
        read(9, *) (msh%xc(j, i),j=5,8)
        read(9, *) (msh%yc(j, i),j=5,8)
        read(9, *) (msh%zc(j, i),j=5,8)
     end if
  end do

  call mesh_free(msh)

  close(9)

end program rea2nbin
