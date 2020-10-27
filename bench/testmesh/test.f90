program mfile
  use neko
  implicit none
  type(map_t) :: m
  type(file_t) :: mf, rf
  type(mesh_t) :: msh
  type(mesh_fld_t) :: mff
  character(len=NEKO_FNAME_LEN) :: fname
  type(tuple4_i4_t) :: t4
  type(tuple_i4_t) :: t2, facet_map
  type(space_t ):: Xh
  type(field_t) :: fld
  integer :: i, j
  type(dofmap_t) :: d

  call neko_init


  call get_command_argument(1, fname)
  !fname = '/home/usr/martin/neko/bench/hemi.rea'
  rf = file_t(fname)
  call rf%read(msh)
  call mesh_generate_conn(msh)
  call space_init(Xh, GLL, 8, 8, 8)
  call field_init(fld, msh, Xh, 'types')
  d = dofmap_t(msh, Xh)
  fname = 'dof.vtk'
  mf = file_t(fname)
  print *, 'asdas' 
  call mf%write(d)
  fld%x = 0d0
  do i  = 1, msh%nelv

     if (msh%facet_neigh(1, i) .eq. 0) then
        fld%x(1,1,1,i) = 1
        fld%x(1,8,1,i) = 1

        fld%x(1,1,8,i) = 1
        fld%x(1,8,8,i) = 1
     end if

     if (msh%facet_neigh(2, i) .eq. 0) then
        fld%x(8,1,1,i) = 1
        fld%x(8,8,1,i) = 1

        fld%x(8,1,8,i) = 1
        fld%x(8,8,8,i) = 1
     end if

     if (msh%facet_neigh(3, i) .eq. 0) then
        fld%x(1,1,1,i) = 1
        fld%x(8,1,1,i) = 1

        fld%x(8,1,8,i) = 1
        fld%x(1,1,8,i) = 1
     end if

     if (msh%facet_neigh(4, i) .eq. 0) then
        fld%x(1,8,1,i) = 1
        fld%x(8,8,1,i) = 1

        fld%x(8,8,8,i) = 1
        fld%x(1,8,8,i) = 1
     end if

     if (msh%facet_neigh(5, i) .eq. 0) then
        fld%x(1,1,1,i) = 1
        fld%x(8,1,1,i) = 1

        fld%x(8,8,1,i) = 1
        fld%x(1,8,1,i) = 1
     end if

     if (msh%facet_neigh(6, i) .eq. 0) then
        fld%x(1,1,8,i) = 1
        fld%x(8,1,8,i) = 1

        fld%x(8,8,8,i) = 1
        fld%x(1,8,8,i) = 1
     end if
  end do

  fname = 'fld.vtk'
  mf = file_t(fname)
  call mf%write(fld)  
   
  call neko_finalize
  
end program mfile
