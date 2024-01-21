!> Program to interpolate the field from GLL points onto equidistant points in elements in a homogenous direction
!> The number of points inside an element of the input and output field file is supposed to be the same
!> Shiyu Du 16/01-24 note: the program is designed for Cartesian Box meshes with equidistant elements on the hom-dir,
!                          the version for rotational coordinates is not yet been thought about.
!> Shiyu Du 16/01-24 note: the program is currently tested on CPUs, GPU is not tested yet.
!! Shiyu Du 16/01-24
program oneD_equidistant_field
  use neko
  use fast3d
  use tensor
  implicit none
  
  character(len=NEKO_FNAME_LEN) :: inputchar, mesh_fname, field_fname, hom_dir, output_fname
  type(file_t) :: field_file, output_file, mesh_file
  real(kind=rp) :: x_equid
  real(kind=rp), allocatable :: wt(:,:), wtt(:,:), ident(:,:)
  type(fld_file_data_t) :: field_data
  type(space_t) :: Xh
  type(mesh_t) :: msh !only needed to count the elements
  type(vector_ptr_t), allocatable :: fields(:)
  integer :: argc, i, lx, j
  
  argc = command_argument_count()

  if ((argc .lt. 4) .or. (argc .gt. 4)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./oneD_equidistant_field mesh.nmsh field.fld dir(x, y, z) outfield.fld' 
        write(*,*) 'Example command: ./oneD_equidistant_field mesh.nmsh fieldblabla.fld x outfield.fld'
        write(*,*) 'Redistribute the points to be equidistant in elements '
        write(*,*) 'in x of the field fieldblabla.nek5000 and stores in outfield.fld'
     end if
     stop
  end if
  
  call neko_init 

  call get_command_argument(1, inputchar) 
  read(inputchar, *) mesh_fname
  mesh_file = file_t(trim(mesh_fname))
  call get_command_argument(2, inputchar) 
  read(inputchar, *) field_fname
  field_file = file_t(trim(field_fname))
  call get_command_argument(3, inputchar) 
  read(inputchar, *) hom_dir
  call get_command_argument(4, inputchar) 
  read(inputchar, *) output_fname

  if (trim(hom_dir) .ne. 'x' .and. trim(hom_dir) .ne. 'y' .and. trim(hom_dir) .ne. 'z') then
     call neko_error('The homogenous direction should be "x", "y" or "z"')
  end if
  
  call mesh_file%read(msh)

  ! read field at t=0
  call field_data%init(msh%nelv,msh%offset_el)
  if (pe_rank .eq. 0) write(*,*) 'Reading file:', 1
  call field_file%read(field_data)

  lx = field_data%lx

  call Xh%init(GLL, field_data%lx, field_data%ly, field_data%lz)

  ! construct the weight matrix
  allocate(wt(lx,lx))
  allocate(wtt(lx,lx))
  allocate(ident(lx,lx))
  do i = 1, lx
    x_equid = -1.0_rp + (i-1) * 2.0_rp/(lx-1)
    call fd_weights_full(x_equid, Xh%zg(:,1), lx-1, 0, wtt(:,i))
    wt(i,:) = wtt(:,i)
    ident(i,i) = 1.0_rp
  end do

  ! redistribute the coordinates to be equidistant within elements
  if (trim(hom_dir) .eq. 'x') then
    call tnsr1_3d(field_data%x%x, lx, lx, wt, ident, ident, msh%nelv)
  else if (trim(hom_dir) .eq. 'y') then
    call tnsr1_3d(field_data%y%x, lx, lx, ident, wtt, ident, msh%nelv)
  else if (trim(hom_dir) .eq. 'z') then
    call tnsr1_3d(field_data%z%x, lx, lx, ident, ident, wtt, msh%nelv)
  end if

  ! interpolate the field at t=0
  allocate(fields(field_data%size()))
  call field_data%get_list(fields,field_data%size())
  do i = 1, field_data%size()
     if (trim(hom_dir) .eq. 'x') then
        call tnsr1_3d(fields(i)%v%x, lx, lx, wt, ident, ident, msh%nelv)
     else if (trim(hom_dir) .eq. 'y') then
        call tnsr1_3d(fields(i)%v%x, lx, lx, ident, wtt, ident, msh%nelv)
     else if (trim(hom_dir) .eq. 'z') then
        call tnsr1_3d(fields(i)%v%x, lx, lx, ident, ident, wtt, msh%nelv)
     end if
  end do 

  ! output at t=0  
  output_file = file_t(trim(output_fname))
  call output_file%write(field_data, field_data%time)

  ! interpolate field for t>0
  do i = 1, field_data%meta_nsamples-1
     if (pe_rank .eq. 0) write(*,*) 'Reading file:', i+1
     call field_file%read(field_data)
     call field_data%get_list(fields,field_data%size())
     do j = 1, field_data%size()
        if (trim(hom_dir) .eq. 'x') then
           call tnsr1_3d(fields(j)%v%x, lx, lx, wt, ident, ident, msh%nelv)
        else if (trim(hom_dir) .eq. 'y') then
           call tnsr1_3d(fields(j)%v%x, lx, lx, ident, wtt, ident, msh%nelv)
        else if (trim(hom_dir) .eq. 'z') then
           call tnsr1_3d(fields(j)%v%x, lx, lx, ident, ident, wtt, msh%nelv)
        end if
     end do 
  end do
  ! output for t>0
  call output_file%write(field_data, field_data%time)
  
  if (pe_rank .eq. 0) write(*,*) 'Done'
  
  call neko_finalize

end program oneD_equidistant_field
