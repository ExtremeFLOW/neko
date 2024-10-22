!> Program to sum up a field in space in a homogenous direction
!! Martin Karp 17/02-23
program average_field_in_space
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, mesh_fname, field_fname, hom_dir, output_fname
  type(file_t) :: field_file, output_file, mesh_file
  real(kind=rp) :: start_time, el_h, el_dim(3,3), domain_height
  real(kind=rp), allocatable :: temp_el(:,:,:)
  type(fld_file_data_t) :: field_data
  type(fld_file_data_t) :: output_data
  type(coef_t) :: coef
  type(dofmap_t), target :: dof
  type(space_t) :: Xh
  type(mesh_t) :: msh
  type(gs_t) :: gs_h
  type(map_1d_t) :: map_1d
  type(map_2d_t) :: map_2d
  type(vector_ptr_t), allocatable :: fields(:), fields2d(:)
  type(matrix_t) :: avg_matrix
  type(vector_t) :: volume_per_gll_lvl
  integer :: argc, i, n, lx, j, e, n_levels, dir, ierr, n_1d, tstep, k
  logical :: avg_to_1d = .false.
  integer :: nelv_2d, glb_nelv_2d, offset_el_2d, lxy, n_2d
  integer, allocatable :: idx_2d(:)
  real(kind=rp), pointer, dimension(:,:,:,:) :: x_ptr, y_ptr
  real(kind=rp) :: coord

  argc = command_argument_count()

  if ((argc .lt. 4) .or. (argc .gt. 4)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./average_field_in_space mesh.nmsh field.fld dir(x, y, z, xz, xy, yz)  outfield.(fld,csv)'
        write(*,*) '----'
        write(*,*) 'Example command for avg in 1 direction: ./average_field_in_space mesh.nmsh fieldblabla.fld x  outfield.fld'
        write(*,*) 'Computes spatial average in 1 direction and saves it in outfield.fld'
        write(*,*) '----'
        write(*,*) 'Example command: ./average_field_in_space mesh.nmsh fieldblabla.fld xy  out.csv'
        write(*,*) 'Computes the spatial average in 2 directions directly of the field fieldblabla.nek5000 and stores in out.csv'
        write(*,*) '----'
        write(*,*) 'In out.csv the first col are the coords of the GLL points'
        write(*,*) 'In columns 2-n_fields are the averages for all fields in fieldblabla.fld'
        write(*,*) 'If averaging in 2 directions output must be .csv and for 1 direction .fld'
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

  call mesh_file%read(msh)

  call field_data%init(msh%nelv,msh%offset_el)
  call field_file%read(field_data)

  lx = field_data%lx
  !To make sure any deformation made in the user file is passed onto here as well
  do i = 1,msh%nelv
     msh%elements(i)%e%pts(1)%p%x(1) = field_data%x%x(linear_index(1,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(2)%p%x(1) = field_data%x%x(linear_index(lx,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(3)%p%x(1) = field_data%x%x(linear_index(1,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(4)%p%x(1) = field_data%x%x(linear_index(lx,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(5)%p%x(1) = field_data%x%x(linear_index(1,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(6)%p%x(1) = field_data%x%x(linear_index(lx,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(7)%p%x(1) = field_data%x%x(linear_index(1,lx,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(8)%p%x(1) = field_data%x%x(linear_index(lx,lx,lx,i,lx,lx,lx))

     msh%elements(i)%e%pts(1)%p%x(2) = field_data%y%x(linear_index(1,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(2)%p%x(2) = field_data%y%x(linear_index(lx,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(3)%p%x(2) = field_data%y%x(linear_index(1,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(4)%p%x(2) = field_data%y%x(linear_index(lx,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(5)%p%x(2) = field_data%y%x(linear_index(1,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(6)%p%x(2) = field_data%y%x(linear_index(lx,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(7)%p%x(2) = field_data%y%x(linear_index(1,lx,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(8)%p%x(2) = field_data%y%x(linear_index(lx,lx,lx,i,lx,lx,lx))

     msh%elements(i)%e%pts(1)%p%x(3) = field_data%z%x(linear_index(1,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(2)%p%x(3) = field_data%z%x(linear_index(lx,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(3)%p%x(3) = field_data%z%x(linear_index(1,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(4)%p%x(3) = field_data%z%x(linear_index(lx,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(5)%p%x(3) = field_data%z%x(linear_index(1,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(6)%p%x(3) = field_data%z%x(linear_index(lx,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(7)%p%x(3) = field_data%z%x(linear_index(1,lx,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(8)%p%x(3) = field_data%z%x(linear_index(lx,lx,lx,i,lx,lx,lx))
  end do

  call Xh%init(GLL, field_data%lx, field_data%ly, field_data%lz)

  call dof%init(msh, Xh)
  call gs_h%init(dof)
  call coef%init(gs_h)

  ! 1 corresponds to x, 2 to y, 3 to z
  if (trim(hom_dir) .eq. 'x') then
     dir = 1
     avg_to_1d = .false.
  else if (trim(hom_dir) .eq. 'y') then
     dir = 2
     avg_to_1d = .false.
  else if (trim(hom_dir) .eq. 'z') then
     dir = 3
     avg_to_1d = .false.
  else if (trim(hom_dir) .eq. 'yz') then
     dir = 1
     avg_to_1d = .true.
  else if (trim(hom_dir) .eq. 'xz') then
     dir = 2
     avg_to_1d = .true.
  else if (trim(hom_dir) .eq. 'xy') then
     dir = 3
     avg_to_1d = .true.
  else
     call neko_error('homogenous direction not supported')
  end if
  if (avg_to_1d) then
     call map_1d%init(coef, dir, 1e-7_rp)
  else
     call map_2d%init(coef, dir, 1e-7_rp)
  end if
  
  !allocate array with pointers to all vectors in the file

  output_file = file_t(trim(output_fname))
  do tstep = 0, field_data%meta_nsamples-1
    if (pe_rank .eq. 0) write(*,*) 'Averaging field:', tstep
     if (tstep .gt. 0) call field_file%read(field_data)
     if (avg_to_1d) then
        call map_1d%average_planes(avg_matrix, fields)
        call output_file%write(avg_matrix,field_data%time)
        ! Compute averages in 1 direction and store in a 3d field (lots of redundant data, sorry)
        ! Should output a 2d field in principle
     else
        call map_2d%average(output_data,field_data)
        call output_file%write(output_data,field_data%time)
     end if
  end do
  if (pe_rank .eq. 0) write(*,*) 'Done'
  call neko_finalize

end program average_field_in_space

