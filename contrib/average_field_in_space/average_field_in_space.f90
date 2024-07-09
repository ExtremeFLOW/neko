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
  type(field_t), pointer :: u, avg_u, old_u, el_heights
  type(vector_ptr_t), allocatable :: fields(:)
  type(matrix_t) :: avg_matrix
  type(vector_t) :: volume_per_gll_lvl
  integer :: argc, i, n, lx, j, e, n_levels, dir, ierr, n_1d, tstep, k
  logical :: avg_to_1d = .false.
  integer :: nelv_2d, glb_nelv_2d, offset_el_2d, lxy, n_2d
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

  dof = dofmap_t(msh, Xh)
  call gs_h%init(dof)
  call coef%init(gs_h)

  call neko_field_registry%add_field(dof, 'u')
  u => neko_field_registry%get_field('u')
  call neko_field_registry%add_field(dof, 'avg_u')
  avg_u => neko_field_registry%get_field('avg_u')
  call neko_field_registry%add_field(dof, 'old_u')
  old_u => neko_field_registry%get_field('old_u')
  call neko_field_registry%add_field(dof, 'el_heights')
  el_heights => neko_field_registry%get_field('el_heights')
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
  call map_1d%init(coef,  dir, 1e-7_rp)
  n_levels = map_1d%n_el_lvls
  n = u%dof%size()

  if ( .not. avg_to_1d) then
     nelv_2d = 0
     do i = 1, msh%nelv
        if (map_1d%el_lvl(i) .eq. 1) then
           nelv_2d = nelv_2d + 1
        end if
     end do
     print *, nelv_2d
     glb_nelv_2d = 0
     call MPI_Allreduce(nelv_2d, glb_nelv_2d, 1, &
          MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

     offset_el_2d = 0
     call MPI_Exscan(nelv_2d, offset_el_2d, 1, &
          MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

     print *, map_1d%n_el_lvls*glb_nelv_2d
     call output_data%init(nelv_2d,offset_el_2d)
     print *, 'bro1'
     allocate(output_data%idx(nelv_2d))
     print *, 'bro2'
     output_data%gdim = 2
     output_data%lx = dof%Xh%lx
     output_data%ly = dof%Xh%ly
     output_data%lz = 1
     output_Data%glb_nelv = glb_nelv_2d
     do i = 1, nelv_2d
        output_data%idx(i) = offset_el_2d + i
     end do
     print *, 'bro3'
     lxy = output_data%lx*output_data%ly
     n_2d = lxy*nelv_2d
     call output_data%x%init(n_2d)
     call output_data%y%init(n_2d)
     j = 0
     if(dir .eq. 1) then
        x_ptr => dof%z
        y_ptr => dof%y
     end if
     if(dir .eq. 2) then
        x_ptr => dof%x
        y_ptr => dof%z
     end if
     if(dir .eq. 3) then
        x_ptr => dof%x
        y_ptr => dof%y
     end if

     do e = 1, msh%nelv
        if (map_1d%el_lvl(e) .eq. 1) then
           if (map_1d%dir_el(e) .eq. 1) then
              do i = 1, lxy
                 output_data%x%x(j*lxy+i) = x_ptr(1,i,1,e)
                 output_data%y%x(j*lxy+i) = y_ptr(1,i,1,e)
              end do
           end if
           if (map_1d%dir_el(e) .eq. 2) then
              do i = 1, lx
                 do k = 1, lx
                    output_data%x%x(j*lxy+k+lx*(i-1)) = x_ptr(k,1,i,e)
                    output_data%y%x(j*lxy+k+lx*(i-1)) = y_ptr(k,1,i,e)
                 end do
              end do
           end if
           if (map_1d%dir_el(e) .eq. 3) then
              do i = 1, lxy
                 output_data%x%x(j*lxy+i) = x_ptr(i,1,1,e)
                 output_data%y%x(j*lxy+i) = y_ptr(i,1,1,e)
              end do
           end if
           j = j +1
        end if
     end do
     call output_data%p%init(n_2d)
     output_data%p = 1.3_rp
     output_file = file_t(trim(output_fname))
     call output_file%write(output_data,13.2_rp)
  end if
  stop
  call field_data%init(msh%nelv,msh%offset_el)
  !allocate array with pointers to all vectors in the file
  allocate(fields(field_data%size()))
  ! Compute average in two direction directly and store in a csv file
  output_file = file_t(trim(output_fname))
  do tstep = 0, field_data%meta_nsamples-1
     if (tstep .gt. 0) call field_file%read(field_data)
     call field_data%get_list(fields,field_data%size())
     if (pe_rank .eq. 0) write(*,*) 'Averaging field:', tstep
     if (avg_to_1d) then
        call map_1d%average_planes(avg_matrix, fields)
        call output_file%write(avg_matrix,field_data%time)
        ! Compute averages in 1 direction and store in a 3d field (lots of redundant data, sorry)
        ! Should output a 2d field in principle
     else
        do i = 1, msh%nelv
           !find height in hom-dir
           !direction in local coords (r,s,t) that is hom is stored in map_1d%dir_el
           !set element to height
           !we assume elements are stacked on eachother...
           el_dim(1,:) = abs(msh%elements(i)%e%pts(1)%p%x-msh%elements(i)%e%pts(2)%p%x)
           el_dim(2,:) = abs(msh%elements(i)%e%pts(1)%p%x-msh%elements(i)%e%pts(3)%p%x)
           el_dim(3,:) = abs(msh%elements(i)%e%pts(1)%p%x-msh%elements(i)%e%pts(5)%p%x)
           ! 1 corresponds to x, 2 to y, 3 to z
           el_h = el_dim(map_1d%dir_el(i),dir)
           el_heights%x(:,:,:,i) = el_h
        end do
        !Need to compute mapping for 3d to 2d
        !Does order matter? Think its ok as long as values written in same order

        call copy(u%x,el_heights%x,n)
        call copy(old_u%x,el_heights%x,n)
        call copy(avg_u%x,el_heights%x,n)
        call perform_global_summation(u, avg_u, old_u, n_levels, &
             map_1d%dir_el,gs_h, coef%mult, msh%nelv, lx)
        domain_height = u%x(1,1,1,1)


        do i = 1, field_data%size()
           call copy(old_u%x,fields(i)%ptr%x,n)
           call perform_local_summation(u,old_u, el_heights, domain_height, &
                map_1d%dir_el, coef, msh%nelv, lx)
           call copy(old_u%x,u%x,n)
           call copy(avg_u%x,u%x,n)
           call perform_global_summation(u, avg_u, old_u, n_levels, &
                map_1d%dir_el,gs_h, coef%mult, msh%nelv, lx)
           call copy(fields(i)%ptr%x,u%x,n)
        end do

        call output_file%write(field_data,field_data%time)
     end if
  end do
  if (pe_rank .eq. 0) write(*,*) 'Done'
  call neko_finalize

end program average_field_in_space

subroutine perform_global_summation(u, avg_u, old_u, n_levels, hom_dir_el, gs_h, mult, nelv,lx)
  use neko
  implicit none
  type(field_t), intent(inout) :: u, avg_u, old_u
  type(gs_t), intent(inout) :: gs_h
  integer, intent(in) :: n_levels, nelv, lx
  integer, intent(in) :: hom_dir_el(nelv)
  real(kind=rp), intent(in) :: mult(nelv*lx**3)
  real(kind=rp) :: temp_el(lx,lx,lx)
  integer :: n, i, j, e, k

  n = u%dof%size()

  do i = 1, n_levels-1
     !compute average
     if (NEKO_BCKND_DEVICE .eq. 1) &
          call device_memcpy(u%x, u%x_d, n, &
                             HOST_TO_DEVICE, sync=.false.)
     call gs_h%op(u,GS_OP_ADD)
     if (NEKO_BCKND_DEVICE .eq. 1) &
          call device_memcpy(u%x, u%x_d, n, DEVICE_TO_HOST, sync=.true.)
     call col2(u%x,mult,n)
     do e = 1, nelv
        temp_el = 2.0*u%x(:,:,:,e)-old_u%x(:,:,:,e)
        if (hom_dir_el(e) .eq. 1) then
           u%x(1,:,:,e) = temp_el(lx,:,:)
           avg_u%x(1,:,:,e) = avg_u%x(1,:,:,e)+temp_el(1,:,:)
           u%x(lx,:,:,e) = temp_el(1,:,:)
        else if (hom_dir_el(e) .eq. 2) then
           u%x(:,1,:,e) = temp_el(:,lx,:)
           avg_u%x(:,1,:,e) = avg_u%x(:,1,:,e)+temp_el(:,1,:)
           u%x(:,lx,:,e) = temp_el(:,1,:)
        else if (hom_dir_el(e) .eq. 3) then
           u%x(:,:,1,e) = temp_el(:,:,lx)
           avg_u%x(:,:,1,e) = avg_u%x(:,:,1,e)+temp_el(:,:,1)
           u%x(:,:,lx,e) = temp_el(:,:,1)
        end if
        old_u%x(:,:,:,e) = u%x(:,:,:,e)
     end do
  end do
  do e = 1, nelv
     do i = 1,lx
        do j = 1, lx
           if (hom_dir_el(e) .eq. 1) then
              u%x(:,i,j,e) = avg_u%x(1,i,j,e)
           else if (hom_dir_el(e) .eq. 2) then
              u%x(i,:,j,e) = avg_u%x(i,1,j,e)
           else if (hom_dir_el(e) .eq. 3) then
              u%x(i,j,:,e) = avg_u%x(i,j,1,e)
           end if
        end do
     end do
  end do
end subroutine perform_global_summation

subroutine perform_local_summation(u_out, u, el_heights,domain_height, hom_dir_el, coef, nelv,lx)
  use neko
  implicit none
  type(field_t), intent(inout) :: u, u_out, el_heights
  type(coef_t), intent(inout) :: coef
  integer, intent(in) :: nelv, lx
  integer, intent(in) :: hom_dir_el(nelv)
  real(kind=rp), intent(in) :: domain_height
  real(kind=rp) :: wt
  integer :: n, i, j, e, k

  n = nelv*lx**3

  call col2(u%x,el_heights%x,n)
  call cmult(u%x, 1.0_rp/(2.0*domain_height),n)
  call rzero(u_out%x,n)


  do e = 1, nelv
     do i = 1,lx
        do j = 1, lx
           do k = 1, lx
              wt = coef%Xh%wx(k)
              if (hom_dir_el(e) .eq. 1) then
                 u_out%x(1,i,j,e) = u_out%x(1,i,j,e)+wt*u%x(k,i,j,e)
              else if (hom_dir_el(e) .eq. 2) then
                 u_out%x(i,1,j,e) = u_out%x(i,1,j,e)+wt*u%x(i,k,j,e)
              else if (hom_dir_el(e) .eq. 3) then
                 u_out%x(i,j,1,e) = u_out%x(i,j,1,e)+wt*u%x(i,j,k,e)
              end if
           end do
        end do
     end do
  end do

  do e = 1, nelv
     do i = 1,lx
        do j = 1, lx
           if (hom_dir_el(e) .eq. 1) then
              u_out%x(:,i,j,e) = u_out%x(1,i,j,e)
           else if (hom_dir_el(e) .eq. 2) then
              u_out%x(i,:,j,e) = u_out%x(i,1,j,e)
           else if (hom_dir_el(e) .eq. 3) then
              u_out%x(i,j,:,e) = u_out%x(i,j,1,e)
           end if
        end do
     end do
  end do
end subroutine perform_local_summation
