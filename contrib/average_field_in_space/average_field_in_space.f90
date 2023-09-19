!> Program to sum up a field in space in a homogenous direction
!! Martin Karp 17/02-23
program average_field_in_space
  use neko
  use mean_flow
  implicit none
  
  character(len=NEKO_FNAME_LEN) :: inputchar, mesh_fname, field_fname, hom_dir, output_fname
  type(file_t) :: field_file, output_file, mesh_file
  real(kind=rp) :: start_time, el_h, el_dim(3,3), domain_height
  real(kind=rp), allocatable :: temp_el(:,:,:)
  type(fld_file_data_t) :: field_data
  type(coef_t) :: coef
  type(dofmap_t) :: dof
  type(space_t) :: Xh
  type(mesh_t) :: msh
  type(gs_t) :: gs_h
  type(field_t), pointer :: u, avg_u, old_u, el_heights
  type(vector_ptr_t), allocatable :: fields(:)
  integer, allocatable :: hom_dir_el(:)
  integer :: argc, i, n, lx, j, e, n_levels
  
  argc = command_argument_count()

  if ((argc .lt. 5) .or. (argc .gt. 5)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./average_field_in_space mesh.nmsh field.fld dir(x, y, z) n_levels  outfield.fld' 
        write(*,*) 'Example command: ./average_field_in_space mesh.nmsh fieldblabla.fld x n_levels outfield.fld'
        write(*,*) 'Computes the spatial average in x of the field fieldblabla.nek5000 and stores in outfield.fld'
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
  read(inputchar, *) n_levels
  call get_command_argument(5, inputchar) 
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

  call space_init(Xh, GLL, field_data%lx, field_data%ly, field_data%lz)

  dof = dofmap_t(msh, Xh)
  call gs_init(gs_h, dof)
  call coef%init(gs_h)

  call neko_field_registry%add_field(dof, 'u')
  u => neko_field_registry%get_field('u')
  call neko_field_registry%add_field(dof, 'avg_u')
  avg_u => neko_field_registry%get_field('avg_u')
  call neko_field_registry%add_field(dof, 'old_u')
  old_u => neko_field_registry%get_field('old_u')
  call neko_field_registry%add_field(dof, 'el_heights')
  el_heights => neko_field_registry%get_field('el_heights')

  !test
  allocate(hom_dir_el(msh%nelv))

  do i = 1, msh%nelv
     !find height in hom-dir
     !store direction which is hom
     !set element to height
     !we assume elements are stacked on eachother...
     el_dim(1,:) = abs(msh%elements(i)%e%pts(1)%p%x-msh%elements(i)%e%pts(2)%p%x)
     el_dim(2,:) = abs(msh%elements(i)%e%pts(1)%p%x-msh%elements(i)%e%pts(3)%p%x)
     el_dim(3,:) = abs(msh%elements(i)%e%pts(1)%p%x-msh%elements(i)%e%pts(5)%p%x)
     ! 1 corresponds to r, 2 to s, 3 to t
     if (trim(hom_dir) .eq. 'x') then
        hom_dir_el(i) = maxloc(el_dim(:,1),dim=1)
        el_h = el_dim(1,hom_dir_el(i))
     else if (trim(hom_dir) .eq. 'y') then
        hom_dir_el(i) = maxloc(el_dim(:,2),dim=1)
        el_h = el_dim(2,hom_dir_el(i))
     else if (trim(hom_dir) .eq. 'z') then
        hom_dir_el(i) = maxloc(el_dim(:,3),dim=1)
        el_h = el_dim(3,hom_dir_el(i))
     else 
        call neko_error('homogenous direction not supported')
     end if
     el_heights%x(:,:,:,i) = el_h
  end do
  n = u%dof%size()
   
  call copy(u%x,el_heights%x,n)
  call copy(old_u%x,el_heights%x,n)
  call copy(avg_u%x,el_heights%x,n)
  call perform_global_summation(u, avg_u, old_u, n_levels, &
       hom_dir_el,gs_h, coef%mult, msh%nelv, lx)
  domain_height = u%x(1,1,1,1)

  allocate(fields(field_data%size()))
  
  call field_data%get_list(fields,field_data%size())

  do i = 1, field_data%size()
     call copy(old_u%x,fields(i)%v%x,n)
     call perform_local_summation(u,old_u, el_heights, domain_height, &
          hom_dir_el, coef, msh%nelv, lx)
     call copy(old_u%x,u%x,n)
     call copy(avg_u%x,u%x,n)
     call perform_global_summation(u, avg_u, old_u, n_levels, &
          hom_dir_el,gs_h, coef%mult, msh%nelv, lx)
     call copy(fields(i)%v%x,u%x,n)
  end do 
  output_file = file_t(trim(output_fname))
  call output_file%write(field_data)
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
        call device_memcpy(u%x, u%x_d, n, HOST_TO_DEVICE)
     call gs_op(gs_h,u,GS_OP_ADD)
     if (NEKO_BCKND_DEVICE .eq. 1) &
        call device_memcpy(u%x, u%x_d, n, DEVICE_TO_HOST)
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
