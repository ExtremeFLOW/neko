!> Somewhat preliminary program to sum up averaged fields computed for statistics and mean field
!! Martin Karp 27/01-23
program postprocess_fluid_stats
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, mesh_fname, stats_fname
  type(file_t) :: stats_file, output_file, mesh_file
  real(kind=rp) :: start_time
  type(fld_file_data_t) :: stats_data
  type(fluid_stats_t) :: fld_stats
  type(coef_t) :: coef
  type(dofmap_t) :: dof
  type(space_t) :: Xh
  type(mesh_t) :: msh
  type(gs_t) :: gs_h
  type(field_t), pointer :: u, v, w, p
  type(field_t), target :: pp, uu, vv, ww, uv, uw, vw, tmp1, tmp2
  type(field_list_t) :: reynolds, mean_vel_grad
  integer :: argc, i, n, lx

  argc = command_argument_count()

  if ((argc .lt. 3) .or. (argc .gt. 3)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./postprocess_fluid_stats mesh.nmsh stats.fld'
        write(*,*) 'Example command: ./postprocess_fluid_stats mesh.nmsh statsblabla.fld'
        write(*,*) 'Computes the statstics from the (3D) fld files described in statsblabla.nek5000'
        write(*,*) 'Recommended to use PyNekTools for more advanced postprocessing.'
        write(*,*) 'If this does not work, switch to that.'
        write(*,*) 'Currently we output two new fld files reynolds and mean_vel_grad'
        write(*,*) 'In Reynolds the fields are ordered as:'
        write(*,*) 'x-velocity=<u`u`>'
        write(*,*) 'y-velocity=<v`v`>'
        write(*,*) 'z-velocity=<w`w`>'
        write(*,*) 'pressure=<p`p`>'
        write(*,*) 'temperature=<u`v`>'
        write(*,*) 's1=<u`w`>'
        write(*,*) 's2=<v`w`>'
        write(*,*) 'In mean_vel_grad:'
        write(*,*) 'x-velocity=dudx'
        write(*,*) 'y-velocity=dudy'
        write(*,*) 'z-velocity=dudz'
        write(*,*) 'pressure=dvdx'
        write(*,*) 'temperature=dvdy'
        write(*,*) 's1=dvdz'
        write(*,*) 's2=dwdx'
        write(*,*) 's2=dwdy'
        write(*,*) 's2=dwdz'
     end if
     stop
  end if

  call neko_init

  call get_command_argument(1, inputchar)
  read(inputchar, *) mesh_fname
  mesh_file = file_t(trim(mesh_fname))
  call get_command_argument(2, inputchar)
  read(inputchar, *) stats_fname
  stats_file = file_t(trim(stats_fname))

  call mesh_file%read(msh)

  call stats_data%init(msh%nelv,msh%offset_el)
  call stats_file%read(stats_data)

  do i = 1,msh%nelv
     lx = stats_data%lx
     msh%elements(i)%e%pts(1)%p%x(1) = stats_data%x%x(linear_index(1,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(2)%p%x(1) = stats_data%x%x(linear_index(lx,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(3)%p%x(1) = stats_data%x%x(linear_index(1,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(4)%p%x(1) = stats_data%x%x(linear_index(lx,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(5)%p%x(1) = stats_data%x%x(linear_index(1,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(6)%p%x(1) = stats_data%x%x(linear_index(lx,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(7)%p%x(1) = stats_data%x%x(linear_index(1,lx,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(8)%p%x(1) = stats_data%x%x(linear_index(lx,lx,lx,i,lx,lx,lx))

     msh%elements(i)%e%pts(1)%p%x(2) = stats_data%y%x(linear_index(1,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(2)%p%x(2) = stats_data%y%x(linear_index(lx,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(3)%p%x(2) = stats_data%y%x(linear_index(1,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(4)%p%x(2) = stats_data%y%x(linear_index(lx,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(5)%p%x(2) = stats_data%y%x(linear_index(1,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(6)%p%x(2) = stats_data%y%x(linear_index(lx,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(7)%p%x(2) = stats_data%y%x(linear_index(1,lx,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(8)%p%x(2) = stats_data%y%x(linear_index(lx,lx,lx,i,lx,lx,lx))

     msh%elements(i)%e%pts(1)%p%x(3) = stats_data%z%x(linear_index(1,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(2)%p%x(3) = stats_data%z%x(linear_index(lx,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(3)%p%x(3) = stats_data%z%x(linear_index(1,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(4)%p%x(3) = stats_data%z%x(linear_index(lx,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(5)%p%x(3) = stats_data%z%x(linear_index(1,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(6)%p%x(3) = stats_data%z%x(linear_index(lx,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(7)%p%x(3) = stats_data%z%x(linear_index(1,lx,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(8)%p%x(3) = stats_data%z%x(linear_index(lx,lx,lx,i,lx,lx,lx))
  end do

  call Xh%init(GLL, stats_data%lx, stats_data%ly, stats_data%lz)

  call dof%init(msh, Xh)
  call gs_h%init(dof)
  call coef%init(gs_h)

  call neko_field_registry%add_field(dof, 'u')
  call neko_field_registry%add_field(dof, 'v')
  call neko_field_registry%add_field(dof, 'w')
  call neko_field_registry%add_field(dof, 'p')

  u => neko_field_registry%get_field('u')
  v => neko_field_registry%get_field('v')
  w => neko_field_registry%get_field('w')
  p => neko_field_registry%get_field('p')

  call fld_stats%init(coef,u,v,w,p)
  n = stats_data%u%n

  call copy(fld_stats%stat_fields%items(1)%ptr%x, stats_data%p%x,n)
  call copy(fld_stats%stat_fields%items(2)%ptr%x, stats_data%u%x,n)
  call copy(fld_stats%stat_fields%items(3)%ptr%x, stats_data%v%x,n)
  call copy(fld_stats%stat_fields%items(4)%ptr%x, stats_data%w%x,n)
  call copy(fld_stats%stat_fields%items(5)%ptr%x, stats_data%t%x,n)
  do i = 6, fld_stats%stat_fields%size()
     call copy(fld_stats%stat_fields%items(i)%ptr%x, stats_data%s(i-5)%x,n)
  end do

  call reynolds%init(7)
  !Temp fields used for the computations to come
  call uu%init(dof)
  call vv%init(dof)
  call ww%init(dof)
  call uv%init(dof)
  call uw%init(dof)
  call vw%init(dof)
  call pp%init(dof)
  call tmp1%init(dof)
  call tmp2%init(dof)

  call reynolds%assign_to_field(1, pp)
  call reynolds%assign_to_field(2, uu)
  call reynolds%assign_to_field(3, vv)
  call reynolds%assign_to_field(4, ww)
  call reynolds%assign_to_field(5, uv)
  call reynolds%assign_to_field(6, uw)
  call reynolds%assign_to_field(7, vw)

  call fld_stats%post_process(reynolds=reynolds)
  output_file = file_t('reynolds.fld')
  if (pe_rank .eq. 0) write(*,*) 'Wrtiting Reynolds stresses into reynolds'
  call output_file%write(reynolds, stats_data%time)
  if (pe_rank .eq. 0) write(*,*) 'Done'

  call mean_vel_grad%init(9)
  call mean_vel_grad%assign_to_field(1, pp)
  call mean_vel_grad%assign_to_field(2, uu)
  call mean_vel_grad%assign_to_field(3, vv)
  call mean_vel_grad%assign_to_field(4, ww)
  call mean_vel_grad%assign_to_field(5, uv)
  call mean_vel_grad%assign_to_field(6, uw)
  call mean_vel_grad%assign_to_field(7, vw)
  call mean_vel_grad%assign_to_field(8, tmp1)
  call mean_vel_grad%assign_to_field(9, tmp2)

  call fld_stats%post_process(mean_vel_grad=mean_vel_grad)
  !Fix order of gradients
  call mean_vel_grad%assign_to_field(1, pp)
  call mean_vel_grad%assign_to_field(2, uu)
  call mean_vel_grad%assign_to_field(3, vv)
  call mean_vel_grad%assign_to_field(4, ww)
  call mean_vel_grad%assign_to_field(5, uv)
  call mean_vel_grad%assign_to_field(6, uw)
  call mean_vel_grad%assign_to_field(7, vw)
  call mean_vel_grad%assign_to_field(8, tmp1)
  call mean_vel_grad%assign_to_field(9, tmp2)


  if (pe_rank .eq. 0) write(*,*) 'Writing mean velocity gradient into mean_vel_grad'
  output_file = file_t('mean_vel_grad.fld')
  call output_file%write(mean_vel_grad, stats_data%time)
  if (pe_rank .eq. 0) write(*,*) 'Done'

  call neko_finalize

end program postprocess_fluid_stats
