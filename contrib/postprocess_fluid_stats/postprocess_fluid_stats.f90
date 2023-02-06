!> Program to sum up averaged fields computed for statistics and mean field
!! Martin Karp 27/01-23
program postprocess_fluid_stats
  use neko
  use mean_flow
  implicit none
  
  character(len=NEKO_FNAME_LEN) :: inputchar, mesh_fname, stats_fname, mean_fname
  type(file_t) :: mean_file, stats_file, output_file, mesh_file
  real(kind=rp) :: start_time
  type(fld_file_data_t) :: stats_data, mean_data
  type(mean_flow_t) :: avg_flow
  type(fluid_stats_t) :: fld_stats
  type(coef_t) :: coef
  type(dofmap_t) :: dof
  type(space_t) :: Xh
  type(mesh_t) :: msh
  type(gs_t) :: gs_h
  type(field_t), pointer :: u, v, w, p
  type(field_t), target :: pp, uu, vv, ww, uv, uw, vw, tmp1, tmp2
  type(field_list_t) :: reynolds, mean_vel_grad
  integer :: argc, i, n
  
  argc = command_argument_count()

  if ((argc .lt. 3) .or. (argc .gt. 3)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./postprocess_fluid_stats mesh.nmsh mean_field.fld stats.fld' 
        write(*,*) 'Example command: ./average_fields mesh.nmsh mean_fieldblabla.fld statsblabla.fld'
        write(*,*) 'Computes the statstics from the fld files described in mean_field104.nek5000 stats104.nek5000'
        write(*,*) 'Currently we output two new fld files reynolds and mean_vei_grad'
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
  read(inputchar, *) mean_fname
  mean_file = file_t(trim(mean_fname))
  call get_command_argument(3, inputchar) 
  read(inputchar, *) stats_fname
  stats_file = file_t(trim(stats_fname))
  
  call mesh_file%read(msh)
   
  call mean_data%init(msh%nelv,msh%offset_el)
  call stats_data%init(msh%nelv,msh%offset_el)
  call mean_file%read(mean_data)
  call stats_file%read(stats_data)

  call space_init(Xh, GLL, mean_data%lx, mean_data%ly, mean_data%lz)
  dof = dofmap_t(msh, Xh)
  call gs_init(gs_h, dof)
  call coef_init(coef, gs_h)

  call neko_field_registry%add_field(dof, 'u')
  call neko_field_registry%add_field(dof, 'v')
  call neko_field_registry%add_field(dof, 'w')
  call neko_field_registry%add_field(dof, 'p')

  u => neko_field_registry%get_field('u')
  v => neko_field_registry%get_field('v')
  w => neko_field_registry%get_field('w')
  p => neko_field_registry%get_field('p')

  call avg_flow%init(u, v, w, p)
  call fld_stats%init(coef)
  n = mean_data%u%n
  call copy(avg_flow%u%mf%x,mean_data%u%x,n)
  call copy(avg_flow%v%mf%x,mean_data%v%x,n)
  call copy(avg_flow%w%mf%x,mean_data%w%x,n)
  call copy(avg_flow%p%mf%x,mean_data%p%x,n)
  
  call copy(fld_stats%stat_fields%fields(1)%f%x,stats_data%p%x,n)
  call copy(fld_stats%stat_fields%fields(2)%f%x,stats_data%u%x,n)
  call copy(fld_stats%stat_fields%fields(3)%f%x,stats_data%v%x,n)
  call copy(fld_stats%stat_fields%fields(4)%f%x,stats_data%w%x,n)
  call copy(fld_stats%stat_fields%fields(5)%f%x,stats_data%t%x,n)
  do i = 6, size(fld_stats%stat_fields%fields)
     call copy(fld_stats%stat_fields%fields(i)%f%x,stats_data%s(i-5)%x,n)
  end do

  allocate(reynolds%fields(7))

  call field_init(uu,dof)
  call field_init(vv,dof)
  call field_init(ww,dof)
  call field_init(uv,dof)
  call field_init(uw,dof)
  call field_init(vw,dof)
  call field_init(pp,dof)
  call field_init(tmp1,dof)
  call field_init(tmp2,dof)

  reynolds%fields(1)%f => pp
  reynolds%fields(2)%f => uu
  reynolds%fields(3)%f => vv
  reynolds%fields(4)%f => ww
  reynolds%fields(5)%f => uv
  reynolds%fields(6)%f => uw
  reynolds%fields(7)%f => vw

  call fld_stats%post_process(reynolds=reynolds)
  output_file = file_t('reynolds.fld')
  if (pe_rank .eq. 0) write(*,*) 'Wrtiting Reynolds stresses into reynolds'
  call output_file%write(reynolds, stats_data%time)
  if (pe_rank .eq. 0) write(*,*) 'Done'

  allocate(mean_vel_grad%fields(9))
  mean_vel_grad%fields(1)%f => pp
  mean_vel_grad%fields(2)%f => uu
  mean_vel_grad%fields(3)%f => vv
  mean_vel_grad%fields(4)%f => ww
  mean_vel_grad%fields(5)%f => uv
  mean_vel_grad%fields(6)%f => uw
  mean_vel_grad%fields(7)%f => vw
  mean_vel_grad%fields(8)%f => tmp1
  mean_vel_grad%fields(9)%f => tmp2

  call fld_stats%post_process(mean_vel_grad=mean_vel_grad)
  !Fix order of gradients
  mean_vel_grad%fields(2)%f => pp
  mean_vel_grad%fields(3)%f => uu
  mean_vel_grad%fields(4)%f => vv
  mean_vel_grad%fields(1)%f => ww
  mean_vel_grad%fields(5)%f => uv
  mean_vel_grad%fields(6)%f => uw
  mean_vel_grad%fields(7)%f => vw
  mean_vel_grad%fields(8)%f => tmp1
  mean_vel_grad%fields(9)%f => tmp2


  if (pe_rank .eq. 0) write(*,*) 'Writing mean velocity gradient into mean_vel_grad'
  output_file = file_t('mean_vel_grad.fld')
  call output_file%write(mean_vel_grad, stats_data%time)
  if (pe_rank .eq. 0) write(*,*) 'Done'
  
  call neko_finalize

end program postprocess_fluid_stats
