!> Program to calculate the force and acting on a single boundary zone as well as
!! the torque around a point if one changes the value of center (defaults to 0,0,0).
!! Outputs the x,y,z prjections of the pressure and viscous forces and torques and
!! additionally saves the distribution of these quantities along a selected homogenous
!! direction to a csv file.
!! Martin Karp 17/01-24
program calc_lift_from_field
  use neko
  use mean_flow
  use matrix
  implicit none
  
  character(len=NEKO_FNAME_LEN) :: inputchar, mesh_fname, field_fname, hom_dir, output_fname
  type(file_t) :: field_file, mesh_file, output_file
  real(kind=rp) :: start_time
  type(fld_file_data_t) :: field_data
  type(coef_t) :: coef
  type(dofmap_t), target :: dof
  type(space_t) :: Xh
  type(mesh_t) :: msh
  type(gs_t) :: gs_h
  type(map_1d_t) :: map_1d
  type(field_t) :: u, v, w, p
  type(field_t) :: s11, s22, s33, s12, s13, s23
  real(kind=rp) :: s11_, s22_, s33_, s12_, s13_, s23_, center(3)
  type(matrix_t) :: drag_torq
  real(kind=rp), allocatable :: lvl_coords(:)
  real(kind=rp), pointer :: line(:,:,:,:)
  integer :: argc, i, j, k, n, lx, e, zone_id, dir, f, glb_n_gll_pts, ierr, mem, t
  real(kind=rp) :: visc, nv(3), dgtq(12)
  
  argc = command_argument_count()

  if ((argc .lt. 6) .or. (argc .gt. 6)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./calc_lift_from_field mesh.nmsh field.fld zone_number viscosity function_of_coord output.csv' 
        write(*,*) 'Example command: ./calc_lift_from_field mesh.nmsh fieldblabla.fld 5 0.04 y out.csv'
        write(*,*) 'Outputs the total force and torque on zone 5 using velocity values from fieldblabla.fld'
        write(*,*)  'as well as writes the distribution of the force and torque across y to output.csv'
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
  read(inputchar, *) zone_id
  call get_command_argument(4, inputchar) 
  read(inputchar, *) visc
  call get_command_argument(5, inputchar) 
  read(inputchar, *) hom_dir
  call get_command_argument(6, inputchar) 
  read(inputchar, *) output_fname
  output_file = file_t(trim(output_fname))

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
  ! Center around which we calculate the torque
  center = 0.0_rp

  ! 1 corresponds to r, 2 to s, 3 to t
  if (trim(hom_dir) .eq. 'x') then
     dir = 1
     line => dof%x
  else if (trim(hom_dir) .eq. 'y') then
     dir = 2
     line => dof%y
  else if (trim(hom_dir) .eq. 'z') then
     dir = 3
     line => dof%z
  else 
     call neko_error('The homogeneous direction should be "x", "y"or "z"')
  end if
  call map_1d%init(coef, dir, 1e-7_rp)


  call u%init(dof)
  call v%init(dof)
  call w%init(dof)
  call p%init(dof)
  call s11%init(dof)
  call s22%init(dof)
  call s33%init(dof)
  call s12%init(dof)
  call s13%init(dof)
  call s23%init(dof)
  glb_n_gll_pts = map_1d%n_el_lvls*lx
  call drag_torq%init(glb_n_gll_pts, 14)
  if (pe_rank .eq. 0) call output_file%set_header('time, coord, forcepx, forcepy, &
             & forcepz, forcevx, forcevy, forcevz, torqpx, torqpy, &
             & torqpz, torqvx, torqvy, torqvz')

  do t = 1, field_data%meta_nsamples
     if (t .ne. 1) call field_file%read(field_data)
     if(pe_rank .eq. 0) write(*,*) 'Calc drag/lift at t= ', field_data%time

     call copy(u%x,field_data%u%x,dof%size())
     call copy(v%x,field_data%v%x,dof%size())
     call copy(w%x,field_data%w%x,dof%size())
     call copy(p%x,field_data%p%x,dof%size())
     drag_torq = 0.0_rp
     !set coords to somoething big
     drag_torq%x(:,2) = 1e15
   
   
     if(pe_rank .eq. 0) write(*,*) 'Total drag'
   
     call strain_rate(s11%x, s22%x, s33%x, s12%x, s13%x, s23%x, u, v, w, coef) 
     call drag_torque_zone(dgtq,field_data%t_counter, msh%labeled_zones(zone_id), center,&
                           s11%x, s22%x, s33%x, s12%x, s13%x, s23%x,&
                           p, coef, visc)
     if (pe_rank .eq. 0) then
         write(*,*) field_data%t_counter,dgtq(1)+dgtq(4),dgtq(1),dgtq(4),'dragx'
         write(*,*) field_data%t_counter,dgtq(2)+dgtq(5),dgtq(2),dgtq(5),'dragy'
         write(*,*) field_data%t_counter,dgtq(3)+dgtq(6),dgtq(3),dgtq(6),'dragz'
      end if


     do mem  = 1,msh%labeled_zones(zone_id)%size
        e = msh%labeled_zones(zone_id)%facet_el(mem)%x(2)
        f = msh%labeled_zones(zone_id)%facet_el(mem)%x(1)
        do k = 1, lx
           do j = 1, lx
              do i = 1, lx
                 if (index_is_on_facet(i,j,k,lx,lx,lx,f)) then
                    nv = coef%get_normal(i,j,k,e,f)*coef%get_area(i,j,k,e,f)
                    s11_ = s11%x(i,j,k,e)
                    s12_ = s12%x(i,j,k,e)
                    s22_ = s22%x(i,j,k,e)
                    s13_ = s13%x(i,j,k,e)
                    s23_ = s23%x(i,j,k,e)
                    s33_ = s33%x(i,j,k,e)
                    call drag_torque_pt(dgtq,dof%x(i,j,k,e), dof%y(i,j,k,e),dof%z(i,j,k,e),&
                                        center, &
                                        s11_, s22_, s33_, s12_, s13_, s23_,&
                                        p%x(i,j,k,e), nv(1), nv(2), nv(3), visc)
                    drag_torq%x(map_1d%pt_lvl(i,j,k,e),3:14) = &
                    drag_torq%x(map_1d%pt_lvl(i,j,k,e),3:14) + dgtq
                    drag_torq%x(map_1d%pt_lvl(i,j,k,e),2) =  line(i,j,k,e)
                 end if
              end do
           end do
        end do
     end do
   
     call MPI_Allreduce(MPI_IN_PLACE,drag_torq%x(1,3), 12*glb_n_gll_pts, &
        MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
     call MPI_Allreduce(MPI_IN_PLACE,drag_torq%x(1,2), glb_n_gll_pts, &
        MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
     drag_torq%x(:,1) = field_data%time
     if (pe_rank .eq. 0) then
        call output_file%write(drag_torq)
     end if
  end do
     
  call neko_finalize

end program calc_lift_from_field
