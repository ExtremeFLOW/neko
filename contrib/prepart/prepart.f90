program prepart
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: fname, nprtschr, output_
  type(mesh_t) :: msh, new_msh
  type(file_t) :: new_msh_file, nmsh_file
  integer :: argc, nprts, i, j, tmp, sum, idx, p_idx, rank
  integer, allocatable :: new_el(:), idx_cntr(:), idx_map(:)
  type(mesh_fld_t) :: parts

  argc = command_argument_count()

  if ((argc .lt. 2) .or. (argc .gt. 3)) then
     write(*,*) 'Usage: ./prepart <neko mesh> <nparts>'
     stop
  end if

  call neko_init 

  call get_command_argument(1, fname)
  call get_command_argument(2, nprtschr)
  read(nprtschr, *) nprts

  nmsh_file = file_t(fname)

  call nmsh_file%read(msh)

  ! Reset possible periodic ids
  call mesh_reset_periodic_ids(msh)

  ! Compute new partitions
  call parmetis_partmeshkway(msh, parts, nprts=nprts)


  !
  ! Compute offsets for new distribution
  !
  allocate(new_el(0:nprts - 1))
  new_el = 0
  do i = 1, msh%nelv
     new_el(parts%data(i)) = new_el(parts%data(i)) + 1
  end do

  sum = 0
  do i = 0, nprts - 1
     tmp = new_el(i)
     new_el(i) = sum
     sum = sum + tmp
  end do

  allocate(idx_cntr(0:nprts - 1))
  idx_cntr = 1

  allocate(idx_map(msh%nelv))

  !
  ! Create redistributed mesh
  !

  call mesh_init(new_msh, msh%gdim, msh%nelv)  
  do i = 1, msh%nelv
     rank = parts%data(i)     
     idx = idx_cntr(rank) + new_el(rank)
     idx_cntr(rank) = idx_cntr(rank) + 1
     idx_map(i) = idx
     call mesh_add_element(new_msh, idx, &
          msh%elements(i)%e%pts(1)%p, &
          msh%elements(i)%e%pts(2)%p, &
          msh%elements(i)%e%pts(3)%p, &
          msh%elements(i)%e%pts(4)%p, &
          msh%elements(i)%e%pts(5)%p, &
          msh%elements(i)%e%pts(6)%p, &
          msh%elements(i)%e%pts(7)%p, &
          msh%elements(i)%e%pts(8)%p)

  end do

  deallocate(new_el)
  deallocate(idx_cntr)

  !
  ! Add zones
  ! 
  do i = 1, msh%wall%size
     idx = idx_map(msh%wall%facet_el(i)%x(2))
     call mesh_mark_wall_facet(new_msh, msh%wall%facet_el(i)%x(1), idx)
  end do

  do i = 1, msh%inlet%size
     idx = idx_map(msh%inlet%facet_el(i)%x(2))
     call mesh_mark_inlet_facet(new_msh, msh%inlet%facet_el(i)%x(1), idx)
  end do

  do i = 1, msh%outlet%size
     idx = idx_map(msh%outlet%facet_el(i)%x(2))
     call mesh_mark_outlet_facet(new_msh, msh%outlet%facet_el(i)%x(1), idx)
  end do

  do i = 1, msh%sympln%size
     idx = idx_map(msh%sympln%facet_el(i)%x(2))
     call mesh_mark_sympln_facet(new_msh, msh%sympln%facet_el(i)%x(1), idx)
  end do

  do i =1, msh%periodic%size
     idx = idx_map(msh%periodic%facet_el(i)%x(2))
     p_idx = idx_map(msh%periodic%p_facet_el(i)%x(2))
     call mesh_apply_periodic_facet(new_msh, msh%periodic%facet_el(i)%x(1), idx, &
          msh%periodic%p_facet_el(i)%x(1), p_idx, msh%periodic%p_ids(i)%x)
  end do

  call mesh_finalize(new_msh)

  deallocate(idx_map)
  call mesh_free(msh)

  output_ = trim(fname(1:scan(trim(fname), &
       '.', back=.true.) - 1))//'_'//trim(nprtschr)//'.nmsh' 

  new_msh_file = file_t(output_)
  call new_msh_file%write(new_msh)
  call mesh_free(new_msh)

  call neko_finalize

end program prepart
