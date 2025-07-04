program prepart
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: fname, nprtschr, output_
  type(mesh_t) :: msh, new_msh
  type(file_t) :: new_msh_file, nmsh_file
  integer :: argc, nprts, i, j, tmp, sum, idx, p_idx, rank, label
  integer, allocatable :: new_el(:), idx_cntr(:), idx_map(:)
  type(mesh_fld_t) :: parts

  argc = command_argument_count()

  if ((argc .lt. 2) .or. (argc .gt. 3)) then
     write(*,*) 'Usage: ./prepart <neko mesh> <nparts>'
     stop
  end if

  call neko_init()

  call get_command_argument(1, fname)
  call get_command_argument(2, nprtschr)
  read(nprtschr, *) nprts

  call nmsh_file%init(fname)
  call nmsh_file%read(msh)

  ! Reset possible periodic ids
  call msh%reset_periodic_ids()
  call neko_log%message('Splitting mesh')
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

  call neko_log%message('Generating new mesh')
  !
  ! Create redistributed mesh
  !

  new_msh%lgenc = .false.
  call new_msh%init(msh%gdim, msh%nelv)
  do i = 1, msh%nelv
     rank = parts%data(i)
     idx = idx_cntr(rank) + new_el(rank)
     idx_cntr(rank) = idx_cntr(rank) + 1
     idx_map(i) = idx
     call new_msh%add_element(idx, idx, &
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

  do i = 1, msh%periodic%size
     idx = idx_map(msh%periodic%facet_el(i)%x(2))
     p_idx = idx_map(msh%periodic%p_facet_el(i)%x(2))
     call new_msh%mark_periodic_facet(msh%periodic%facet_el(i)%x(1), idx, &
          msh%periodic%p_facet_el(i)%x(1), p_idx, msh%periodic%p_ids(i)%x)
  end do

  do j = 1, NEKO_MSH_MAX_ZLBLS
     do i = 1, msh%labeled_zones(j)%size
        idx = idx_map(msh%labeled_zones(j)%facet_el(i)%x(2))
        label = j ! adhere to standards...
        call new_msh%mark_labeled_facet(msh%labeled_zones(j)%facet_el(i)%x(1), &
             idx, label)
     end do
  end do

  do i = 1, msh%periodic%size
     idx = idx_map(msh%periodic%facet_el(i)%x(2))
     p_idx = idx_map(msh%periodic%p_facet_el(i)%x(2))
     call new_msh%apply_periodic_facet(msh%periodic%facet_el(i)%x(1), idx, &
          msh%periodic%p_facet_el(i)%x(1), p_idx, msh%periodic%p_ids(i)%x)
  end do

  do i = 1, msh%curve%size
     idx = idx_map(msh%curve%curve_el(i)%el_idx)
     call new_msh%mark_curve_element(idx, msh%curve%curve_el(i)%curve_data, &
          msh%curve%curve_el(i)%curve_type)
  end do

  call new_msh%finalize()

  deallocate(idx_map)
  call msh%free()

  output_ = trim(fname(1:scan(trim(fname), '.', back = .true.) - 1)) // &
       '_' // trim(nprtschr) // '.nmsh'

  call new_msh_file%init(output_)
  call new_msh_file%write(new_msh)
  call new_msh%free()

  call neko_finalize

end program prepart
