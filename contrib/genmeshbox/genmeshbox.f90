!> Simple program to generate a mesh for a box
!! Gives some insight into how one can use built in functions in Neko
!! Martin Karp 26/10-22
program genmeshbox
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, fname = 'box.nmsh'
  real(kind=dp) :: x0, x1
  real(kind=dp) :: y0, y1
  real(kind=dp) :: z0, z1
  real(kind=dp) :: coord(3)
  integer :: nelx, nely, nelz, nel
  type(mesh_t) :: msh
  type(file_t) :: nmsh_file
  type(htable_pt_t) :: htable_pts
  logical :: period_x, period_y, period_z
  type(point_t) :: p(2,2,2)
  type(tuple4_i4_t) :: facet_ord
  integer :: argc, gdim = 3
  integer :: el_idx, p_el_idx, pt_idx
  integer :: i, zone_id, e_x, e_y, e_z, ix, iy, iz, e_id
  real(kind=dp), allocatable :: el_len_x(:), el_len_y(:), el_len_z(:)
  real(kind=dp), allocatable :: cumm_x(:), cumm_y(:), cumm_z(:)
  type(vector_t) :: dist_x, dist_y, dist_z
  type(file_t) :: dist_x_file, dist_y_file, dist_z_file
  character(len = 80) :: dist_x_fname, dist_y_fname, dist_z_fname
  character(len=80) :: log_fname = "genmeshbox.log"
  logical :: file_exists

  argc = command_argument_count()

  if ((argc .lt. 9) .or. (argc .gt. 15)) then
     write(*,*) 'Usage: ./genmeshbox x0 x1 y0 y1 z0 z1 nel_x nel_y nel_z'
     write(*,*) '**optional** periodic in x? (.true./.false.) y? ', &
          '(.true./.false.) z? (.true./.false.)'
     write(*,*) '**optional** Point distribution in x? (fname/uniform) y? ', &
          '(fname/uniform) z? (fname/uniform)'
     write(*,*)
     write(*,*) 'Example command: ./genmeshbox 0 1 0 1 0 1 8 8 8 ', &
          '.true. .true. .false.'
     write(*,*) 'This examples generates a cube (box.nmsh) with side length ', &
          '1 and with 8 elements in each spatial direction and periodic ', &
          'boundaries in x-y.'
     write(*,*) 'BCs for face 5,6 (z zones) can then be set by setting ', &
          'bc_labels(5), bc_labels(6) in the parameter file'
     write(*,*) 'If you want a specific distribution of vertices in the ', &
          'directions, give the filename where it is stored'
     write(*,*) 'Example command: ./genmeshbox 0 1 0 1 0 1 8 8 8 ', &
          '.false. .false. .false. uniform uniform dist.csv'
     write(*,*) 'This example uses the vertex distribution in the z ', &
          'direction given in dist.csv and keep the other'
     write(*,*) 'directions uniform'
     stop
  end if

  call neko_init
  msh%lgenc = .false.
  period_x = .false.
  period_y = .false.
  period_z = .false.
  dist_x_fname = 'uniform'
  dist_y_fname = 'uniform'
  dist_z_fname = 'uniform'

  call get_command_argument(1, inputchar)
  read(inputchar, *) x0
  call get_command_argument(2, inputchar)
  read(inputchar, *) x1
  call get_command_argument(3, inputchar)
  read(inputchar, *) y0
  call get_command_argument(4, inputchar)
  read(inputchar, *) y1
  call get_command_argument(5, inputchar)
  read(inputchar, *) z0
  call get_command_argument(6, inputchar)
  read(inputchar, *) z1
  call get_command_argument(7, inputchar)
  read(inputchar, *) nelx
  call get_command_argument(8, inputchar)
  read(inputchar, *) nely
  call get_command_argument(9, inputchar)
  read(inputchar, *) nelz
  if (argc .gt. 9) then
     call get_command_argument(10, inputchar)
     read(inputchar, *) period_x
     call get_command_argument(11, inputchar)
     read(inputchar, *) period_y
     call get_command_argument(12, inputchar)
     read(inputchar, *) period_z
  end if
  if (argc .gt. 12) then
     call get_command_argument(13, inputchar)
     read(inputchar, *) dist_x_fname
     call get_command_argument(14, inputchar)
     read(inputchar, *) dist_y_fname
     call get_command_argument(15, inputchar)
     read(inputchar, *) dist_z_fname
  end if

  ! Write a log of what parameters we used with genmeshbox
  if (pe_rank .eq. 0) then

     do i = 1, 1000

        inquire(file = trim(log_fname), exist = file_exists)
        if (.not. file_exists) then

           open(unit=10, file=trim(log_fname), status = 'new', action = 'write')
           write (10, '(A,2(F10.6," "),I4,L2)') "xmin, xmax, Nel, periodic:", &
                x0, x1, nelx, period_x
           write (10, '(A,2(F10.6," "),I4,L2)') "ymin, ymax, Nel, periodic:", &
                y0, y1, nely, period_y
           write (10, '(A,2(F10.6," "),I4,L2)') "zmin, zmax, Nel, periodic:", &
                z0, z1, nelz, period_z
           close(10)
           exit

        end if

        ! if the original genmeshbox.log does not exist, we create new
        ! files with genmeshbox.log.1, .2, .3, etc
        if (i .lt. 10) then
           write(log_fname, '(A,I1,A)') "genmeshbox_", i, ".log"
        else if (i .lt. 100) then
           write(log_fname, '(A,I2,A)') "genmeshbox_", i, ".log"
        else if (i .lt. 1000) then
           write(log_fname, '(A,I3,A)') "genmeshbox_", i, ".log"
        else
           write(log_fname, '(A,I4,A)') "genmeshbox_", i, ".log"
        end if

     end do
  end if

  nel = nelx*nely*nelz
  call msh%init(gdim, nel)
  call htable_pts%init( nel, gdim)

  allocate(el_len_x(nelx), el_len_y(nely), el_len_z(nelz))
  allocate(cumm_x(nelx+1), cumm_y(nely+1), cumm_z(nelz+1))

  if (trim(dist_x_fname) .eq. 'uniform') then
     el_len_x(:) = (x1 - x0)/nelx
  else
     dist_x_file = file_t(trim(dist_x_fname))
     call dist_x%init(nelx+1)
     call dist_x_file%read(dist_x)
     do i = 1, nelx
        el_len_x(i) = dist_x%x(i+1) - dist_x%x(i)
     end do
     call file_free(dist_x_file)
     call dist_x%free()
  end if

  if (trim(dist_y_fname) .eq. 'uniform') then
     el_len_y(:) = (y1 - y0)/nely
  else
     dist_y_file = file_t(trim(dist_y_fname))
     call dist_y%init(nely+1)
     call dist_y_file%read(dist_y)
     do i = 1, nely
        el_len_y(i) = dist_y%x(i+1) - dist_y%x(i)
     end do
     call file_free(dist_y_file)
     call dist_y%free()
  end if

  if (trim(dist_z_fname) .eq. 'uniform') then
     el_len_z(:) = (z1 - z0)/nelz
  else
     dist_z_file = file_t(trim(dist_z_fname))
     call dist_z%init(nelz+1)
     call dist_z_file%read(dist_z)
     do i = 1, nelz
        el_len_z(i) = dist_z%x(i+1) - dist_z%x(i)
     end do
     call file_free(dist_z_file)
     call dist_z%free()
  end if

  cumm_x(1) = x0
  do i = 1, nelx
     cumm_x(i+1) = cumm_x(i) + el_len_x(i)
  end do

  cumm_y(1) = y0
  do i = 1, nely
     cumm_y(i+1) = cumm_y(i) + el_len_y(i)
  end do

  cumm_z(1) = z0
  do i = 1, nelz
     cumm_z(i+1) = cumm_z(i) + el_len_z(i)
  end do

  i = 1
  do e_z = 0, (nelz-1)
     do e_y = 0, (nely-1)
        do e_x = 0, (nelx-1)
           do iz = 0,1
              do iy = 0,1
                 do ix = 0,1
                    coord(1) = cumm_x(e_x + 1 + ix)
                    coord(2) = cumm_y(e_y + 1 + iy)
                    coord(3) = cumm_z(e_z + 1 + iz)
                    pt_idx = 1 + (ix + e_x) + (iy + e_y)*(nelx + 1) + &
                         (iz + e_z)*(nelx + 1)*(nely + 1)
                    p(ix+1, iy+1, iz+1) = point_t(coord, pt_idx)
                 end do
              end do
           end do
           call msh%add_element(i, p(1,1,1), p(2,1,1),p(1,2,1),p(2,2,1),&
                p(1,1,2), p(2,1,2), p(1,2,2), p(2,2,2))
           i = i + 1
        end do
     end do
  end do

  !Compute x zones
  do zone_id = 1,2
     do e_y = 0, (nely-1)
        do e_z = 0, (nelz-1)
           e_id = 1+e_y*nelx+e_z*nely*nelx
           if (zone_id .eq. 1) then
              el_idx = 1+e_y*nelx+e_z*nely*nelx
              p_el_idx = 1+(nelx-1)+e_y*nelx+e_z*nely*nelx
           else
              p_el_idx = 1+e_y*nelx+e_z*nely*nelx
              el_idx = 1+(nelx-1)+e_y*nelx+e_z*nely*nelx
           end if
           if (period_x) then
              call msh%elements(e_id)%e%facet_order(facet_ord, zone_id)
              call msh%mark_periodic_facet(zone_id, el_idx, &
                   1+mod(zone_id,2), p_el_idx, facet_ord%x)
           else
              call msh%mark_labeled_facet(zone_id, el_idx, zone_id)
           end if
        end do
     end do
  end do

  !Compute y zones
  do zone_id = 3,4
     do e_x = 0, (nelx-1)
        do e_z = 0, (nelz-1)
           e_id = 1+e_x+e_z*nely*nelx
           if (zone_id .eq. 3) then
              el_idx = 1+e_x+e_z*nely*nelx
              p_el_idx = 1+e_x+(nely-1)*nelx+e_z*nely*nelx
           else
              el_idx = 1+e_x+(nely-1)*nelx+e_z*nely*nelx
              p_el_idx = 1+e_x+e_z*nely*nelx
           end if
           if (period_y) then
              call msh%elements(e_id)%e%facet_order(facet_ord,3)
              call msh%mark_periodic_facet(zone_id, el_idx, &
                   3+mod(zone_id,2), p_el_idx, facet_ord%x)
           else
              call msh%mark_labeled_facet(zone_id, el_idx, zone_id)
           end if
        end do
     end do
  end do

  !Compute z zones
  do zone_id = 5,6
     do e_x = 0, (nelx-1)
        do e_y = 0, (nely-1)
           e_id = 1+e_x+e_y*nelx
           if (zone_id .eq. 5) then
              el_idx = 1+e_x+e_y*nelx
              p_el_idx = 1+e_x+e_y*nelx+(nelz-1)*nely*nelx
           else
              p_el_idx = 1+e_x+e_y*nelx
              el_idx = 1+e_x+e_y*nelx+(nelz-1)*nely*nelx
           end if
           if (period_z) then
              call msh%elements(e_id)%e%facet_order(facet_ord,5)
              call msh%mark_periodic_facet(zone_id, el_idx, &
                   5+mod(zone_id,2), p_el_idx, facet_ord%x)
           else
              call msh%mark_labeled_facet(zone_id, el_idx, zone_id)
           end if
        end do
     end do
  end do
  ! No elegance just brute force...
  ! To setup the correct global ids for the periodic points
  if (period_x) then
     do zone_id = 1,2
        do e_y = 0, (nely-1)
           do e_z = 0, (nelz-1)
              e_id = 1+e_y*nelx+e_z*nely*nelx
              if (zone_id .eq. 1) then
                 el_idx = 1 + e_y*nelx + e_z*nely*nelx
                 p_el_idx = 1 + (nelx-1) + e_y*nelx + e_z*nely*nelx
              else
                 p_el_idx = 1 + e_y*nelx + e_z*nely*nelx
                 el_idx = 1 + (nelx-1) + e_y*nelx + e_z*nely*nelx
              end if
              call msh%create_periodic_ids(zone_id, el_idx, &
                   1+mod(zone_id,2), p_el_idx)
           end do
        end do
     end do
  end if

  if (period_y) then
     !Compute y zones
     do zone_id = 3,4
        do e_x = 0, (nelx-1)
           do e_z = 0, (nelz-1)
              e_id = 1+e_x+e_z*nely*nelx
              if (zone_id .eq. 3) then
                 el_idx = 1+e_x+e_z*nely*nelx
                 p_el_idx = 1+e_x+(nely-1)*nelx+e_z*nely*nelx
              else
                 el_idx = 1+e_x+(nely-1)*nelx+e_z*nely*nelx
                 p_el_idx = 1+e_x+e_z*nely*nelx
              end if
              call msh%create_periodic_ids(zone_id, el_idx, &
                   3+mod(zone_id,2), p_el_idx)
           end do
        end do
     end do
  end if

  if (period_z) then
     do zone_id = 5,6
        do e_x = 0, (nelx-1)
           do e_y = 0, (nely-1)
              e_id = 1+e_x+e_y*nelx
              if (zone_id .eq. 5) then
                 el_idx = 1+e_x+e_y*nelx
                 p_el_idx = 1+e_x+e_y*nelx+(nelz-1)*nely*nelx
              else
                 p_el_idx = 1+e_x+e_y*nelx
                 el_idx = 1+e_x+e_y*nelx+(nelz-1)*nely*nelx
              end if
              call msh%create_periodic_ids(zone_id, el_idx, &
                   5+mod(zone_id,2), p_el_idx)
           end do
        end do
     end do
  end if

  call msh%finalize()


  nmsh_file = file_t(fname)

  call nmsh_file%write(msh)

  call msh%free()

  call neko_finalize

end program genmeshbox
