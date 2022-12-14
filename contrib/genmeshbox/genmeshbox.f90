!> Simple program to generate a mesh for a box
!! Gives some insight into how one can use built in functions in Neko
!! Martin Karp 26/10-22
program genmeshbox
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, fname='box.nmsh'
  real(kind=dp) :: x0, x1, el_len_x
  real(kind=dp) :: y0, y1, el_len_y
  real(kind=dp) :: z0, z1, el_len_z
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
  
  argc = command_argument_count()

  if ((argc .lt. 9) .or. (argc .gt. 12)) then
     write(*,*) 'Usage: ./genmeshbox x0 x1 y0 y1 z0 z1 nel_x nel_y nel_z'
     write(*,*) '**optional** periodic in x? (.true./.false.) y? (.true./.false.) z? (.true./.false.)'
     write(*,*)
     write(*,*) 'Example command: ./genmeshbox 0 1 0 1 0 1 8 8 8 .true. .true. .false.'
write(*,*) 'This examples generates a cube (box.nmsh) with side length 1 and with',&
           ' 8 elements in each spatial direction and periodic  boundaries in x-y.'
       write(*,*) 'BCs for face 5,6 (z zones) can then be set by setting bc_labels(5), bc_labels(6) in the parameter file'
     stop
  end if
  
  call neko_init 
  msh%lgenc = .false.
  period_x = .false.
  period_y = .false.
  period_z = .false.

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

  nel = nelx*nely*nelz
  call mesh_init(msh, gdim, nel)
  call htable_pts%init( nel, gdim)

  el_len_x = (x1 - x0)/nelx
  el_len_y = (y1 - y0)/nely
  el_len_z = (z1 - z0)/nelz
  i = 1
  do e_z = 0, (nelz-1)
     do e_y = 0, (nely-1)
        do e_x = 0, (nelx-1)
           do iz = 0,1
              do iy = 0,1
                 do ix = 0,1
                    coord(1) = x0 + (ix+e_x)*el_len_x
                    coord(2) = y0 + (iy+e_y)*el_len_y
                    coord(3) = z0 + (iz+e_z)*el_len_z
                    pt_idx = 1 + (ix+e_x) + (iy+e_y)*(nelx+1) + (iz+e_z)*(nelx+1)*(nely+1)
                    p(ix+1, iy+1, iz+1) = point_t(coord, pt_idx)
                 end do
              end do
           end do
           call mesh_add_element(msh, i, p(1,1,1), p(2,1,1),p(1,2,1),p(2,2,1),&
                            p(1,1,2), p(2,1,2),p(1,2,2),p(2,2,2))
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
              call msh%elements(e_id)%e%facet_order(facet_ord,zone_id)
              call mesh_mark_periodic_facet(msh, zone_id, el_idx, &
              1+mod(zone_id,2), p_el_idx, facet_ord%x)
           else
              call mesh_mark_labeled_facet(msh, zone_id, el_idx,zone_id)
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
              call mesh_mark_periodic_facet(msh, zone_id, el_idx, &
              3+mod(zone_id,2), p_el_idx, facet_ord%x)
           else
              call mesh_mark_labeled_facet(msh, zone_id, el_idx,zone_id)
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
              call mesh_mark_periodic_facet(msh, zone_id, el_idx, &
              5+mod(zone_id,2), p_el_idx, facet_ord%x)
           else
              call mesh_mark_labeled_facet(msh, zone_id, el_idx,zone_id)
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
                 el_idx = 1+e_y*nelx+e_z*nely*nelx
                 p_el_idx = 1+(nelx-1)+e_y*nelx+e_z*nely*nelx
              else 
                 p_el_idx = 1+e_y*nelx+e_z*nely*nelx
                 el_idx = 1+(nelx-1)+e_y*nelx+e_z*nely*nelx
              end if 
                 call mesh_create_periodic_ids(msh, zone_id, el_idx, &
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
              call mesh_create_periodic_ids(msh, zone_id, el_idx, &
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
              call mesh_create_periodic_ids(msh, zone_id, el_idx, &
              5+mod(zone_id,2), p_el_idx)
           end do
        end do
     end do
  end if

  call mesh_finalize(msh) 
  
  
  nmsh_file = file_t(fname)
  
  call nmsh_file%write(msh)
  
  call mesh_free(msh)
  
  call neko_finalize

end program genmeshbox
