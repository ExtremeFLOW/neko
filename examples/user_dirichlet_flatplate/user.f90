module user
  use neko
  implicit none

  ! ----- For the boundary condition stuff ------------------------------------
  type(tuple_i4_t), allocatable :: map_outlet(:), map_inlet(:) ! Lists of tuples 
  type(matrix_t) :: csv_in_o, csv_in_v
  ! ---------------------------------------------------------------------------

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules     => initialize ! Initialize parameters
    u%user_check            => check
    u%user_finalize_modules => finalize   ! Finalize
    u%user_dirichlet_update => user_bc    ! Update boundary conditions
  end subroutine user_setup

  subroutine initialize(t, u, v, w, p, coef, params)
    implicit none

    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    type(file_t) :: csv_file_in_o, csv_file_in_v
    integer :: lx, ierr

    lx = coef%Xh%lx

    !
    ! Initialize matrices to read inlet/outlet values
    ! the sizes of the arrays were obtained by counting the lines in
    ! each file using the command "wc -l data/o.csv".
    ! The columns of each file are respectively:
    ! x coordinate | y coordinate | pressure | u velocity | v velocity
    !
    call csv_in_o%init(560, 5)
    call csv_in_v%init(992, 5)

    csv_file_in_o = file_t("data/o.csv")
    csv_file_in_v = file_t("data/v.csv")

    call csv_file_in_o%read(csv_in_o)
    call csv_file_in_v%read(csv_in_v)

    call MPI_Bcast(csv_in_o%x , csv_in_o%n , MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(csv_in_v%x , csv_in_v%n , MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)

    call file_free(csv_file_in_o)
    call file_free(csv_file_in_v)

    ! 
    ! Allocate necessary arrays to map the csv data
    ! These mapping arrays are lists of tuples (i,j) where
    ! i is a linear index of the gll point on which to apply 
    ! the value j in the csv file (see the subroutine map_boundary below
    ! for more info)
    !
    allocate(map_inlet(u%msh%labeled_zones(1)%size*lx*lx))
    allocate(map_outlet(u%msh%labeled_zones(2)%size*lx*lx))

    !
    ! Map the csv values with the inlet boundary points
    !
    if (u%msh%labeled_zones(1)%size .ne. 0) then
      
       call neko_log%message("mapping inlet...")
       call map_boundary(u%dof%x, u%dof%y, u%msh%labeled_zones(1)%facet_el, & 
           csv_in_v, lx, map_inlet) 
   
       ! Make sure to sync all csv_in values to the device
       ! NOTE: This is actually not necessary since we do the copying in user_bc
       if ( NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(csv_in_v%x, csv_in_v%x_d, csv_in_v%n, HOST_TO_DEVICE, sync=.true. )
       end if
    end if 

    !
    ! Map the csv values with the outlet boundary points
    !
    if (u%msh%labeled_zones(2)%size .ne. 0) then
       
       call neko_log%message("mapping outlet...")
       call map_boundary(u%dof%x, u%dof%y, u%msh%labeled_zones(2)%facet_el, & 
               csv_in_o, lx, map_outlet)
       
       ! Make sure to sync all csv_in values to the device
       ! NOTE: This is actually not necessary since we do the copying in user_bc
       if ( NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(csv_in_o%x, csv_in_o%x_d, csv_in_o%n, HOST_TO_DEVICE, sync=.true.)
       end if
    end if

  end subroutine initialize

  !> Set boundary conditions
  subroutine user_bc(field_bc_list, t, tstep)
    type(field_list_t), intent(inout) :: field_bc_list
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    integer :: nlidx(4), i, n, lx

    lx = field_bc_list%fields(1)%f%Xh%lx
    
    associate( u => field_bc_list%fields(1)%f, & 
               v => field_bc_list%fields(2)%f, &
               w => field_bc_list%fields(3)%f, &
               p => field_bc_list%fields(4)%f)

    n = max(size(map_outlet), size(map_inlet))

    call neko_log%message("applying BC")
    do i = 1, n

       ! Apply inlet
       if (i .le. size(map_inlet)) then
          nlidx = nonlinear_index(map_inlet(i)%x(1), lx, lx, lx)
          u%x(nlidx(1), nlidx(2), nlidx(3), nlidx(4)) = csv_in_v%x(map_inlet(i)%x(2), 4)
          v%x(nlidx(1), nlidx(2), nlidx(3), nlidx(4)) = csv_in_v%x(map_inlet(i)%x(2), 5)
          w%x(nlidx(1), nlidx(2), nlidx(3), nlidx(4)) = 0.0_rp
       end if
       
       ! Apply outlet
       if (i .le. size(map_outlet)) then
          nlidx = nonlinear_index(map_outlet(i)%x(1), lx, lx, lx)
          p%x(nlidx(1), nlidx(2), nlidx(3), nlidx(4)) = csv_in_o%x(map_outlet(i)%x(2), 3)
       end if

    end do
    
    call neko_log%message("end applying BC")
    end associate

    !
    ! Force copy the field_bcs to the device
    !
    call neko_log%message("copying BC")
    do i = 1, 4
       call device_memcpy(field_bc_list%fields(i)%f%x, field_bc_list%fields(i)%f%x_D, &
            field_bc_list%fields(i)%f%dof%size(), HOST_TO_DEVICE, sync = .true.)
    end do
    call neko_log%message("end copying BC")
  end subroutine user_bc

  ! usrcheck, this is called at the end of every time step
  subroutine check(t, tstep,u, v, w, p, coef, param)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: param

  end subroutine check

  ! Free relevant objects 
  subroutine finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    call csv_in_o%free
    call csv_in_v%free
    deallocate(map_outlet)
    deallocate(map_inlet)

  end subroutine finalize

  !> Map points from target GLL points to those given in a
  !! csv_file
  subroutine map_boundary(x_GLL, y_GLL, facet_el, csv_data, lx, map_array)
    real(kind=rp), intent(in) :: x_GLL(:,:,:,:)
    real(kind=rp), intent(in) :: y_GLL(:,:,:,:)
    type(tuple_i4_t), intent(in) :: facet_el(:)
    type(matrix_t), intent(in) :: csv_data
    integer, intent(in) :: lx
    type(tuple_i4_t), intent(inout) :: map_array(:)

    ! -- indices
    integer :: i,j,k, ii,jj,kk, idx, nl_index(4), l_index
    integer :: start_ii, start_jj, start_kk, end_ii, end_jj, end_kk
    integer :: map_index

    ! -- sizes
    integer :: n_outlets

    ! -- misc
    integer :: f, el
    real(kind=rp) :: uvp(3), x_target, y_target

    !
    ! Initialize stuff
    ! 
    n_outlets = size(facet_el)

    ! 
    ! Perform the search and mapping
    !
    map_index = 1
    do i = 1, n_outlets

       f = facet_el(i)%x(1)
       el = facet_el(i)%x(2)

       ! Set boundaries for the loop through the GLL points depending on the facet
       select case (f)
       case (1)
          start_ii = 1
          start_jj = 1
          start_kk = 1
          end_ii   = 1
          end_jj   = lx
          end_kk   = lx
       case (2)
          start_ii = lx
          start_jj = 1
          start_kk = 1
          end_ii   = lx
          end_jj   = lx
          end_kk   = lx
       case (3)
          start_ii = 1
          start_jj = 1
          start_kk = 1
          end_ii   = lx
          end_jj   = 1
          end_kk   = lx
       case (4)
          start_ii = 1
          start_jj = lx
          start_kk = 1
          end_ii   = lx
          end_jj   = lx
          end_kk   = lx
       case (5)
          start_ii = 1
          start_jj = 1
          start_kk = 1
          end_ii   = lx
          end_jj   = lx
          end_kk   = 1
       case (6)
          start_ii = 1
          start_jj = 1
          start_kk = lx
          end_ii   = lx
          end_jj   = lx
          end_kk   = lx
       end select

       do kk = start_kk, end_kk
          do jj = start_jj, end_jj
             do ii = start_ii, end_ii

                x_target =  x_GLL(ii,jj,kk,el)
                y_target =  y_GLL(ii,jj,kk,el)

                ! Find the (x,y) coordinates in the csv file that are closest
                ! to our target x,y coordinates and return the corresponding
                ! u,v,p values + its index in the csv array.
                uvp = find_coordinates_and_values(x_target, y_target, &
                     csv_data, idx)
                
                ! Build the mask such that every linear index k in the mesh
                ! is mapped to an entry in o.csv
                k = linear_index(ii,jj,kk,el,lx,lx,lx) ! index of the gll point
                map_array(map_index)%x(1) = k
                map_array(map_index)%x(2) = idx ! index of the entry in the csv 
                map_index = map_index + 1

             end do
          end do
       end do
    end do

  end subroutine map_boundary

  ! Search an array for matching coordinates and output u,v,p when found
  function find_coordinates_and_values(x_target, y_target, search_array, idx) result(res)
    real(kind=rp), intent(in) :: x_target
    real(kind=rp), intent(in) :: y_target
    type(matrix_t), intent(in) :: search_array
    integer, intent(out) :: idx

    real(kind=rp) :: res(3) ! will store u,v,p
    integer :: i
    real(kind=rp) :: dist, dist_ref, x, y

    res = 0.0_rp
    dist = 1d10
    dist_ref = 1d10

    do i = 1, search_array%nrows

       x = search_array%x(i,1)
       y = search_array%x(i,2)
       dist = sqrt((x-x_target)**2 + (y-y_target)**2)

       if (dist .lt. dist_ref) then
          dist_ref = dist
          res(1) = search_array%x(i,4)  ! u
          res(2) = search_array%x(i,5)  ! v
          res(3) = search_array%x(i,3)  ! p
          idx = i ! This will help us build a mapping with each GLL point
       end if

    end do

    if (dist_ref .gt. 1d-3) call neko_error("Found distance less than &
         1d-3, point may not be mapped properly!")

  end function find_coordinates_and_values

  ! Read field file from filename and load it into data
  ! NOTE: data has to be initialized before it is used
  subroutine read_field_file(fld_filename, data)
    character(len=*), intent(in) :: fld_filename
    type(fld_file_data_t), intent(inout) :: data

    type(file_t) :: fld_file

    fld_file = file_t(trim(fld_filename))
    call fld_file%read(data)
    call file_free(fld_file)

  end subroutine read_field_file

end module user
