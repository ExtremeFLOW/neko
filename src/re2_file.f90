!> NEKTON mesh data in re2 format
!! @details This module is used to read/write binary NEKTION mesh data
module re2_file
  use generic_file
  use num_types
  use utils
  use mesh
  use point
  use comm
  use mpi
  use mpi_types
  use datadist
  use re2
  use map
  use map_file
  use htable
  implicit none
  private
  

  !> Interface for NEKTON re2 files
  type, public, extends(generic_file_t) :: re2_file_t
   contains
     procedure :: read => re2_file_read
     procedure :: write => re2_file_write
  end type re2_file_t

contains

  !> Load a binary NEKTON mesh from a re2 file
  subroutine re2_file_read(this, data)
    class(re2_file_t) :: this
    class(*), target, intent(inout) :: data
    type(re2_xy_t), allocatable :: re2_data_xy(:)
    type(re2_xyz_t), allocatable :: re2_data_xyz(:)
    type(re2_bc_t), allocatable :: re2_data_bc(:)
    type(mesh_t), pointer :: msh
    character(len=5) :: hdr_ver
    character(len=54) :: hdr_str
    integer :: i, j, k, fh, nel, ndim, nelv, ierr, pt_idx, el_idx
    integer :: status(MPI_STATUS_SIZE)
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset
    real(kind=sp) :: test
    real(kind=dp) :: t2
    type(point_t) :: p(8)
    integer :: ncurv, nbcs
    type(linear_dist_t) :: dist
    type(map_t) :: nm
    type(map_file_t) :: map_file
    character(len=80) :: map_fname
    logical :: read_map
    integer :: element_offset
    integer :: re2_data_xy_size
    integer :: re2_data_xyz_size
    integer :: re2_data_cv_size
    integer :: re2_data_bc_size
    integer :: pids(4)
    type(htable_pt_t) :: htp
    integer :: sym_facet
    integer, parameter, dimension(6) :: facet_map = (/3, 2, 4, 1, 5, 6/)
    integer :: p_el_idx, p_facet

    
    select type(data)
    type is (mesh_t)
       msh => data
    class default
       call neko_error('Invalid output data')
    end select


    call MPI_Type_size(MPI_RE2_DATA_XY, re2_data_xy_size, ierr)
    call MPI_Type_size(MPI_RE2_DATA_XYZ, re2_data_xyz_size, ierr)
    call MPI_Type_size(MPI_RE2_DATA_CV, re2_data_cv_size, ierr)
    call MPI_Type_size(MPI_RE2_DATA_BC, re2_data_bc_size, ierr)
    
    open(unit=9,file=trim(this%fname), status='old', iostat=ierr)
    if (pe_rank .eq. 0) then
       write(*, '(A,A)') " Reading binary NEKTON file ", this%fname
    end if
    read(9, '(a5,i9,i3,i9,a54)') hdr_ver, nel, ndim, nelv, hdr_str
    if (hdr_ver .eq. '#v002') then
       call neko_error("Double precision re2 files not supported yet")
    end if

    if (pe_rank .eq. 0) write(*,1) ndim, nelv
1   format(1x,'ndim = ', i1, ', nelements =', i7)
    close(9)

    call filename_chsuffix(this%fname, map_fname,'map')

    inquire(file=map_fname, exist=read_map)
    if (read_map) then
       call map_init(nm, nelv, 2**ndim)
       call map_file%init(map_fname)
       call map_file%read(nm)
    else
       if (pe_rank .eq. 0) call neko_warning('No NEKTON map file found')
    end if

    call MPI_File_open(NEKO_COMM, trim(this%fname), &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    
    if (ierr .ne. 0) then
       call neko_error("Can't open binary NEKTON file ")
    end if
    dist = linear_dist_t(nelv, pe_rank, pe_size, NEKO_COMM)
    nelv = dist%num_local()
    element_offset = dist%start_idx()

    call mesh_init(msh, ndim, dist)

    call htp%init((2*ndim) * nel, ndim)
   
    ! Set offset (header)
    mpi_offset = RE2_HDR_SIZE * MPI_CHARACTER_SIZE

    call MPI_File_read_at_all(fh, mpi_offset, test, 1, MPI_REAL, status, ierr)
    mpi_offset = mpi_offset + MPI_REAL_SIZE
    
    if (abs(RE2_ENDIAN_TEST - test) .gt. 1e-4) then
       call neko_error('Invalid endian of re2 file, byte swap not implemented yet')
    end if
    
    pt_idx = 0
    if (ndim .eq. 2) then
       allocate(re2_data_xy(nelv))
       mpi_offset = mpi_offset + element_offset * re2_data_xy_size
       call MPI_File_read_at_all(fh, mpi_offset, &
            re2_data_xy, nelv, MPI_RE2_DATA_XY, status, ierr)
       do i = 1, nelv
          do j = 1, 4             
             p(j) = point_t(dble(re2_data_xy(i)%x(j)), &
                  dble(re2_data_xy(i)%y(j)), 0d0)
             call re2_file_add_point(htp, p(j), pt_idx)
          end do

          call mesh_add_element(msh, i, p(1), p(2), p(3), p(4))
       end do
       deallocate(re2_data_xy)
    else if (ndim .eq. 3) then
       allocate(re2_data_xyz(nelv))
       mpi_offset = mpi_offset + element_offset * re2_data_xyz_size
       call MPI_File_read_at_all(fh, mpi_offset, &
            re2_data_xyz, nelv, MPI_RE2_DATA_XYZ, status, ierr)
       do i = 1, nelv
          do j = 1, 8             
             p(j) = point_t(dble(re2_data_xyz(i)%x(j)), &
                  dble(re2_data_xyz(i)%y(j)),&
                  dble(re2_data_xyz(i)%z(j)))
             call re2_file_add_point(htp, p(j), pt_idx)
          end do

          call mesh_add_element(msh, i, &
               p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8))          
       end do
       deallocate(re2_data_xyz)
    end if

    call htp%free()


    ! Set offset to start of curved side data
    mpi_offset = RE2_HDR_SIZE * MPI_CHARACTER_SIZE + MPI_REAL_SIZE
    if (ndim .eq. 2) then
       mpi_offset = mpi_offset + dist%num_global() * re2_data_xy_size
    else
       mpi_offset = mpi_offset + dist%num_global() * re2_data_xyz_size
    end if

    !> @todo Add support for curved side data
    !! Skip curved side data
    call MPI_File_read_at_all(fh, mpi_offset, ncurv, 1, MPI_INTEGER, status, ierr)
    mpi_offset = mpi_offset + MPI_INTEGER_SIZE + ncurv * RE2_DATA_CV_SIZE

    call MPI_File_read_at_all(fh, mpi_offset, nbcs, 1, MPI_INTEGER, status, ierr)
    mpi_offset = mpi_offset + MPI_INTEGER_SIZE

    allocate(re2_data_bc(nbcs))

    call MPI_File_read_at_all(fh, mpi_offset, re2_data_bc, nbcs, &
         MPI_RE2_DATA_BC, status, ierr)

    !> @todo Use element offset in parallel
    do i = 1, nbcs
       el_idx = re2_data_bc(i)%elem - dist%start_idx()
       sym_facet = facet_map(re2_data_bc(i)%face)
       select case(trim(re2_data_bc(i)%type))
       case ('W')
          call mesh_mark_wall_facet(msh, sym_facet, el_idx)
       case ('v', 'V')
          call mesh_mark_inlet_facet(msh, sym_facet, el_idx)
       case ('O', 'o')
          call mesh_mark_outlet_facet(msh, sym_facet, el_idx)
       case ('SYM')
          call mesh_mark_sympln_facet(msh, sym_facet, el_idx)
       case ('P')
          p_el_idx = int(re2_data_bc(i)%bc_data(1))
          p_facet = facet_map(int(re2_data_bc(i)%bc_data(2)))
          call mesh_get_periodic_ids(msh, sym_facet, el_idx, &
                                     p_facet, p_el_idx, pids)
          call mesh_mark_periodic_facet(msh, sym_facet, el_idx, &
                                        p_facet, p_el_idx, pids)
       end select
    end do
    !> @todo Use element offset in parallel
    do i = 1, nbcs
       el_idx = re2_data_bc(i)%elem - dist%start_idx()
       sym_facet = facet_map(re2_data_bc(i)%face)
       select case(trim(re2_data_bc(i)%type))
       case ('P')
          p_el_idx = int(re2_data_bc(i)%bc_data(1))
          p_facet = facet_map(int(re2_data_bc(i)%bc_data(2)))
          call mesh_create_periodic_ids(msh, sym_facet, el_idx, &
              p_facet, p_el_idx) 
       end select
    end do
    !> @todo Use element offset in parallel
    do i = 1, nbcs
       el_idx = re2_data_bc(i)%elem - dist%start_idx()
       sym_facet = facet_map(re2_data_bc(i)%face)
       select case(trim(re2_data_bc(i)%type))
       case ('P')
          p_el_idx = int(re2_data_bc(i)%bc_data(1))
          p_facet = facet_map(int(re2_data_bc(i)%bc_data(2)))
          call mesh_create_periodic_ids(msh, sym_facet, el_idx, &
              p_facet, p_el_idx) 
       end select
    end do   
    !> @todo Use element offset in parallel
    do i = 1, nbcs
       el_idx = re2_data_bc(i)%elem - dist%start_idx()
       sym_facet = facet_map(re2_data_bc(i)%face)
       select case(trim(re2_data_bc(i)%type))
       case ('P')
          p_el_idx = int(re2_data_bc(i)%bc_data(1))
          p_facet = facet_map(int(re2_data_bc(i)%bc_data(2)))
          call mesh_create_periodic_ids(msh, sym_facet, el_idx, &
              p_facet, p_el_idx) 
       end select
    end do   
 
    !>@ todo process bc data here 
    deallocate(re2_data_bc)

    call MPI_FILE_close(fh, ierr)

    call mesh_finalize(msh)

    if (pe_rank .eq. 0) write(*,*) 'Done'

    
  end subroutine re2_file_read

  subroutine re2_file_write(this, data, t)
    class(re2_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=dp), intent(in), optional :: t
    type(re2_xy_t), allocatable :: re2_data_xy(:)
    type(re2_xyz_t), allocatable :: re2_data_xyz(:)
    type(mesh_t), pointer :: msh
    character(len=5), parameter :: RE2_HDR_VER = '#v001'
    character(len=54), parameter :: RE2_HDR_STR = 'RE2 exported by NEKO'
    integer :: i, j, k, fh, ierr, pt_idx, nelgv
    integer :: status(MPI_STATUS_SIZE)
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset
    integer :: element_offset
    integer :: re2_data_xy_size
    integer :: re2_data_xyz_size
        
    select type(data)
    type is (mesh_t)
       msh => data
    class default
       call neko_error('Invalid output data')
    end select

    call MPI_Type_size(MPI_RE2_DATA_XY, re2_data_xy_size, ierr)
    call MPI_Type_size(MPI_RE2_DATA_XYZ, re2_data_xyz_size, ierr)
    call MPI_Reduce(msh%nelv, nelgv, 1, &
         MPI_INTEGER, MPI_SUM, 0, NEKO_COMM, ierr)
    element_offset = 0
    call MPI_Exscan(msh%nelv, element_offset, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    if (pe_rank .eq. 0) then
       open(unit=9,file=trim(this%fname), status='new', iostat=ierr)
       write(*, '(A,A)') " Writing data as a binary NEKTON file ", this%fname
       write(9, '(a5,i9,i3,i9,a54)') RE2_HDR_VER, nelgv, msh%gdim,&
            nelgv, RE2_HDR_STR
       close(9)
    end if

    call MPI_Barrier(NEKO_COMM, ierr)
    call MPI_File_open(NEKO_COMM, trim(this%fname), &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)
    mpi_offset = RE2_HDR_SIZE * MPI_CHARACTER_SIZE
    
    call MPI_File_write_at(fh, mpi_offset, RE2_ENDIAN_TEST, 1, &
         MPI_REAL, status, ierr)
    mpi_offset = mpi_offset + MPI_REAL_SIZE

    if (msh%gdim .eq. 2) then
       allocate(re2_data_xy(msh%nelv))
       do i = 1, msh%nelv
          re2_data_xy(i)%rgroup = 1.0 ! Not used
          do j = 1, 4
             re2_data_xy(i)%x(j) = real(msh%elements(i)%e%pts(j)%p%x(1))
             re2_data_xy(i)%y(j) = real(msh%elements(i)%e%pts(j)%p%x(2))
          end do
       end do
       mpi_offset = mpi_offset + element_offset * re2_data_xy_size
       call MPI_File_write_at(fh, mpi_offset, &
            re2_data_xy, msh%nelv, MPI_RE2_DATA_XY, status, ierr)

       deallocate(re2_data_xy)

    else if (msh%gdim .eq. 3) then
       allocate(re2_data_xyz(msh%nelv))
       do i = 1, msh%nelv
          re2_data_xyz(i)%rgroup = 1.0 ! Not used
          do j = 1, 8 
             re2_data_xyz(i)%x(j) = real(msh%elements(i)%e%pts(j)%p%x(1))
             re2_data_xyz(i)%y(j) = real(msh%elements(i)%e%pts(j)%p%x(2))
             re2_data_xyz(i)%z(j) = real(msh%elements(i)%e%pts(j)%p%x(3))
          end do
       end do
       mpi_offset = mpi_offset + element_offset * re2_data_xyz_size
       call MPI_File_write_at(fh, mpi_offset, &
            re2_data_xyz, msh%nelv, MPI_RE2_DATA_XYZ, status, ierr)
       
       deallocate(re2_data_xyz)
    else
       call neko_error("Invalid dimension of mesh")
    end if
        
    call MPI_FILE_close(fh, ierr)
    write(*,*) 'Done'

    !> @todo Add support for curved side data
    
  end subroutine re2_file_write

  subroutine re2_file_add_point(htp, p, idx)
    type(htable_pt_t), intent(inout) :: htp
    type(point_t), intent(inout) :: p
    integer, intent(inout) :: idx
    integer :: tmp
    
    if (htp%get(p, tmp) .gt. 0) then
       idx = idx + 1
       call htp%set(p, idx)
       call p%set_id(idx)
    else
       call p%set_id(tmp)
    end if
    
  end subroutine re2_file_add_point
  
end module re2_file
