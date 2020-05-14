!> Neko binary mesh data
module nmsh_file
  use generic_file
  use comm
  use mesh
  use utils
  use nmsh
  use datadist
  use mpi_types
  use mpi
  implicit none
  
  private

  !> Interface for Neko nmsh files
  type, public, extends(generic_file_t) :: nmsh_file_t
   contains
     procedure :: read => nmsh_file_read
     procedure :: write => nmsh_file_write
  end type nmsh_file_t

contains

  !> Load a mesh from a binary Neko nmsh file
  subroutine nmsh_file_read(this, data)
    class(nmsh_file_t) :: this
    class(*), target, intent(inout) :: data    
    type(nmsh_hex_t), allocatable :: nmsh_hex(:)
    type(nmsh_quad_t), allocatable :: nmsh_quad(:)
    type(mesh_t), pointer :: msh
    integer :: status(MPI_STATUS_SIZE)
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset
    integer :: i, j, ierr, fh, nelgv, element_offset
    integer :: nmsh_quad_size, nmsh_hex_size
    class(element_t), pointer :: ep
    integer nelv, gdim
    type(point_t) :: p(8)
    type(linear_dist_t) :: dist

    select type(data)
    type is(mesh_t)
       msh => data
    class default
       call neko_error('Invalid output data')
    end select

    if (pe_rank .eq. 0) then
       write(*, '(A,A)') " Reading a binary Neko file ", this%fname
    end if

    call MPI_Type_size(MPI_NMSH_HEX, nmsh_hex_size, ierr)
    call MPI_Type_size(MPI_NMSH_QUAD, nmsh_quad_size, ierr)

    call MPI_File_open(NEKO_COMM, trim(this%fname), &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_File_read_all(fh, nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, gdim, 1, MPI_INTEGER, status, ierr)

    if (pe_rank .eq. 0) then
       write(*,1) gdim, nelv
1      format(1x,'gdim = ', i1, ', nelements =', i7)
    end if
       
    dist = linear_dist_t(nelv, pe_rank, pe_size, NEKO_COMM)
    nelv = dist%num_local()
    element_offset = dist%start_idx()
    
    call mesh_init(msh, gdim, nelv)
   

    if (msh%gdim .eq. 2) then
       allocate(nmsh_quad(msh%nelv))
       mpi_offset = 2 * MPI_INTEGER_SIZE + element_offset * nmsh_quad_size
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_quad, msh%nelv, MPI_NMSH_QUAD, status, ierr)
       do i = 1, nelv
          do j = 1, 4
             p(j) = point_t(nmsh_quad(i)%v(j)%v_xyz, nmsh_quad(i)%v(j)%v_idx)
          end do
          call mesh_add_element(msh, i, p(1), p(2), p(3), p(4))
       end do
       deallocate(nmsh_quad)
    else if (msh%gdim .eq. 3) then
       allocate(nmsh_hex(msh%nelv))
       mpi_offset = 2 * MPI_INTEGER_SIZE + element_offset * nmsh_hex_size
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_hex, msh%nelv, MPI_NMSH_HEX, status, ierr)
       do i = 1, nelv
          do j = 1, 8
             p(j) = point_t(nmsh_hex(i)%v(j)%v_xyz, nmsh_hex(i)%v(j)%v_idx)
          end do
          call mesh_add_element(msh, i, &
               p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8))
       end do
       deallocate(nmsh_hex)
    else        
       if (pe_rank .eq. 0) call neko_error('Invalid dimension of mesh')
    end if

    call MPI_File_close(fh, ierr)
    if (pe_rank .eq. 0) write(*,*) 'Done'
       
  end subroutine nmsh_file_read

    !> Load a mesh from a binary Neko nmsh file
  subroutine nmsh_file_write(this, data)
    class(nmsh_file_t), intent(in) :: this
    class(*), target, intent(in) :: data  
    type(nmsh_quad_t), allocatable :: nmsh_quad(:)
    type(nmsh_hex_t), allocatable :: nmsh_hex(:)
    type(mesh_t), pointer :: msh
    integer :: status(MPI_STATUS_SIZE)
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset
    integer :: i, j, ierr, fh, nelgv, element_offset
    integer :: nmsh_quad_size, nmsh_hex_size
    class(element_t), pointer :: ep

    select type(data)
    type is (mesh_t)
       msh => data
    class default
       call neko_error('Invalid output data')
    end select

    call MPI_Type_size(MPI_NMSH_QUAD, nmsh_quad_size, ierr)
    call MPI_Type_size(MPI_NMSH_HEX, nmsh_hex_size, ierr)

    call MPI_Reduce(msh%nelv, nelgv, 1, MPI_INTEGER, &
         MPI_SUM, 0, NEKO_COMM, ierr)
    element_offset = 0
    call MPI_Exscan(msh%nelv, element_offset, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    if (pe_rank .eq. 0) then
       write(*, '(A,A)') " Writing data as a binary Neko file ", this%fname
    end if

    call MPI_File_open(NEKO_COMM, trim(this%fname), &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

    call MPI_File_write_all(fh, msh%nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_write_all(fh, msh%gdim, 1, MPI_INTEGER, status, ierr)

    if (msh%gdim .eq. 2) then
       allocate(nmsh_quad(msh%nelv))       
       do i = 1, msh%nelv
          ep => msh%elements(i)%e
          nmsh_quad(i)%el_idx = ep%id()
          do j = 1, 4
             nmsh_quad(i)%v(j)%v_idx = ep%pts(j)%p%id()
             nmsh_quad(i)%v(j)%v_xyz = ep%pts(j)%p%x
          end do
       end do
       mpi_offset = 2 * MPI_INTEGER_SIZE + element_offset * nmsh_quad_size
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_quad, msh%nelv, MPI_NMSH_QUAD, status, ierr)
       deallocate(nmsh_quad)
    else if (msh%gdim .eq. 3) then
       allocate(nmsh_hex(msh%nelv))       
       do i = 1, msh%nelv
          ep => msh%elements(i)%e
          nmsh_hex(i)%el_idx = ep%id()
          do j = 1, 8
             nmsh_hex(i)%v(j)%v_idx = ep%pts(j)%p%id()
             nmsh_hex(i)%v(j)%v_xyz = ep%pts(j)%p%x
          end do
       end do
       mpi_offset = 2 * MPI_INTEGER_SIZE + element_offset * nmsh_hex_size
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_HEX, msh%nelv, MPI_NMSH_HEX, status, ierr)

       deallocate(nmsh_hex)
    else 
       call neko_error('Invalid dimension of mesh')
    end if


    call MPI_File_close(fh, ierr)
    if (pe_rank .eq. 0) write(*,*) 'Done'

  end subroutine nmsh_file_write
  
end module nmsh_file
  
