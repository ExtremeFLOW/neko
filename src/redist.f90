!> Redistribution routines
module redist
  use mesh_field
  use mpi_types
  use point
  use stack
  use comm
  use mesh
  use nmsh
  use mpi  
  implicit none

contains

  !> Redistribute a mesh @a msh according to new partitions
  subroutine redist_mesh(msh, parts)
    type(mesh_t), intent(inout), target :: msh    !< Mesh
    type(mesh_fld_t), intent(in) :: parts         !< Partitions
    type(stack_nh_t), allocatable :: new_dist(:)
    type(nmsh_hex_t) :: el
    type(nmsh_hex_t), pointer :: np(:)
    type(nmsh_hex_t), allocatable :: recv_buf(:)
    class(element_t), pointer :: ep
    integer :: status(MPI_STATUS_SIZE)
    integer :: i, j, ierr, max_recv, src, dst, recv_size
    integer :: gdim
    type(point_t) :: p(8)

    allocate(new_dist(0:pe_size - 1))
    do i = 0, pe_size - 1
       call new_dist(i)%init()
    end do

    do i = 1, msh%nelv
       ep => msh%elements(i)%e
       el%el_idx = ep%id()
       do j = 1, 8
          el%v(j)%v_idx = ep%pts(j)%p%id()
          el%v(j)%v_xyz = ep%pts(j)%p%x
       end do       
       call new_dist(parts%data(i))%push(el)
    end do

    gdim = msh%gdim    
    call mesh_free(msh)

    max_recv = 0
    do i = 0, pe_size - 1
       max_recv = max(max_recv, new_dist(i)%size())
    end do
    
    call MPI_Allreduce(MPI_IN_PLACE, max_recv, 1, MPI_INTEGER, &
         MPI_MAX, NEKO_COMM, ierr)
    allocate(recv_buf(max_recv))

    do i = 1, pe_size - 1
       src = modulo(pe_rank - i + pe_size, pe_size)
       dst = modulo(pe_rank + i, pe_size)

       call MPI_Sendrecv(new_dist(dst)%array(), new_dist(dst)%size(), &
            MPI_NMSH_HEX, dst, 0, recv_buf, max_recv, MPI_NMSH_HEX, src, 0,&
            NEKO_COMM, status, ierr)
       call MPI_Get_count(status, MPI_NMSH_HEX, recv_size, ierr)

       do j = 1, recv_size
          call new_dist(pe_rank)%push(recv_buf(j))
       end do       
    end do

    call mesh_init(msh, gdim, new_dist(pe_rank)%size())

    np => new_dist(pe_rank)%array()
    do i = 1, new_dist(pe_rank)%size()
       do j = 1, 8
          p(j) = point_t(np(i)%v(j)%v_xyz, np(i)%v(j)%v_idx)
       end do
       call mesh_add_element(msh, i, &
            p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8))
    end do

    call mesh_finalize(msh)
    
  end subroutine redist_mesh

end module redist
