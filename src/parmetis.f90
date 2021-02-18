!> Interface to ParMETIS
module parmetis
  use comm
  use point
  use mesh_field 
  use mesh, only : mesh_t
  use, intrinsic :: iso_c_binding
  implicit none
  private

  integer, private, parameter :: METIS_OK = 1, METIS_ERROR = -4

  interface
     integer (c_int) function parmetis_v3_partgeom &
          (vtxdist, ndims, xyz, part) &
          bind(c, name='ParMETIS_V3_PartGeom_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
       ! idx_t variables
       type(c_ptr), value :: vtxdist, ndims, part
       ! real_t variables
       type(c_ptr), value :: xyz
     end function parmetis_v3_partgeom
  end interface

contains

  !> Compute a k-way partitioning of a mesh @a msh using
  !! a coordinated-based space-filing curves method
  subroutine parmetis_partgeom(msh, p)
    type(mesh_t), intent(inout) :: msh   !< Mesh
    type(mesh_fld_t), intent(inout) :: p !< Partitions
    integer(kind=c_int), target :: ndims, rcode
    real(kind=c_float), allocatable, target, dimension(:) :: xyz
    integer(kind=c_int), allocatable, target, dimension(:) :: vtxdist, part
    type(point_t) :: c
    integer :: i, j, ierr

    ndims = msh%gdim

    allocate(part(msh%nelv), xyz(ndims * msh%nelv))
    allocate(vtxdist(0:pe_size))

    call parmetis_vtxdist(vtxdist, msh%nelv)
    
    i = 1
    do j = 1, msh%nelv
       c = msh%elements(j)%e%centroid()
       xyz(i) = c%x(1)
       xyz(i + 1) = c%x(2)
       xyz(i + 2) = c%x(3)
       i = i + 3
    end do
    
    rcode = parmetis_v3_partgeom(c_loc(vtxdist), c_loc(ndims), &
         c_loc(xyz), c_loc(part))
    
    if (rcode .eq. METIS_OK) then
       ! Mark partitions
       call mesh_field_init(p, msh, 'partitions')
       do i = 1, msh%nelv
          p%data(i) = part(i)
       end do
    else
       call neko_error(rcode)
    end if

    deallocate(part, xyz)
    
  end subroutine parmetis_partgeom

  !> Compute the (parallel) vertex distribution of the dual graph
  subroutine parmetis_vtxdist(vtxdist, nelv)
    integer(kind=c_int), intent(inout) :: vtxdist(0:pe_size)
    integer, intent(in) :: nelv
    integer :: i, ierr, tmp, sum

    vtxdist(pe_rank) = nelv

    call MPI_Allgather(nelv, 1, MPI_INTEGER, vtxdist, 1, &
         MPI_INTEGER, NEKO_COMM, ierr)

    sum = vtxdist(0)
    do i = 1, pe_size
       tmp = vtxdist(i)
       vtxdist(i) = sum
       sum = tmp + sum
    end do
    
  end subroutine parmetis_vtxdist
  
end module parmetis
