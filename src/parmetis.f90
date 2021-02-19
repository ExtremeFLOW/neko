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

  public :: parmetis_partgeom, parmetis_partmeshkway

  interface
     integer (c_int) function parmetis_v3_partmeshkway &
          (elmdist, eptr, eind, elmwgt, wgtflag, numflag, ncon, &
          ncommonnodes, nparts, tpwgts, ubvec, options, edgecut, part) &
          bind(c, name='ParMETIS_V3_PartMeshKway_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
       ! idx_t variables
       type(c_ptr), value :: elmdist, eptr, eind, elmwgt, wgtflag, &
            numflag, ncon, ncommonnodes, nparts, options, edgecut, part
       ! real_t variables
       type(c_ptr), value :: tpwgts, ubvec
     end function parmetis_v3_partmeshkway
  end interface
  
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

  !> Compute a k-way partitioning of a mesh @a msh 
  subroutine parmetis_partmeshkway(msh, parts, weights)
    type(mesh_t), intent(inout) :: msh                !< Mesh
    type(mesh_fld_t), intent(inout) :: parts          !< Partitions
    type(mesh_fld_t), intent(in), optional :: weights !< Weights
    integer(kind=c_int), target :: wgtflag, numflag, ncon, ncommonnodes
    integer(kind=c_int), target :: nparts, options(3), edgecut, rcode
    real(kind=c_float), allocatable, target, dimension(:) :: tpwgts, ubvec
    integer(kind=c_int), allocatable, target, dimension(:) :: &
         elmdist, eptr, eind, elmwgt, part
    integer :: i, j, k, ierr

    !> @note We use C-style numbering in the element distribution.
    !! Otherwise, the partitions return by Metis will not be numbered
    !! in the same way as the MPI ranks
    numflag = 0
    ncon = 1
    nparts = pe_size
    ncommonnodes = 2**(msh%gdim - 1)
    options(1) = 1
    options(2) = 1
    options(3) = 15 * pe_rank
    wgtflag = 2
    
    allocate(elmdist(0:pe_size), eptr(0:msh%nelv))
    allocate(eind(0:(msh%nelv * msh%npts)), part(msh%nelv))
    allocate(elmwgt(msh%nelv), tpwgts(ncon * pe_size), ubvec(ncon)) 

    call parmetis_dist(elmdist, msh%nelv)

    if (present(weights)) then
       call parmetis_wgt(msh, elmwgt, tpwgts, ubvec, nparts, ncon, weights)
    else
       call parmetis_wgt(msh, elmwgt, tpwgts, ubvec, nparts, ncon)
    end if

    eptr(0) = 0
    do i = 1, msh%nelv
       eptr(i) = eptr(i - 1) + msh%npts
    end do

    k = 0
    do i = 1, msh%nelv
       do j = 1, msh%npts
          eind(k) = msh%elements(i)%e%pts(j)%p%id() - 1
          k = k + 1
       end do
    end do

    rcode = parmetis_v3_partmeshkway(c_loc(elmdist), c_loc(eptr), c_loc(eind), &
         c_loc(elmwgt), c_loc(wgtflag), c_loc(numflag), c_loc(ncon), &
         c_loc(ncommonnodes), c_loc(nparts), c_loc(tpwgts), c_loc(ubvec),&
         c_loc(options), c_loc(edgecut), c_loc(part))

    if (rcode .eq. METIS_OK) then
       call parmetis_mark_parts(parts, msh, part)
    else
       call neko_error(rcode)
    end if
    
    deallocate(elmdist, eptr, eind, part, elmwgt, tpwgts, ubvec)

  end subroutine parmetis_partmeshkway
  
  !> Compute a k-way partitioning of a mesh @a msh using
  !! a coordinated-based space-filing curves method
  subroutine parmetis_partgeom(msh, parts)
    type(mesh_t), intent(inout) :: msh       !< Mesh
    type(mesh_fld_t), intent(inout) :: parts !< Partitions
    integer(kind=c_int), target :: ndims, rcode
    real(kind=c_float), allocatable, target, dimension(:) :: xyz
    integer(kind=c_int), allocatable, target, dimension(:) :: vtxdist, part
    type(point_t) :: c
    integer :: i, j, ierr

    ndims = msh%gdim

    allocate(part(msh%nelv), xyz(ndims * msh%nelv))
    allocate(vtxdist(0:pe_size))

    call parmetis_dist(vtxdist, msh%nelv)
    
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
       call parmetis_mark_parts(parts, msh, part)
    else
       call neko_error(rcode)
    end if

    deallocate(part, xyz, vtxdist)    
    
  end subroutine parmetis_partgeom

  !> Fill mesh field according to new partitions
  subroutine parmetis_mark_parts(parts, msh, part)
    type(mesh_fld_t), intent(inout) :: parts
    type(mesh_t), intent(in) :: msh
    integer(kind=c_int), allocatable, intent(in) :: part(:)
    integer :: i

    call mesh_field_init(parts, msh, 'partitions')

    do i = 1, msh%nelv
       parts%data(i) = part(i)
    end do
    
  end subroutine parmetis_mark_parts

  !> Setup weights and balance constraints for the dual graph
  subroutine parmetis_wgt(msh, wgt, tpwgts, ubvec, nparts, ncon, weight)
    type(mesh_t), intent(in) :: msh
    integer(kind=c_int), allocatable, intent(inout) :: wgt(:)
    real(kind=c_float), allocatable, intent(inout) :: tpwgts(:)
    real(kind=c_float), allocatable, intent(inout) :: ubvec(:)
    integer, intent(in) :: nparts, ncon
    type(mesh_fld_t), intent(in), optional :: weight
    integer :: i
    
    if (present(weight)) then
       do i = 1, msh%nelv
          wgt(i) = weight%data(i)
       end do
    else
       wgt = 1.0
    end if
    
    do i = 1, nparts
       tpwgts(i) = 1.0e0 / real(nparts)
    end do

    do i = 1, ncon
       ubvec(i) = 1.05e0
    end do

  end subroutine parmetis_wgt
  
  !> Compute the (parallel) vertex distribution of the dual graph
  subroutine parmetis_dist(dist, nelv)
    integer(kind=c_int), intent(inout) :: dist(0:pe_size)
    integer, intent(in) :: nelv
    integer :: i, ierr, tmp, sum

    dist(pe_rank) = nelv

    call MPI_Allgather(nelv, 1, MPI_INTEGER, dist, 1, &
         MPI_INTEGER, NEKO_COMM, ierr)

    sum = dist(0)
    dist(0) = 0
    do i = 1, pe_size
       tmp = dist(i)
       dist(i) = sum
       sum = tmp + sum
    end do

  end subroutine parmetis_dist
  
end module parmetis
