! Copyright (c) 2021, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Interface to ParMETIS
module parmetis
  use comm
  use point
  use utils
  use num_types
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

!
! Define data types depending on ParMETIS' configuration
!
! parmetis_real converts between NEKO and ParMETIS' real representation
!
! parmetis_idx converts between NEKO and ParMETIS' integer representation
! neko_idx converts between ParMETIS and NEKO's integer representation
!

#ifdef HAVE_PARMETIS
#ifdef HAVE_PARMETIS_REAL64
#define M_REAL c_double
#define parmetis_real(i) real((i), 8)
#else
#define M_REAL c_float
#define parmetis_real(i) real((i), 4)
#endif
#ifdef HAVE_PARMETIS_INT64
#define M_INT c_int64_t
#define parmetis_idx(i) int((i), 8)
#define neko_idx(i) int((i), 4)
#else
#define M_INT c_int32_t
#define parmetis_idx(i) (i)
#define neko_idx(i) (i)
#endif
#endif
  
contains

#ifdef HAVE_PARMETIS

  !> Compute a k-way partitioning of a mesh @a msh 
  subroutine parmetis_partmeshkway(msh, parts, weights, nprts)
    type(mesh_t), intent(inout) :: msh                !< Mesh
    type(mesh_fld_t), intent(inout) :: parts          !< Partitions
    type(mesh_fld_t), intent(in), optional :: weights !< Weights
    integer, intent(in), optional :: nprts            !< Number of partitions
    integer(kind=M_INT), target :: wgtflag, numflag, ncon, ncommonnodes
    integer(kind=M_INT), target :: nparts, options(3), edgecut, rcode
    real(kind=M_REAL), allocatable, target, dimension(:) :: tpwgts, ubvec
    integer(kind=M_INT), allocatable, target, dimension(:) :: &
         elmdist, eptr, eind, elmwgt, part
    integer :: i, j, k, ierr

    !> @note We use C-style numbering in the element distribution.
    !! Otherwise, the partitions return by Metis will not be numbered
    !! in the same way as the MPI ranks
    numflag = 0
    ncon = 1
    ncommonnodes = 2**(msh%gdim - 1)
    options(1) = 1
    options(2) = 1
    options(3) = 15
    wgtflag = 2

    if (present(nprts)) then
       nparts = nprts
    else          
       nparts = pe_size
    end if
    
    allocate(elmdist(0:pe_size), eptr(0:msh%nelv))
    allocate(eind(0:(msh%nelv * msh%npts)), part(msh%nelv))
    allocate(elmwgt(msh%nelv), tpwgts(ncon * nparts), ubvec(ncon)) 

    call parmetis_dist(elmdist, msh%nelv)

    if (present(weights)) then
       call parmetis_wgt(msh, elmwgt, tpwgts, ubvec, nparts, ncon, weights)
    else
       call parmetis_wgt(msh, elmwgt, tpwgts, ubvec, nparts, ncon)
    end if

    eptr(0) = 0
    do i = 1, msh%nelv
       eptr(i) = parmetis_idx(eptr(i - 1) + msh%npts)
    end do

    k = 0
    do i = 1, msh%nelv
       do j = 1, msh%npts
          eind(k) = parmetis_idx(msh%elements(i)%e%pts(j)%p%id() - 1)
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
    integer(kind=M_INT), target :: ndims
    real(kind=M_REAL), allocatable, target, dimension(:) :: xyz
    integer(kind=M_INT), allocatable, target, dimension(:) :: vtxdist, part
    type(point_t) :: c
    integer :: i, j, ierr, rcode

    ndims = msh%gdim

    allocate(part(msh%nelv), xyz(ndims * msh%nelv))
    allocate(vtxdist(0:pe_size))

    call parmetis_dist(vtxdist, msh%nelv)
    
    i = 1
    do j = 1, msh%nelv
       c = msh%elements(j)%e%centroid()
       xyz(i) = parmetis_real(c%x(1))
       xyz(i + 1) = parmetis_real(c%x(2))
       xyz(i + 2) = parmetis_real(c%x(3))
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
    integer(kind=M_INT), allocatable, intent(in) :: part(:)
    integer :: i

    call mesh_field_init(parts, msh, 'partitions')

    do i = 1, msh%nelv
       parts%data(i) = neko_idx(part(i))
    end do
    
  end subroutine parmetis_mark_parts

  !> Setup weights and balance constraints for the dual graph
  subroutine parmetis_wgt(msh, wgt, tpwgts, ubvec, nparts, ncon, weight)
    type(mesh_t), intent(in) :: msh
    integer(kind=M_INT), allocatable, intent(inout) :: wgt(:)
    real(kind=M_REAL), allocatable, intent(inout) :: tpwgts(:)
    real(kind=M_REAL), allocatable, intent(inout) :: ubvec(:)
    integer, intent(in) :: nparts, ncon
    type(mesh_fld_t), intent(in), optional :: weight
    integer :: i
    
    if (present(weight)) then
       do i = 1, msh%nelv
          wgt(i) = parmetis_idx(weight%data(i))
       end do
    else
       wgt = parmetis_idx(1)
    end if
    
    do i = 1, (ncon * nparts)
       tpwgts(i) = parmetis_real(1) / parmetis_real(nparts)
    end do

    do i = 1, ncon
       ubvec(i) = parmetis_real(1.05d0)
    end do

  end subroutine parmetis_wgt
  
  !> Compute the (parallel) vertex distribution of the dual graph
  subroutine parmetis_dist(dist, nelv)
    integer(kind=M_INT), intent(inout) :: dist(0:pe_size)
    integer, intent(in) :: nelv
    integer(kind=M_INT) :: tmp, sum
    integer :: i, ierr

    dist(pe_rank) = parmetis_idx(nelv)

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

#else

  !> Compute a k-way partitioning of a mesh @a msh 
  subroutine parmetis_partmeshkway(msh, parts, weights, nprts)
    type(mesh_t), intent(inout) :: msh                !< Mesh
    type(mesh_fld_t), intent(inout) :: parts          !< Partitions
    type(mesh_fld_t), intent(in), optional :: weights !< Weights
    integer, intent(in), optional :: nprts            !< Number of partitions
    call neko_error('NEKO needs to be built with ParMETIS support')
  end subroutine parmetis_partmeshkway

  !> Compute a k-way partitioning of a mesh @a msh using
  !! a coordinated-based space-filing curves method
  subroutine parmetis_partgeom(msh, parts)
    type(mesh_t), intent(inout) :: msh       !< Mesh
    type(mesh_fld_t), intent(inout) :: parts !< Partitions
    call neko_error('NEKO needs to be built with ParMETIS support')
  end subroutine parmetis_partgeom
  
#endif
  
end module parmetis
