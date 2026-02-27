! Copyright (c) 2026, The Neko Authors
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
!> HDF5 file format
module vtkhdf_file
  use num_types, only : rp, dp, sp
  use generic_file, only : generic_file_t
  use checkpoint, only : chkp_t
  use utils, only : neko_error, filename_suffix_pos
  use mesh, only : mesh_t
  use field, only : field_t, field_ptr_t
  use field_list, only : field_list_t
  use field_series, only : field_series_t, field_series_ptr_t
  use dofmap, only : dofmap_t
  use logger, only : neko_log
  use comm, only : pe_rank, NEKO_COMM
  use device, only : DEVICE_TO_HOST
  use mpi_f08, only : MPI_INFO_NULL, MPI_Allreduce, MPI_Allgather, MPI_IN_PLACE, &
       MPI_INTEGER, MPI_INTEGER8, MPI_SUM, MPI_Comm_size
#ifdef HAVE_HDF5
  use hdf5
#endif
  implicit none
  private

  !> Interface for HDF5 files
  type, public, extends(generic_file_t) :: vtkhdf_file_t
   contains
     procedure :: read => vtkhdf_file_read
     procedure :: write => vtkhdf_file_write
     procedure :: set_overwrite => vtkhdf_file_set_overwrite
  end type vtkhdf_file_t

  integer, dimension(2), parameter :: vtkhdf_version = [2, 4]

contains

#ifdef HAVE_HDF5

  !> Write data in HDF5 format following official VTKHDF UnstructuredGrid specification
  subroutine vtkhdf_file_write(this, data, t)
    class(vtkhdf_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dof
    type(field_t), pointer :: fld, u, v, w
    type(field_ptr_t), allocatable :: fp(:)
    logical, allocatable :: field_written(:)
    integer :: ierr, info, drank, i, j, ie, n_fields
    integer(hid_t) :: plist_id, file_id, dset_id, grp_id, attr_id, vtkhdf_grp
    integer(hid_t) :: dcpl_id
    integer(hid_t) :: pointdata_grp, celldata_grp, fielddata_grp
    integer(hid_t) :: filespace, memspace
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer(hsize_t), dimension(2) :: dcount2, doffset2
    integer(hsize_t), dimension(:), allocatable :: dims
    integer(hsize_t), dimension(2) :: vdims, maxdims
    integer(hsize_t), dimension(1) :: chunkdims
    integer :: suffix_pos, lx, ly, lz, nelv, mpts, npts
    integer :: num_partitions
    integer :: local_points, local_cells, local_conn
    integer :: total_points, total_cells, total_conn, total_offsets
    integer :: point_offset, cell_offset, conn_offset, offsets_offset
    integer, allocatable :: part_points(:), part_cells(:), part_conns(:)
    character(len=5) :: id_str
    character(len=1024) :: fname
    character(len=16), dimension(1) :: type_str
    character(len=128) :: field_name
    integer :: name_idx, clean_idx
    real(kind=rp), allocatable :: coords(:,:)
    integer, allocatable :: connectivity(:), offsets(:)
    integer, allocatable :: cell_types(:)
    integer(kind=1), allocatable :: cell_types_byte(:)
    integer, dimension(8), parameter :: vcyc_to_sym = [1, 2, 4, 3, 5, 6, 8, 7]
    integer, dimension(8) :: id
    real(kind=rp), allocatable :: point_data(:,:)
    integer :: npts_per_cell, nodes_per_cell, subcells_per_el
    integer :: ii, jj, kk, local_idx, conn_idx
    integer :: base

    ! Determine mesh and field data
    select type(data)
    type is (mesh_t)
       msh => data
       fld => null()
       n_fields = 0
    type is (field_t)
       msh => data%msh
       fld => data
       dof => data%dof
       n_fields = 1
       allocate(fp(1))
       allocate(field_written(1), source=.false.)
       fp(1)%ptr => data
    type is (field_list_t)
       msh => data%msh(1)
       dof => data%dof(1)
       fld => null()
       n_fields = data%size()
       allocate(fp(n_fields))
       allocate(field_written(n_fields), source=.false.)
       do i = 1, n_fields
          fp(i)%ptr => data%items(i)%ptr
       end do
    type is (dofmap_t)
       dof => data
       msh => data%msh
       fld => null()
       n_fields = 0
    class default
       call neko_error('Invalid data type for vtkhdf_file_write')
    end select

    ! Syn all the fields
    do i = 1, n_fields
       if (associated(fp(i)%ptr)) then
          call fp(i)%ptr%copy_from(DEVICE_TO_HOST, sync = i .eq. n_fields)
       end if
    end do

    if (.not. this%overwrite) call this%increment_counter()
    fname = trim(this%get_fname())

    call h5open_f(ierr)
    call vtkhdf_file_determine_real(H5T_NEKO_REAL)

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(plist_id, NEKO_COMM%mpi_val, info, ierr)

    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, &
         file_id, ierr, access_prp = plist_id)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Create VTKHDF root group with vtkhdf_version and type attributes
    call h5gcreate_f(file_id, "VTKHDF", vtkhdf_grp, ierr, lcpl_id=h5p_default_f, &
         gcpl_id=h5p_default_f, gapl_id=h5p_default_f)

    ! Write Version attribute [2, 5] as an array of 2 64bit integers
    vdims(1) = 2
    call h5screate_simple_f(1, vdims(1:1), filespace, ierr)
    call h5acreate_f(vtkhdf_grp, "Version", H5T_NATIVE_INTEGER, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, vtkhdf_version, vdims(1:1), ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(filespace, ierr)

    ! Write Type attribute "UnstructuredGrid" as a fixed-length string, with
    ! NULL padding to the full length of the string (16 chars)
    type_str(1) = "UnstructuredGrid"
    vdims(1) = 1
    call h5screate_f(H5S_SCALAR_F, filespace, ierr)

    call h5tcopy_f(H5T_FORTRAN_S1, memspace, ierr)
    call h5tset_size_f(memspace, int(len_trim(type_str(1)), kind=8), ierr)
    call h5tset_strpad_f(memspace, H5T_STR_NULLTERM_F, ierr)

    call h5acreate_f(vtkhdf_grp, "Type", memspace, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, memspace, type_str, vdims(1:1), ierr)
    call h5aclose_f(attr_id, ierr)

    call h5tclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)

    ! Write mesh information if present
    if (associated(msh)) then
       call MPI_Comm_size(NEKO_COMM, num_partitions, ierr)
       local_cells = msh%nelv
       local_points = msh%mpts
       local_conn = msh%npts * msh%nelv

       if (dof%Xh%lx < 2 .or. dof%Xh%ly < 2) then
          call neko_error('VTKHDF linear output requires lx, ly >= 2')
       end if
       if (msh%gdim .eq. 3 .and. dof%Xh%lz < 2) then
          call neko_error('VTKHDF linear output requires lz >= 2 in 3D')
       end if

       npts_per_cell = dof%Xh%lx * dof%Xh%ly * dof%Xh%lz
       if (msh%gdim .eq. 3) then
          nodes_per_cell = 8
          subcells_per_el = (dof%Xh%lx - 1) * (dof%Xh%ly - 1) * (dof%Xh%lz - 1)
       else
          nodes_per_cell = 4
          subcells_per_el = (dof%Xh%lx - 1) * (dof%Xh%ly - 1)
       end if

       local_points = msh%nelv * npts_per_cell
       local_cells = msh%nelv * subcells_per_el
       local_conn = local_cells * nodes_per_cell

       allocate(part_points(num_partitions))
       allocate(part_cells(num_partitions))
       allocate(part_conns(num_partitions))

       call MPI_Allgather(local_points, 1, MPI_INTEGER8, part_points, 1, MPI_INTEGER8, NEKO_COMM, ierr)
       call MPI_Allgather(local_cells, 1, MPI_INTEGER8, part_cells, 1, MPI_INTEGER8, NEKO_COMM, ierr)
       call MPI_Allgather(local_conn, 1, MPI_INTEGER8, part_conns, 1, MPI_INTEGER8, NEKO_COMM, ierr)

       total_points = sum(part_points)
       total_cells = sum(part_cells)
       total_conn = sum(part_conns)
       total_offsets = total_cells + num_partitions

       point_offset = 0
       cell_offset = 0
       conn_offset = 0
       offsets_offset = 0
       do i = 1, num_partitions
          if (i >= pe_rank + 1) exit
          point_offset = point_offset + part_points(i)
          cell_offset = cell_offset + part_cells(i)
          conn_offset = conn_offset + part_conns(i)
          offsets_offset = offsets_offset + part_cells(i) + 1
       end do

       ! Write NumberOfPoints dataset
       vdims(1) = 1_hsize_t
       maxdims(1) = H5S_UNLIMITED_F
       call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       chunkdims(1) = 1_hsize_t
       call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
       call h5dcreate_f(vtkhdf_grp, "NumberOfPoints", H5T_NATIVE_INTEGER, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       if (pe_rank .eq. 0) then
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, total_points, vdims(1:1), ierr)
       end if
       call h5dclose_f(dset_id, ierr)
       call h5pclose_f(dcpl_id, ierr)
       call h5sclose_f(filespace, ierr)

       ! Write NumberOfCells dataset
       vdims(1) = 1_hsize_t
       maxdims(1) = H5S_UNLIMITED_F
       call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       chunkdims(1) = 1_hsize_t
       call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
       call h5dcreate_f(vtkhdf_grp, "NumberOfCells", H5T_NATIVE_INTEGER, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       if (pe_rank .eq. 0) then
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, total_cells, vdims(1:1), ierr)
       end if
       call h5dclose_f(dset_id, ierr)
       call h5pclose_f(dcpl_id, ierr)
       call h5sclose_f(filespace, ierr)

       ! Write NumberOfConnectivityIds dataset
       vdims(1) = 1_hsize_t
       maxdims(1) = H5S_UNLIMITED_F
       call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       chunkdims(1) = 1_hsize_t
       call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
       call h5dcreate_f(vtkhdf_grp, "NumberOfConnectivityIds", H5T_NATIVE_INTEGER, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       if (pe_rank .eq. 0) then
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, total_conn, vdims(1:1), ierr)
       end if
       call h5dclose_f(dset_id, ierr)
       call h5pclose_f(dcpl_id, ierr)
       call h5sclose_f(filespace, ierr)

       ! Write Points dataset (global coordinates)
       allocate(coords(3, local_points))
       ! Write points directly from dofmap in (i,j,k) tensor-product order
       do i = 1, msh%nelv
          local_idx = (i - 1) * npts_per_cell
          do kk = 1, dof%Xh%lz
             do jj = 1, dof%Xh%ly
                do ii = 1, dof%Xh%lx
                   local_idx = local_idx + 1
                   coords(1, local_idx) = dof%x(ii, jj, kk, i)
                   coords(2, local_idx) = dof%y(ii, jj, kk, i)
                   coords(3, local_idx) = dof%z(ii, jj, kk, i)
                end do
             end do
          end do
       end do

       vdims = [3_hsize_t, int(total_points, hsize_t)]
       maxdims = [3_hsize_t, H5S_UNLIMITED_F]
       call h5screate_simple_f(2, vdims, filespace, ierr, maxdims)
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       chunkdims(1) = max(1_hsize_t, min(int(local_points, hsize_t), vdims(2)))
       call h5pset_chunk_f(dcpl_id, 2, [3_hsize_t, chunkdims(1)], ierr)
       call h5dcreate_f(vtkhdf_grp, "Points", H5T_NEKO_REAL, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5dget_space_f(dset_id, filespace, ierr)
       dcount2 = [3_hsize_t, int(local_points, hsize_t)]
       doffset2 = [0_hsize_t, int(point_offset, hsize_t)]
       call h5screate_simple_f(2, dcount2, memspace, ierr)
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
            doffset2, dcount2, ierr)
       call h5dwrite_f(dset_id, H5T_NEKO_REAL, coords, dcount2, ierr, &
            file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
       call h5sclose_f(memspace, ierr)
       call h5dclose_f(dset_id, ierr)
       call h5pclose_f(dcpl_id, ierr)
       call h5sclose_f(filespace, ierr)
       deallocate(coords)

       ! Write Connectivity dataset (linear cells from GLL point subdivision)
       allocate(connectivity(local_conn))
       conn_idx = 0
       do i = 1, msh%nelv
          base = (i - 1) * npts_per_cell
          if (msh%gdim .eq. 3) then
             do ii = 1, dof%Xh%lx - 1
                do jj = 1, dof%Xh%ly - 1
                   do kk = 1, dof%Xh%lz - 1
                      connectivity(conn_idx + 1) = point_offset + base + &
                           (kk - 1) * dof%Xh%lx * dof%Xh%ly + &
                           (jj - 1) * dof%Xh%lx + ii - 1
                      connectivity(conn_idx + 2) = point_offset + base + &
                           (kk - 1) * dof%Xh%lx * dof%Xh%ly + &
                           (jj - 1) * dof%Xh%lx + (ii + 1) - 1
                      connectivity(conn_idx + 3) = point_offset + base + &
                           (kk - 1) * dof%Xh%lx * dof%Xh%ly + &
                           (jj + 1 - 1) * dof%Xh%lx + (ii + 1) - 1
                      connectivity(conn_idx + 4) = point_offset + base + &
                           (kk - 1) * dof%Xh%lx * dof%Xh%ly + &
                           (jj + 1 - 1) * dof%Xh%lx + ii - 1
                      connectivity(conn_idx + 5) = point_offset + base + &
                           (kk + 1 - 1) * dof%Xh%lx * dof%Xh%ly + &
                           (jj - 1) * dof%Xh%lx + ii - 1
                      connectivity(conn_idx + 6) = point_offset + base + &
                           (kk + 1 - 1) * dof%Xh%lx * dof%Xh%ly + &
                           (jj - 1) * dof%Xh%lx + (ii + 1) - 1
                      connectivity(conn_idx + 7) = point_offset + base + &
                           (kk + 1 - 1) * dof%Xh%lx * dof%Xh%ly + &
                           (jj + 1 - 1) * dof%Xh%lx + (ii + 1) - 1
                      connectivity(conn_idx + 8) = point_offset + base + &
                           (kk + 1 - 1) * dof%Xh%lx * dof%Xh%ly + &
                           (jj + 1 - 1) * dof%Xh%lx + ii - 1
                      conn_idx = conn_idx + 8
                   end do
                end do
             end do
          else
             do jj = 1, dof%Xh%ly - 1
                do ii = 1, dof%Xh%lx - 1
                   connectivity(conn_idx + 1) = point_offset + base + &
                        (jj - 1) * dof%Xh%lx + ii - 1
                   connectivity(conn_idx + 2) = point_offset + base + &
                        (jj - 1) * dof%Xh%lx + (ii + 1) - 1
                   connectivity(conn_idx + 3) = point_offset + base + &
                        (jj + 1 - 1) * dof%Xh%lx + (ii + 1) - 1
                   connectivity(conn_idx + 4) = point_offset + base + &
                        (jj + 1 - 1) * dof%Xh%lx + ii - 1
                   conn_idx = conn_idx + 4
                end do
             end do
          end if
       end do

       vdims(1) = int(total_conn, hsize_t)
       maxdims(1) = H5S_UNLIMITED_F
       call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       chunkdims(1) = max(1_hsize_t, min(int(local_conn, hsize_t), vdims(1)))
       call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
       call h5dcreate_f(vtkhdf_grp, "Connectivity", H5T_NATIVE_INTEGER, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5dget_space_f(dset_id, filespace, ierr)
       dcount(1) = int(local_conn, hsize_t)
       doffset(1) = int(conn_offset, hsize_t)
       call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, connectivity, dcount(1:1), ierr, &
            file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
       call h5sclose_f(memspace, ierr)
       call h5dclose_f(dset_id, ierr)
       call h5pclose_f(dcpl_id, ierr)
       call h5sclose_f(filespace, ierr)
       deallocate(connectivity)

       ! Write Offsets dataset (cell offsets into connectivity)
       allocate(offsets(local_cells + 1))
       do i = 1, local_cells
          offsets(i) = (i - 1) * nodes_per_cell + conn_offset
       end do

       offsets(local_cells + 1) = local_conn + conn_offset

       vdims(1) = int(total_offsets, hsize_t)
       maxdims(1) = H5S_UNLIMITED_F
       call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       chunkdims(1) = max(1_hsize_t, min(int(local_cells + 1, hsize_t), vdims(1)))
       call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
       call h5dcreate_f(vtkhdf_grp, "Offsets", H5T_NATIVE_INTEGER, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5dget_space_f(dset_id, filespace, ierr)
       dcount(1) = int(local_cells + 1, hsize_t)
       doffset(1) = int(offsets_offset, hsize_t)
       call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, offsets, dcount(1:1), ierr, &
            file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
       call h5sclose_f(memspace, ierr)
       call h5dclose_f(dset_id, ierr)
       call h5pclose_f(dcpl_id, ierr)
       call h5sclose_f(filespace, ierr)
       deallocate(offsets)

       ! Write Types dataset (VTK cell types)
       allocate(cell_types(local_cells))
       allocate(cell_types_byte(local_cells))
       if (msh%gdim .eq. 3) then
          cell_types = 12 ! VTK_HEXAHEDRON
       else
          cell_types = 9 ! VTK_QUAD
       end if
       ! Convert integer array to byte array
       cell_types_byte = int(cell_types, kind=1)

       vdims(1) = int(total_cells, hsize_t)
       maxdims(1) = H5S_UNLIMITED_F
       call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       chunkdims(1) = max(1_hsize_t, min(int(local_cells, hsize_t), vdims(1)))
       call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
       call h5dcreate_f(vtkhdf_grp, "Types", H5T_STD_U8LE, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5dget_space_f(dset_id, filespace, ierr)
       dcount(1) = int(local_cells, hsize_t)
       doffset(1) = int(cell_offset, hsize_t)
       call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dwrite_f(dset_id, H5T_STD_U8LE, cell_types_byte, dcount(1:1), ierr, &
            file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
       call h5sclose_f(memspace, ierr)
       call h5dclose_f(dset_id, ierr)
       call h5pclose_f(dcpl_id, ierr)
       call h5sclose_f(filespace, ierr)
       deallocate(cell_types, cell_types_byte)
       deallocate(part_points, part_cells, part_conns)
    end if

    ! Write field data in PointData group
    if (associated(msh) .and. n_fields > 0) then
       call h5gcreate_f(vtkhdf_grp, "PointData", pointdata_grp, ierr, lcpl_id=h5p_default_f, &
            gcpl_id=h5p_default_f, gapl_id=h5p_default_f)

       ! Get field dimensions to create common filespace and memspace
       fld => fp(1)%ptr
       lx = fld%Xh%lx
       ly = fld%Xh%ly
       lz = fld%Xh%lz
       nelv = fld%msh%nelv
       mpts = fld%msh%mpts

       dcount(1) = int(local_points, hsize_t)
       doffset(1) = int(point_offset, hsize_t)
       ddim(1) = int(total_points, hsize_t)
       drank = 1

       ! Loop through all fields and write them
       do i = 1, n_fields
          if (field_written(i)) cycle
          fld => fp(i)%ptr

          ! Remove underscores from field name
          field_name = fld%name

          if (field_name .eq. 'p') field_name = 'Pressure'

          if (field_name .eq. 'u' .or. &
               field_name .eq. 'v' .or. &
               field_name .eq. 'w') then
             u => null()
             v => null()
             w => null()
             do j = 1, n_fields
                select case (fp(j)%ptr%name)
                case ('u')
                   u => fp(j)%ptr
                   field_written(j) = .true.
                case ('v')
                   v => fp(j)%ptr
                   field_written(j) = .true.
                case ('w')
                   w => fp(j)%ptr
                   field_written(j) = .true.
                end select
             end do
             field_name = 'Velocity'

             ! Allocate temporary array to hold velocity components (3, local_points)
             allocate(point_data(3, local_points))

             ! Pack u, v, w components into point_data array using same indexing as Points
             local_idx = 0
             do ie = 1, msh%nelv
                local_idx = (ie - 1) * npts_per_cell
                do kk = 1, dof%Xh%lz
                   do jj = 1, dof%Xh%ly
                      do ii = 1, dof%Xh%lx
                         local_idx = local_idx + 1
                         point_data(1, local_idx) = u%x(ii, jj, kk, ie) ! u component
                         point_data(2, local_idx) = v%x(ii, jj, kk, ie) ! v component
                         if (associated(w)) then
                            point_data(3, local_idx) = w%x(ii, jj, kk, ie) ! w component
                         else
                            point_data(3, local_idx) = 0.0_rp ! zero for 2D
                         end if
                      end do
                   end do
                end do
             end do

             ! Create filespace with unlimited dimension for 2D array (3, total_points)
             vdims = [3_hsize_t, int(total_points, hsize_t)]
             maxdims = [3_hsize_t, H5S_UNLIMITED_F]
             call h5screate_simple_f(2, vdims, filespace, ierr, maxdims)
             call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
             chunkdims(1) = max(1_hsize_t, min(int(local_points, hsize_t), vdims(2)))
             call h5pset_chunk_f(dcpl_id, 2, [3_hsize_t, chunkdims(1)], ierr)
             call h5dcreate_f(pointdata_grp, trim(field_name), H5T_NEKO_REAL, &
                  filespace, dset_id, ierr, dcpl_id = dcpl_id)
             call h5dget_space_f(dset_id, filespace, ierr)
             dcount2 = [3_hsize_t, int(local_points, hsize_t)]
             doffset2 = [0_hsize_t, int(point_offset, hsize_t)]
             call h5screate_simple_f(2, dcount2, memspace, ierr)
             call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                  doffset2, dcount2, ierr)
             call h5dwrite_f(dset_id, H5T_NEKO_REAL, point_data, dcount2, ierr, &
                  file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
             call h5sclose_f(filespace, ierr)
             call h5sclose_f(memspace, ierr)
             call h5dclose_f(dset_id, ierr)
             call h5pclose_f(dcpl_id, ierr)
             deallocate(point_data)

          else

             ! Create filespace with unlimited dimension
             vdims(1) = int(total_points, hsize_t)
             maxdims(1) = H5S_UNLIMITED_F
             call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
             call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
             chunkdims(1) = max(1_hsize_t, min(int(local_points, hsize_t), vdims(1)))
             call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
             call h5dcreate_f(pointdata_grp, trim(field_name), H5T_NEKO_REAL, &
                  filespace, dset_id, ierr, dcpl_id = dcpl_id)
             call h5dget_space_f(dset_id, filespace, ierr)
             call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
             call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                  doffset(1:1), dcount(1:1), ierr)
             call h5dwrite_f(dset_id, H5T_NEKO_REAL, &
                  fld%x(1,1,1,1), &
                  dcount(1:1), ierr, file_space_id = filespace, &
                  mem_space_id = memspace, xfer_prp = plist_id)
             call h5sclose_f(filespace, ierr)
             call h5sclose_f(memspace, ierr)
             call h5dclose_f(dset_id, ierr)
             call h5pclose_f(dcpl_id, ierr)
          end if
       end do

       call h5gclose_f(pointdata_grp, ierr)
    end if

    ! Write Time attribute if present
    if (present(t)) then
       ddim(1) = 1
       call h5screate_simple_f(1, ddim(1:1), filespace, ierr)
       call h5acreate_f(vtkhdf_grp, "Time", H5T_NEKO_REAL, filespace, attr_id, &
            ierr, h5p_default_f, h5p_default_f)
       call h5awrite_f(attr_id, H5T_NEKO_REAL, t, ddim(1:1), ierr)
       call h5aclose_f(attr_id, ierr)
       call h5sclose_f(filespace, ierr)
    end if

    call h5gclose_f(vtkhdf_grp, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)

    if (allocated(fp)) deallocate(fp)
    if (allocated(field_written)) deallocate(field_written)

  end subroutine vtkhdf_file_write

  !> Read data in HDF5 format following official VTKHDF specification
  subroutine vtkhdf_file_read(this, data)
    class(vtkhdf_file_t) :: this
    class(*), target, intent(inout) :: data
    integer(hid_t) :: plist_id, file_id, dset_id, grp_id, attr_id, vtkhdf_grp
    integer(hid_t) :: pointdata_grp
    integer(hid_t) :: filespace, memspace
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: i, j, ierr, info, drank, n_fields
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dof
    type(field_t), pointer :: fld
    type(field_ptr_t), allocatable :: fp(:)
    real(kind=rp) :: t
    character(len=1024) :: fname
    character(len=16), dimension(1) :: type_str
    character(len=128) :: field_name
    integer :: name_idx, clean_idx
    integer :: vtkhdf_version(2)
    integer :: num_partitions, npts_per_cell
    integer :: local_points, total_points, point_offset, part_index
    integer, allocatable :: part_points(:)

    select type(data)
    type is (mesh_t)
       msh => data
       fld => null()
       n_fields = 0
    type is (field_t)
       msh => data%msh
       fld => data
       dof => data%dof
       n_fields = 1
       allocate(fp(1))
       fp(1)%ptr => data
    type is (field_list_t)
       msh => data%msh(1)
       dof => data%dof(1)
       fld => null()
       n_fields = data%size()
       allocate(fp(n_fields))
       do i = 1, n_fields
          fp(i)%ptr => data%items(i)%ptr
       end do
    type is (dofmap_t)
       dof => data
       msh => data%msh
       fld => null()
       n_fields = 0
    class default
       call neko_error('Invalid data type for vtkhdf_file_read')
    end select

    fname = trim(this%get_fname())

    call h5open_f(ierr)
    call vtkhdf_file_determine_real(H5T_NEKO_REAL)

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(plist_id, NEKO_COMM%mpi_val, info, ierr)

    call h5fopen_f(fname, H5F_ACC_RDONLY_F, &
         file_id, ierr, access_prp = plist_id)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Open VTKHDF root group
    call h5gopen_f(file_id, 'VTKHDF', vtkhdf_grp, ierr, gapl_id=h5p_default_f)

    ! Read Version attribute
    ddim(1) = 2
    call h5aopen_name_f(vtkhdf_grp, 'Version', attr_id, ierr)
    if (ierr .eq. 0) then
       call h5aread_f(attr_id, H5T_NATIVE_INTEGER, vtkhdf_version, ddim, ierr)
       call h5aclose_f(attr_id, ierr)
    end if

    ! Read Type attribute
    ddim(1) = 1
    call h5aopen_name_f(vtkhdf_grp, 'Type', attr_id, ierr)
    if (ierr .eq. 0) then
       call h5aread_f(attr_id, H5T_FORTRAN_S1, type_str, ddim, ierr)
       call h5aclose_f(attr_id, ierr)
    end if

    ! Read time attribute if present
    ddim(1) = 1
    call h5aopen_name_f(vtkhdf_grp, 'Time', attr_id, ierr)
    if (ierr .eq. 0) then
       call h5aread_f(attr_id, H5T_NEKO_REAL, t, ddim, ierr)
       call h5aclose_f(attr_id, ierr)
    end if

    ! Read field data from PointData group if present
    if (n_fields > 0 .and. associated(dof)) then
       call h5gopen_f(vtkhdf_grp, 'PointData', pointdata_grp, ierr, gapl_id=h5p_default_f)

       if (ierr .eq. 0) then
          ! Calculate point offsets and counts (same as writer)
          call MPI_Comm_size(NEKO_COMM, num_partitions, ierr)
          npts_per_cell = dof%Xh%lx * dof%Xh%ly * dof%Xh%lz
          local_points = msh%nelv * npts_per_cell

          allocate(part_points(num_partitions))
          call MPI_Allgather(local_points, 1, MPI_INTEGER, part_points, &
               1, MPI_INTEGER, NEKO_COMM, ierr)

          total_points = sum(part_points)
          part_index = pe_rank + 1
          point_offset = 0
          do j = 1, part_index - 1
             point_offset = point_offset + part_points(j)
          end do

          ! Read each field - data is in Neko's natural (i,j,k,e) order
          do i = 1, n_fields
             fld => fp(i)%ptr
             drank = 1
             dcount(1) = int(local_points, hsize_t)
             doffset(1) = int(point_offset, hsize_t)
             ddim(1) = int(total_points, hsize_t)

             ! Remove underscores from field name
             field_name = ''
             clean_idx = 1
             do name_idx = 1, len_trim(fld%name)
                if (fld%name(name_idx:name_idx) /= '_') then
                   field_name(clean_idx:clean_idx) = fld%name(name_idx:name_idx)
                   clean_idx = clean_idx + 1
                end if
             end do

             call h5dopen_f(pointdata_grp, trim(field_name), dset_id, ierr)
             if (ierr .eq. 0) then
                call h5dget_space_f(dset_id, filespace, ierr)
                call h5screate_simple_f(drank, dcount, memspace, ierr)
                call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                     doffset, dcount, ierr)

                ! Read directly into field array (no reordering needed)
                call h5dread_f(dset_id, H5T_NEKO_REAL, fld%x(1,1,1,1), dcount, &
                     ierr, file_space_id = filespace, mem_space_id = memspace, &
                     xfer_prp = plist_id)

                call h5dclose_f(dset_id, ierr)
                call h5sclose_f(filespace, ierr)
                call h5sclose_f(memspace, ierr)
             end if
          end do

          deallocate(part_points)
          call h5gclose_f(pointdata_grp, ierr)
       end if
    end if

    call h5gclose_f(vtkhdf_grp, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)

    if (allocated(fp)) deallocate(fp)

  end subroutine vtkhdf_file_read



  !> Determine hdf5 real type corresponding to NEKO_REAL
  !! @note This must be called after h5open_f, otherwise
  !! the H5T_NATIVE_XYZ types has a value of 0
  subroutine vtkhdf_file_determine_real(H5T_NEKO_REAL)
    integer(hid_t), intent(inout) :: H5T_NEKO_REAL
    select case(rp)
    case(dp)
       H5T_NEKO_REAL = H5T_NATIVE_DOUBLE
    case(sp)
       H5T_NEKO_REAL = H5T_NATIVE_REAL
    case default
       call neko_error("Unsupported real type")
    end select
  end subroutine vtkhdf_file_determine_real

  !> Set the overwrite flag for HDF5 files
  subroutine vtkhdf_file_set_overwrite(this, overwrite)
    class(vtkhdf_file_t), intent(inout) :: this
    logical, intent(in) :: overwrite
    this%overwrite = overwrite
  end subroutine vtkhdf_file_set_overwrite

#else

  !> Write data in HDF5 format (no HDF5 support)
  subroutine vtkhdf_file_write(this, data, t)
    class(vtkhdf_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine vtkhdf_file_write

  !> Read data in HDF5 format (no HDF5 support)
  subroutine vtkhdf_file_read(this, data)
    class(vtkhdf_file_t) :: this
    class(*), target, intent(inout) :: data
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine vtkhdf_file_read

  !> Set the overwrite flag for HDF5 files (no HDF5 support)
  subroutine vtkhdf_file_set_overwrite(this, overwrite)
    class(vtkhdf_file_t), intent(inout) :: this
    logical, intent(in) :: overwrite
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine vtkhdf_file_set_overwrite

#endif

end module vtkhdf_file
