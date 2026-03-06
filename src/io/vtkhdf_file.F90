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
       MPI_INTEGER, MPI_SUM, MPI_Comm_size
#ifdef HAVE_HDF5
  use hdf5
#endif
  implicit none
  private

  !> Interface for HDF5 files
  type, public, extends(generic_file_t) :: vtkhdf_file_t
     logical :: amr_enabled = .false.
   contains
     procedure :: read => vtkhdf_file_read
     procedure :: write => vtkhdf_file_write
     procedure :: set_overwrite => vtkhdf_file_set_overwrite
     procedure :: enable_amr => vtkhdf_file_enable_amr
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
    type(field_ptr_t), allocatable :: fp(:)
    logical, allocatable :: field_written(:)
    integer :: ierr, info, i, n_fields
    integer(hid_t) :: plist_id, plist_coll, file_id, attr_id, vtkhdf_grp
    integer(hid_t) :: filespace, memspace
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(2) :: vdims
    integer :: lx, ly, lz
    integer :: num_partitions
    integer :: local_points, local_cells, local_conn
    integer :: total_points, total_cells, total_conn
    integer :: point_offset
    integer :: max_local_points
    integer, allocatable :: part_points(:), part_cells(:), part_conns(:)
    character(len=1024) :: fname
    character(len=16), dimension(1) :: type_str
    integer(hsize_t) :: pointdata_time_offset
    logical :: link_exists, file_exists

    ! Determine mesh and field data
    select type(data)
    type is (field_t)
       msh => data%msh
       dof => data%dof
       n_fields = 1
       allocate(fp(1))
       allocate(field_written(1), source=.false.)
       fp(1)%ptr => data
    type is (field_list_t)
       msh => data%msh(1)
       dof => data%dof(1)
       n_fields = data%size()
       allocate(fp(n_fields))
       allocate(field_written(n_fields), source=.false.)
       do i = 1, n_fields
          fp(i)%ptr => data%items(i)%ptr
       end do
    class default
       call neko_error('Invalid data type for vtkhdf_file_write')
    end select

    ! Check conditions to ensure the input data is supported.
    if (.not. associated(msh)) then
       call neko_error('Mesh must be associated for vtkhdf_file_write')
    end if
    if (dof%Xh%lx < 2 .or. dof%Xh%ly < 2) then
       call neko_error('VTKHDF linear output requires lx, ly >= 2')
    end if
    if (msh%gdim .eq. 3 .and. dof%Xh%lz < 2) then
       call neko_error('VTKHDF linear output requires lz >= 2 in 3D')
    end if

    ! Assign commonly used values
    lx = dof%Xh%lx
    ly = dof%Xh%ly
    lz = dof%Xh%lz

    call this%increment_counter()
    fname = trim(this%get_base_fname())

    call h5open_f(ierr)
    call vtkhdf_file_determine_real(H5T_NEKO_REAL)

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(plist_id, NEKO_COMM%mpi_val, info, ierr)

    inquire(file = fname, exist = file_exists)
    if (file_exists) then
       call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, ierr, access_prp = plist_id)
    else
       call h5fcreate_f(fname, H5F_ACC_TRUNC_F, &
            file_id, ierr, access_prp = plist_id)
    end if

    ! Create collective transfer property list
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_coll, ierr)
    call h5pset_dxpl_mpio_f(plist_coll, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Create/open VTKHDF root group with vtkhdf_version and type attributes
    call h5lexists_f(file_id, "VTKHDF", link_exists, ierr)
    if (link_exists) then
       call h5gopen_f(file_id, "VTKHDF", vtkhdf_grp, ierr)
    else
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

       ! Write Type attribute "UnstructuredGrid" as a fixed-length string
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
    end if

    ! Write mesh information if present

    call MPI_Comm_size(NEKO_COMM, num_partitions, ierr)

    if (dof%Xh%lx < 2 .or. dof%Xh%ly < 2) then
       call neko_error('VTKHDF linear output requires lx, ly >= 2')
    end if
    if (msh%gdim .eq. 3 .and. dof%Xh%lz < 2) then
       call neko_error('VTKHDF linear output requires lz >= 2 in 3D')
    end if

    local_points = msh%nelv * lx * ly * lz
    if (msh%gdim .eq. 3) then
       local_cells = msh%nelv * (lx - 1) * (ly - 1) * (lz - 1)
       local_conn = local_cells * 8
    else
       local_cells = msh%nelv * (lx - 1) * (ly - 1)
       local_conn = local_cells * 4
    end if

    allocate(part_points(num_partitions))
    allocate(part_cells(num_partitions))
    allocate(part_conns(num_partitions))

    call MPI_Allgather(local_points, 1, MPI_INTEGER, part_points, 1, MPI_INTEGER, NEKO_COMM, ierr)
    call MPI_Allgather(local_cells, 1, MPI_INTEGER, part_cells, 1, MPI_INTEGER, NEKO_COMM, ierr)
    call MPI_Allgather(local_conn, 1, MPI_INTEGER, part_conns, 1, MPI_INTEGER, NEKO_COMM, ierr)

    total_points = sum(part_points)
    total_cells = sum(part_cells)
    total_conn = sum(part_conns)

    max_local_points = maxval(part_points)

    point_offset = sum(part_points(1:pe_rank))

    ! For static mesh, only write geometry on the first call
    call h5lexists_f(vtkhdf_grp, "Points", link_exists, ierr)
    if (this%amr_enabled .or. .not. link_exists) then
       call vtkhdf_write_mesh(vtkhdf_grp, dof, msh, plist_coll, &
            H5T_NEKO_REAL, num_partitions, &
            local_points, total_points, total_cells, total_conn, &
            point_offset, max_local_points, part_cells, part_conns)
    end if ! write mesh conditional
    deallocate(part_points, part_cells, part_conns)

    pointdata_time_offset = 0_hsize_t

    if (present(t)) then
       call vtkhdf_write_steps(vtkhdf_grp, plist_coll, H5T_NEKO_REAL, t, &
            this%amr_enabled, num_partitions, &
            total_points, total_cells, total_conn, &
            pointdata_time_offset)
    end if

    ! Write field data in PointData group
    if (n_fields > 0) then
       call vtkhdf_write_pointdata(vtkhdf_grp, plist_coll, H5T_NEKO_REAL, &
            fp, field_written, msh%nelv, &
            local_points, point_offset, max_local_points, total_points, &
            pointdata_time_offset, lx, ly, lz, present(t))
    end if

    call h5gclose_f(vtkhdf_grp, ierr)
    call h5pclose_f(plist_coll, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)

    if (allocated(fp)) deallocate(fp)
    if (allocated(field_written)) deallocate(field_written)

  end subroutine vtkhdf_file_write

  !> Build local connectivity for VTK sub-cells from a spectral element
  !! tensor-product grid. Subdivides each spectral element into linear
  !! sub-cells based on the GLL node positions.
  !! @param conn Output connectivity array (pre-allocated)
  !! @param vtk_type VTK cell type: 12 = VTK_HEXAHEDRON, 9 = VTK_QUAD
  !! @param lx Number of GLL points in x-direction
  !! @param ly Number of GLL points in y-direction
  !! @param lz Number of GLL points in z-direction
  !! @param nelv Number of elements
  subroutine vtkhdf_build_connectivity(conn, vtk_type, lx, ly, lz, nelv)
    integer, intent(out) :: conn(:)
    integer, intent(in) :: vtk_type, lx, ly, lz, nelv
    integer :: ie, ii, jj, kk, base, idx, npts_per_cell

    npts_per_cell = lx * ly * lz
    idx = 0

    do ie = 1, nelv
       base = (ie - 1) * npts_per_cell

       select case (vtk_type)
       case (12) ! VTK_HEXAHEDRON
          do ii = 1, lx - 1
             do jj = 1, ly - 1
                do kk = 1, lz - 1
                   conn(idx + 1) = base + &
                        (kk - 1) * lx * ly + (jj - 1) * lx + ii - 1
                   conn(idx + 2) = base + &
                        (kk - 1) * lx * ly + (jj - 1) * lx + (ii + 1) - 1
                   conn(idx + 3) = base + &
                        (kk - 1) * lx * ly + (jj + 1 - 1) * lx + (ii + 1) - 1
                   conn(idx + 4) = base + &
                        (kk - 1) * lx * ly + (jj + 1 - 1) * lx + ii - 1
                   conn(idx + 5) = base + &
                        (kk + 1 - 1) * lx * ly + (jj - 1) * lx + ii - 1
                   conn(idx + 6) = base + &
                        (kk + 1 - 1) * lx * ly + (jj - 1) * lx + (ii + 1) - 1
                   conn(idx + 7) = base + &
                        (kk + 1 - 1) * lx * ly + (jj + 1 - 1) * lx + (ii + 1) - 1
                   conn(idx + 8) = base + &
                        (kk + 1 - 1) * lx * ly + (jj + 1 - 1) * lx + ii - 1
                   idx = idx + 8
                end do
             end do
          end do

       case (9) ! VTK_QUAD
          do jj = 1, ly - 1
             do ii = 1, lx - 1
                conn(idx + 1) = base + (jj - 1) * lx + ii - 1
                conn(idx + 2) = base + (jj - 1) * lx + (ii + 1) - 1
                conn(idx + 3) = base + (jj + 1 - 1) * lx + (ii + 1) - 1
                conn(idx + 4) = base + (jj + 1 - 1) * lx + ii - 1
                idx = idx + 4
             end do
          end do

       case default
          call neko_error('Unsupported VTK cell type in vtkhdf_build_connectivity')
       end select
    end do

  end subroutine vtkhdf_build_connectivity

  !> Write mesh geometry datasets to the VTKHDF group.
  !! Writes NumberOfPoints, NumberOfCells, NumberOfConnectivityIds,
  !! Points, Connectivity, Offsets, and Types datasets.
  !! @param vtkhdf_grp HDF5 group ID for VTKHDF root group
  !! @param dof Dofmap for coordinate data
  !! @param msh Mesh object
  !! @param plist_coll Collective transfer property list
  !! @param H5T_NEKO_REAL HDF5 real type matching NEKO precision
  !! @param num_partitions Number of MPI partitions
  !! @param local_points Number of points on this rank
  !! @param local_cells Number of sub-cells on this rank
  !! @param local_conn Connectivity array size on this rank
  !! @param total_points Global total number of points
  !! @param total_cells Global total number of sub-cells
  !! @param total_conn Global total connectivity size
  !! @param point_offset This rank's offset into the Points dataset
  !! @param max_local_points Maximum local_points across all ranks
  !! @param part_cells Per-partition cell counts from MPI_Allgather
  !! @param part_conns Per-partition connectivity sizes from MPI_Allgather
  subroutine vtkhdf_write_mesh(vtkhdf_grp, dof, msh, plist_coll, &
       H5T_NEKO_REAL, num_partitions, &
       local_points, total_points, total_cells, total_conn, &
       point_offset, max_local_points, part_cells, part_conns)
    type(dofmap_t), intent(in) :: dof
    type(mesh_t), intent(in) :: msh
    integer(hid_t), intent(in) :: vtkhdf_grp, plist_coll, H5T_NEKO_REAL
    integer, intent(in) :: num_partitions
    integer, intent(in) :: local_points
    integer, intent(in) :: total_points, total_cells, total_conn
    integer, intent(in) :: point_offset, max_local_points
    integer, intent(in) :: part_cells(:), part_conns(:)

    integer :: ierr, i, ii, jj, kk, local_idx
    integer :: lx, ly, lz, npts_per_cell, nodes_per_cell
    integer :: local_cells, local_conn
    integer :: total_offsets, cell_offset, conn_offset, offsets_offset
    integer :: max_local_cells, max_local_conn
    integer(hid_t) :: dset_id, dcpl_id, filespace, memspace
    integer(hsize_t), dimension(1) :: dcount, doffset_1d, chunkdims
    integer(hsize_t), dimension(2) :: vdims, maxdims, dcount2, doffset2
    real(kind=rp), allocatable :: coords(:,:)
    integer, allocatable :: connectivity(:), offsets(:)
    integer, allocatable :: cell_types(:)
    integer(kind=1), allocatable :: cell_types_byte(:)
    logical :: link_exists

    lx = dof%Xh%lx
    ly = dof%Xh%ly
    lz = dof%Xh%lz
    npts_per_cell = lx * ly * lz
    if (msh%gdim .eq. 3) then
       nodes_per_cell = 8
    else
       nodes_per_cell = 4
    end if

    local_cells = part_cells(pe_rank + 1)
    local_conn = part_conns(pe_rank + 1)
    cell_offset = sum(part_cells(1:pe_rank))
    conn_offset = sum(part_conns(1:pe_rank))
    offsets_offset = cell_offset + pe_rank
    max_local_cells = maxval(part_cells)
    max_local_conn = maxval(part_conns)
    total_offsets = total_cells + num_partitions

    ! --- NumberOfPoints dataset (per-partition) ---
    call h5lexists_f(vtkhdf_grp, "NumberOfPoints", link_exists, ierr)
    if (link_exists) call h5ldelete_f(vtkhdf_grp, "NumberOfPoints", ierr)
    vdims(1) = int(num_partitions, hsize_t)
    call h5screate_simple_f(1, vdims(1:1), filespace, ierr)
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    chunkdims(1) = int(num_partitions, hsize_t)
    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
    call h5dcreate_f(vtkhdf_grp, "NumberOfPoints", H5T_NATIVE_INTEGER, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5sclose_f(filespace, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)
    dcount(1) = 1_hsize_t
    doffset_1d(1) = int(pe_rank, hsize_t)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset_1d, dcount, ierr)
    call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, local_points, &
         dcount(1:1), ierr, file_space_id = filespace, &
         mem_space_id = memspace, xfer_prp = plist_coll)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(dcpl_id, ierr)

    ! --- NumberOfCells dataset (per-partition) ---
    call h5lexists_f(vtkhdf_grp, "NumberOfCells", link_exists, ierr)
    if (link_exists) call h5ldelete_f(vtkhdf_grp, "NumberOfCells", ierr)
    vdims(1) = int(num_partitions, hsize_t)
    call h5screate_simple_f(1, vdims(1:1), filespace, ierr)
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    chunkdims(1) = int(num_partitions, hsize_t)
    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
    call h5dcreate_f(vtkhdf_grp, "NumberOfCells", H5T_NATIVE_INTEGER, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5sclose_f(filespace, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)
    dcount(1) = 1_hsize_t
    doffset_1d(1) = int(pe_rank, hsize_t)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset_1d, dcount, ierr)
    call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, local_cells, &
         dcount(1:1), ierr, file_space_id = filespace, &
         mem_space_id = memspace, xfer_prp = plist_coll)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(dcpl_id, ierr)

    ! --- NumberOfConnectivityIds dataset (per-partition) ---
    call h5lexists_f(vtkhdf_grp, "NumberOfConnectivityIds", link_exists, ierr)
    if (link_exists) call h5ldelete_f(vtkhdf_grp, "NumberOfConnectivityIds", ierr)
    vdims(1) = int(num_partitions, hsize_t)
    call h5screate_simple_f(1, vdims(1:1), filespace, ierr)
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    chunkdims(1) = int(num_partitions, hsize_t)
    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
    call h5dcreate_f(vtkhdf_grp, "NumberOfConnectivityIds", H5T_NATIVE_INTEGER, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5sclose_f(filespace, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)
    dcount(1) = 1_hsize_t
    doffset_1d(1) = int(pe_rank, hsize_t)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset_1d, dcount, ierr)
    call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, local_conn, &
         dcount(1:1), ierr, file_space_id = filespace, &
         mem_space_id = memspace, xfer_prp = plist_coll)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(dcpl_id, ierr)

    ! --- Points dataset (global coordinates) ---
    call h5lexists_f(vtkhdf_grp, "Points", link_exists, ierr)
    if (link_exists) call h5ldelete_f(vtkhdf_grp, "Points", ierr)
    allocate(coords(3, local_points))
    do i = 1, msh%nelv
       local_idx = (i - 1) * npts_per_cell
       do kk = 1, lz
          do jj = 1, ly
             do ii = 1, lx
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
    chunkdims(1) = max(1_hsize_t, min(int(max_local_points, hsize_t), vdims(2)))
    call h5pset_chunk_f(dcpl_id, 2, [3_hsize_t, chunkdims(1)], ierr)
    call h5dcreate_f(vtkhdf_grp, "Points", H5T_NEKO_REAL, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5sclose_f(filespace, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)
    dcount2 = [3_hsize_t, int(local_points, hsize_t)]
    doffset2 = [0_hsize_t, int(point_offset, hsize_t)]
    call h5screate_simple_f(2, dcount2, memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset2, dcount2, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, coords, dcount2, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_coll)
    call h5sclose_f(memspace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(dcpl_id, ierr)
    call h5sclose_f(filespace, ierr)
    deallocate(coords)

    ! --- Connectivity dataset ---
    call h5lexists_f(vtkhdf_grp, "Connectivity", link_exists, ierr)
    if (link_exists) call h5ldelete_f(vtkhdf_grp, "Connectivity", ierr)
    allocate(connectivity(local_conn))
    if (msh%gdim .eq. 3) then
       call vtkhdf_build_connectivity(connectivity, 12, lx, ly, lz, msh%nelv)
    else
       call vtkhdf_build_connectivity(connectivity, 9, lx, ly, lz, msh%nelv)
    end if

    vdims(1) = int(total_conn, hsize_t)
    maxdims(1) = H5S_UNLIMITED_F
    call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    chunkdims(1) = max(1_hsize_t, min(int(max_local_conn, hsize_t), vdims(1)))
    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
    call h5dcreate_f(vtkhdf_grp, "Connectivity", H5T_NATIVE_INTEGER, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5sclose_f(filespace, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)
    dcount(1) = int(local_conn, hsize_t)
    doffset_1d(1) = int(conn_offset, hsize_t)
    call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset_1d, dcount, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, connectivity, dcount(1:1), ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_coll)
    call h5sclose_f(memspace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(dcpl_id, ierr)
    call h5sclose_f(filespace, ierr)
    deallocate(connectivity)

    ! --- Offsets dataset ---
    call h5lexists_f(vtkhdf_grp, "Offsets", link_exists, ierr)
    if (link_exists) call h5ldelete_f(vtkhdf_grp, "Offsets", ierr)
    allocate(offsets(local_cells + 1))
    do i = 1, local_cells
       offsets(i) = (i - 1) * nodes_per_cell
    end do
    offsets(local_cells + 1) = local_conn

    vdims(1) = int(total_offsets, hsize_t)
    maxdims(1) = H5S_UNLIMITED_F
    call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    chunkdims(1) = max(1_hsize_t, min(int(max_local_cells + 1, hsize_t), vdims(1)))
    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
    call h5dcreate_f(vtkhdf_grp, "Offsets", H5T_NATIVE_INTEGER, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5sclose_f(filespace, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)
    dcount(1) = int(local_cells + 1, hsize_t)
    doffset_1d(1) = int(offsets_offset, hsize_t)
    call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset_1d, dcount, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, offsets, dcount(1:1), ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_coll)
    call h5sclose_f(memspace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(dcpl_id, ierr)
    call h5sclose_f(filespace, ierr)
    deallocate(offsets)

    ! --- Types dataset (VTK cell types) ---
    call h5lexists_f(vtkhdf_grp, "Types", link_exists, ierr)
    if (link_exists) call h5ldelete_f(vtkhdf_grp, "Types", ierr)
    allocate(cell_types(local_cells))
    allocate(cell_types_byte(local_cells))
    if (msh%gdim .eq. 3) then
       cell_types = 12 ! VTK_HEXAHEDRON
    else
       cell_types = 9 ! VTK_QUAD
    end if
    cell_types_byte = int(cell_types, kind=1)

    vdims(1) = int(total_cells, hsize_t)
    maxdims(1) = H5S_UNLIMITED_F
    call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    chunkdims(1) = max(1_hsize_t, min(int(max_local_cells, hsize_t), vdims(1)))
    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
    call h5dcreate_f(vtkhdf_grp, "Types", H5T_STD_U8LE, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5sclose_f(filespace, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)
    dcount(1) = int(local_cells, hsize_t)
    doffset_1d(1) = int(cell_offset, hsize_t)
    call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset_1d, dcount, ierr)
    call h5dwrite_f(dset_id, H5T_STD_U8LE, cell_types_byte, dcount(1:1), ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_coll)
    call h5sclose_f(memspace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(dcpl_id, ierr)
    call h5sclose_f(filespace, ierr)
    deallocate(cell_types, cell_types_byte)

  end subroutine vtkhdf_write_mesh

  !> Write temporal Steps group metadata to the VTKHDF group.
  !! Writes Values, NumberOfParts, PartOffsets, PointOffsets,
  !! CellOffsets, ConnectivityIdOffsets datasets, and NSteps attribute.
  !! @param vtkhdf_grp HDF5 group ID for VTKHDF root group
  !! @param plist_coll Collective transfer property list
  !! @param H5T_NEKO_REAL HDF5 real type matching NEKO precision
  !! @param t Current simulation time
  !! @param amr_enabled Whether AMR is enabled
  !! @param num_partitions Number of MPI partitions
  !! @param total_points Global total number of points
  !! @param total_cells Global total number of sub-cells
  !! @param total_conn Global total connectivity size
  !! @param pointdata_time_offset Output: offset for PointData at this timestep
  subroutine vtkhdf_write_steps(vtkhdf_grp, plist_coll, H5T_NEKO_REAL, t, &
       amr_enabled, num_partitions, &
       total_points, total_cells, total_conn, &
       pointdata_time_offset)
    integer(hid_t), intent(in) :: vtkhdf_grp, plist_coll, H5T_NEKO_REAL
    real(kind=rp), intent(in) :: t
    logical, intent(in) :: amr_enabled
    integer, intent(in) :: num_partitions
    integer, intent(in) :: total_points, total_cells, total_conn
    integer(hsize_t), intent(out) :: pointdata_time_offset

    integer :: ierr, N_Steps
    integer(hid_t) :: grp_id, dset_id, dcpl_id, filespace, memspace, attr_id
    integer(hsize_t), dimension(1) :: step_dims, step_maxdims
    integer(hsize_t), dimension(1) :: step_count, step_offset, chunkdims, ddim
    real(kind=rp), dimension(1) :: time_value
    integer(kind=8), dimension(1) :: i8_value
    logical :: link_exists, attr_exists

    ! Create or open Steps group
    call h5lexists_f(vtkhdf_grp, "Steps", link_exists, ierr)
    if (link_exists) then
       call h5gopen_f(vtkhdf_grp, "Steps", grp_id, ierr)
    else
       call h5gcreate_f(vtkhdf_grp, "Steps", grp_id, ierr)
    end if

    ! --- Values dataset (time values, real type) ---
    call h5lexists_f(grp_id, "Values", link_exists, ierr)
    if (link_exists) then
       call h5dopen_f(grp_id, "Values", dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sget_simple_extent_dims_f(filespace, step_dims, step_maxdims, ierr)
       call h5sclose_f(filespace, ierr)
    else
       step_dims(1) = 0_hsize_t
       step_maxdims(1) = H5S_UNLIMITED_F
       call h5screate_simple_f(1, step_dims, filespace, ierr, step_maxdims)
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       chunkdims(1) = 1_hsize_t
       call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
       call h5dcreate_f(grp_id, "Values", H5T_NEKO_REAL, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5pclose_f(dcpl_id, ierr)
       call h5sclose_f(filespace, ierr)
    end if

    step_dims(1) = step_dims(1) + 1_hsize_t
    call h5dset_extent_f(dset_id, step_dims, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)
    step_count(1) = 1_hsize_t
    step_offset(1) = step_dims(1) - 1_hsize_t
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         step_offset, step_count, ierr)
    call h5screate_simple_f(1, step_count, memspace, ierr)
    time_value(1) = t
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, time_value, step_count, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_coll)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)

    N_Steps = int(step_dims(1))
    pointdata_time_offset = int(N_Steps - 1, hsize_t) * int(total_points, hsize_t)

    ! --- NumberOfParts ---
    i8_value(1) = int(num_partitions, kind=8)
    call vtkhdf_append_step_i8(grp_id, "NumberOfParts", i8_value(1), plist_coll)

    ! --- PartOffsets ---
    if (amr_enabled) then
       i8_value(1) = int(N_Steps - 1, kind=8) * int(num_partitions, kind=8)
    else
       i8_value(1) = 0
    end if
    call vtkhdf_append_step_i8(grp_id, "PartOffsets", i8_value(1), plist_coll)

    ! --- PointOffsets ---
    if (amr_enabled) then
       i8_value(1) = int(N_Steps - 1, kind=8) * int(total_points, kind=8)
    else
       i8_value(1) = 0
    end if
    call vtkhdf_append_step_i8(grp_id, "PointOffsets", i8_value(1), plist_coll)

    ! --- CellOffsets ---
    if (amr_enabled) then
       i8_value(1) = int(N_Steps - 1, kind=8) * int(total_cells, kind=8)
    else
       i8_value(1) = 0
    end if
    call vtkhdf_append_step_i8(grp_id, "CellOffsets", i8_value(1), plist_coll)

    ! --- ConnectivityIdOffsets ---
    if (amr_enabled) then
       i8_value(1) = int(N_Steps - 1, kind=8) * int(total_conn, kind=8)
    else
       i8_value(1) = 0
    end if
    call vtkhdf_append_step_i8(grp_id, "ConnectivityIdOffsets", i8_value(1), plist_coll)

    ! --- NSteps attribute ---
    ddim(1) = 1_hsize_t
    call h5aexists_f(grp_id, "NSteps", attr_exists, ierr)
    if (attr_exists) then
       call h5aopen_f(grp_id, "NSteps", attr_id, ierr)
    else
       call h5screate_f(H5S_SCALAR_F, filespace, ierr)
       call h5acreate_f(grp_id, "NSteps", H5T_NATIVE_INTEGER, filespace, &
            attr_id, ierr, h5p_default_f, h5p_default_f)
       call h5sclose_f(filespace, ierr)
    end if
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, N_Steps, ddim(1:1), ierr)
    call h5aclose_f(attr_id, ierr)

    call h5gclose_f(grp_id, ierr)

  end subroutine vtkhdf_write_steps

  !> Append an integer(kind=8) scalar to a 1D chunked HDF5 dataset.
  !! Opens the dataset if it exists, or creates a new empty chunked dataset.
  !! Extends by 1 element and writes the value at the end.
  !! @param grp_id HDF5 group containing the dataset
  !! @param name Dataset name
  !! @param value Value to append
  !! @param plist_coll Collective transfer property list
  subroutine vtkhdf_append_step_i8(grp_id, name, value, plist_coll)
    integer(hid_t), intent(in) :: grp_id, plist_coll
    character(len=*), intent(in) :: name
    integer(kind=8), intent(in) :: value

    integer :: ierr
    integer(hid_t) :: dset_id, dcpl_id, filespace, memspace
    integer(hsize_t), dimension(1) :: dims, maxdims, cnt, off, chunkdims
    integer(kind=8), dimension(1) :: buf
    logical :: link_exists

    call h5lexists_f(grp_id, name, link_exists, ierr)
    if (link_exists) then
       call h5dopen_f(grp_id, name, dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sget_simple_extent_dims_f(filespace, dims, maxdims, ierr)
       call h5sclose_f(filespace, ierr)
    else
       dims(1) = 0_hsize_t
       maxdims(1) = H5S_UNLIMITED_F
       call h5screate_simple_f(1, dims, filespace, ierr, maxdims)
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       chunkdims(1) = 1_hsize_t
       call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
       call h5dcreate_f(grp_id, name, H5T_STD_I64LE, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5pclose_f(dcpl_id, ierr)
       call h5sclose_f(filespace, ierr)
    end if

    dims(1) = dims(1) + 1_hsize_t
    call h5dset_extent_f(dset_id, dims, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)
    cnt(1) = 1_hsize_t
    off(1) = dims(1) - 1_hsize_t
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, off, cnt, ierr)
    call h5screate_simple_f(1, cnt, memspace, ierr)
    buf(1) = value
    call h5dwrite_f(dset_id, H5T_STD_I64LE, buf, cnt, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_coll)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)

  end subroutine vtkhdf_append_step_i8

  !> Write field data into the VTKHDF PointData group.
  !! Groups u, v, w fields into a 3-component Velocity vector dataset.
  !! All other fields are written as scalar datasets.
  !! For temporal output, PointDataOffsets are appended under Steps.
  !! @param vtkhdf_grp Root VTKHDF group
  !! @param plist_coll Collective transfer property list
  !! @param H5T_NEKO_REAL HDF5 real type matching Neko precision
  !! @param fp Array of field pointers to write
  !! @param field_written Tracking array for already-written fields
  !! @param n_fields Number of fields
  !! @param nelv Number of local elements
  !! @param local_points Number of local points
  !! @param point_offset Point offset for this rank
  !! @param max_local_points Maximum local points across all ranks
  !! @param total_points Total points across all ranks
  !! @param pointdata_time_offset Temporal offset for PointData datasets
  !! @param lx Polynomial order in x
  !! @param ly Polynomial order in y
  !! @param lz Polynomial order in z
  !! @param has_time Whether temporal data (Steps) are being written
  subroutine vtkhdf_write_pointdata(vtkhdf_grp, plist_coll, H5T_NEKO_REAL, &
       fp, field_written, nelv, &
       local_points, point_offset, max_local_points, total_points, &
       pointdata_time_offset, lx, ly, lz, has_time)
    integer(hid_t), intent(in) :: vtkhdf_grp, plist_coll, H5T_NEKO_REAL
    type(field_ptr_t), intent(in) :: fp(:)
    logical, intent(inout) :: field_written(:)
    integer, intent(in) :: nelv
    integer, intent(in) :: local_points, point_offset
    integer, intent(in) :: max_local_points, total_points
    integer(hsize_t), intent(in) :: pointdata_time_offset
    integer, intent(in) :: lx, ly, lz
    logical, intent(in) :: has_time

    integer :: ierr, i, j, ie, ii, jj, kk, local_idx
    integer :: n_fields, npts_per_cell
    integer(hid_t) :: pointdata_grp, pointdata_offsets_grp, grp_id
    integer(hid_t) :: dset_id, dcpl_id, filespace, memspace
    integer(hsize_t), dimension(1) :: dcount, doffset, chunkdims
    integer(hsize_t), dimension(2) :: dcount2, doffset2
    integer(hsize_t), dimension(1) :: pd_dims1, pd_maxdims1
    integer(hsize_t), dimension(2) :: pd_dims2, pd_maxdims2
    type(field_t), pointer :: fld, u, v, w
    real(kind=rp), allocatable :: point_data(:,:)
    character(len=128) :: field_name
    logical :: link_exists, is_vector

    n_fields = size(fp)
    npts_per_cell = lx * ly * lz

    ! Sync all the fields
    do i = 1, n_fields
       if (associated(fp(i)%ptr)) then
          call fp(i)%ptr%copy_from(DEVICE_TO_HOST, sync = i .eq. n_fields)
       end if
    end do

    ! Create or open PointData group
    call h5lexists_f(vtkhdf_grp, "PointData", link_exists, ierr)
    if (link_exists) then
       call h5gopen_f(vtkhdf_grp, "PointData", pointdata_grp, ierr)
    else
       call h5gcreate_f(vtkhdf_grp, "PointData", pointdata_grp, ierr, &
            lcpl_id=h5p_default_f, gcpl_id=h5p_default_f, gapl_id=h5p_default_f)
    end if

    dcount(1) = int(local_points, hsize_t)
    doffset(1) = pointdata_time_offset + int(point_offset, hsize_t)

    do i = 1, n_fields
       if (field_written(i)) cycle
       fld => fp(i)%ptr

       field_name = fld%name
       if (field_name .eq. 'p') field_name = 'Pressure'

       ! Determine if this is a velocity component to group as a vector
       is_vector = (field_name .eq. 'u' .or. &
            field_name .eq. 'v' .or. &
            field_name .eq. 'w')

       if (is_vector) then
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
       end if

       ! Write PointDataOffsets under Steps for temporal output
       if (has_time) then
          call h5gopen_f(vtkhdf_grp, "Steps", grp_id, ierr)
          call h5lexists_f(grp_id, "PointDataOffsets", link_exists, ierr)
          if (link_exists) then
             call h5gopen_f(grp_id, "PointDataOffsets", &
                  pointdata_offsets_grp, ierr)
             if (ierr .ne. 0) then
                call h5ldelete_f(grp_id, "PointDataOffsets", ierr)
                call h5gcreate_f(grp_id, "PointDataOffsets", &
                     pointdata_offsets_grp, ierr)
             end if
          else
             call h5gcreate_f(grp_id, "PointDataOffsets", &
                  pointdata_offsets_grp, ierr)
          end if
          call vtkhdf_append_step_i8(pointdata_offsets_grp, &
               trim(field_name), int(pointdata_time_offset, kind=8), &
               plist_coll)
          call h5gclose_f(pointdata_offsets_grp, ierr)
          call h5gclose_f(grp_id, ierr)
       end if

       if (is_vector) then
          ! Assemble 3-component vector from u, v, w fields
          allocate(point_data(3, local_points))
          do ie = 1, nelv
             local_idx = (ie - 1) * npts_per_cell
             do kk = 1, lz
                do jj = 1, ly
                   do ii = 1, lx
                      local_idx = local_idx + 1
                      point_data(1, local_idx) = u%x(ii, jj, kk, ie)
                      point_data(2, local_idx) = v%x(ii, jj, kk, ie)
                      if (associated(w)) then
                         point_data(3, local_idx) = w%x(ii, jj, kk, ie)
                      else
                         point_data(3, local_idx) = 0.0_rp
                      end if
                   end do
                end do
             end do
          end do

          ! Create or extend 2D dataset (3 x total_points)
          call h5lexists_f(pointdata_grp, trim(field_name), link_exists, ierr)
          if (link_exists) then
             call h5dopen_f(pointdata_grp, trim(field_name), dset_id, ierr)
             call h5dget_space_f(dset_id, filespace, ierr)
             call h5sget_simple_extent_dims_f(filespace, pd_dims2, &
                  pd_maxdims2, ierr)
             call h5sclose_f(filespace, ierr)
          else
             pd_dims2 = [3_hsize_t, 0_hsize_t]
             pd_maxdims2 = [3_hsize_t, H5S_UNLIMITED_F]
             call h5screate_simple_f(2, pd_dims2, filespace, ierr, pd_maxdims2)
             call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
             chunkdims(1) = max(1_hsize_t, int(max_local_points, hsize_t))
             call h5pset_chunk_f(dcpl_id, 2, [3_hsize_t, chunkdims(1)], ierr)
             call h5dcreate_f(pointdata_grp, trim(field_name), H5T_NEKO_REAL, &
                  filespace, dset_id, ierr, dcpl_id = dcpl_id)
             call h5pclose_f(dcpl_id, ierr)
             call h5sclose_f(filespace, ierr)
             pd_dims2 = [3_hsize_t, 0_hsize_t]
          end if

          pd_dims2(2) = pd_dims2(2) + int(total_points, hsize_t)
          call h5dset_extent_f(dset_id, pd_dims2, ierr)
          call h5dget_space_f(dset_id, filespace, ierr)

          dcount2 = [3_hsize_t, int(local_points, hsize_t)]
          doffset2 = [0_hsize_t, doffset(1)]
          call h5screate_simple_f(2, dcount2, memspace, ierr)
          call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
               doffset2, dcount2, ierr)
          call h5dwrite_f(dset_id, H5T_NEKO_REAL, point_data, dcount2, ierr, &
               file_space_id = filespace, mem_space_id = memspace, &
               xfer_prp = plist_coll)
          call h5sclose_f(filespace, ierr)
          call h5sclose_f(memspace, ierr)
          call h5dclose_f(dset_id, ierr)
          deallocate(point_data)

       else
          ! Write scalar field as 1D dataset
          call h5lexists_f(pointdata_grp, trim(field_name), link_exists, ierr)
          if (link_exists) then
             call h5dopen_f(pointdata_grp, trim(field_name), dset_id, ierr)
             call h5dget_space_f(dset_id, filespace, ierr)
             call h5sget_simple_extent_dims_f(filespace, pd_dims1, &
                  pd_maxdims1, ierr)
             call h5sclose_f(filespace, ierr)
          else
             pd_dims1(1) = 0_hsize_t
             pd_maxdims1(1) = H5S_UNLIMITED_F
             call h5screate_simple_f(1, pd_dims1, filespace, ierr, pd_maxdims1)
             call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
             chunkdims(1) = max(1_hsize_t, int(max_local_points, hsize_t))
             call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
             call h5dcreate_f(pointdata_grp, trim(field_name), H5T_NEKO_REAL, &
                  filespace, dset_id, ierr, dcpl_id = dcpl_id)
             call h5pclose_f(dcpl_id, ierr)
             call h5sclose_f(filespace, ierr)
             pd_dims1(1) = 0_hsize_t
          end if

          pd_dims1(1) = pd_dims1(1) + int(total_points, hsize_t)
          call h5dset_extent_f(dset_id, pd_dims1, ierr)
          call h5dget_space_f(dset_id, filespace, ierr)
          call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
          call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
               doffset(1:1), dcount(1:1), ierr)
          call h5dwrite_f(dset_id, H5T_NEKO_REAL, fld%x(1,1,1,1), &
               dcount(1:1), ierr, file_space_id = filespace, &
               mem_space_id = memspace, xfer_prp = plist_coll)
          call h5sclose_f(filespace, ierr)
          call h5sclose_f(memspace, ierr)
          call h5dclose_f(dset_id, ierr)
       end if
    end do

    call h5gclose_f(pointdata_grp, ierr)

  end subroutine vtkhdf_write_pointdata

  !> Read data in HDF5 format following official VTKHDF specification
  subroutine vtkhdf_file_read(this, data)
    class(vtkhdf_file_t) :: this
    class(*), target, intent(inout) :: data

    call neko_error('VTKHDF file reading is not yet implemented')

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

#endif

  !> Set the overwrite flag for HDF5 files
  subroutine vtkhdf_file_set_overwrite(this, overwrite)
    class(vtkhdf_file_t), intent(inout) :: this
    logical, intent(in) :: overwrite
    this%overwrite = overwrite
  end subroutine vtkhdf_file_set_overwrite

  !> Enable support for Adaptive Mesh Refinement
  subroutine vtkhdf_file_enable_amr(this)
    class(vtkhdf_file_t), intent(inout) :: this
    this%amr_enabled = .false.
  end subroutine vtkhdf_file_enable_amr

end module vtkhdf_file
