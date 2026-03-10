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
!> VTKHDF file format
module vtkhdf_file
  use num_types, only : rp, sp, dp
  use generic_file, only : generic_file_t
  use checkpoint, only : chkp_t
  use utils, only : neko_error, neko_warning, filename_suffix_pos, &
       nonlinear_index
  use mesh, only : mesh_t
  use field, only : field_t, field_ptr_t
  use field_list, only : field_list_t
  use field_series, only : field_series_t, field_series_ptr_t
  use dofmap, only : dofmap_t
  use logger, only : neko_log
  use comm, only : pe_rank, pe_size, NEKO_COMM
  use device, only : DEVICE_TO_HOST
  use mpi_f08, only : MPI_INFO_NULL, MPI_Allreduce, MPI_Allgather, &
       MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, MPI_MAX, MPI_Comm_size, MPI_Exscan, &
       MPI_Barrier
#ifdef HAVE_HDF5
  use hdf5
#endif
  implicit none
  private

  !> Interface for HDF5 files
  type, public, extends(generic_file_t) :: vtkhdf_file_t
     logical :: amr_enabled = .false.
     integer :: precision = 0
   contains
     procedure :: read => vtkhdf_file_read
     procedure :: write => vtkhdf_file_write
     procedure :: set_overwrite => vtkhdf_file_set_overwrite
     procedure :: enable_amr => vtkhdf_file_enable_amr
     procedure :: set_precision => vtkhdf_file_set_precision
  end type vtkhdf_file_t

  integer, dimension(2), parameter :: vtkhdf_version = [2, 6]

contains

  ! -------------------------------------------------------------------------- !
  ! Well defined subroutines

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

#ifdef HAVE_HDF5
  ! -------------------------------------------------------------------------- !
  ! HDF5 Required subroutines

  !> Write data in HDF5 format following official VTKHDF UnstructuredGrid
  !! specification
  subroutine vtkhdf_file_write(this, data, t)
    class(vtkhdf_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dof
    type(field_ptr_t), allocatable :: fp(:)
    integer :: ierr, mpi_info, mpi_comm, i, n_fields
    integer(hid_t) :: plist_id, file_id, attr_id, vtkhdf_grp
    integer(hid_t) :: filespace, memspace
    integer(hsize_t), dimension(1) :: vdims
    integer :: lx, ly, lz
    integer :: local_points, local_cells, local_conn
    integer :: total_points, total_cells, total_conn
    integer :: point_offset
    integer :: max_local_points
    integer, allocatable :: part_points(:), part_cells(:), part_conns(:)
    character(len=1024) :: fname
    character(len=16) :: type_str
    logical :: link_exists, file_exists
    integer(kind=1) :: VTK_cell_type
    integer :: counter

    ! Determine mesh and field data
    select type(data)
    type is (field_t)
       msh => data%msh
       dof => data%dof
       n_fields = 1
       allocate(fp(1))
       fp(1)%ptr => data
    type is (field_list_t)
       msh => data%msh(1)
       dof => data%dof(1)
       n_fields = data%size()
       allocate(fp(n_fields))
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

    ! Ensure precision is set and are valid.
    if (this%precision .gt. rp) then
       this%precision = rp
       call neko_warning('Requested precision is higher than working precision')
    else if (this%precision .eq. 0) then
       this%precision = rp
    end if

    if (msh%gdim .eq. 3) then
       VTK_cell_type = 12 ! VTK_HEXAHEDRON
    else
       VTK_cell_type = 9 ! VTK_QUAD
    end if

    call this%increment_counter()
    fname = trim(this%get_base_fname())
    counter = this%get_counter()

    mpi_info = MPI_INFO_NULL%mpi_val
    mpi_comm = NEKO_COMM%mpi_val

    call h5open_f(ierr)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    call h5pset_fapl_mpio_f(plist_id, mpi_comm, mpi_info, ierr)

    inquire(file = fname, exist = file_exists)
    if (file_exists) then
       call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, ierr, &
            access_prp = plist_id)
    else
       call h5fcreate_f(fname, H5F_ACC_TRUNC_F, &
            file_id, ierr, access_prp = plist_id)
    end if

    ! Create/open VTKHDF root group with vtkhdf_version and type attributes
    call h5lexists_f(file_id, "VTKHDF", link_exists, ierr)
    if (link_exists) then
       call h5gopen_f(file_id, "VTKHDF", vtkhdf_grp, ierr)
    else
       call h5gcreate_f(file_id, "VTKHDF", vtkhdf_grp, ierr)

       ! Write Version attribute
       vdims = 2
       call h5screate_simple_f(1, vdims, filespace, ierr)
       call h5acreate_f(vtkhdf_grp, "Version", H5T_NATIVE_INTEGER, filespace, &
            attr_id, ierr)
       call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, vtkhdf_version, vdims, ierr)
       call h5aclose_f(attr_id, ierr)
       call h5sclose_f(filespace, ierr)

       ! Write Type attribute "UnstructuredGrid" as a fixed-length string
       type_str = "UnstructuredGrid"
       vdims = 1
       call h5screate_f(H5S_SCALAR_F, filespace, ierr)

       call h5tcopy_f(H5T_FORTRAN_S1, memspace, ierr)
       call h5tset_size_f(memspace, int(len_trim(type_str), kind=hsize_t), ierr)
       call h5tset_strpad_f(memspace, H5T_STR_NULLTERM_F, ierr)

       call h5acreate_f(vtkhdf_grp, "Type", memspace, filespace, attr_id, ierr)
       call h5awrite_f(attr_id, memspace, [type_str], vdims, ierr)
       call h5aclose_f(attr_id, ierr)

       call h5tclose_f(memspace, ierr)
       call h5sclose_f(filespace, ierr)
    end if

    if (present(t)) then
       call vtkhdf_write_steps(vtkhdf_grp, counter, t)
    end if

    if (associated(msh)) then
       call vtkhdf_write_mesh(vtkhdf_grp, dof, msh, VTK_cell_type, &
            this%amr_enabled, counter, t)
    end if

    ! Write field data in PointData group
    if (n_fields > 0) then
       call vtkhdf_write_pointdata(vtkhdf_grp, fp, this%precision, counter, t)
    end if

    call h5gclose_f(vtkhdf_grp, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)

    if (allocated(fp)) deallocate(fp)

  end subroutine vtkhdf_file_write

  !> Build local connectivity for VTK sub-cells from a spectral element
  !! tensor-product grid. Subdivides each spectral element into linear
  !! sub-cells based on the GLL node positions.
  !! @param conn Output connectivity array (pre-allocated)
  !! @param vtk_type VTK cell type: 12 = VTK_HEXAHEDRON, 9 = VTK_QUAD
  !! @param msh Mesh object containing element information
  !! @param dof Dofmap containing the lx, ly, lz dimensions of the spectral
  !!            element grid
  subroutine vtkhdf_build_connectivity(conn, vtk_type, msh, dof)
    integer, intent(out) :: conn(:)
    integer(kind=1), intent(in) :: vtk_type
    type(mesh_t) :: msh
    type(dofmap_t) :: dof
    integer :: lx, ly, lz, nelv
    integer :: ie, ii, jj, kk, base, idx, npts_per_cell

    nelv = msh%nelv
    lx = dof%Xh%lx
    ly = dof%Xh%ly
    lz = dof%Xh%lz
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
                        (kk - 1) * lx * ly + jj * lx + (ii + 1) - 1
                   conn(idx + 4) = base + &
                        (kk - 1) * lx * ly + jj * lx + ii - 1
                   conn(idx + 5) = base + &
                        kk * lx * ly + (jj - 1) * lx + ii - 1
                   conn(idx + 6) = base + &
                        kk * lx * ly + (jj - 1) * lx + (ii + 1) - 1
                   conn(idx + 7) = base + &
                        kk * lx * ly + jj * lx + (ii + 1) - 1
                   conn(idx + 8) = base + &
                        kk * lx * ly + jj * lx + ii - 1
                   idx = idx + 8
                end do
             end do
          end do

       case (9) ! VTK_QUAD
          do jj = 1, ly - 1
             do ii = 1, lx - 1
                conn(idx + 1) = base + (jj - 1) * lx + ii - 1
                conn(idx + 2) = base + (jj - 1) * lx + (ii + 1) - 1
                conn(idx + 3) = base + jj * lx + (ii + 1) - 1
                conn(idx + 4) = base + jj * lx + ii - 1
                idx = idx + 4
             end do
          end do

       case default
          call neko_error('Unsupported VTK cell type')
       end select
    end do

  end subroutine vtkhdf_build_connectivity

  !> Write mesh geometry datasets to the VTKHDF group.
  !! Writes NumberOfPoints, NumberOfCells, NumberOfConnectivityIds,
  !! Points, Connectivity, Offsets, and Types datasets.
  !! @param vtkhdf_grp HDF5 group ID for VTKHDF root group
  !! @param dof Dofmap for coordinate data
  !! @param msh Mesh object
  !! @param VTK_cell_type VTK cell type (e.g. 12 for hexahedra, 9 for quads)
  !! @param amr AMR flag to determine if mesh should be rewritten at every time
  !!            step
  !! @param t Optional time value for time-dependent mesh output (e.g. for AMR)
  subroutine vtkhdf_write_mesh(vtkhdf_grp, dof, msh, VTK_cell_type, amr, &
       counter, t)
    type(dofmap_t), intent(in) :: dof
    type(mesh_t), intent(in) :: msh
    integer(hid_t), intent(in) :: vtkhdf_grp
    integer(kind=1), intent(in) :: VTK_cell_type
    logical, intent(in) :: amr
    integer, intent(in) :: counter
    real(kind=rp), intent(in), optional :: t

    integer :: ierr, i, ii, jj, kk, el, local_idx
    integer :: lx, ly, lz, npts_per_cell, nodes_per_cell, cells_per_element
    integer :: local_points, local_cells, local_conn
    integer :: total_points, total_cells, total_conn
    integer :: point_offset, max_local_points
    integer :: total_offsets, cell_offset, conn_offset, offsets_offset
    integer :: max_local_cells, max_local_conn
    integer(hid_t) :: xf_id, dset_id, dcpl_id, grp_id, attr_id
    integer(hid_t) :: filespace, memspace
    integer(hsize_t), dimension(1) :: dcount, doffset_1d, chunkdims
    integer(hsize_t), dimension(2) :: vdims, maxdims, dcount2, doffset2
    integer(kind=8) :: i8_value
    logical :: link_exists
    integer, dimension(3) :: component_sizes
    integer, dimension(3) :: component_offsets
    integer, dimension(3) :: component_max_sizes


    lx = dof%Xh%lx
    ly = dof%Xh%ly
    lz = dof%Xh%lz
    npts_per_cell = lx * ly * lz

    select case(VTK_cell_type)
    case(12)
       cells_per_element = (lx - 1) * (ly - 1) * (lz - 1)
       nodes_per_cell = 8
    case(9)
       cells_per_element = (lx - 1) * (ly - 1)
       nodes_per_cell = 4
    case default
       call neko_error('Unsupported VTK cell type')
    end select

    ! --- Build the number of cells and the connectivity
    local_points = dof%size()
    local_cells = msh%nelv * cells_per_element
    local_conn = local_cells * nodes_per_cell

    total_points = dof%global_size()
    total_cells = msh%glb_nelv * cells_per_element
    total_conn = total_cells * nodes_per_cell

    component_sizes = [local_points, local_cells, local_conn]
    component_offsets = 0
    component_max_sizes = 0

    call MPI_Exscan(component_sizes, component_offsets, 3, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)
    call MPI_Allreduce(component_sizes, component_max_sizes, 3, MPI_INTEGER, &
         MPI_MAX, NEKO_COMM, ierr)

    point_offset = component_offsets(1)
    cell_offset = component_offsets(2)
    conn_offset = component_offsets(3)
    max_local_points = component_max_sizes(1)
    max_local_cells = component_max_sizes(2)
    max_local_conn = component_max_sizes(3)

    offsets_offset = cell_offset + pe_rank
    total_offsets = total_cells + pe_size

    ! Create collective transfer property list
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! --- Shared sizes for 1D datasets ---
    dcount(1) = 1_hsize_t
    vdims(1) = int(pe_size, hsize_t)
    chunkdims(1) = int(pe_size, hsize_t)
    doffset_1d(1) = int(pe_rank, hsize_t)

    ! --- Create filespaces and memspaces for 1D datasets ---
    call h5screate_simple_f(1, vdims(1:1), filespace, ierr)
    call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset_1d, dcount, ierr)
    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)

    ! --- NumberOfPoints dataset (per-partition) ---
    call h5lexists_f(vtkhdf_grp, "NumberOfPoints", link_exists, ierr)
    if (.not. link_exists) then
       call h5dcreate_f(vtkhdf_grp, "NumberOfPoints", H5T_NATIVE_INTEGER, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, local_points, &
            dcount(1:1), ierr, file_space_id = filespace, &
            mem_space_id = memspace, xfer_prp = xf_id)
       call h5dclose_f(dset_id, ierr)
    end if

    ! --- NumberOfCells dataset (per-partition) ---
    call h5lexists_f(vtkhdf_grp, "NumberOfCells", link_exists, ierr)
    if (.not. link_exists) then
       call h5dcreate_f(vtkhdf_grp, "NumberOfCells", H5T_NATIVE_INTEGER, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, local_cells, &
            dcount(1:1), ierr, file_space_id = filespace, &
            mem_space_id = memspace, xfer_prp = xf_id)
       call h5dclose_f(dset_id, ierr)
    end if

    ! --- NumberOfConnectivityIds dataset (per-partition) ---
    call h5lexists_f(vtkhdf_grp, "NumberOfConnectivityIds", link_exists, ierr)
    if (.not. link_exists) then
       call h5dcreate_f(vtkhdf_grp, "NumberOfConnectivityIds", &
            H5T_NATIVE_INTEGER, filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, local_conn, &
            dcount(1:1), ierr, file_space_id = filespace, &
            mem_space_id = memspace, xfer_prp = xf_id)
       call h5dclose_f(dset_id, ierr)
    end if

    ! --- Close the mempory and file spaces ---
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5pclose_f(dcpl_id, ierr)

    ! --- Points dataset (global coordinates) ---
    call h5lexists_f(vtkhdf_grp, "Points", link_exists, ierr)
    if (.not. link_exists) then

       vdims = [3_hsize_t, int(total_points, hsize_t)]
       maxdims = [3_hsize_t, H5S_UNLIMITED_F]
       chunkdims(1) = max(1_hsize_t, min(int(max_local_points, hsize_t), &
            vdims(2)))
       dcount2 = [3_hsize_t, int(local_points, hsize_t)]
       doffset2 = [0_hsize_t, int(point_offset, hsize_t)]

       call h5screate_simple_f(2, vdims, filespace, ierr, maxdims)
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       call h5pset_chunk_f(dcpl_id, 2, [3_hsize_t, chunkdims(1)], ierr)
       call h5dcreate_f(vtkhdf_grp, "Points", H5T_NATIVE_DOUBLE, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5screate_simple_f(2, dcount2, memspace, ierr)
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
            doffset2, dcount2, ierr)

       block
         real(kind=dp), allocatable :: coords(:,:)

         allocate(coords(3, local_points))
         do concurrent (local_idx = 1:local_points)
            block
              integer :: idx(4)
              idx = nonlinear_index(local_idx, lx, ly, lz)

              coords(1, local_idx) = dof%x(idx(1), idx(2), idx(3), idx(4))
              coords(2, local_idx) = dof%y(idx(1), idx(2), idx(3), idx(4))
              coords(3, local_idx) = dof%z(idx(1), idx(2), idx(3), idx(4))
            end block
         end do
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, coords, dcount2, ierr, &
              file_space_id = filespace, mem_space_id = memspace, &
              xfer_prp = xf_id)
         deallocate(coords)
       end block

       call h5sclose_f(memspace, ierr)
       call h5dclose_f(dset_id, ierr)
       call h5pclose_f(dcpl_id, ierr)
       call h5sclose_f(filespace, ierr)
    end if

    ! --- Connectivity dataset ---
    call h5lexists_f(vtkhdf_grp, "Connectivity", link_exists, ierr)
    if (link_exists) call h5ldelete_f(vtkhdf_grp, "Connectivity", ierr)

    vdims(1) = int(total_conn, hsize_t)
    maxdims(1) = H5S_UNLIMITED_F
    chunkdims(1) = max(1_hsize_t, min(int(max_local_conn, hsize_t), vdims(1)))
    dcount(1) = int(local_conn, hsize_t)
    doffset_1d(1) = int(conn_offset, hsize_t)

    call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
    call h5dcreate_f(vtkhdf_grp, "Connectivity", H5T_NATIVE_INTEGER, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset_1d, dcount, ierr)

    block
      integer, allocatable :: connectivity(:)

      allocate(connectivity(local_conn))
      call vtkhdf_build_connectivity(connectivity, VTK_cell_type, msh, dof)
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, connectivity, dcount(1:1), &
           ierr, file_space_id = filespace, mem_space_id = memspace, &
           xfer_prp = xf_id)
      deallocate(connectivity)
    end block

    call h5sclose_f(memspace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(dcpl_id, ierr)
    call h5sclose_f(filespace, ierr)

    ! --- Offsets dataset ---
    call h5lexists_f(vtkhdf_grp, "Offsets", link_exists, ierr)
    if (link_exists) call h5ldelete_f(vtkhdf_grp, "Offsets", ierr)

    vdims(1) = int(total_offsets, hsize_t)
    maxdims(1) = H5S_UNLIMITED_F
    chunkdims(1) = max(1_hsize_t, min(int(max_local_cells + 1, hsize_t), &
         vdims(1)))
    dcount(1) = int(local_cells + 1, hsize_t)
    doffset_1d(1) = int(offsets_offset, hsize_t)

    call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
    call h5dcreate_f(vtkhdf_grp, "Offsets", H5T_NATIVE_INTEGER, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset_1d, dcount, ierr)

    block
      integer, allocatable :: offsets(:)

      allocate(offsets(local_cells + 1))
      do concurrent (i = 1:local_cells)
         offsets(i) = (i - 1) * nodes_per_cell
      end do
      offsets(local_cells + 1) = local_conn
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, offsets, dcount(1:1), ierr, &
           file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
      deallocate(offsets)
    end block

    call h5sclose_f(memspace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(dcpl_id, ierr)
    call h5sclose_f(filespace, ierr)

    ! --- Types dataset (VTK cell types) ---
    call h5lexists_f(vtkhdf_grp, "Types", link_exists, ierr)
    if (link_exists) call h5ldelete_f(vtkhdf_grp, "Types", ierr)

    vdims(1) = int(total_cells, hsize_t)
    maxdims(1) = H5S_UNLIMITED_F
    chunkdims(1) = max(1_hsize_t, min(int(max_local_cells, hsize_t), vdims(1)))
    dcount(1) = int(local_cells, hsize_t)
    doffset_1d(1) = int(cell_offset, hsize_t)

    call h5screate_simple_f(1, vdims(1:1), filespace, ierr, maxdims(1:1))
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
    call h5dcreate_f(vtkhdf_grp, "Types", H5T_STD_U8LE, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5screate_simple_f(1, dcount(1:1), memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset_1d, dcount, ierr)

    block
      integer(kind=1), allocatable :: cell_types(:)
      allocate(cell_types(local_cells), source=VTK_cell_type)
      call h5dwrite_f(dset_id, H5T_STD_U8LE, cell_types, dcount(1:1), ierr, &
           file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
      deallocate(cell_types)
    end block

    call h5sclose_f(memspace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(dcpl_id, ierr)
    call h5sclose_f(filespace, ierr)

    if (present(t)) then
       ! Open Steps group
       call h5gopen_f(vtkhdf_grp, "Steps", grp_id, ierr)

       ! --- NumberOfParts ---
       call vtkhdf_write_i8_at(grp_id, "NumberOfParts", int(pe_size, kind=8), &
            counter)

       ! --- PartOffsets ---
       i8_value = 0
       if (amr) then
          i8_value = int(counter - 1, kind=8) * int(pe_size, kind=8)
          call vtkhdf_write_i8_at(grp_id, "PartOffsets", i8_value, counter)

          i8_value = int(counter - 1, kind=8) * int(total_points, kind=8)
          call vtkhdf_write_i8_at(grp_id, "PointOffsets", i8_value, counter)

          i8_value = int(counter - 1, kind=8) * int(total_cells, kind=8)
          call vtkhdf_write_i8_at(grp_id, "CellOffsets", i8_value, counter)

          i8_value = int(counter - 1, kind=8) * int(total_conn, kind=8)
          call vtkhdf_write_i8_at(grp_id, "ConnectivityIdOffsets", i8_value, &
               counter)

       else
          i8_value = int(0, kind=8)
          call vtkhdf_write_i8_at(grp_id, "PartOffsets", i8_value, counter)
          call vtkhdf_write_i8_at(grp_id, "PointOffsets", i8_value, counter)
          call vtkhdf_write_i8_at(grp_id, "CellOffsets", i8_value, counter)
          call vtkhdf_write_i8_at(grp_id, "ConnectivityIdOffsets", i8_value, &
               counter)
       end if

       call h5gclose_f(grp_id, ierr)
    end if

    call h5pclose_f(xf_id, ierr)

  end subroutine vtkhdf_write_mesh

  !> Write temporal Steps group metadata to the VTKHDF group.
  !! Writes Values, NumberOfParts, PartOffsets, PointOffsets,
  !! CellOffsets, ConnectivityIdOffsets datasets, and NSteps attribute.
  !! @param vtkhdf_grp HDF5 group ID for VTKHDF root group
  !! @param t Current simulation time
  !! @param counter Current counter for how many steps have been written.
  subroutine vtkhdf_write_steps(vtkhdf_grp, counter, t)
    integer(hid_t), intent(in) :: vtkhdf_grp
    integer, intent(in) :: counter
    real(kind=rp), intent(in) :: t

    integer(hid_t) :: xf_id
    integer :: ierr
    integer(hid_t) :: grp_id, dset_id, dcpl_id, filespace, memspace, attr_id
    integer(hsize_t), dimension(1) :: step_dims, step_maxdims
    integer(hsize_t), dimension(1) :: step_count, step_offset, chunkdims, ddim
    real(kind=dp), dimension(1) :: time_value
    integer(kind=8) :: i8_value
    logical :: link_exists, attr_exists

    ! Create collective transfer property list
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

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
       call h5sget_simple_extent_dims_f(filespace, step_dims, step_maxdims, &
            ierr)

       ! We have not written this timestep yet, expand the array
       if (step_dims(1) .eq. counter) then
          step_dims(1) = int(counter + 1, hsize_t)
          call h5dset_extent_f(dset_id, step_dims, ierr)
       else if (step_dims(1) .lt. counter) then
          call neko_error("VTKHDF: Time steps written out of order.")
       end if
    else
       step_dims(1) = 1_hsize_t
       step_maxdims(1) = H5S_UNLIMITED_F
       chunkdims(1) = 1_hsize_t

       call h5screate_simple_f(1, step_dims, filespace, ierr, step_maxdims)
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
       call h5dcreate_f(grp_id, "Values", H5T_NATIVE_DOUBLE, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5pclose_f(dcpl_id, ierr)
    end if

    step_count(1) = 1_hsize_t
    step_offset(1) = int(counter, hsize_t)

    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         step_offset, step_count, ierr)
    call h5screate_simple_f(1, step_count, memspace, ierr)

    time_value(1) = real(t, kind=dp)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time_value, step_count, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)

    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)

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

    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, counter, ddim, ierr)

    call h5aclose_f(attr_id, ierr)
    call h5gclose_f(grp_id, ierr)
    call h5pclose_f(xf_id, ierr)

  end subroutine vtkhdf_write_steps

  !> Append an integer(kind=8) scalar to a 1D chunked HDF5 dataset.
  !! Opens the dataset if it exists, or creates a new empty chunked dataset.
  !! Extends by 1 element and writes the value at the end.
  !! @param grp_id HDF5 group containing the dataset
  !! @param name Dataset name
  !! @param value Value to append
  !! @param counter Position to write to
  subroutine vtkhdf_write_i8_at(grp_id, name, value, counter)
    integer(hid_t), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer(kind=8), intent(in) :: value
    integer, intent(in) :: counter

    integer :: ierr
    integer(hid_t) :: dset_id, dcpl_id, filespace, memspace
    integer(hsize_t), dimension(1) :: dims, maxdims, cnt, off, chunkdims
    integer(kind=8), dimension(1) :: buf
    integer(hid_t) :: xf_id
    logical :: link_exists

    ! Create collective transfer property list
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    call h5lexists_f(grp_id, name, link_exists, ierr)
    if (link_exists) then
       call h5dopen_f(grp_id, name, dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sget_simple_extent_dims_f(filespace, dims, maxdims, ierr)

       if (counter .eq. dims(1)) then
          dims(1) = int(counter + 1, hsize_t)
          call h5dset_extent_f(dset_id, dims, ierr)
       else if (counter .gt. dims(1)) then
          call neko_error("VTKHDF: Values written out of order.")
       end if
    else
       dims(1) = 1_hsize_t
       maxdims(1) = H5S_UNLIMITED_F
       chunkdims(1) = 1_hsize_t

       call h5screate_simple_f(1, dims, filespace, ierr, maxdims)
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
       call h5dcreate_f(grp_id, name, H5T_STD_I64LE, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5pclose_f(dcpl_id, ierr)
    end if

    cnt(1) = 1_hsize_t
    off(1) = int(counter, hsize_t)
    buf(1) = value

    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, off, cnt, ierr)
    call h5screate_simple_f(1, cnt, memspace, ierr)
    call h5dwrite_f(dset_id, H5T_STD_I64LE, buf, cnt, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)

    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(xf_id, ierr)

  end subroutine vtkhdf_write_i8_at

  !> Write field data into the VTKHDF PointData group.
  !! Groups u, v, w fields into a 3-component Velocity vector dataset.
  !! All other fields are written as scalar datasets.
  !! For temporal output, PointDataOffsets are appended under Steps.
  !! @param vtkhdf_grp Root VTKHDF group
  !! @param fp Array of field pointers to write
  !! @param precision Output precision (optional, for VTKHDF files)
  !! @param t Current simulation time (optional, for temporal output)
  subroutine vtkhdf_write_pointdata(vtkhdf_grp, fp, precision, counter, t)
    integer(hid_t), intent(in) :: vtkhdf_grp
    type(field_ptr_t), intent(in) :: fp(:)
    integer, intent(in), optional :: precision
    integer, intent(in) :: counter
    real(kind=rp), intent(in), optional :: t

    logical, allocatable :: field_written(:)
    integer(kind=8) :: time_offset
    integer :: nelv
    integer :: local_points, point_offset
    integer :: lx, ly, lz
    integer :: max_local_points, total_points
    integer(hid_t) :: precision_hdf
    integer(hid_t) :: xf_id
    integer :: ierr, i, j, ie, ii, jj, kk, local_idx
    integer :: n_fields, npts_per_cell
    integer(hid_t) :: pointdata_grp, grp_id, step_grp_id
    integer(hid_t) :: dset_id, dcpl_id, attr_id, filespace, memspace
    integer(hsize_t), dimension(1) :: dcount, doffset, chunkdims
    integer(hsize_t), dimension(2) :: dcount2, doffset2
    integer(hsize_t), dimension(1) :: pd_dims1, pd_maxdims1
    integer(hsize_t), dimension(2) :: pd_dims2, pd_maxdims2
    type(field_t), pointer :: fld, u, v, w
    character(len=128) :: field_name
    logical :: link_exists, is_vector

    ! Create collective transfer property list
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    precision_hdf = vtkhdf_file_determine_real(precision)

    n_fields = size(fp)

    ! --- Build the number of cells and the connectivity
    local_points = fp(1)%ptr%dof%size()
    total_points = fp(1)%ptr%dof%global_size()
    nelv = fp(1)%ptr%msh%nelv
    lx = fp(1)%ptr%dof%Xh%lx
    ly = fp(1)%ptr%dof%Xh%ly
    lz = fp(1)%ptr%dof%Xh%lz

    point_offset = 0
    max_local_points = 0
    call MPI_Exscan(local_points, point_offset, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)
    call MPI_Allreduce(local_points, max_local_points, 1, MPI_INTEGER, &
         MPI_MAX, NEKO_COMM, ierr)
    npts_per_cell = lx * ly * lz

    ! Sync all the fields
    do i = 1, n_fields
       if (associated(fp(i)%ptr)) then
          call fp(i)%ptr%copy_from(DEVICE_TO_HOST, sync = i .eq. n_fields)
       end if
    end do

    allocate(field_written(n_fields), source=.false.)

    ! Create or open PointData group
    call h5lexists_f(vtkhdf_grp, "PointData", link_exists, ierr)
    if (link_exists) then
       call h5gopen_f(vtkhdf_grp, "PointData", pointdata_grp, ierr)
    else
       call h5gcreate_f(vtkhdf_grp, "PointData", pointdata_grp, ierr)
    end if

    ! Read how many steps have been written so far
    if (present(t)) then
       call h5gopen_f(vtkhdf_grp, "Steps", step_grp_id, ierr)
       time_offset = int(counter, kind=8) * int(total_points, kind=8)
    else
       time_offset = 0_8
    end if

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
       if (present(t)) then
          call h5lexists_f(step_grp_id, "PointDataOffsets", link_exists, ierr)
          if (link_exists) then
             call h5gopen_f(step_grp_id, "PointDataOffsets", grp_id, ierr)
          else
             call h5gcreate_f(step_grp_id, "PointDataOffsets", grp_id, ierr)
          end if
          call vtkhdf_write_i8_at(grp_id, trim(field_name), time_offset, &
               counter)
          call h5gclose_f(grp_id, ierr)
       end if

       if (is_vector) then
          ! Create or extend 2D dataset (3 x total_points)
          call h5lexists_f(pointdata_grp, trim(field_name), link_exists, ierr)
          if (link_exists .and. present(t)) then
             call h5dopen_f(pointdata_grp, trim(field_name), dset_id, ierr)
             call h5dget_space_f(dset_id, filespace, ierr)
             call h5sget_simple_extent_dims_f(filespace, pd_dims2, &
                  pd_maxdims2, ierr)

             if (pd_dims2(2) .eq. counter * total_points) then
                pd_dims2(2) = int((counter + 1) * total_points, hsize_t)
                call h5dset_extent_f(dset_id, pd_dims2, ierr)
             else if (pd_dims2(2) .lt. counter * total_points) then
                call neko_error("VTKHDF: Data written out of order.")
             end if

          else if (.not. link_exists) then
             pd_dims2 = [3_hsize_t, int(total_points, hsize_t)]
             pd_maxdims2 = [3_hsize_t, H5S_UNLIMITED_F]
             chunkdims(1) = max(1_hsize_t, int(max_local_points, hsize_t))

             call h5screate_simple_f(2, pd_dims2, filespace, ierr, pd_maxdims2)
             call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
             call h5pset_chunk_f(dcpl_id, 2, [3_hsize_t, chunkdims(1)], ierr)
             call h5dcreate_f(pointdata_grp, trim(field_name), precision_hdf, &
                  filespace, dset_id, ierr, dcpl_id = dcpl_id)
             call h5pclose_f(dcpl_id, ierr)
          end if

          dcount2 = [3_hsize_t, int(local_points, hsize_t)]
          doffset2 = [0_hsize_t, int(time_offset + point_offset, hsize_t)]

          call h5screate_simple_f(2, dcount2, memspace, ierr)
          call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
               doffset2, dcount2, ierr)

          if (precision .eq. sp) then
             call write_vector_single(dset_id, u, v, w, dcount2, ierr, &
                  filespace, memspace, xf_id)
          else if (precision .eq. dp) then
             call write_vector_double(dset_id, u, v, w, dcount2, ierr, &
                  filespace, memspace, xf_id)
          end if

          call h5sclose_f(filespace, ierr)
          call h5sclose_f(memspace, ierr)
          call h5dclose_f(dset_id, ierr)

       else
          ! Write scalar field as 1D dataset
          call h5lexists_f(pointdata_grp, trim(field_name), link_exists, ierr)
          if (link_exists) then
             call h5dopen_f(pointdata_grp, trim(field_name), dset_id, ierr)
             call h5dget_space_f(dset_id, filespace, ierr)
             call h5sget_simple_extent_dims_f(filespace, pd_dims1, &
                  pd_maxdims1, ierr)

             if (pd_dims1(1) .eq. counter * total_points) then
                pd_dims1(1) = int((counter + 1) * total_points, hsize_t)
                call h5dset_extent_f(dset_id, pd_dims1, ierr)
             else if (pd_dims1(1) .lt. counter * total_points) then
                call neko_error("VTKHDF: Data written out of order.")
             end if
          else
             pd_dims1(1) = int(total_points, hsize_t)
             pd_maxdims1(1) = H5S_UNLIMITED_F
             chunkdims(1) = max(1_hsize_t, int(max_local_points, hsize_t))

             call h5screate_simple_f(1, pd_dims1, filespace, ierr, pd_maxdims1)
             call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
             call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
             call h5dcreate_f(pointdata_grp, trim(field_name), precision_hdf, &
                  filespace, dset_id, ierr, dcpl_id = dcpl_id)
             call h5pclose_f(dcpl_id, ierr)
          end if

          dcount(1) = int(local_points, hsize_t)
          doffset(1) = int(time_offset + point_offset, hsize_t)

          call h5screate_simple_f(1, dcount, memspace, ierr)
          call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
               doffset, dcount, ierr)

          if (precision .eq. rp) then
             call h5dwrite_f(dset_id, precision_hdf, fld%x, dcount2, ierr, &
                  file_space_id = filespace, mem_space_id = memspace, &
                  xfer_prp = xf_id)

          else if (precision .eq. sp) then
             block
               real(kind=rp), pointer :: x(:)
               real(kind=sp), allocatable :: data_1d(:)

               x(1:local_points) => fld%x
               allocate(data_1d(local_points))
               do concurrent (local_idx = 1:local_points)
                  data_1d(local_idx) = real(x(local_idx), sp)
               end do

               call h5dwrite_f(dset_id, precision_hdf, data_1d, dcount(1:1), &
                    ierr, file_space_id = filespace, mem_space_id = memspace, &
                    xfer_prp = xf_id)
               deallocate(data_1d)
             end block

          else if (precision .eq. dp) then
             block
               real(kind=rp), pointer :: x(:)
               real(kind=dp), allocatable :: data_1d(:)

               x(1:local_points) => fld%x
               allocate(data_1d(local_points))
               do concurrent (local_idx = 1:local_points)
                  data_1d(local_idx) = real(x(local_idx), dp)
               end do

               call h5dwrite_f(dset_id, precision_hdf, data_1d, dcount(1:1), &
                    ierr, file_space_id = filespace, mem_space_id = memspace, &
                    xfer_prp = xf_id)
               deallocate(data_1d)
             end block

          end if

          call h5sclose_f(filespace, ierr)
          call h5sclose_f(memspace, ierr)
          call h5dclose_f(dset_id, ierr)
       end if
    end do

    if (present(t)) call h5gclose_f(step_grp_id, ierr)
    call h5gclose_f(pointdata_grp, ierr)

  end subroutine vtkhdf_write_pointdata

  ! -------------------------------------------------------------------------- !
  ! Helper functions and routines

  !> Write a 3-component vector field as a 2D dataset with single precision,
  !! converting from double if necessary.
  !! The dataset is organized as (3, total_points) for VTK compatibility.
  subroutine write_vector_single(dset_id, u, v, w, dcount2, ierr, filespace, &
       memspace, xf_id)
    integer(hid_t), intent(in) :: dset_id
    type(field_t), intent(in), pointer :: u, v, w
    integer(hsize_t), intent(in) :: dcount2(:)
    integer, intent(inout) :: ierr
    integer(hid_t), intent(inout) :: filespace, memspace, xf_id
    real(kind=sp), allocatable :: data_2d(:,:)
    integer :: ie, kk, ii, jj, local_idx
    integer :: lx, ly, lz, nelv, npts_per_cell, local_points

    lx = u%dof%Xh%lx
    ly = u%dof%Xh%ly
    lz = u%dof%Xh%lz
    nelv = u%msh%nelv
    local_points = lx * ly * lz * nelv
    npts_per_cell = lx * ly * lz

    ! Assemble 3-component vector from u, v, w fields
    allocate(data_2d(3, local_points))
    do ie = 1, nelv
       local_idx = (ie - 1) * npts_per_cell
       do kk = 1, lz
          do jj = 1, ly
             do ii = 1, lx
                local_idx = local_idx + 1
                data_2d(1, local_idx) = real(u%x(ii, jj, kk, ie), sp)
                data_2d(2, local_idx) = real(v%x(ii, jj, kk, ie), sp)
                if (associated(w)) then
                   data_2d(3, local_idx) = real(w%x(ii, jj, kk, ie), sp)
                else
                   data_2d(3, local_idx) = 0.0_sp
                end if
             end do
          end do
       end do
    end do

    call h5dwrite_f(dset_id, H5T_NATIVE_REAL, data_2d, dcount2, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = xf_id)

    deallocate(data_2d)
  end subroutine write_vector_single

  !> Write a 3-component vector field as a 2D dataset with double precision,
  !! converting from single if necessary.
  !! The dataset is organized as (3, total_points) for VTK compatibility.
  subroutine write_vector_double(dset_id, u, v, w, dcount2, ierr, filespace, &
       memspace, xf_id)
    integer(hid_t), intent(in) :: dset_id
    type(field_t), intent(in), pointer :: u, v, w
    integer(hsize_t), intent(in) :: dcount2(:)
    integer, intent(inout) :: ierr
    integer(hid_t), intent(inout) :: filespace, memspace, xf_id
    real(kind=dp), allocatable :: data_2d(:,:)
    integer :: ie, kk, ii, jj, local_idx
    integer :: lx, ly, lz, nelv, npts_per_cell, local_points

    lx = u%dof%Xh%lx
    ly = u%dof%Xh%ly
    lz = u%dof%Xh%lz
    nelv = u%msh%nelv
    local_points = lx * ly * lz * nelv
    npts_per_cell = lx * ly * lz

    ! Assemble 3-component vector from u, v, w fields
    allocate(data_2d(3, local_points))
    do ie = 1, nelv
       local_idx = (ie - 1) * npts_per_cell
       do kk = 1, lz
          do jj = 1, ly
             do ii = 1, lx
                local_idx = local_idx + 1
                data_2d(1, local_idx) = real(u%x(ii, jj, kk, ie), dp)
                data_2d(2, local_idx) = real(v%x(ii, jj, kk, ie), dp)
                if (associated(w)) then
                   data_2d(3, local_idx) = real(w%x(ii, jj, kk, ie), dp)
                else
                   data_2d(3, local_idx) = 0.0_dp
                end if
             end do
          end do
       end do
    end do

    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_2d, dcount2, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = xf_id)

    deallocate(data_2d)
  end subroutine write_vector_double

  !> Read data in HDF5 format following official VTKHDF specification
  subroutine vtkhdf_file_read(this, data)
    class(vtkhdf_file_t) :: this
    class(*), target, intent(inout) :: data

    call neko_error('VTKHDF file reading is not yet implemented')

  end subroutine vtkhdf_file_read

  !> Determine hdf5 real type corresponding to NEKO_REAL
  !! @note This must be called after h5open_f, otherwise
  !! the H5T_NATIVE_XYZ types has a value of 0
  function vtkhdf_file_determine_real(precision) result(H5T_NEKO_REAL)
    integer, intent(in) :: precision
    integer(hid_t) :: H5T_NEKO_REAL

    select case(precision)
    case(sp)
       H5T_NEKO_REAL = H5T_NATIVE_REAL
    case(dp)
       H5T_NEKO_REAL = H5T_NATIVE_DOUBLE
    case default
       call neko_error("Unsupported real type")
    end select
  end function vtkhdf_file_determine_real

  !> Set the precision for VTKHDF output (single or double)
  subroutine vtkhdf_file_set_precision(this, precision)
    class(vtkhdf_file_t), intent(inout) :: this
    integer, intent(in) :: precision
    this%precision = precision
  end subroutine vtkhdf_file_set_precision


#else
  ! -------------------------------------------------------------------------- !
  ! Dummy functions and subroutines

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

end module vtkhdf_file
