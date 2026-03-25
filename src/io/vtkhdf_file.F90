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
  use num_types, only : rp, sp, dp, qp, i8
  use generic_file, only : generic_file_t
  use checkpoint, only : chkp_t
  use utils, only : neko_error, neko_warning, filename_split, &
       nonlinear_index, linear_index
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
  use vtk, only : vtk_ordering
#ifdef HAVE_HDF5
  use hdf5, only : &
       hid_t, hsize_t, &
       h5open_f, h5close_f, &
       h5fcreate_f, h5fopen_f, h5fclose_f, &
       h5gcreate_f, h5gopen_f, h5gclose_f, &
       h5acreate_f, h5awrite_f, h5aclose_f, &
       h5dcreate_f, h5dopen_f, h5dwrite_f, h5dclose_f, &
       h5screate_f, h5screate_simple_f, h5sclose_f, h5sselect_hyperslab_f, &
       h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f, h5pset_dxpl_mpio_f, &
       h5lexists_f, h5tset_strpad_f, h5tset_size_f, h5tcopy_f, &
       H5P_FILE_ACCESS_F, H5P_DATASET_XFER_F, H5P_DATASET_CREATE_F, &
       H5F_ACC_TRUNC_F, H5F_ACC_RDWR_F, &
       H5T_STD_U8LE, H5T_NATIVE_INTEGER, H5T_FORTRAN_S1, H5T_STR_NULLTERM_F, &
       h5kind_to_type, H5_REAL_KIND, H5_INTEGER_KIND, &
       H5S_SCALAR_F, H5S_SELECT_SET_F, H5FD_MPIO_COLLECTIVE_F, h5p_default_f, &
       H5S_UNLIMITED_F
#endif
  implicit none
  private

  !> Interface for HDF5 files
  type, public, extends(generic_file_t) :: vtkhdf_file_t
     logical :: amr_enabled = .false.
     logical :: subdivide = .false.
     integer :: precision = -1
   contains
     procedure :: get_vtkhdf_fname => vtkhdf_file_get_fname
     procedure :: read => vtkhdf_file_read
     procedure :: write => vtkhdf_file_write
     procedure :: set_overwrite => vtkhdf_file_set_overwrite
     procedure :: enable_amr => vtkhdf_file_enable_amr
     procedure :: set_precision => vtkhdf_file_set_precision
     procedure :: set_subdivide => vtkhdf_file_set_subdivide
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

  !> Set the precision for VTKHDF output (single or double)
  subroutine vtkhdf_file_set_precision(this, precision)
    class(vtkhdf_file_t), intent(inout) :: this
    integer, intent(in) :: precision
    this%precision = precision
  end subroutine vtkhdf_file_set_precision

  !> Return the file name with the start counter.
  function vtkhdf_file_get_fname(this) result(base_fname)
    class(vtkhdf_file_t), intent(in) :: this
    character(len=1024) :: base_fname
    character(len=1024) :: fname
    character(len=1024) :: path, name, suffix

    fname = trim(this%get_base_fname())
    call filename_split(fname, path, name, suffix)

    write(base_fname, '(A,A,"_",I0,A)') &
         trim(path), trim(name), this%get_start_counter(), trim(suffix)

  end function vtkhdf_file_get_fname

  !> Enable or disable subdivision of spectral elements into linear sub-cells.
  !! When subdivision is enabled, each spectral element is written as multiple
  !! linear VTK cells (VTK_HEXAHEDRON in 3D, VTK_QUAD in 2D) with connectivity
  !! corresponding to the tensor-product grid of the spectral element.
  !! @param subdivide Whether to subdivide into linear sub-cells.
  subroutine vtkhdf_file_set_subdivide(this, subdivide)
    class(vtkhdf_file_t), intent(inout) :: this
    logical, intent(in) :: subdivide
    this%subdivide = subdivide
  end subroutine vtkhdf_file_set_subdivide

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
    logical :: exists
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
    if (msh%gdim .lt. 2 .or. msh%gdim .gt. 3) then
       call neko_error('VTKHDF output only supports 2D and 3D meshes')
    end if

    ! Ensure precision is set and are valid.
    if (this%precision .gt. rp) then
       this%precision = rp
       call neko_warning('Requested precision is higher than working precision')
    else if (this%precision .eq. -1) then
       this%precision = rp
    end if

    call this%increment_counter()
    fname = trim(this%get_vtkhdf_fname())
    counter = this%get_counter() - this%get_start_counter()

    mpi_info = MPI_INFO_NULL%mpi_val
    mpi_comm = NEKO_COMM%mpi_val

    call h5open_f(ierr)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    call h5pset_fapl_mpio_f(plist_id, mpi_comm, mpi_info, ierr)

    if (counter .eq. 0) then
       ! First write: always create a fresh file to avoid stale data
       call h5fcreate_f(fname, H5F_ACC_TRUNC_F, &
            file_id, ierr, access_prp = plist_id)
    else
       call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, ierr, &
            access_prp = plist_id)
    end if

    ! Create/open VTKHDF root group with vtkhdf_version and type attributes
    call h5lexists_f(file_id, "VTKHDF", exists, ierr)
    if (exists) then
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
       call vtkhdf_write_mesh(vtkhdf_grp, dof, msh, &
            this%amr_enabled, counter, this%subdivide, t)
    end if

    ! Write field data in PointData group
    if (n_fields > 0) then
       call vtkhdf_write_pointdata(vtkhdf_grp, fp, this%precision, counter, &
            fname, t)
    end if

    call h5gclose_f(vtkhdf_grp, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)

    if (allocated(fp)) deallocate(fp)

  end subroutine vtkhdf_file_write

  !> Write mesh geometry datasets to the VTKHDF group.
  !! Writes NumberOfPoints, NumberOfCells, NumberOfConnectivityIds,
  !! Points, Connectivity, Offsets, and Types datasets.
  !! @param vtkhdf_grp HDF5 group ID for VTKHDF root group
  !! @param dof Dofmap for coordinate data
  !! @param msh Mesh object
  !! @param amr AMR flag to determine if mesh should be rewritten at every time
  !!            step
  !! @param t Optional time value for time-dependent mesh output (e.g. for AMR)
  subroutine vtkhdf_write_mesh(vtkhdf_grp, dof, msh, amr, &
       counter, subdivide, t)
    type(dofmap_t), intent(in) :: dof
    type(mesh_t), intent(in) :: msh
    integer(hid_t), intent(in) :: vtkhdf_grp
    logical, intent(in) :: amr
    integer, intent(in) :: counter
    logical, intent(in) :: subdivide
    real(kind=rp), intent(in), optional :: t

    integer(kind=1) :: VTK_cell_type
    integer :: ierr, i, ii, jj, kk, el, local_idx
    integer :: lx, ly, lz, npts_per_cell, nodes_per_cell, cells_per_element
    integer :: local_points, local_cells, local_conn
    integer :: total_points, total_cells, total_conn
    integer :: point_offset, max_local_points
    integer :: total_offsets, cell_offset, conn_offset, offsets_offset
    integer :: max_local_cells, max_local_conn
    integer(hid_t) :: xf_id, dset_id, dcpl_id, grp_id, attr_id
    integer(hid_t) :: filespace, memspace, H5T_NEKO_DOUBLE, H5T_NEKO_INT8
    integer(hsize_t), dimension(1) :: dcount, vdims, maxdims, doffset, chunkdims
    integer(hsize_t), dimension(2) :: dcount2, vdims2, maxdims2, doffset2
    integer(kind=i8) :: i8_value
    logical :: exists
    integer, dimension(3) :: component_sizes
    integer, dimension(3) :: component_offsets
    integer, dimension(3) :: component_max_sizes

    lx = dof%Xh%lx
    ly = dof%Xh%ly
    lz = dof%Xh%lz

    if (subdivide .and. msh%gdim .eq. 3) then
       VTK_cell_type = 12 ! VTK_HEXAHEDRON
       cells_per_element = (lx - 1) * (ly - 1) * (lz - 1)
       nodes_per_cell = 8
    else if (subdivide .and. msh%gdim .eq. 2) then
       VTK_cell_type = 9 ! VTK_QUAD
       cells_per_element = (lx - 1) * (ly - 1)
       nodes_per_cell = 4
    else if (msh%gdim .eq. 3) then
       VTK_cell_type = 72 ! VTK_LAGRANGE_HEXAHEDRON
       cells_per_element = 1
       nodes_per_cell = lx * ly * lz
    else if (msh%gdim .eq. 2) then
       VTK_cell_type = 70 ! VTK_LAGRANGE_QUADRILATERAL
       cells_per_element = 1
       nodes_per_cell = lx * ly
    end if

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

    ! --- NumberOfPoints, NumberOfCells, NumberOfConnectivityIds ---
    ! These datasets must accumulate nPieces entries per timestep,
    ! giving a total size of nSteps * nPieces. VTK's reader computes
    ! numberOfPieces = dims[0] / nSteps, so missing entries cause
    ! garbage reads.
    block
      integer(hsize_t), dimension(1) :: nof_dims, nof_maxdims
      integer(hsize_t), dimension(1) :: nof_count, nof_offset, nof_chunk
      integer(hid_t) :: nof_filespace, nof_memspace, nof_dcpl
      integer(kind=i8) :: nof_values(3)

      nof_count(1) = 1_hsize_t
      nof_offset(1) = int(counter, hsize_t) * int(pe_size, hsize_t) &
           + int(pe_rank, hsize_t)
      nof_chunk(1) = max(1_hsize_t, int(pe_size, hsize_t))
      nof_values = [int(local_points, kind=i8), int(local_cells, kind=i8), &
           int(local_conn, kind=i8)]

      call h5pcreate_f(H5P_DATASET_CREATE_F, nof_dcpl, ierr)
      call h5pset_chunk_f(nof_dcpl, 1, nof_chunk, ierr)

      call vtkhdf_write_numberof(vtkhdf_grp, "NumberOfPoints", &
           nof_values(1), nof_offset, nof_count, nof_dcpl, &
           counter, xf_id, ierr)
      call vtkhdf_write_numberof(vtkhdf_grp, "NumberOfCells", &
           nof_values(2), nof_offset, nof_count, nof_dcpl, &
           counter, xf_id, ierr)
      call vtkhdf_write_numberof(vtkhdf_grp, "NumberOfConnectivityIds", &
           nof_values(3), nof_offset, nof_count, nof_dcpl, &
           counter, xf_id, ierr)

      call h5pclose_f(nof_dcpl, ierr)
    end block

    ! --- Points dataset (global coordinates) ---
    call h5lexists_f(vtkhdf_grp, "Points", exists, ierr)
    if (.not. exists) then

       vdims2 = [3_hsize_t, int(total_points, hsize_t)]
       maxdims2 = [3_hsize_t, H5S_UNLIMITED_F]
       chunkdims(1) = int(max(1, min(max_local_points, total_points)), hsize_t)
       dcount2 = [3_hsize_t, int(local_points, hsize_t)]
       doffset2 = [0_hsize_t, int(point_offset, hsize_t)]
       H5T_NEKO_DOUBLE = h5kind_to_type(dp, H5_REAL_KIND)

       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       call h5screate_simple_f(2, dcount2, memspace, ierr)
       call h5screate_simple_f(2, vdims2, filespace, ierr, maxdims2)

       call h5pset_chunk_f(dcpl_id, 2, [3_hsize_t, chunkdims(1)], ierr)
       call h5dcreate_f(vtkhdf_grp, "Points", H5T_NEKO_DOUBLE, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
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
         call h5dwrite_f(dset_id, H5T_NEKO_DOUBLE, coords, dcount2, ierr, &
              file_space_id = filespace, mem_space_id = memspace, &
              xfer_prp = xf_id)
         deallocate(coords)
       end block

       call h5dclose_f(dset_id, ierr)
       call h5sclose_f(filespace, ierr)
       call h5sclose_f(memspace, ierr)
       call h5pclose_f(dcpl_id, ierr)
    end if

    ! --- Connectivity dataset ---
    call h5lexists_f(vtkhdf_grp, "Connectivity", exists, ierr)
    if (exists) call h5ldelete_f(vtkhdf_grp, "Connectivity", ierr)

    vdims = int(total_conn, hsize_t)
    maxdims = H5S_UNLIMITED_F
    chunkdims = int(max(1, min(max_local_conn, total_conn)), hsize_t)
    dcount = int(local_conn, hsize_t)
    doffset = int(conn_offset, hsize_t)

    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    call h5screate_simple_f(1, dcount, memspace, ierr)
    call h5screate_simple_f(1, vdims, filespace, ierr, maxdims)

    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
    call h5dcreate_f(vtkhdf_grp, "Connectivity", H5T_NATIVE_INTEGER, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset, dcount, ierr)

    block
      integer, allocatable :: connectivity(:)

      allocate(connectivity(local_conn))
      call vtkhdf_build_connectivity(connectivity, VTK_cell_type, msh, dof, &
           subdivide)
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, connectivity, dcount, &
           ierr, file_space_id = filespace, mem_space_id = memspace, &
           xfer_prp = xf_id)
      deallocate(connectivity)
    end block

    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)
    call h5pclose_f(dcpl_id, ierr)

    ! --- Offsets dataset ---
    call h5lexists_f(vtkhdf_grp, "Offsets", exists, ierr)
    if (exists) call h5ldelete_f(vtkhdf_grp, "Offsets", ierr)

    vdims = int(total_offsets, hsize_t)
    maxdims = H5S_UNLIMITED_F
    chunkdims = int(max(1, min(max_local_cells + 1, total_offsets)), hsize_t)
    dcount = int(local_cells + 1, hsize_t)
    doffset = int(offsets_offset, hsize_t)
    H5T_NEKO_INT8 = h5kind_to_type(i8, H5_INTEGER_KIND)

    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    call h5screate_simple_f(1, dcount, memspace, ierr)
    call h5screate_simple_f(1, vdims, filespace, ierr, maxdims)

    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
    call h5dcreate_f(vtkhdf_grp, "Offsets", H5T_NEKO_INT8, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset, dcount, ierr)

    block
      integer(kind=i8), allocatable :: offsets(:)

      allocate(offsets(local_cells + 1))
      do concurrent (i = 1:local_cells)
         offsets(i) = int((i - 1) * nodes_per_cell, kind=i8)
      end do
      offsets(local_cells + 1) = int(local_conn, kind=i8)
      call h5dwrite_f(dset_id, H5T_NEKO_INT8, offsets, dcount, ierr, &
           file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
      deallocate(offsets)
    end block

    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)
    call h5pclose_f(dcpl_id, ierr)

    ! --- Types dataset (VTK cell types) ---
    call h5lexists_f(vtkhdf_grp, "Types", exists, ierr)
    if (exists) call h5ldelete_f(vtkhdf_grp, "Types", ierr)

    vdims = int(total_cells, hsize_t)
    maxdims = H5S_UNLIMITED_F
    chunkdims = int(max(1, min(max_local_cells, total_cells)), hsize_t)
    dcount = int(local_cells, hsize_t)
    doffset = int(cell_offset, hsize_t)

    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    call h5screate_simple_f(1, dcount, memspace, ierr)
    call h5screate_simple_f(1, vdims, filespace, ierr, maxdims)

    call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
    call h5dcreate_f(vtkhdf_grp, "Types", H5T_STD_U8LE, &
         filespace, dset_id, ierr, dcpl_id = dcpl_id)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset, dcount, ierr)

    block
      integer(kind=1), allocatable :: cell_types(:)
      allocate(cell_types(local_cells), source=VTK_cell_type)
      call h5dwrite_f(dset_id, H5T_STD_U8LE, cell_types, dcount, ierr, &
           file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
      deallocate(cell_types)
    end block

    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)
    call h5pclose_f(dcpl_id, ierr)

    if (present(t)) then
       ! Open Steps group
       call h5gopen_f(vtkhdf_grp, "Steps", grp_id, ierr)

       ! --- NumberOfParts ---
       call vtkhdf_write_i8_at(grp_id, "NumberOfParts", int(pe_size, kind=i8), &
            counter)

       ! --- PartOffsets ---
       i8_value = 0
       if (amr) then
          i8_value = int(counter - 1, kind=i8) * int(pe_size, kind=i8)
          call vtkhdf_write_i8_at(grp_id, "PartOffsets", i8_value, counter)

          i8_value = int(counter - 1, kind=i8) * int(total_points, kind=i8)
          call vtkhdf_write_i8_at(grp_id, "PointOffsets", i8_value, counter)

          i8_value = int(counter - 1, kind=i8) * int(total_cells, kind=i8)
          call vtkhdf_write_i8_at(grp_id, "CellOffsets", i8_value, counter)

          i8_value = int(counter - 1, kind=i8) * int(total_conn, kind=i8)
          call vtkhdf_write_i8_at(grp_id, "ConnectivityIdOffsets", i8_value, &
               counter)

       else
          i8_value = 0_i8
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

    integer(hid_t) :: xf_id, H5T_NEKO_DOUBLE
    integer :: ierr
    integer(hid_t) :: grp_id, dset_id, dcpl_id, filespace, memspace, attr_id
    integer(hsize_t), dimension(1) :: step_dims, step_maxdims
    integer(hsize_t), dimension(1) :: step_count, step_offset, chunkdims, ddim
    real(kind=dp), dimension(1) :: time_value
    integer(kind=i8) :: i8_value
    logical :: exists, attr_exists

    ! Create collective transfer property list
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    H5T_NEKO_DOUBLE = h5kind_to_type(dp, H5_REAL_KIND)

    ! Create or open Steps group
    call h5lexists_f(vtkhdf_grp, "Steps", exists, ierr)
    if (exists) then
       call h5gopen_f(vtkhdf_grp, "Steps", grp_id, ierr)
    else
       call h5gcreate_f(vtkhdf_grp, "Steps", grp_id, ierr)
    end if

    ! --- Values dataset (time values, real type) ---
    call h5lexists_f(grp_id, "Values", exists, ierr)
    if (exists) then
       call h5dopen_f(grp_id, "Values", dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sget_simple_extent_dims_f(filespace, step_dims, step_maxdims, &
            ierr)
       call h5sclose_f(filespace, ierr)

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
       call h5dcreate_f(grp_id, "Values", H5T_NEKO_DOUBLE, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5sclose_f(filespace, ierr)
       call h5pclose_f(dcpl_id, ierr)
    end if

    step_count(1) = 1_hsize_t
    step_offset(1) = int(counter, hsize_t)

    call h5dget_space_f(dset_id, filespace, ierr)
    call h5screate_simple_f(1, step_count, memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         step_offset, step_count, ierr)

    time_value(1) = real(t, kind=dp)
    call h5dwrite_f(dset_id, H5T_NEKO_DOUBLE, time_value, step_count, ierr, &
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

    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, counter + 1, ddim, ierr)

    call h5aclose_f(attr_id, ierr)
    call h5gclose_f(grp_id, ierr)
    call h5pclose_f(xf_id, ierr)

  end subroutine vtkhdf_write_steps

  !> Write field data into the VTKHDF PointData group.
  !! Groups u, v, w fields into a 3-component Velocity vector dataset.
  !! All other fields are written as scalar datasets.
  !! For temporal output (when t is present), each timestep's data is stored
  !! in a separate external HDF5 file and the main file uses HDF5 Virtual
  !! Datasets (VDS) with printf-style source file patterns to map them.
  !! For non-temporal output (when t is absent), fields are written directly
  !! into the main file's PointData group as regular datasets.
  !! @param vtkhdf_grp Root VTKHDF group
  !! @param fp Array of field pointers to write
  !! @param precision Output precision
  !! @param counter Current timestep counter
  !! @param fname Main VTKHDF file path (used to derive external file names)
  !! @param t Current simulation time (optional, for temporal output)
  subroutine vtkhdf_write_pointdata(vtkhdf_grp, fp, precision, counter, &
       fname, t)
    integer(hid_t), intent(in) :: vtkhdf_grp
    type(field_ptr_t), intent(in) :: fp(:)
    integer, intent(in) :: precision
    integer, intent(in) :: counter
    character(len=*), intent(in) :: fname
    real(kind=rp), intent(in), optional :: t

    integer(kind=i8) :: time_offset
    integer :: local_points, point_offset, total_points
    integer(hid_t) :: precision_hdf
    integer :: ierr, i, j
    integer :: n_fields
    integer(hid_t) :: pointdata_grp, grp_id, step_grp_id
    integer(hid_t) :: dset_id, dcpl_id, filespace
    integer(hsize_t), dimension(1) :: pd_dims1, pd_maxdims1
    integer(hsize_t), dimension(2) :: pd_dims2, pd_maxdims2
    type(field_t), pointer :: fld, u, v, w
    character(len=128) :: field_name
    logical :: exists, is_vector

    ! VDS and per-timestep external file variables
    character(len=1024) :: ext_fname, ext_path, src_pattern
    character(len=1024) :: main_path, main_name, main_suffix
    integer(hid_t) :: ext_file_id, ext_plist_id, vds_src_space
    integer(hid_t) :: write_target
    integer :: mpi_info, mpi_comm

    ! Collected field info for VDS phase
    integer :: fields_written
    character(len=128), allocatable :: name_list(:)
    logical, allocatable :: vector_list(:)

    mpi_info = MPI_INFO_NULL%mpi_val
    mpi_comm = NEKO_COMM%mpi_val

    n_fields = size(fp)

    ! Compute local/global point counts and MPI offsets
    local_points = fp(1)%ptr%dof%size()
    total_points = fp(1)%ptr%dof%global_size()
    point_offset = 0
    call MPI_Exscan(local_points, point_offset, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)

    ! Sync all the fields
    do i = 1, n_fields
       if (associated(fp(i)%ptr)) then
          call fp(i)%ptr%copy_from(DEVICE_TO_HOST, sync = i .eq. n_fields)
       end if
    end do

    fields_written = 0
    allocate(name_list(n_fields))
    allocate(vector_list(n_fields))

    ! Create PointData group if missing
    call h5lexists_f(vtkhdf_grp, "PointData", exists, ierr)
    if (.not. exists) then
       call h5gcreate_f(vtkhdf_grp, "PointData", pointdata_grp, ierr)
       call h5gclose_f(pointdata_grp, ierr)
    end if

    ! ------------------------------------------------------------------------ !
    ! Construct the target where data is written

    if (present(t)) then

       ! Derive base path from main filename for external files
       call filename_split(fname, main_path, main_name, main_suffix)
       write(ext_path, '(A,A,".data/")') trim(main_path), trim(main_name)
       write(ext_fname, '(A,I0,".h5")') trim(ext_path), counter
       write(src_pattern, '(A,".data/%b.h5")') trim(main_name)

       if (pe_rank == 0) then
          inquire(file = trim(ext_path), exist = exists)
          if (.not. exists) then
             call execute_command_line("mkdir -p '" // trim(ext_path) // "'")
          end if
       end if
       call MPI_Barrier(NEKO_COMM, ierr)

       call h5pcreate_f(H5P_FILE_ACCESS_F, ext_plist_id, ierr)
       call h5pset_fapl_mpio_f(ext_plist_id, mpi_comm, mpi_info, ierr)
       call h5fcreate_f(trim(ext_fname), H5F_ACC_TRUNC_F, write_target, ierr, &
            access_prp = ext_plist_id)
       call h5pclose_f(ext_plist_id, ierr)

    else
       ! Non-temporal: write directly into the main file's PointData group
       call h5gopen_f(vtkhdf_grp, "PointData", write_target, ierr)
    end if

    ! ------------------------------------------------------------------------ !
    ! Write field data

    do i = 1, n_fields
       fld => fp(i)%ptr
       field_name = fld%name
       if (field_name .eq. 'p') field_name = 'Pressure'

       ! Determine if this is a velocity component to group as a vector
       is_vector = .false.
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
             case ('v')
                v => fp(j)%ptr
             case ('w')
                w => fp(j)%ptr
             end select
          end do

          if (associated(u) .and. associated(v) .and. associated(w)) then
             is_vector = .true.
             field_name = 'Velocity'
          end if
       end if

       ! Skip duplicate fields (e.g. fluid_rho added by both fluid output
       ! and field_writer when sharing the same output file)
       exists = .false.
       do j = 1, fields_written
          if (trim(name_list(j)) .eq. trim(field_name)) then
             exists = .true.
             exit
          end if
       end do

       ! Skip duplicate fields
       if (exists) cycle

       ! Track unique field names for VDS phase
       fields_written = fields_written + 1
       name_list(fields_written) = field_name
       vector_list(fields_written) = is_vector

       ! Write field data to the target
       if (is_vector) then
          call write_vector_field(write_target, field_name, u%x, v%x, w%x, &
               local_points, precision, total_points, point_offset)
       else
          call write_scalar_field(write_target, field_name, fld%x, &
               local_points, precision, total_points, point_offset)
       end if
    end do

    ! Close write target
    if (present(t)) then
       call h5fclose_f(write_target, ierr)
    else
       call h5gclose_f(write_target, ierr)
    end if

    ! ------------------------------------------------------------------------ !
    ! Manage temporal datasets through VDS

    if (present(t)) then

       ! Temporal: set up per-timestep external file as write target
       call h5gopen_f(vtkhdf_grp, "Steps", step_grp_id, ierr)
       time_offset = int(counter, kind=i8) * int(total_points, kind=i8)

       ! Write PointDataOffsets under Steps for each unique field
       call h5lexists_f(step_grp_id, "PointDataOffsets", exists, ierr)
       if (exists) then
          call h5gopen_f(step_grp_id, "PointDataOffsets", grp_id, ierr)
       else
          call h5gcreate_f(step_grp_id, "PointDataOffsets", grp_id, ierr)
       end if
       do i = 1, fields_written
          call vtkhdf_write_i8_at(grp_id, trim(name_list(i)), &
               time_offset, counter)
       end do
       call h5gclose_f(grp_id, ierr)
       call h5gclose_f(step_grp_id, ierr)

       ! ===== Create or extend VDS in main file =====
       call h5gopen_f(vtkhdf_grp, "PointData", pointdata_grp, ierr)

       do i = 1, fields_written
          field_name = name_list(i)
          is_vector = vector_list(i)

          if (counter .eq. 0) then
             ! First write: create VDS with pattern-based mapping
             call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
             precision_hdf = h5kind_to_type(precision, H5_REAL_KIND)

             if (is_vector) then
                pd_dims2 = [3_hsize_t, int(total_points, hsize_t)]
                call h5screate_simple_f(2, pd_dims2, vds_src_space, ierr)
                call h5sselect_all_f(vds_src_space, ierr)

                pd_maxdims2 = [3_hsize_t, H5S_UNLIMITED_F]
                call h5screate_simple_f(2, pd_dims2, filespace, ierr, &
                     pd_maxdims2)

                call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                     [0_hsize_t, 0_hsize_t], &
                     [1_hsize_t, H5S_UNLIMITED_F], &
                     ierr, &
                     stride = [3_hsize_t, int(total_points, hsize_t)], &
                     block = [3_hsize_t, int(total_points, hsize_t)])

                call h5pset_virtual_f(dcpl_id, filespace, trim(src_pattern), &
                     trim(field_name), vds_src_space, ierr)
                call h5sclose_f(vds_src_space, ierr)

                call h5dcreate_f(pointdata_grp, trim(field_name), &
                     precision_hdf, filespace, dset_id, ierr, &
                     dcpl_id = dcpl_id)
                call h5sclose_f(filespace, ierr)
             else
                pd_dims1 = int(total_points, hsize_t)
                call h5screate_simple_f(1, pd_dims1, vds_src_space, ierr)
                call h5sselect_all_f(vds_src_space, ierr)

                pd_maxdims1(1) = H5S_UNLIMITED_F
                call h5screate_simple_f(1, pd_dims1, filespace, ierr, &
                     pd_maxdims1)

                call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                     [0_hsize_t], &
                     [H5S_UNLIMITED_F], &
                     ierr, &
                     stride = [int(total_points, hsize_t)], &
                     block = [int(total_points, hsize_t)])

                call h5pset_virtual_f(dcpl_id, filespace, trim(src_pattern), &
                     trim(field_name), vds_src_space, ierr)
                call h5sclose_f(vds_src_space, ierr)

                call h5dcreate_f(pointdata_grp, trim(field_name), &
                     precision_hdf, filespace, dset_id, ierr, &
                     dcpl_id = dcpl_id)
                call h5sclose_f(filespace, ierr)
             end if

             call h5pclose_f(dcpl_id, ierr)
             call h5dclose_f(dset_id, ierr)

          else
             ! Subsequent write: extend VDS to include new timestep
             call h5dopen_f(pointdata_grp, trim(field_name), dset_id, ierr)

             if (is_vector) then
                pd_dims2 = [3_hsize_t, &
                     int((counter + 1) * total_points, hsize_t)]
                call h5dset_extent_f(dset_id, pd_dims2, ierr)
             else
                pd_dims1 = int((counter + 1) * total_points, hsize_t)
                call h5dset_extent_f(dset_id, pd_dims1, ierr)
             end if

             call h5dclose_f(dset_id, ierr)
          end if
       end do

       call h5gclose_f(pointdata_grp, ierr)
    end if

    ! ------------------------------------------------------------------------ !
    ! Cleanup before returning

    deallocate(name_list)
    deallocate(vector_list)

  end subroutine vtkhdf_write_pointdata

  ! -------------------------------------------------------------------------- !
  ! Helper functions and routines

  !> Build local connectivity for VTK cells from a spectral element
  !! tensor-product grid. For linear types (12, 9) subdivides each spectral
  !! element into linear sub-cells. For Lagrange types (72, 70) writes one
  !! high-order cell per element with VTK node ordering.
  !! @param conn Output connectivity array (pre-allocated)
  !! @param vtk_type VTK cell type: 12 = VTK_HEXAHEDRON, 9 = VTK_QUAD,
  !!        72 = VTK_LAGRANGE_HEXAHEDRON, 70 = VTK_LAGRANGE_QUADRILATERAL
  !! @param msh Mesh object containing element information
  !! @param dof Dofmap containing the lx, ly, lz dimensions of the spectral
  !!            element grid
  !! @param subdivide Logical flag indicating whether to subdivide to linear
  !!        elements
  subroutine vtkhdf_build_connectivity(conn, vtk_type, msh, dof, subdivide)
    integer, intent(inout) :: conn(:)
    integer(kind=1), intent(in) :: vtk_type
    type(mesh_t), intent(in) :: msh
    type(dofmap_t), intent(in) :: dof
    logical, intent(in) :: subdivide
    integer :: lx, ly, lz, nelv
    integer :: ie, ii, n_pts_per_elem, n_conn_per_elem
    integer, allocatable :: node_order(:)

    nelv = msh%nelv
    lx = dof%Xh%lx
    ly = dof%Xh%ly
    lz = dof%Xh%lz
    n_pts_per_elem = lx * ly * lz

    if (subdivide .and. vtk_type .eq. 12) then
       node_order = subdivide_to_hex_ordering(lx, ly, lz)
    else if (subdivide .and. vtk_type .eq. 9) then
       node_order = subdivide_to_quad_ordering(lx, ly)
    else
       node_order = vtk_ordering(vtk_type, lx, ly, lz)
    end if

    n_conn_per_elem = size(node_order)

    do concurrent (ie = 1:nelv, ii = 1:n_conn_per_elem)
       block
         integer :: idx, base
         idx = (ie - 1) * n_conn_per_elem
         base = (ie - 1) * n_pts_per_elem
         conn(idx + ii) = base + node_order(ii)
       end block
    end do

    deallocate(node_order)

  end subroutine vtkhdf_build_connectivity

  !> Create-or-extend a 1D NumberOf* dataset, then write one entry.
  !! On step 0 the dataset is created with unlimited max extent.
  !! On subsequent steps it is extended to (counter+1)*nPieces.
  !! @param grp        Parent HDF5 group (VTKHDF root)
  !! @param dset_name  Dataset name, e.g. "NumberOfPoints"
  !! @param value      The integer value to write for this rank
  !! @param offset     1-element array: counter*pe_size + pe_rank
  !! @param cnt        1-element array, always [1]
  !! @param dcpl       Dataset creation property list (chunked)
  !! @param index      Current timestep index (0-based)
  !! @param xf_id      Collective transfer property list
  !! @param ierr       HDF5 error code (output)
  subroutine vtkhdf_write_numberof(grp, dset_name, value, offset, cnt, &
       dcpl, index, xf_id, ierr)
    integer(hid_t), intent(in) :: grp, dcpl, xf_id
    character(len=*), intent(in) :: dset_name
    integer(kind=i8), intent(in) :: value
    integer, intent(in):: index
    integer(hsize_t), dimension(1), intent(in) :: offset, cnt
    integer, intent(out) :: ierr

    integer(hid_t) :: dset_id, fspace, mspace
    integer(hsize_t), dimension(1) :: dims, maxdims
    integer(hid_t) :: H5T_NEKO_INTEGER
    logical :: exists

    H5T_NEKO_INTEGER = h5kind_to_type(i8, H5_INTEGER_KIND)

    call h5lexists_f(grp, dset_name, exists, ierr)
    if (exists) then
       call h5dopen_f(grp, dset_name, dset_id, ierr)
       dims(1) = int(index + 1, hsize_t) * int(pe_size, hsize_t)
       call h5dset_extent_f(dset_id, dims, ierr)
    else
       dims(1) = int(pe_size, hsize_t)
       maxdims(1) = H5S_UNLIMITED_F
       call h5screate_simple_f(1, dims, fspace, ierr, maxdims)
       call h5dcreate_f(grp, dset_name, H5T_NEKO_INTEGER, &
            fspace, dset_id, ierr, dcpl_id = dcpl)
       call h5sclose_f(fspace, ierr)
    end if

    call h5dget_space_f(dset_id, fspace, ierr)
    call h5screate_simple_f(1, cnt, mspace, ierr)
    call h5sselect_hyperslab_f(fspace, H5S_SELECT_SET_F, offset, cnt, ierr)

    call h5dwrite_f(dset_id, H5T_NEKO_INTEGER, value, cnt, ierr, &
         file_space_id = fspace, mem_space_id = mspace, xfer_prp = xf_id)

    call h5sclose_f(mspace, ierr)
    call h5sclose_f(fspace, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine vtkhdf_write_numberof

  !> Append an integer scalar to a 1D chunked HDF5 dataset.
  !! Opens the dataset if it exists, or creates a new empty chunked dataset.
  !! Extends by 1 element and writes the value at the end.
  !! @param grp_id HDF5 group containing the dataset
  !! @param name Dataset name
  !! @param value Value to append
  !! @param index Position to write to
  subroutine vtkhdf_write_i8_at(grp_id, name, value, index)
    integer(hid_t), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer(kind=i8), intent(in) :: value
    integer, intent(in) :: index

    integer :: ierr
    integer(hid_t) :: dset_id, dcpl_id, xf_id, filespace, memspace
    integer(hsize_t), dimension(1) :: dims, maxdims, count, offset, chunkdims
    integer(hid_t) :: H5T_NEKO_INTEGER
    logical :: exists

    H5T_NEKO_INTEGER = h5kind_to_type(i8, H5_INTEGER_KIND)

    ! Create collective transfer property list
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    call h5lexists_f(grp_id, name, exists, ierr)
    if (exists) then
       call h5dopen_f(grp_id, name, dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sget_simple_extent_dims_f(filespace, dims, maxdims, ierr)
       call h5sclose_f(filespace, ierr)

       if (index .eq. dims(1)) then
          dims(1) = int(index + 1, hsize_t)
          call h5dset_extent_f(dset_id, dims, ierr)
       else if (index .gt. dims(1)) then
          call neko_error("VTKHDF: Values written out of order.")
       end if
    else
       dims = 1_hsize_t
       maxdims = H5S_UNLIMITED_F
       chunkdims = 1_hsize_t

       call h5screate_simple_f(1, dims, filespace, ierr, maxdims)
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       call h5pset_chunk_f(dcpl_id, 1, chunkdims, ierr)
       call h5dcreate_f(grp_id, name, H5T_NEKO_INTEGER, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5sclose_f(filespace, ierr)
       call h5pclose_f(dcpl_id, ierr)
    end if

    count = 1_hsize_t
    offset = int(index, hsize_t)

    call h5dget_space_f(dset_id, filespace, ierr)
    call h5screate_simple_f(1, count, memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_INTEGER, value, count, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)

    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(xf_id, ierr)

  end subroutine vtkhdf_write_i8_at

  !> @brief Write a scalar field as a 1D dataset.
  !! @details
  !! This function simplifies the interfacing with HDF5 for the purpose of
  !! writing an array.
  !! @param hdf_root HDF5 group or file identifier to write under
  !! @param name Dataset name
  !! @param x Local array of field data to write
  !! @param n_local Number of local points in x
  !! @param precision Desired output precision (sp, dp, or rp)
  !! @param n_total Total number of points across all MPI ranks (optional)
  !! @param offset Starting index offset for this rank (optional)
  subroutine write_scalar_field(hdf_root, name, x, n_local, &
       precision, n_total, offset)
    integer, intent(in) :: n_local
    integer(hid_t), intent(in) :: hdf_root
    character(len=*), intent(in) :: name
    real(kind=rp), dimension(n_local), intent(in) :: x
    integer, intent(in), optional :: precision
    integer, intent(in), optional :: n_total, offset

    integer(hsize_t), dimension(1) :: dims, dcount, doffset
    integer(hid_t) :: xf_id, dset_id, filespace, memspace, precision_hdf
    integer :: i, ierr, precision_local, n_tot, off

    ! Setup data sizes, offsets and precision
    dcount = int(n_local, hsize_t)

    if (present(n_total)) then
       dims = int(n_total, hsize_t)
    else
       call MPI_Allreduce(n_local, n_tot, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, &
            ierr)
       dims = int(n_tot, hsize_t)
    end if

    if (present(offset)) then
       doffset = int(offset, hsize_t)
    else
       call MPI_Exscan(n_local, off, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, &
            ierr)
       doffset = int(off, hsize_t)
    end if

    if (present(precision)) then
       precision_local = precision
    else
       precision_local = rp
    end if
    precision_hdf = h5kind_to_type(precision_local, H5_REAL_KIND)

    ! Prepare memory and filespaces
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    call h5screate_simple_f(1, dims, filespace, ierr)
    call h5screate_simple_f(1, dcount, memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset, dcount, ierr)

    call h5dcreate_f(hdf_root, trim(name), precision_hdf, &
         filespace, dset_id, ierr)

    if (precision_local .eq. rp) then
       call h5dwrite_f(dset_id, precision_hdf, x, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = xf_id)

    else if (precision_local .eq. sp) then
       block
         real(kind=sp), allocatable :: x_sp(:)

         allocate(x_sp(n_local))
         do concurrent (i = 1:n_local)
            x_sp(i) = real(x(i), sp)
         end do

         call h5dwrite_f(dset_id, precision_hdf, x_sp, dcount, ierr, &
              file_space_id = filespace, mem_space_id = memspace, &
              xfer_prp = xf_id)
         deallocate(x_sp)
       end block

    else if (precision_local .eq. dp) then
       block
         real(kind=dp), allocatable :: x_dp(:)

         allocate(x_dp(n_local))
         do concurrent (i = 1:n_local)
            x_dp(i) = real(x(i), dp)
         end do

         call h5dwrite_f(dset_id, precision_hdf, x_dp, dcount, ierr, &
              file_space_id = filespace, mem_space_id = memspace, &
              xfer_prp = xf_id)
         deallocate(x_dp)
       end block

    else if (precision_local .eq. qp) then
       block
         real(kind=qp), allocatable :: x_qp(:)

         allocate(x_qp(n_local))
         do concurrent (i = 1:n_local)
            x_qp(i) = real(x(i), qp)
         end do

         call h5dwrite_f(dset_id, precision_hdf, x_qp, dcount, ierr, &
              file_space_id = filespace, mem_space_id = memspace, &
              xfer_prp = xf_id)
         deallocate(x_qp)
       end block

    else
       call neko_error("Unsupported precision in HDF5 write_scalar_field")
    end if

    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(xf_id, ierr)
  end subroutine write_scalar_field

  !> @brief Write 3 vector components as a single 2D dataset with shape (3, n_local).
  !! @details
  !! This function simplifies the interfacing with HDF5 for the purpose of
  !! writing an array.
  !! @param hdf_root HDF5 group or file identifier to write under
  !! @param name Dataset name
  !! @param u Local array of u component to write
  !! @param v Local array of v component to write
  !! @param w Local array of w component to write
  !! @param n_local Number of entries in each component
  !! @param precision Desired output precision (sp, dp, or rp)
  !! @param n_total Total number of points across all MPI ranks (optional)
  !! @param offset Starting index offset for this rank (optional)
  subroutine write_vector_field(hdf_root, name, u, v, w, n_local, &
       precision, n_total, offset)
    integer, intent(in) :: n_local
    integer(hid_t), intent(in) :: hdf_root
    character(len=*), intent(in) :: name
    real(kind=rp), dimension(n_local), intent(in) :: u, v, w
    integer, intent(in), optional :: precision
    integer, intent(in), optional :: n_total, offset

    integer(hsize_t), dimension(2) :: dims, dcount, doffset
    integer(hid_t) :: xf_id, dset_id, filespace, memspace, precision_hdf
    integer :: i, ierr, precision_local, n_tot, off

    ! Setup data sizes, offsets and precision
    dcount = [3_hsize_t, int(n_local, hsize_t)]

    if (present(n_total)) then
       dims = [3_hsize_t, int(n_total, hsize_t)]
    else
       call MPI_Allreduce(n_local, n_tot, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, &
            ierr)
       dims = [3_hsize_t, int(n_tot, hsize_t)]
    end if

    if (present(offset)) then
       doffset = [0_hsize_t, int(offset, hsize_t)]
    else
       call MPI_Exscan(n_local, off, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, &
            ierr)
       doffset = [0_hsize_t, int(off, hsize_t)]
    end if

    if (present(precision)) then
       precision_local = precision
    else
       precision_local = rp
    end if
    precision_hdf = h5kind_to_type(precision_local, H5_REAL_KIND)

    ! Prepare memory and filespaces
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    call h5screate_simple_f(2, dims, filespace, ierr)
    call h5screate_simple_f(2, dcount, memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset, dcount, ierr)

    call h5dcreate_f(hdf_root, trim(name), precision_hdf, &
         filespace, dset_id, ierr)

    if (precision_local .eq. sp) then
       block
         real(kind=sp), allocatable :: f(:,:)

         allocate(f(3, n_local))
         do concurrent (i = 1:n_local)
            f(1, i) = real(u(i), sp)
            f(2, i) = real(v(i), sp)
            f(3, i) = real(w(i), sp)
         end do

         call h5dwrite_f(dset_id, precision_hdf, f, dcount, ierr, &
              file_space_id = filespace, mem_space_id = memspace, &
              xfer_prp = xf_id)
         deallocate(f)
       end block

    else if (precision_local .eq. dp) then
       block
         real(kind=dp), allocatable :: f(:,:)

         allocate(f(3, n_local))
         do concurrent (i = 1:n_local)
            f(1, i) = real(u(i), dp)
            f(2, i) = real(v(i), dp)
            f(3, i) = real(w(i), dp)
         end do

         call h5dwrite_f(dset_id, precision_hdf, f, dcount, ierr, &
              file_space_id = filespace, mem_space_id = memspace, &
              xfer_prp = xf_id)
         deallocate(f)
       end block

    else if (precision_local .eq. qp) then
       block
         real(kind=qp), allocatable :: f(:,:)

         allocate(f(3, n_local))
         do concurrent (i = 1:n_local)
            f(1, i) = real(u(i), qp)
            f(2, i) = real(v(i), qp)
            f(3, i) = real(w(i), qp)
         end do

         call h5dwrite_f(dset_id, precision_hdf, f, dcount, ierr, &
              file_space_id = filespace, mem_space_id = memspace, &
              xfer_prp = xf_id)
         deallocate(f)
       end block

    else
       call neko_error("Unsupported precision in HDF5 write_vector_field")
    end if

    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(xf_id, ierr)
  end subroutine write_vector_field

  !> Read data in HDF5 format following official VTKHDF specification
  subroutine vtkhdf_file_read(this, data)
    class(vtkhdf_file_t) :: this
    class(*), target, intent(inout) :: data

    call neko_error('VTKHDF file reading is not yet implemented')

  end subroutine vtkhdf_file_read

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

  ! -------------------------------------------------------------------------- !
  ! Sub-cell node ordering functions for VTK compatibility

  !> Build linear hexahedron sub-cell node ordering for a spectral element.
  !! Returns an array of 0-based tensor-product indices that subdivides
  !! the lx*ly*lz grid into (lx-1)*(ly-1)*(lz-1) linear hexahedra
  !! (VTK_HEXAHEDRON, type 12), with 8 nodes per sub-cell.
  !! @param lx Number of points in x-direction
  !! @param ly Number of points in y-direction
  !! @param lz Number of points in z-direction
  !! @return Array of size 8*(lx-1)*(ly-1)*(lz-1) with 0-based indices
  pure function subdivide_to_hex_ordering(lx, ly, lz) result(node_order)
    integer, intent(in) :: lx, ly, lz
    integer :: node_order(8 * (lx - 1) * (ly - 1) * (lz - 1))
    integer :: ii, jj, kk, idx

    idx = 0

    do ii = 1, lx - 1
       do jj = 1, ly - 1
          do kk = 1, lz - 1
             node_order(idx + 1) = (kk - 1) * lx * ly + (jj - 1) * lx + ii - 1
             node_order(idx + 2) = (kk - 1) * lx * ly + (jj - 1) * lx + ii
             node_order(idx + 3) = (kk - 1) * lx * ly + jj * lx + ii
             node_order(idx + 4) = (kk - 1) * lx * ly + jj * lx + ii - 1
             node_order(idx + 5) = kk * lx * ly + (jj - 1) * lx + ii - 1
             node_order(idx + 6) = kk * lx * ly + (jj - 1) * lx + ii
             node_order(idx + 7) = kk * lx * ly + jj * lx + ii
             node_order(idx + 8) = kk * lx * ly + jj * lx + ii - 1
             idx = idx + 8
          end do
       end do
    end do

  end function subdivide_to_hex_ordering

  !> Build linear quadrilateral sub-cell node ordering for a spectral element.
  !! Returns an array of 0-based tensor-product indices that subdivides
  !! the lx*ly grid into (lx-1)*(ly-1) linear quadrilaterals
  !! (VTK_QUAD, type 9), with 4 nodes per sub-cell.
  !! @param lx Number of points in x-direction
  !! @param ly Number of points in y-direction
  !! @return Array of size 4*(lx-1)*(ly-1) with 0-based indices
  pure function subdivide_to_quad_ordering(lx, ly) result(node_order)
    integer, intent(in) :: lx, ly
    integer :: node_order(4 * (lx - 1) * (ly - 1))
    integer :: ii, jj, idx

    idx = 0

    do jj = 1, ly - 1
       do ii = 1, lx - 1
          node_order(idx + 1) = (jj - 1) * lx + ii - 1
          node_order(idx + 2) = (jj - 1) * lx + ii
          node_order(idx + 3) = jj * lx + ii
          node_order(idx + 4) = jj * lx + ii - 1
          idx = idx + 4
       end do
    end do

  end function subdivide_to_quad_ordering

end module vtkhdf_file
