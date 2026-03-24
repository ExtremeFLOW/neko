! Copyright (c) 2024-2025, The Neko Authors
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
module hdf5_file
  use num_types, only : rp, dp, sp
  use generic_file, only : generic_file_t
  use checkpoint, only : chkp_t
  use utils, only : neko_error, neko_warning, filename_suffix_pos, filename_split
  use mesh, only : mesh_t
  use field, only : field_t, field_ptr_t
  use field_list, only : field_list_t
  use field_series, only : field_series_t, field_series_ptr_t
  use dofmap, only : dofmap_t
  use logger, only : neko_log
  use vector, only : vector_t
  use matrix, only : matrix_t
  use field, only : field_t
  use device, only : DEVICE_TO_HOST
  use datadist, only : linear_dist_t
  use comm, only : pe_rank, pe_size, NEKO_COMM
  use mpi_f08, only : MPI_INFO_NULL, MPI_Allreduce, MPI_Allgather, &
       MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, MPI_MAX, MPI_Comm_size, MPI_Exscan, &
       MPI_Barrier, MPI_INTEGER8, MPI_Scan
#ifdef HAVE_HDF5
  use hdf5
#endif
  implicit none
  private

  !> Interface for HDF5 files
  type, public, extends(generic_file_t) :: hdf5_file_t

     ! HDF5 members
#ifdef HAVE_HDF5
     integer(hid_t) :: file_id = -1_hid_t
     integer(hid_t) :: active_group_id = -1_hid_t
     integer(hid_t) :: plist_id = -1_hid_t
#endif
     character(len=1) :: mode
     integer :: precision = -1
     integer :: offset = 0
     integer :: count = 0

   contains
     ! General methods for reading/writing HDF5 files
     procedure :: read => hdf5_file_read
     procedure :: write => hdf5_file_write
     procedure :: set_overwrite => hdf5_file_set_overwrite
     ! Granular methods for dealing with HDF5 files
     procedure :: open => hdf5_file_open
     procedure :: close => hdf5_file_close
     procedure :: set_active_group => hdf5_file_set_group
     procedure :: set_precision => hdf5_file_set_precision
     procedure :: get_fname => file_get_fname
     procedure, pass(this) :: write_vector => hdf5_file_write_vector
     procedure, pass(this) :: write_matrix => hdf5_file_write_matrix
     procedure, pass(this) :: write_field => hdf5_file_write_field
     procedure, pass(this) :: write_int_attribute => hdf5_file_write_int_attribute
     procedure, pass(this) :: write_rp_attribute => hdf5_file_write_rp_attribute
     procedure, pass(this) :: read_vector => hdf5_file_read_vector
     procedure, pass(this) :: read_matrix => hdf5_file_read_matrix
     procedure :: write_dataset => hdf5_file_write_dataset
     procedure :: read_dataset => hdf5_file_read_dataset
     procedure :: write_attribute => hdf5_file_write_attribute
  end type hdf5_file_t

contains

  !> Set the overwrite flag for HDF5 files
  subroutine hdf5_file_set_overwrite(this, overwrite)
    class(hdf5_file_t), intent(inout) :: this
    logical, intent(in) :: overwrite
    this%overwrite = overwrite
  end subroutine hdf5_file_set_overwrite

  !> Return the file name with the start counter.
  function file_get_fname(this) result(base_fname)
    class(hdf5_file_t), intent(in) :: this
    character(len=1024) :: base_fname
    character(len=1024) :: fname
    character(len=1024) :: path, name, suffix

    fname = trim(this%get_base_fname())
    call filename_split(fname, path, name, suffix)

    ! Append a counter
    !write(base_fname, '(A,A,"_",I0,A)') &
    !     trim(path), trim(name), this%get_start_counter(), trim(suffix)

    ! Do not append anything
    base_fname = trim(fname)

  end function file_get_fname

  !> Set the precision for the output (single or double)
  subroutine hdf5_file_set_precision(this, precision)
    class(hdf5_file_t), intent(inout) :: this
    integer, intent(in) :: precision
    this%precision = precision
  end subroutine hdf5_file_set_precision


#ifdef HAVE_HDF5

  ! ===============
  ! General methods
  ! ===============

  !> Write data in HDF5 format
  subroutine hdf5_file_write(this, data, t)
    class(hdf5_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dof
    type(field_ptr_t), allocatable :: fp(:)
    type(field_series_ptr_t), allocatable :: fsp(:)
    real(kind=rp), pointer :: dtlag(:)
    real(kind=rp), pointer :: tlag(:)
    integer :: ierr, info, drank, i, j
    integer(hid_t) :: plist_id, file_id, dset_id, grp_id, attr_id
    integer(hid_t) :: filespace, memspace
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: suffix_pos
    character(len=5) :: id_str
    character(len=1024) :: fname

    call hdf5_file_determine_data(data, msh, dof, fp, fsp, dtlag, tlag)

    if (.not. this%overwrite) call this%increment_counter()
    fname = trim(this%get_fname())

    call h5open_f(ierr)

    call hdf5_file_determine_real(H5T_NEKO_REAL)

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(plist_id, NEKO_COMM%mpi_val, info, ierr)

    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, &
         file_id, ierr, access_prp = plist_id)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    call h5screate_f(H5S_SCALAR_F, filespace, ierr)
    ddim = 1

    if (present(t)) then
       call h5acreate_f(file_id, "Time", H5T_NEKO_REAL, filespace, attr_id, &
            ierr, h5p_default_f, h5p_default_f)
       call h5awrite_f(attr_id, H5T_NEKO_REAL, t, ddim, ierr)
       call h5aclose_f(attr_id, ierr)
    end if

    if (associated(dof)) then
       call h5acreate_f(file_id, "Lx", H5T_NATIVE_INTEGER, filespace, attr_id, &
            ierr, h5p_default_f, h5p_default_f)
       call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, dof%Xh%lx, ddim, ierr)
       call h5aclose_f(attr_id, ierr)
    end if

    if (associated(msh)) then
       call h5gcreate_f(file_id, "Mesh", grp_id, ierr, lcpl_id=h5p_default_f, &
            gcpl_id=h5p_default_f, gapl_id=h5p_default_f)

       call h5acreate_f(grp_id, "Elements", H5T_NATIVE_INTEGER, filespace, attr_id, &
            ierr, h5p_default_f, h5p_default_f)
       call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, msh%glb_nelv, ddim, ierr)
       call h5aclose_f(attr_id, ierr)

       call h5acreate_f(grp_id, "Dimension", H5T_NATIVE_INTEGER, filespace, attr_id, &
            ierr, h5p_default_f, h5p_default_f)
       call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, msh%gdim, ddim, ierr)
       call h5aclose_f(attr_id, ierr)

       call h5gclose_f(grp_id, ierr)
    end if


    call h5sclose_f(filespace, ierr)

    !
    ! Write restart group (tlag, dtlag)
    !
    if (associated(tlag) .and. associated(dtlag)) then
       call h5gcreate_f(file_id, "Restart", grp_id, ierr, lcpl_id=h5p_default_f, &
            gcpl_id=h5p_default_f, gapl_id=h5p_default_f)

       drank = 1
       ddim = size(tlag)
       doffset(1) = 0
       if (pe_rank .eq. 0) then
          dcount = size(tlag)
       else
          dcount = 0
       end if

       call h5screate_simple_f(drank, ddim, filespace, ierr)

       call h5dcreate_f(grp_id,'tlag', H5T_NEKO_REAL, &
            filespace, dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dwrite_f(dset_id, H5T_NEKO_REAL, tlag, &
            ddim, ierr, xfer_prp = plist_id)
       call h5dclose_f(dset_id, ierr)

       call h5dcreate_f(grp_id,'dtlag', H5T_NEKO_REAL, &
            filespace, dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dwrite_f(dset_id, H5T_NEKO_REAL, dtlag, &
            ddim, ierr, xfer_prp = plist_id)
       call h5dclose_f(dset_id, ierr)

       call h5sclose_f(filespace, ierr)
       call h5gclose_f(grp_id, ierr)

    end if


    !
    ! Write fields group
    !
    if (allocated(fp) .or. allocated(fsp)) then
       call h5gcreate_f(file_id, "Fields", grp_id, ierr, lcpl_id=h5p_default_f, &
            gcpl_id=h5p_default_f, gapl_id=h5p_default_f)

       dcount(1) = int(dof%size(), 8)
       doffset(1) = int(msh%offset_el, 8) * int((dof%Xh%lx**3),8)
       ddim = int(dof%size(), 8)
       drank = 1
       call MPI_Allreduce(MPI_IN_PLACE, ddim(1), 1, &
            MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)

       call h5screate_simple_f(drank, ddim, filespace, ierr)
       call h5screate_simple_f(drank, dcount, memspace, ierr)


       if (allocated(fp)) then
          do i = 1, size(fp)
             call h5dcreate_f(grp_id, fp(i)%ptr%name, H5T_NEKO_REAL, &
                  filespace, dset_id, ierr)
             call h5dget_space_f(dset_id, filespace, ierr)
             call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                  doffset, dcount, ierr)
             call h5dwrite_f(dset_id, H5T_NEKO_REAL, &
                  fp(i)%ptr%x(1,1,1,1), &
                  ddim, ierr, file_space_id = filespace, &
                  mem_space_id = memspace, xfer_prp = plist_id)
             call h5dclose_f(dset_id, ierr)
          end do
          deallocate(fp)
       end if

       if (allocated(fsp)) then
          do i = 1, size(fsp)
             do j = 1, fsp(i)%ptr%size()
                call h5dcreate_f(grp_id, fsp(i)%ptr%lf(j)%name, &
                     H5T_NEKO_REAL, filespace, dset_id, ierr)
                call h5dget_space_f(dset_id, filespace, ierr)
                call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                     doffset, dcount, ierr)
                call h5dwrite_f(dset_id, H5T_NEKO_REAL, &
                     fsp(i)%ptr%lf(j)%x(1,1,1,1), &
                     ddim, ierr, file_space_id = filespace, &
                     mem_space_id = memspace, xfer_prp = plist_id)
                call h5dclose_f(dset_id, ierr)
             end do
          end do
          deallocate(fsp)
       end if

       call h5gclose_f(grp_id, ierr)
       call h5sclose_f(filespace, ierr)
       call h5sclose_f(memspace, ierr)
    end if

    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)

    call h5close_f(ierr)

  end subroutine hdf5_file_write

  !> Read data in HDF5 format
  subroutine hdf5_file_read(this, data)
    class(hdf5_file_t) :: this
    class(*), target, intent(inout) :: data
    integer(hid_t) :: plist_id, file_id, dset_id, grp_id, attr_id
    integer(hid_t) :: filespace, memspace
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: i,j, ierr, info, glb_nelv, gdim, lx, drank
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dof
    type(field_ptr_t), allocatable :: fp(:)
    type(field_series_ptr_t), allocatable :: fsp(:)
    real(kind=rp), pointer :: dtlag(:)
    real(kind=rp), pointer :: tlag(:)
    real(kind=rp) :: t
    character(len=1024) :: fname

    fname = trim(this%get_fname())

    call hdf5_file_determine_data(data, msh, dof, fp, fsp, dtlag, tlag)

    call h5open_f(ierr)

    call hdf5_file_determine_real(H5T_NEKO_REAL)

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(plist_id, NEKO_COMM%mpi_val, info, ierr)

    call h5fopen_f(fname, H5F_ACC_RDONLY_F, &
         file_id, ierr, access_prp = plist_id)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ddim = 1
    call h5aopen_name_f(file_id, 'Time', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NEKO_REAL, t, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    select type(data)
    type is(chkp_t)
       data%t = t
    end select

    call h5aopen_name_f(file_id, 'Lx', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, lx, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5gopen_f(file_id, 'Mesh', grp_id, ierr, gapl_id=h5p_default_f)

    call h5aopen_name_f(grp_id, 'Elements', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, glb_nelv, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_name_f(grp_id, 'Dimension', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, gdim, ddim, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5gclose_f(grp_id, ierr)


    if (associated(tlag) .and. associated(dtlag)) then
       drank = 1
       ddim = size(tlag)
       doffset(1) = 0
       if (pe_rank .eq. 0) then
          dcount = size(tlag)
       else
          dcount = 0
       end if

       call h5gopen_f(file_id, 'Restart', grp_id, ierr, gapl_id=h5p_default_f)
       call h5dopen_f(grp_id, 'tlag', dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dread_f(dset_id, H5T_NEKO_REAL, tlag, ddim, ierr, xfer_prp=plist_id)
       call h5dclose_f(dset_id, ierr)
       call h5sclose_f(filespace, ierr)

       call h5dopen_f(grp_id, 'dtlag', dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dread_f(dset_id, H5T_NEKO_REAL, dtlag, ddim, ierr, xfer_prp=plist_id)
       call h5dclose_f(dset_id, ierr)
       call h5sclose_f(filespace, ierr)

       call h5gclose_f(grp_id, ierr)
    end if

    if (allocated(fp) .or. allocated(fsp)) then
       call h5gopen_f(file_id, 'Fields', grp_id, ierr, gapl_id=h5p_default_f)

       dcount(1) = int(dof%size(), 8)
       doffset(1) = int(msh%offset_el, 8) * int((dof%Xh%lx**3),8)
       ddim = int(dof%size(), 8)
       drank = 1

       dcount(1) = int(dof%size(), 8)
       doffset(1) = int(msh%offset_el, 8) * int((dof%Xh%lx**3),8)
       ddim = int(dof%size(), 8)
       drank = 1
       call MPI_Allreduce(MPI_IN_PLACE, ddim(1), 1, &
            MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)

       call h5screate_simple_f(drank, dcount, memspace, ierr)

       if (allocated(fp)) then
          do i = 1, size(fp)
             call h5dopen_f(grp_id, fp(i)%ptr%name, dset_id, ierr)
             call h5dget_space_f(dset_id, filespace, ierr)
             call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
                  doffset, dcount, ierr)
             call h5dread_f(dset_id, H5T_NEKO_REAL, &
                  fp(i)%ptr%x(1,1,1,1), &
                  ddim, ierr, file_space_id = filespace, &
                  mem_space_id = memspace, xfer_prp=plist_id)
             call h5dclose_f(dset_id, ierr)
             call h5sclose_f(filespace, ierr)
          end do
       end if

       if (allocated(fsp)) then
          do i = 1, size(fsp)
             do j = 1, fsp(i)%ptr%size()
                call h5dopen_f(grp_id, fsp(i)%ptr%lf(j)%name, dset_id, ierr)
                call h5dget_space_f(dset_id, filespace, ierr)
                call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
                     doffset, dcount, ierr)
                call h5dread_f(dset_id, H5T_NEKO_REAL, &
                     fsp(i)%ptr%lf(j)%x(1,1,1,1), &
                     ddim, ierr, file_space_id = filespace, &
                     mem_space_id = memspace, xfer_prp=plist_id)
                call h5dclose_f(dset_id, ierr)
                call h5sclose_f(filespace, ierr)
             end do
          end do
       end if
       call h5sclose_f(memspace, ierr)
       call h5gclose_f(grp_id, ierr)
    end if

    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)

    call h5close_f(ierr)

  end subroutine hdf5_file_read


  subroutine hdf5_file_determine_data(data, msh, dof, fp, fsp, dtlag, tlag)
    class(*), target, intent(in) :: data
    type(mesh_t), pointer, intent(inout) :: msh
    type(dofmap_t), pointer, intent(inout) :: dof
    type(field_ptr_t), allocatable, intent(inout) :: fp(:)
    type(field_series_ptr_t), allocatable, intent(inout) :: fsp(:)
    real(kind=rp), pointer, intent(inout) :: dtlag(:)
    real(kind=rp), pointer, intent(inout) :: tlag(:)
    integer :: i, j, fp_size, fp_cur, fsp_size, fsp_cur, scalar_count, ab_count
    character(len=32) :: scalar_name

    select type(data)
    type is (field_t)
       dof => data%dof
       msh => data%msh
       fp_size = 1
       allocate(fp(fp_size))
       fp(1)%ptr => data

       nullify(dtlag)
       nullify(tlag)

    type is (field_list_t)

       if (data%size() .gt. 0) then
          allocate(fp(data%size()))

          dof => data%dof(1)
          msh => data%msh(1)

          do i = 1, data%size()
             fp(i)%ptr => data%items(i)%ptr
          end do
       else
          call neko_error('Empty field list')
       end if

       nullify(dtlag)
       nullify(tlag)

    type is(chkp_t)

       if ( .not. associated(data%u) .or. &
            .not. associated(data%v) .or. &
            .not. associated(data%w) .or. &
            .not. associated(data%p) ) then
          call neko_error('Checkpoint not initialized')
       end if

       fp_size = 4

       if (allocated(data%scalar_lags%items) .and. data%scalar_lags%size() > 0) then
          scalar_count = data%scalar_lags%size()
       else if (associated(data%s)) then
          scalar_count = 1
       else
          scalar_count = 0
       end if

       if (scalar_count .gt. 1) then
          fp_size = fp_size + scalar_count

          ! Add abx1 and abx2 fields for each scalar
          fp_size = fp_size + (scalar_count * 2)
       else if (associated(data%s)) then
          ! Single scalar support
          fp_size = fp_size + 1
          if (associated(data%abs1)) then
             fp_size = fp_size + 2
          end if
       end if

       if (associated(data%abx1)) then
          fp_size = fp_size + 6
       end if

       allocate(fp(fp_size))

       fsp_size = 0
       if (associated(data%ulag)) then
          fsp_size = fsp_size + 3
       end if

       if (scalar_count .gt. 1) then
          if (allocated(data%scalar_lags%items)) then
             fsp_size = fsp_size + data%scalar_lags%size()
          end if
       else if (associated(data%slag)) then
          fsp_size = fsp_size + 1
       end if

       if (fsp_size .gt. 0) then
          allocate(fsp(fsp_size))
          fsp_cur = 1
       end if

       dof => data%u%dof
       msh => data%u%msh

       fp(1)%ptr => data%u
       fp(2)%ptr => data%v
       fp(3)%ptr => data%w
       fp(4)%ptr => data%p

       fp_cur = 5

       if (scalar_count .gt. 1) then
          block
            type(field_series_t), pointer :: slag
            do i = 1, scalar_count
               slag => data%scalar_lags%get(i)
               fp(fp_cur)%ptr => slag%f
               fp_cur = fp_cur + 1
            end do
          end block

          do i = 1, scalar_count
             fp(fp_cur)%ptr => data%scalar_abx1(i)%ptr
             fp_cur = fp_cur + 1
             fp(fp_cur)%ptr => data%scalar_abx2(i)%ptr
             fp_cur = fp_cur + 1
          end do
       else if (associated(data%s)) then
          ! Single scalar support
          fp(fp_cur)%ptr => data%s
          fp_cur = fp_cur + 1

          if (associated(data%abs1)) then
             fp(fp_cur)%ptr => data%abs1
             fp(fp_cur+1)%ptr => data%abs2
             fp_cur = fp_cur + 2
          end if
       end if

       if (associated(data%abx1)) then
          fp(fp_cur)%ptr => data%abx1
          fp(fp_cur+1)%ptr => data%abx2
          fp(fp_cur+2)%ptr => data%aby1
          fp(fp_cur+3)%ptr => data%aby2
          fp(fp_cur+4)%ptr => data%abz1
          fp(fp_cur+5)%ptr => data%abz2
          fp_cur = fp_cur + 6
       end if

       if (associated(data%ulag)) then
          fsp(fsp_cur)%ptr => data%ulag
          fsp(fsp_cur+1)%ptr => data%vlag
          fsp(fsp_cur+2)%ptr => data%wlag
          fsp_cur = fsp_cur + 3
       end if


       if (scalar_count .gt. 1) then
          if (allocated(data%scalar_lags%items)) then
             do j = 1, data%scalar_lags%size()
                fsp(fsp_cur)%ptr => data%scalar_lags%get(j)
                fsp_cur = fsp_cur + 1
             end do
          end if
       else if (associated(data%slag)) then
          fsp(fsp_cur)%ptr => data%slag
          fsp_cur = fsp_cur + 1
       end if

       if (associated(data%tlag)) then
          tlag => data%tlag
          dtlag => data%dtlag
       end if

    class default
       call neko_log%error('Invalid data')
    end select

  end subroutine hdf5_file_determine_data

  !> Determine hdf5 real type corresponding to NEKO_REAL
  !! @note This must be called after h5open_f, otherwise
  !! the H5T_NATIVE_XYZ types has a value of 0
  subroutine hdf5_file_determine_real(H5T_NEKO_REAL)
    integer(hid_t), intent(inout) :: H5T_NEKO_REAL
    select case(rp)
    case(dp)
       H5T_NEKO_REAL = H5T_NATIVE_DOUBLE
    case(sp)
       H5T_NEKO_REAL = H5T_NATIVE_REAL
    case default
       call neko_error("Unsupported real type")
    end select
  end subroutine hdf5_file_determine_real

  ! ================
  ! Granular methods
  ! ================

  !> Open a HDF5 file in a given mode
  subroutine hdf5_file_open(this, mode)
    class(hdf5_file_t), intent(inout) :: this
    character(len=1), intent(in) :: mode
    integer :: ierr, mpi_info, mpi_comm, i, n_fields, counter
    logical :: file_exists
    character(len=1024) :: fname

    ! Set the mode for the file
    this%mode = mode

    ! Ensure precision is set and are valid.
    if (this%precision .gt. rp) then
       this%precision = rp
       call neko_warning('Requested precision is higher than working precision')
    else if (this%precision .eq. -1) then
       this%precision = rp
    end if

    fname = trim(this%get_fname())
    counter = this%get_counter() - this%get_start_counter()

    ! Set the configuration for MPI IO
    mpi_info = MPI_INFO_NULL%mpi_val
    mpi_comm = NEKO_COMM%mpi_val
    call h5open_f(ierr)
    call h5pcreate_f(H5P_FILE_ACCESS_F, this%plist_id, ierr)
    call h5pset_fapl_mpio_f(this%plist_id, mpi_comm, mpi_info, ierr)

    ! Open the file
    inquire(file = fname, exist = file_exists)
    if (file_exists) then
       call h5fopen_f(fname, H5F_ACC_RDWR_F, this%file_id, ierr, &
            access_prp = this%plist_id)
    else
       call h5fcreate_f(fname, H5F_ACC_TRUNC_F, &
            this%file_id, ierr, access_prp = this%plist_id)
    end if

    ! Set the active group to the root of the file
    call this%set_active_group()

    if (pe_rank .eq.0) then
       write(*,*) "Opened HDF5 file: ", trim(fname), " with counter: ", counter
    end if

  end subroutine hdf5_file_open

  !> Close the file
  subroutine hdf5_file_close(this)
    class(hdf5_file_t), intent(inout) :: this
    integer :: ierr

    if (this%active_group_id .ne. -1_hid_t .and. this%active_group_id .ne. this%file_id) then
       call h5gclose_f(this%active_group_id, ierr)
    end if
    this%active_group_id = -1_hid_t

    call h5pclose_f(this%plist_id, ierr)
    this%plist_id = -1_hid_t
    call h5fclose_f(this%file_id, ierr)
    this%file_id = -1_hid_t
    call h5close_f(ierr)

    if (pe_rank .eq.0) then
       write(*,*) "Closed HDF5 file: ", trim(this%get_fname())
    end if

  end subroutine hdf5_file_close

  !> Set the active group for HDF5 files
  !! @param this The HDF5 file object
  !! @param An array of strings that show the path to the group to create or open.
  subroutine hdf5_file_set_group(this, group_name)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in), optional :: group_name(:)

    integer(hid_t) :: current_id, group_id
    integer :: ierr, i, num_groups
    logical :: group_exists


    ! Close previous active group if one is open
    if (this%active_group_id .ne. -1_hid_t .and. this%active_group_id .ne. this%file_id) then
       call h5gclose_f(this%active_group_id, ierr)
    end if
    this%active_group_id = -1_hid_t

    ! Start from root location = file
    current_id = this%file_id
    ! Return the root directory if no group name is given
    if (.not. present(group_name)) then
       this%active_group_id = current_id
       return
    end if

    ! Iterate through the group names
    num_groups = size(group_name)

    do i = 1, num_groups
       call h5lexists_f(current_id, trim(group_name(i)), group_exists, ierr)

       ! Only create groups if they dont exist and we are in write mode "w"
       if (group_exists) then
          call h5gopen_f(current_id, trim(group_name(i)), group_id, ierr)
       else
          if (this%mode == "r") then
             call neko_error("Group " // trim(group_name(i)) // " does not exist in file " // trim(this%get_fname()))
          end if
          call h5gcreate_f(current_id, trim(group_name(i)), group_id, ierr)
       end if

       ! Close previous location only if it was an opened group, not the file
       if (i > 1) then
          call h5gclose_f(current_id, ierr)
       end if

       current_id = group_id
    end do

    this%active_group_id = current_id
  end subroutine hdf5_file_set_group


  subroutine hdf5_file_write_dataset(this, data)
    class(hdf5_file_t), intent(inout) :: this
    class(*), intent(inout) :: data

    select type (d => data)
    type is (vector_t)
       call this%write_vector(d)
    type is (matrix_t)
       call this%write_matrix(d)
    type is (field_t)
       call this%write_field(d)
    class default
       call neko_error("write_dataset not implemented for this data type")
    end select
  end subroutine hdf5_file_write_dataset

  subroutine hdf5_file_read_dataset(this, keyword, data, strategy)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: keyword
    class(*), intent(inout) :: data
    character(len=*), intent(in), optional :: strategy

    select type (d => data)
    type is (vector_t)
       call this%read_vector(keyword, d, strategy)
    type is (matrix_t)
       call this%read_matrix(keyword, d, strategy)
    type is (field_t)
       call neko_error("Reading a field_t is not supported yet")
    class default
       call neko_error("read_dataset not implemented for this data type")
    end select
  end subroutine hdf5_file_read_dataset

  subroutine hdf5_file_write_attribute(this, data, data_name)
    class(hdf5_file_t), intent(inout) :: this
    class(*), intent(inout) :: data
    character(len=*), intent(in) :: data_name

    select type (d => data)
    type is (integer)
       call this%write_int_attribute(d, data_name)
    type is (real(kind=rp))
       call this%write_rp_attribute(d, data_name)
    class default
       call neko_error("write_attribute not implemented for this data type")
    end select
  end subroutine hdf5_file_write_attribute


  subroutine hdf5_file_write_vector(this, vec)
    class(hdf5_file_t), intent(inout) :: this
    type(vector_t), intent(inout) :: vec
    integer :: ierr, counts, offset, total_count, dset_rank, max_count
    integer(hsize_t) :: append_offset
    integer(hid_t) :: precision_hdf
    integer(hid_t) :: xf_id, filespace, dset_id, memspace, dcpl_id
    integer(hsize_t), dimension(1) :: dcount, doffset
    integer(hsize_t), dimension(1) :: ddims, ddims_max, chunkdims
    integer(hsize_t), dimension(1) :: tempddims, tempmaxddims
    logical :: dset_exists
    real(kind=sp), allocatable :: write_buffer_sp(:) ! Write buffer single
    real(kind=dp), allocatable :: write_buffer_dp(:) ! Write buffer double

    ! ===============
    ! Get vector info
    ! ===============
    counts = vec%size()
    append_offset = 0_hsize_t
    offset = 0
    total_count = 0
    max_count = 0
    call MPI_Scan(counts, offset, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)
    offset = offset - counts ! Not using exclusive scan
    call MPI_Allreduce(counts, total_count, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)
    call MPI_Allreduce(counts, max_count, 1, MPI_INTEGER, &
         MPI_MAX, NEKO_COMM, ierr)

    ! ===============
    ! Configure MPIIO
    ! ===============
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    precision_hdf = h5kind_to_type(this%precision, H5_REAL_KIND)

    ! ===================
    ! Create the data set
    ! ===================
    dset_rank = 1 ! rank 1 array, i.e. a vector
    ddims = [int(total_count, hsize_t)] ! global size of the vector
    chunkdims = [max(int(max_count, hsize_t), 1_hsize_t)] ! Enable chunking to be able to append
    ddims_max = [H5S_UNLIMITED_F] ! allow unlimited size for appending
    call h5lexists_f(this%active_group_id, trim(vec%name), dset_exists, ierr)
    if (dset_exists) then
       if (this%overwrite) then
          ! retrieve the dset id for the existing data set
          call h5dopen_f(this%active_group_id, trim(vec%name), dset_id, ierr)
       else
          ! Retreive the existing data set
          call h5dopen_f(this%active_group_id, trim(vec%name), dset_id, ierr)
          ! Retrieve the current filespace (shape space)
          call h5dget_space_f(dset_id, filespace, ierr)
          ! Get the current shape
          call h5sget_simple_extent_dims_f(filespace, tempddims, tempmaxddims, ierr)
          ! Clean up the opened file space
          call h5sclose_f(filespace, ierr)
          ! Overwrite the new full shape
          ddims(1) = ddims(1) + tempddims(1) ! New size
          append_offset = tempddims(1) ! current size which is the offset
          ! Extend the data set to the new shape
          call h5dset_extent_f(dset_id, ddims, ierr)
       end if
    else
       ! create file space of this shape
       call h5screate_simple_f(dset_rank, ddims, filespace, ierr, ddims_max)
       ! Create chunk property list (needed to be able to append)
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       call h5pset_chunk_f(dcpl_id, dset_rank, chunkdims, ierr)
       ! create the data set with the given shape
       call h5dcreate_f(this%active_group_id, trim(vec%name), precision_hdf, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       ! clean opened ids
       call h5sclose_f(filespace, ierr)
       call h5pclose_f(dcpl_id, ierr)
    end if

    ! ===========================
    ! Set up writing the data set
    ! ===========================
    dcount = [int(counts, hsize_t)] ! local size of the vector
    doffset = [int(offset, hsize_t) + append_offset] ! offset for this rank in the global vector
    ! Get the total file space (shape) of the data set
    call h5dget_space_f(dset_id, filespace, ierr)
    ! Get only the slice where my rank writes
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, ierr)
    ! Create the corresponding memory space (buffer) for my local data
    call h5screate_simple_f(dset_rank, dcount, memspace, ierr)


    ! =======================
    ! Cast and write the data
    ! =======================
    if (this%precision == sp) then
       allocate(write_buffer_sp(vec%size()))
       if (vec%size() > 0) write_buffer_sp = real(vec%x, kind=sp)
       ! Write the data
       call h5dwrite_f(dset_id, precision_hdf, write_buffer_sp, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = xf_id)
       deallocate(write_buffer_sp)
    else if (this%precision == dp) then
       allocate(write_buffer_dp(vec%size()))
       if (vec%size() > 0) write_buffer_dp = real(vec%x, kind=dp)
       ! Write the data
       call h5dwrite_f(dset_id, precision_hdf, write_buffer_dp, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = xf_id)
       deallocate(write_buffer_dp)
    else
       call neko_error("Unsupported precision")
    end if

    ! =======================
    ! Clean up
    ! =======================
    call h5pclose_f(xf_id, ierr)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)

  end subroutine hdf5_file_write_vector

  subroutine hdf5_file_write_matrix(this, mat)
    class(hdf5_file_t), intent(inout) :: this
    type(matrix_t), intent(inout) :: mat
    integer :: ierr, counts, offset, total_count, dset_rank, strides, max_count
    integer(hsize_t) :: append_offset
    integer(hid_t) :: precision_hdf
    integer(hid_t) :: xf_id, filespace, dset_id, memspace, dcpl_id
    integer(hsize_t), dimension(2) :: dcount, doffset
    integer(hsize_t), dimension(2) :: ddims, ddims_max, chunkdims
    integer(hsize_t), dimension(2) :: tempddims, tempmaxddims
    logical :: dset_exists
    real(kind=sp), allocatable :: write_buffer_sp(:,:) ! Write buffer single
    real(kind=dp), allocatable :: write_buffer_dp(:,:) ! Write buffer double

    ! ===============
    ! Get Matrix info
    ! ===============
    strides = mat%get_nrows()
    counts = mat%get_ncols()
    append_offset = 0_hsize_t
    total_count = 0
    max_count = 0
    offset = 0
    call MPI_Scan(counts, offset, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)
    offset = offset - counts ! Not using exclusive scan
    call MPI_Allreduce(counts, total_count, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)
    call MPI_Allreduce(counts, max_count, 1, MPI_INTEGER, &
         MPI_MAX, NEKO_COMM, ierr)

    ! ===============
    ! Configure MPIIO
    ! ===============
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    precision_hdf = h5kind_to_type(this%precision, H5_REAL_KIND)

    ! ===================
    ! Create the data set
    ! ===================
    dset_rank = 2 ! rank 2 array, i.e. a matrix
    ddims = [int(strides, hsize_t), int(total_count, hsize_t)] ! global size of the matrix
    chunkdims = [int(strides, hsize_t), max(int(max_count, hsize_t), 1_hsize_t)]
    ddims_max = [int(strides, hsize_t), H5S_UNLIMITED_F]
    call h5lexists_f(this%active_group_id, trim(mat%name), dset_exists, ierr)
    if (dset_exists) then
       if (this%overwrite) then
          ! retrieve the dset id for the existing data set
          if (pe_rank .eq. 0) then
             write(*,*) "Dataset ", trim(mat%name), " already exists in file ", trim(this%get_fname()), " and will be overwritten."
             write(*,*) "This only works if the global shape is the same"
          end if
          call h5dopen_f(this%active_group_id, trim(mat%name), dset_id, ierr)
       else
          call h5dopen_f(this%active_group_id, trim(mat%name), dset_id, ierr)
          call h5dget_space_f(dset_id, filespace, ierr)
          call h5sget_simple_extent_dims_f(filespace, tempddims, tempmaxddims, ierr)
          call h5sclose_f(filespace, ierr)
          ddims(2) = ddims(2) + tempddims(2)
          append_offset = tempddims(2)
          call h5dset_extent_f(dset_id, ddims, ierr)
       end if
    else
       ! create file space of this shape
       call h5screate_simple_f(dset_rank, ddims, filespace, ierr, ddims_max)
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       call h5pset_chunk_f(dcpl_id, dset_rank, chunkdims, ierr)
       ! create the data set with the given shape
       call h5dcreate_f(this%active_group_id, trim(mat%name), precision_hdf, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5sclose_f(filespace, ierr)
       call h5pclose_f(dcpl_id, ierr)
    end if

    ! ===========================
    ! Set up writing the data set
    ! ===========================
    dcount = [int(strides, hsize_t), int(counts, hsize_t)] ! local size of the matrix
    doffset = [0_hsize_t, int(offset, hsize_t) + append_offset] ! offset for this rank in the global matrix
    ! Get the total file space (shape) of the data set
    call h5dget_space_f(dset_id, filespace, ierr)
    ! Get only the slice where my rank writes
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, ierr)
    ! Create the corresponding memory space (buffer) for my local data
    call h5screate_simple_f(dset_rank, dcount, memspace, ierr)

    ! =======================
    ! Cast and write the data
    ! =======================
    if (this%precision == sp) then
       allocate(write_buffer_sp(mat%get_nrows(), mat%get_ncols()))
       if (mat%size() > 0) write_buffer_sp = real(mat%x, kind=sp)
       ! Write the data
       call h5dwrite_f(dset_id, precision_hdf, write_buffer_sp, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = xf_id)
       deallocate(write_buffer_sp)
    else if (this%precision == dp) then
       allocate(write_buffer_dp(mat%get_nrows(), mat%get_ncols()))
       if (mat%size() > 0) write_buffer_dp = real(mat%x, kind=dp)
       ! Write the data
       call h5dwrite_f(dset_id, precision_hdf, write_buffer_dp, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = xf_id)
       deallocate(write_buffer_dp)
    else
       call neko_error("Unsupported precision")
    end if

    ! =======================
    ! Clean up
    ! =======================
    call h5pclose_f(xf_id, ierr)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)

  end subroutine hdf5_file_write_matrix

  subroutine hdf5_file_write_field(this, field)
    class(hdf5_file_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    integer :: ierr, counts, offset, total_count, dset_rank, max_count
    integer :: stride_ax_1, stride_ax_2, stride_ax_3
    integer(hsize_t) :: append_offset
    integer(hid_t) :: precision_hdf
    integer(hid_t) :: xf_id, filespace, dset_id, memspace, dcpl_id
    integer(hsize_t), dimension(4) :: dcount, doffset
    integer(hsize_t), dimension(4) :: ddims, ddims_max, chunkdims
    integer(hsize_t), dimension(4) :: tempddims, tempmaxddims
    logical :: dset_exists
    real(kind=sp), allocatable :: write_buffer_sp(:,:,:,:) ! Write buffer single
    real(kind=dp), allocatable :: write_buffer_dp(:,:,:,:) ! Write buffer double

    ! ==============
    ! Get Field info
    ! ==============
    stride_ax_1 = field%Xh%lx
    stride_ax_2 = field%Xh%ly
    stride_ax_3 = field%Xh%lz
    counts = field%msh%nelv
    append_offset = 0_hsize_t
    total_count = field%msh%glb_nelv
    max_count = 0
    offset = field%msh%offset_el
    call MPI_Allreduce(counts, max_count, 1, MPI_INTEGER, &
         MPI_MAX, NEKO_COMM, ierr)

    ! ===============
    ! Configure MPIIO
    ! ===============
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    precision_hdf = h5kind_to_type(this%precision, H5_REAL_KIND)

    ! ===================
    ! Create the data set
    ! ===================
    dset_rank = 4 ! rank 4 array, i.e. a 4D tensor
    ddims = [int(stride_ax_1, hsize_t), &
         int(stride_ax_2, hsize_t), &
         int(stride_ax_3, hsize_t), &
         int(total_count, hsize_t)] ! global size of the tensor
    chunkdims = [int(stride_ax_1, hsize_t), &
         int(stride_ax_2, hsize_t), &
         int(stride_ax_3, hsize_t), &
         max(int(max_count, hsize_t), 1_hsize_t)]
    ddims_max = [int(stride_ax_1, hsize_t), &
         int(stride_ax_2, hsize_t), &
         int(stride_ax_3, hsize_t), &
         H5S_UNLIMITED_F]
    call h5lexists_f(this%active_group_id, trim(field%name), dset_exists, ierr)
    if (dset_exists) then
       if (this%overwrite) then
          ! retrieve the dset id for the existing data set
          if (pe_rank .eq. 0) then
             write(*,*) "Overwriting Dataset: ", trim(field%name)
             write(*,*) "This only works if the global shape is the same"
          end if
          call h5dopen_f(this%active_group_id, trim(field%name), dset_id, ierr)
       else
          call h5dopen_f(this%active_group_id, trim(field%name), dset_id, ierr)
          call h5dget_space_f(dset_id, filespace, ierr)
          call h5sget_simple_extent_dims_f(filespace, tempddims, tempmaxddims, ierr)
          call h5sclose_f(filespace, ierr)
          ddims(4) = ddims(4) + tempddims(4)
          append_offset = tempddims(4)
          call h5dset_extent_f(dset_id, ddims, ierr)
       end if
    else
       ! create file space of this shape
       call h5screate_simple_f(dset_rank, ddims, filespace, ierr, ddims_max)
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
       call h5pset_chunk_f(dcpl_id, dset_rank, chunkdims, ierr)
       ! create the data set with the given shape
       call h5dcreate_f(this%active_group_id, trim(field%name), precision_hdf, &
            filespace, dset_id, ierr, dcpl_id = dcpl_id)
       call h5sclose_f(filespace, ierr)
       call h5pclose_f(dcpl_id, ierr)
    end if

    ! ===========================
    ! Set up writing the data set
    ! ===========================
    dcount = [int(stride_ax_1, hsize_t), int(stride_ax_2, hsize_t), int(stride_ax_3, hsize_t), int(counts, hsize_t)] ! local size of the tensor
    doffset = [0_hsize_t, 0_hsize_t, 0_hsize_t, int(offset, hsize_t) + append_offset] ! offset for this rank in the global tensor
    ! Get the total file space (shape) of the data set
    call h5dget_space_f(dset_id, filespace, ierr)
    ! Get only the slice where my rank writes
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, ierr)
    ! Create the corresponding memory space (buffer) for my local data
    call h5screate_simple_f(dset_rank, dcount, memspace, ierr)

    ! =======================
    ! Cast and write the data
    ! =======================
    if (this%precision == sp) then
       allocate(write_buffer_sp(field%Xh%lx, field%Xh%ly, field%Xh%lz, field%msh%nelv))
       if (field%msh%nelv > 0) write_buffer_sp = real(field%x, kind=sp)
       ! Write the data
       call h5dwrite_f(dset_id, precision_hdf, write_buffer_sp, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = xf_id)
       deallocate(write_buffer_sp)
    else if (this%precision == dp) then
       allocate(write_buffer_dp(field%Xh%lx, field%Xh%ly, field%Xh%lz, field%msh%nelv))
       if (field%msh%nelv > 0) write_buffer_dp = real(field%x, kind=dp)
       ! Write the data
       call h5dwrite_f(dset_id, precision_hdf, write_buffer_dp, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = xf_id)
       deallocate(write_buffer_dp)
    else
       call neko_error("Unsupported precision")
    end if

    ! =======================
    ! Clean up
    ! =======================
    call h5pclose_f(xf_id, ierr)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)

  end subroutine hdf5_file_write_field

  !> Read a vector
  subroutine hdf5_file_read_vector(this, keyword, vec, strategy)
    class(hdf5_file_t) :: this
    character(len=*), intent(in) :: keyword
    type(vector_t), intent(inout) :: vec
    character(len=*), intent(in), optional :: strategy
    character(len=1000) :: strategy_
    integer :: ierr, counts, offset, total_count, dset_rank
    integer(hid_t) :: precision_hdf
    integer(hid_t) :: xf_id, filespace, dset_id, memspace
    integer(hsize_t), dimension(1) :: dcount, doffset
    integer(hsize_t), dimension(1) :: tempddims, tempmaxddims
    integer :: temprank
    logical :: dset_exists
    type(linear_dist_t) :: dist

    ! Set up strategy
    if (present(strategy)) then
       if (trim(strategy) .eq. "linear" .or. &
            trim(strategy) .eq. "rank_0") then
          strategy_ = strategy
       else
          call neko_error("Unsupported strategy: " // trim(strategy))
       end if
    else
       strategy_ = "linear"
    end if

    ! Free the input
    call vec%free()

    ! ===============
    ! Configure MPIIO
    ! ===============
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    precision_hdf = h5kind_to_type(rp, H5_REAL_KIND)

    ! ===================
    ! Get the data set info
    ! ===================
    call h5lexists_f(this%active_group_id, trim(keyword), dset_exists, ierr)
    if (dset_exists) then
       ! Open the data set
       call h5dopen_f(this%active_group_id, trim(keyword), dset_id, ierr)
       ! Get the current rank of the dataset
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sget_simple_extent_ndims_f(filespace, temprank, ierr)
       if (temprank .ne. 1) then
          call neko_error("Dataset " // trim(keyword) // " is not a rank 1 vector in file " // trim(this%get_fname()))
       end if
       ! Get the current shape and close the filespace
       call h5sget_simple_extent_dims_f(filespace, tempddims, tempmaxddims, ierr)
       call h5sclose_f(filespace, ierr)
    else
       call neko_error("Dataset " // trim(keyword) // " does not exist in current group " // trim(this%get_fname()))
    end if

    ! =============================
    ! Perform the data distribution
    ! =============================
    total_count = int(tempddims(1))
    if (strategy_ .eq. "linear") then
       dist = linear_dist_t(total_count, pe_rank, pe_size, NEKO_COMM)
       counts = dist%num_local()
       offset = 0
    else if (strategy_ .eq. "rank_0") then
       if (pe_rank .eq. 0) then
          counts = total_count
       else
          counts = 0
       end if
       offset = 0
    end if
    call MPI_Scan(counts, offset, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)
    offset = offset - counts ! Not using exclusive scan

    ! ===========================
    ! Set up reading the data set
    ! ===========================
    dset_rank = 1 ! rank 1 array, i.e. a vector
    dcount = [int(counts, hsize_t)] ! local size of the vector
    doffset = [int(offset, hsize_t)] ! offset for this rank in the global vector
    ! Get the total file space (shape) of the data set
    call h5dget_space_f(dset_id, filespace, ierr)
    ! Get only the slice where my rank reads
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, ierr)
    ! Create the corresponding memory space (buffer) for my local data
    call h5screate_simple_f(dset_rank, dcount, memspace, ierr)

    ! =============================
    ! Allocate data. HDF5 will cast
    ! =============================
    call vec%init(counts, trim(keyword)) ! this is rp
    call h5dread_f(dset_id, precision_hdf, vec%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = xf_id)

    ! =======================
    ! Clean up
    ! =======================
    call h5pclose_f(xf_id, ierr)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)

  end subroutine hdf5_file_read_vector

  !> Read a matrix
  subroutine hdf5_file_read_matrix(this, keyword, mat, strategy)
    class(hdf5_file_t) :: this
    character(len=*), intent(in) :: keyword
    type(matrix_t), intent(inout) :: mat
    character(len=*), intent(in), optional :: strategy
    character(len=1000) :: strategy_
    integer :: ierr, counts, offset, total_count, dset_rank
    integer(hid_t) :: precision_hdf
    integer(hid_t) :: xf_id, filespace, dset_id, memspace
    integer(hsize_t), dimension(2) :: dcount, doffset
    integer(hsize_t), dimension(2) :: tempddims, tempmaxddims
    integer :: temprank
    logical :: dset_exists
    type(linear_dist_t) :: dist

    ! Set up strategy
    if (present(strategy)) then
       if (trim(strategy) .eq. "linear" .or. &
            trim(strategy) .eq. "rank_0") then
          strategy_ = strategy
       else
          call neko_error("Unsupported strategy: " // trim(strategy))
       end if
    else
       strategy_ = "linear"
    end if

    ! Free the input
    call mat%free()

    ! ===============
    ! Configure MPIIO
    ! ===============
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    precision_hdf = h5kind_to_type(rp, H5_REAL_KIND)

    ! ===================
    ! Get the data set info
    ! ===================
    call h5lexists_f(this%active_group_id, trim(keyword), dset_exists, ierr)
    if (dset_exists) then
       ! Openr the data set
       call h5dopen_f(this%active_group_id, trim(keyword), dset_id, ierr)
       ! Get the current rank of the of the dataset
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sget_simple_extent_ndims_f(filespace, temprank, ierr)
       if (temprank .ne. 2) then
          call neko_error("Dataset " // trim(keyword) // " is not a rank 2 matrix in file " // trim(this%get_fname()))
       end if
       ! Get the current shape and close the filespace
       call h5sget_simple_extent_dims_f(filespace, tempddims, tempmaxddims, ierr)
       call h5sclose_f(filespace, ierr)
    else
       call neko_error("Dataset " // trim(keyword) // " does not exist in current group " // trim(this%get_fname()))
    end if

    ! =============================
    ! Perform the data distribution
    ! =============================
    total_count = int(tempddims(2))
    if (strategy_ .eq. "linear") then
       dist = linear_dist_t(total_count, pe_rank, pe_size, NEKO_COMM)
       counts = dist%num_local()
       offset = 0
    else if (strategy_ .eq. "rank_0") then
       if (pe_rank .eq. 0) then
          counts = total_count
       else
          counts = 0
       end if
       offset = 0
    end if
    call MPI_Scan(counts, offset, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)
    offset = offset - counts ! Not using exclusive scan

    ! ===========================
    ! Set up reading the data set
    ! ===========================
    dset_rank = 2 ! rank 2 array, i.e. a matrix
    dcount = [int(tempddims(1), hsize_t), int(counts, hsize_t)] ! local size of the matrix
    doffset = [0_hsize_t, int(offset, hsize_t)] ! offset for this rank in the global matrix
    ! Get the total file space (shape) of the data set
    call h5dget_space_f(dset_id, filespace, ierr)
    ! Get only the slice where my rank reads
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, ierr)
    ! Create the corresponding memory space (buffer) for my local data
    call h5screate_simple_f(dset_rank, dcount, memspace, ierr)

    ! =============================
    ! Allocate data. HDF5 will cast
    ! =============================
    call mat%init(int(tempddims(1)), counts, trim(keyword)) ! this is rp
    call h5dread_f(dset_id, precision_hdf, mat%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = xf_id)

    ! =======================
    ! Clean up
    ! =======================
    call h5pclose_f(xf_id, ierr)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)
    call h5dclose_f(dset_id, ierr)


  end subroutine hdf5_file_read_matrix

  !> Write an integer attribute
  subroutine hdf5_file_write_int_attribute(this, attr, attr_name)
    class(hdf5_file_t), intent(inout) :: this
    integer, intent(in) :: attr
    character(len=*), intent(in) :: attr_name
    integer :: ierr
    integer(hid_t) :: filespace, attr_id
    integer(hsize_t), dimension(1) :: dcount
    logical :: attr_exists

    ! ====================
    ! Create the attribute
    ! ====================
    dcount = [int(1, hsize_t)]
    call h5aexists_f(this%active_group_id, trim(attr_name), attr_exists, ierr)
    if (attr_exists) then
       ! retrieve the attr id for the existing attribute
       call h5aopen_f(this%active_group_id, trim(attr_name), attr_id, ierr)
    else
       ! create file space of this shape
       call h5screate_f(H5S_SCALAR_F, filespace, ierr)
       ! create the data set with the given shape
       call h5acreate_f(this%active_group_id, trim(attr_name), H5T_NATIVE_INTEGER, &
            filespace, attr_id, ierr, h5p_default_f, h5p_default_f)
       call h5sclose_f(filespace, ierr)
    end if

    ! ===========================
    ! Set up writing the data set
    ! ===========================
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr, dcount, ierr)

    ! =======================
    ! Clean up
    ! =======================
    call h5aclose_f(attr_id, ierr)

  end subroutine hdf5_file_write_int_attribute

  !> Write a real (kind=rp) attribute
  subroutine hdf5_file_write_rp_attribute(this, attr, attr_name)
    class(hdf5_file_t), intent(inout) :: this
    real(kind=rp), intent(in) :: attr
    character(len=*), intent(in) :: attr_name
    integer :: ierr
    integer(hid_t) :: precision_hdf
    integer(hid_t) :: filespace, attr_id
    integer(hsize_t), dimension(1) :: dcount
    logical :: attr_exists

    ! Get the precision
    precision_hdf = h5kind_to_type(rp, H5_REAL_KIND)

    ! ====================
    ! Create the attribute
    ! ====================
    dcount = [int(1, hsize_t)]
    call h5aexists_f(this%active_group_id, trim(attr_name), attr_exists, ierr)
    if (attr_exists) then
       ! retrieve the attr id for the existing attribute
       call h5aopen_f(this%active_group_id, trim(attr_name), attr_id, ierr)
    else
       ! create file space of this shape
       call h5screate_f(H5S_SCALAR_F, filespace, ierr)
       ! create the data set with the given shape
       call h5acreate_f(this%active_group_id, trim(attr_name), precision_hdf, &
            filespace, attr_id, ierr, h5p_default_f, h5p_default_f)
       call h5sclose_f(filespace, ierr)
    end if

    ! ===========================
    ! Set up writing the data set
    ! ===========================
    call h5awrite_f(attr_id, precision_hdf, attr, dcount, ierr)

    ! =======================
    ! Clean up
    ! =======================
    call h5aclose_f(attr_id, ierr)

  end subroutine hdf5_file_write_rp_attribute

#else

  !> Write data in HDF5 format
  subroutine hdf5_file_write(this, data, t)
    class(hdf5_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_write

  !> Read data in HDF5 format
  subroutine hdf5_file_read(this, data)
    class(hdf5_file_t) :: this
    class(*), target, intent(inout) :: data
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_read

#endif

end module hdf5_file
