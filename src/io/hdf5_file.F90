! Copyright (c) 2024-2026, The Neko Authors
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
!> HDF5 file format.
module hdf5_file
  use num_types, only : rp, dp, sp, i8
  use generic_file, only : generic_file_t
  use checkpoint, only : chkp_t
  use checkpoint_payload, only : checkpoint_payload_t, checkpoint_array_t, &
       checkpoint_mesh_array_t
  use utils, only : neko_error, neko_warning, filename_suffix_pos, &
       filename_split
  use mesh, only : mesh_t
  use field, only : field_t, field_ptr_t
  use field_list, only : field_list_t
  use field_series, only : field_series_t, field_series_ptr_t
  use dofmap, only : dofmap_t
  use space, only : neko_space_t => space_t, GLL
  use interpolation, only : interpolator_t
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_DEBUG
  use vector, only : vector_t
  use matrix, only : matrix_t
  use datadist, only : linear_dist_t
  use comm, only : pe_rank, pe_size, NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_INFO_NULL, MPI_Allreduce, MPI_Allgather, &
       MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, MPI_MAX, MPI_Comm_size, MPI_Exscan, &
       MPI_Barrier, MPI_INTEGER8, MPI_Scan, MPI_Bcast
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
     !> Read data in HDF5 format.
     procedure :: read => hdf5_file_read
     !> Write data in HDF5 format.
     procedure :: write => hdf5_file_write
     procedure :: set_overwrite => hdf5_file_set_overwrite
     ! Granular methods for dealing with HDF5 files
     procedure :: open => hdf5_file_open
     procedure :: close => hdf5_file_close
     procedure :: set_active_group => hdf5_file_set_group
     procedure :: set_precision => hdf5_file_set_precision
     procedure, pass(this) :: write_vector => hdf5_file_write_vector
     procedure, pass(this) :: write_matrix => hdf5_file_write_matrix
     procedure, pass(this) :: write_field => hdf5_file_write_field
     procedure, pass(this) :: write_int_attribute => &
          hdf5_file_write_int_attribute
     procedure, pass(this) :: write_rp_attribute => hdf5_file_write_rp_attribute
     procedure, pass(this) :: read_vector => hdf5_file_read_vector
     procedure, pass(this) :: read_matrix => hdf5_file_read_matrix
     procedure, pass(this) :: read_int_attribute => hdf5_file_read_int_attribute
     procedure, pass(this) :: read_rp_attribute => hdf5_file_read_rp_attribute
     procedure :: write_dataset => hdf5_file_write_dataset
     procedure :: read_dataset => hdf5_file_read_dataset
     procedure :: write_attribute => hdf5_file_write_attribute
     procedure :: read_attribute => hdf5_file_read_attribute
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

  !> Write data in HDF5 format.
  !! @param data Data object to write.
  !! @param[optional] t Simulation time.
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
    integer(hid_t) :: plist_id, access_plist_id
    integer(hid_t) :: file_id, dset_id, grp_id, attr_id
    integer(hid_t) :: filespace, dataset_space, memspace
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: suffix_pos
    character(len=5) :: id_str
    character(len=1024) :: fname
    logical :: checkpoint_data

    call hdf5_file_determine_data(data, msh, dof, fp, fsp, dtlag, tlag)
    checkpoint_data = .false.
    select type (data)
    type is (chkp_t)
       checkpoint_data = .true.
    end select

    if (.not. this%overwrite) call this%increment_counter()
    fname = trim(this%get_fname())

    ! h5open_f is idempotent; the process-wide session owns finalization.
    call h5open_f(ierr)
    call hdf5_file_determine_real(H5T_NEKO_REAL)

    call h5pcreate_f(H5P_FILE_ACCESS_F, access_plist_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(access_plist_id, NEKO_COMM%mpi_val, info, ierr)

    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, &
         file_id, ierr, access_prp = access_plist_id)
    call h5pclose_f(access_plist_id, ierr)

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
       call h5gcreate_f(file_id, "Mesh", grp_id, ierr, &
            lcpl_id = h5p_default_f, gcpl_id = h5p_default_f, &
            gapl_id = h5p_default_f)

       call h5acreate_f(grp_id, "Elements", H5T_NATIVE_INTEGER, filespace, &
            attr_id, ierr, h5p_default_f, h5p_default_f)
       call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, msh%glb_nelv, ddim, ierr)
       call h5aclose_f(attr_id, ierr)

       call h5acreate_f(grp_id, "Dimension", H5T_NATIVE_INTEGER, filespace, &
            attr_id, ierr, h5p_default_f, h5p_default_f)
       call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, msh%gdim, ddim, ierr)
       call h5aclose_f(attr_id, ierr)

       call h5gclose_f(grp_id, ierr)
    end if


    call h5sclose_f(filespace, ierr)

    !
    ! Write restart group (tlag, dtlag)
    !
    if (associated(tlag) .and. associated(dtlag)) then
       call h5gcreate_f(file_id, "Restart", grp_id, ierr, &
            lcpl_id = h5p_default_f, gcpl_id = h5p_default_f, &
            gapl_id = h5p_default_f)

       drank = 1
       ddim = size(tlag)
       doffset(1) = 0
       if (pe_rank .eq. 0) then
          dcount = size(tlag)
       else
          dcount = 0
       end if

       call h5screate_simple_f(drank, ddim, filespace, ierr)

       call h5dcreate_f(grp_id, 'tlag', H5T_NEKO_REAL, &
            filespace, dset_id, ierr)
       call h5dget_space_f(dset_id, dataset_space, ierr)
       call h5sselect_hyperslab_f (dataset_space, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dwrite_f(dset_id, H5T_NEKO_REAL, tlag, &
            ddim, ierr, xfer_prp = plist_id)
       call h5sclose_f(dataset_space, ierr)
       call h5dclose_f(dset_id, ierr)

       call h5dcreate_f(grp_id, 'dtlag', H5T_NEKO_REAL, &
            filespace, dset_id, ierr)
       call h5dget_space_f(dset_id, dataset_space, ierr)
       call h5sselect_hyperslab_f (dataset_space, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dwrite_f(dset_id, H5T_NEKO_REAL, dtlag, &
            ddim, ierr, xfer_prp = plist_id)
       call h5sclose_f(dataset_space, ierr)
       call h5dclose_f(dset_id, ierr)

       call h5sclose_f(filespace, ierr)
       call h5gclose_f(grp_id, ierr)

    end if


    !
    ! Write fields group
    !
    if (checkpoint_data) then
       select type (data)
       type is (chkp_t)
          call hdf5_checkpoint_write_payloads( &
               file_id, plist_id, H5T_NEKO_REAL, data)
       end select
       if (allocated(fp)) deallocate(fp)
       if (allocated(fsp)) deallocate(fsp)
    else if (allocated(fp) .or. allocated(fsp)) then
       call h5gcreate_f(file_id, "Fields", grp_id, ierr, &
            lcpl_id = h5p_default_f, gcpl_id = h5p_default_f, &
            gapl_id = h5p_default_f)

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
             call h5dget_space_f(dset_id, dataset_space, ierr)
             call h5sselect_hyperslab_f(dataset_space, H5S_SELECT_SET_F, &
                  doffset, dcount, ierr)
             call h5dwrite_f(dset_id, H5T_NEKO_REAL, &
                  fp(i)%ptr%x(1,1,1,1), &
                  ddim, ierr, file_space_id = dataset_space, &
                  mem_space_id = memspace, xfer_prp = plist_id)
             call h5sclose_f(dataset_space, ierr)
             call h5dclose_f(dset_id, ierr)
          end do
          deallocate(fp)
       end if

       if (allocated(fsp)) then
          do i = 1, size(fsp)
             do j = 1, fsp(i)%ptr%size()
                call h5dcreate_f(grp_id, fsp(i)%ptr%lf(j)%name, &
                     H5T_NEKO_REAL, filespace, dset_id, ierr)
                call h5dget_space_f(dset_id, dataset_space, ierr)
                call h5sselect_hyperslab_f(dataset_space, H5S_SELECT_SET_F, &
                     doffset, dcount, ierr)
                call h5dwrite_f(dset_id, H5T_NEKO_REAL, &
                     fsp(i)%ptr%lf(j)%x(1,1,1,1), &
                     ddim, ierr, file_space_id = dataset_space, &
                     mem_space_id = memspace, xfer_prp = plist_id)
                call h5sclose_f(dataset_space, ierr)
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

  end subroutine hdf5_file_write

  !> Read data in HDF5 format.
  !! @param data Data object to populate.
  subroutine hdf5_file_read(this, data)
    class(hdf5_file_t) :: this
    class(*), target, intent(inout) :: data
    integer(hid_t) :: plist_id, access_plist_id
    integer(hid_t) :: file_id, dset_id, grp_id, attr_id
    integer(hid_t) :: filespace, memspace
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: i,j, ierr, info, glb_nelv, gdim, lx, drank
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dof
    type(field_ptr_t), allocatable :: fp(:)
    type(field_series_ptr_t), allocatable :: fsp(:)
    type(neko_space_t), target :: source_Xh
    real(kind=rp), pointer :: dtlag(:)
    real(kind=rp), pointer :: tlag(:)
    real(kind=rp) :: t
    character(len=1024) :: fname
    logical :: payloads_exist

    fname = trim(this%get_fname())

    ! h5open_f is idempotent; the process-wide session owns finalization.
    call h5open_f(ierr)
    call hdf5_file_determine_data(data, msh, dof, fp, fsp, dtlag, tlag)
    call hdf5_file_determine_real(H5T_NEKO_REAL)

    call h5pcreate_f(H5P_FILE_ACCESS_F, access_plist_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(access_plist_id, NEKO_COMM%mpi_val, info, ierr)

    call h5fopen_f(fname, H5F_ACC_RDONLY_F, &
         file_id, ierr, access_prp = access_plist_id)
    call h5pclose_f(access_plist_id, ierr)

    call h5lexists_f(file_id, 'Payloads', payloads_exist, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ddim = 1
    call h5aopen_name_f(file_id, 'Time', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NEKO_REAL, t, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    select type (data)
    type is (chkp_t)
       data%t = t
    end select

    call h5aopen_name_f(file_id, 'Lx', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, lx, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5gopen_f(file_id, 'Mesh', grp_id, ierr, gapl_id = h5p_default_f)

    call h5aopen_name_f(grp_id, 'Elements', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, glb_nelv, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_name_f(grp_id, 'Dimension', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, gdim, ddim, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5gclose_f(grp_id, ierr)

    if (glb_nelv .ne. msh%glb_nelv .or. gdim .ne. msh%gdim) then
       call neko_error("HDF5 file mesh does not match the case")
    end if
    call source_Xh%init(GLL, lx, lx, lx)

    select type (data)
    type is (chkp_t)
       call data%previous_Xh%init(GLL, lx, lx)
    end select

    if (associated(tlag) .and. associated(dtlag)) then
       drank = 1
       ddim = size(tlag)
       doffset(1) = 0
       if (pe_rank .eq. 0) then
          dcount = size(tlag)
       else
          dcount = 0
       end if

       call h5gopen_f(file_id, 'Restart', grp_id, ierr, &
            gapl_id = h5p_default_f)
       call h5dopen_f(grp_id, 'tlag', dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dread_f(dset_id, H5T_NEKO_REAL, tlag, ddim, ierr, &
            xfer_prp = plist_id)
       call h5dclose_f(dset_id, ierr)
       call h5sclose_f(filespace, ierr)

       call h5dopen_f(grp_id, 'dtlag', dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dread_f(dset_id, H5T_NEKO_REAL, dtlag, ddim, ierr, &
            xfer_prp = plist_id)
       call h5dclose_f(dset_id, ierr)
       call h5sclose_f(filespace, ierr)

       call h5gclose_f(grp_id, ierr)
    end if

    if (payloads_exist) then
       select type (data)
       type is (chkp_t)
          call hdf5_checkpoint_read_payloads( &
               file_id, plist_id, H5T_NEKO_REAL, data, source_Xh)
       class default
          call neko_error( &
               "HDF5 payload checkpoints require a checkpoint object")
       end select
       if (allocated(fp)) deallocate(fp)
       if (allocated(fsp)) deallocate(fsp)
    else if (allocated(fp) .or. allocated(fsp)) then
       call h5gopen_f(file_id, 'Fields', grp_id, ierr, gapl_id = h5p_default_f)

       if (allocated(fp)) then
          do i = 1, size(fp)
             call hdf5_checkpoint_read_field(grp_id, plist_id, &
                  H5T_NEKO_REAL, fp(i)%ptr, source_Xh)
          end do
       end if

       if (allocated(fsp)) then
          do i = 1, size(fsp)
             do j = 1, fsp(i)%ptr%size()
                call hdf5_checkpoint_read_field(grp_id, plist_id, &
                     H5T_NEKO_REAL, fsp(i)%ptr%lf(j), source_Xh)
             end do
          end do
       end if
       call h5gclose_f(grp_id, ierr)
    end if

    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)
    call source_Xh%free()

  end subroutine hdf5_file_read

  !> Write the format-independent checkpoint payload hierarchy.
  !! @param file_id Open HDF5 checkpoint file.
  !! @param plist_id Collective dataset transfer property list.
  !! @param h5_neko_real HDF5 datatype corresponding to `rp`.
  !! @param chkp Registered checkpoint payloads to write.
  subroutine hdf5_checkpoint_write_payloads(file_id, plist_id, &
       h5_neko_real, chkp)
    integer(hid_t), intent(in) :: file_id, plist_id, h5_neko_real
    type(chkp_t), intent(in) :: chkp
    integer(hid_t) :: payloads_id, payload_id
    integer :: i, j, k, ierr

    call h5gcreate_f(file_id, "Payloads", payloads_id, ierr)

    do i = 1, chkp%payload_count()
       call hdf5_checkpoint_open_group(payloads_id, &
            chkp%payloads(i)%ptr%name, .true., payload_id)

       do j = 1, chkp%payloads(i)%ptr%field_count()
          call hdf5_checkpoint_write_field(payload_id, plist_id, &
               h5_neko_real, chkp%payloads(i)%ptr%fields(j)%ptr)
       end do

       do j = 1, chkp%payloads(i)%ptr%series_count()
          do k = 1, chkp%payloads(i)%ptr%series(j)%ptr%size()
             call hdf5_checkpoint_write_field(payload_id, plist_id, &
                  h5_neko_real, &
                  chkp%payloads(i)%ptr%series(j)%ptr%lf(k))
          end do
       end do

       do j = 1, chkp%payloads(i)%ptr%mesh_array_count()
          call hdf5_checkpoint_write_mesh_array(payload_id, plist_id, &
               h5_neko_real, chkp%payloads(i)%ptr%mesh_arrays(j)%ptr)
       end do

       do j = 1, chkp%payloads(i)%ptr%array_count()
          call hdf5_checkpoint_write_array(payload_id, plist_id, &
               h5_neko_real, chkp%payloads(i)%ptr%arrays(j)%ptr)
       end do

       call h5gclose_f(payload_id, ierr)
    end do

    call h5gclose_f(payloads_id, ierr)

  end subroutine hdf5_checkpoint_write_payloads

  !> Read the format-independent checkpoint payload hierarchy.
  !! @param file_id Open HDF5 checkpoint file.
  !! @param plist_id Collective dataset transfer property list.
  !! @param h5_neko_real HDF5 datatype corresponding to `rp`.
  !! @param chkp Registered checkpoint payloads to populate.
  !! @param source_Xh Function space stored in the HDF5 file.
  subroutine hdf5_checkpoint_read_payloads(file_id, plist_id, &
       h5_neko_real, chkp, source_Xh)
    integer(hid_t), intent(in) :: file_id, plist_id, h5_neko_real
    type(chkp_t), intent(inout) :: chkp
    type(neko_space_t), target, intent(inout) :: source_Xh
    integer(hid_t) :: payloads_id, payload_id
    integer :: i, j, k, ierr

    call h5gopen_f(file_id, "Payloads", payloads_id, ierr)

    do i = 1, chkp%payload_count()
       call hdf5_checkpoint_open_group(payloads_id, &
            chkp%payloads(i)%ptr%name, .false., payload_id)

       do j = 1, chkp%payloads(i)%ptr%field_count()
          call hdf5_checkpoint_read_field(payload_id, plist_id, &
               h5_neko_real, chkp%payloads(i)%ptr%fields(j)%ptr, &
               source_Xh)
       end do

       do j = 1, chkp%payloads(i)%ptr%series_count()
          do k = 1, chkp%payloads(i)%ptr%series(j)%ptr%size()
             call hdf5_checkpoint_read_field(payload_id, plist_id, &
                  h5_neko_real, &
                  chkp%payloads(i)%ptr%series(j)%ptr%lf(k), source_Xh)
          end do
       end do

       do j = 1, chkp%payloads(i)%ptr%mesh_array_count()
          call hdf5_checkpoint_read_mesh_array(payload_id, plist_id, &
               h5_neko_real, chkp%payloads(i)%ptr%mesh_arrays(j)%ptr, &
               source_Xh)
       end do

       do j = 1, chkp%payloads(i)%ptr%array_count()
          call hdf5_checkpoint_read_array(payload_id, plist_id, &
               h5_neko_real, chkp%payloads(i)%ptr%arrays(j)%ptr)
       end do

       call h5gclose_f(payload_id, ierr)
    end do

    call h5gclose_f(payloads_id, ierr)

  end subroutine hdf5_checkpoint_read_payloads

  !> Open a nested group path, optionally creating missing groups.
  !! @param root_id HDF5 group from which the path is resolved.
  !! @param path Slash-separated relative group path.
  !! @param create Whether missing groups should be created.
  !! @param group_id Open HDF5 group at `path`.
  subroutine hdf5_checkpoint_open_group(root_id, path, create, group_id)
    integer(hid_t), intent(in) :: root_id
    character(len=*), intent(in) :: path
    logical, intent(in) :: create
    integer(hid_t), intent(out) :: group_id
    integer(hid_t) :: current_id, next_id
    integer :: first, last, slash, path_len, ierr
    logical :: group_exists
    character(len=:), allocatable :: group_name

    current_id = root_id
    first = 1
    path_len = len_trim(path)

    do
       slash = index(path(first:path_len), "/")
       if (slash .eq. 0) then
          last = path_len
       else
          last = first + slash - 2
       end if

       group_name = path(first:last)
       call h5lexists_f(current_id, group_name, group_exists, ierr)
       if (group_exists) then
          call h5gopen_f(current_id, group_name, next_id, ierr)
       else if (create) then
          call h5gcreate_f(current_id, group_name, next_id, ierr)
       else
          call neko_error("Checkpoint payload group '" // trim(path) // &
               "' is missing")
       end if

       if (current_id .ne. root_id) call h5gclose_f(current_id, ierr)
       current_id = next_id

       if (slash .eq. 0) exit
       first = last + 2
    end do

    group_id = current_id

  end subroutine hdf5_checkpoint_open_group

  !> Write one checkpoint field under its native field name.
  !! @param group_id Open HDF5 payload group.
  !! @param plist_id Collective dataset transfer property list.
  !! @param h5_neko_real HDF5 datatype corresponding to `rp`.
  !! @param fld Field to write.
  subroutine hdf5_checkpoint_write_field(group_id, plist_id, &
       h5_neko_real, fld)
    integer(hid_t), intent(in) :: group_id, plist_id, h5_neko_real
    type(field_t), intent(in) :: fld
    integer(hid_t) :: dset_id, filespace, memspace
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: ierr

    dcount(1) = int(fld%dof%size(), 8)
    doffset(1) = int(fld%msh%offset_el, 8) * int(fld%Xh%lxyz, 8)
    ddim = dcount
    call MPI_Allreduce(MPI_IN_PLACE, ddim(1), 1, &
         MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)

    call h5screate_simple_f(1, ddim, filespace, ierr)
    call h5screate_simple_f(1, dcount, memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset, dcount, ierr)
    call h5dcreate_f(group_id, trim(fld%name), h5_neko_real, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, h5_neko_real, fld%x(1,1,1,1), &
         dcount, ierr, file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)

  end subroutine hdf5_checkpoint_write_field

  !> Read one checkpoint field by its native field name.
  !! @param group_id Open HDF5 payload group.
  !! @param plist_id Collective dataset transfer property list.
  !! @param h5_neko_real HDF5 datatype corresponding to `rp`.
  !! @param fld Field to populate.
  !! @param source_Xh Function space stored in the HDF5 file.
  subroutine hdf5_checkpoint_read_field(group_id, plist_id, &
       h5_neko_real, fld, source_Xh)
    integer(hid_t), intent(in) :: group_id, plist_id, h5_neko_real
    type(field_t), intent(inout) :: fld
    type(neko_space_t), target, intent(inout) :: source_Xh
    integer(hid_t) :: dset_id, filespace, memspace
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    real(kind=rp), allocatable :: checkpoint_data(:)
    type(interpolator_t) :: space_interp
    integer :: ierr
    logical :: dataset_exists

    call h5lexists_f(group_id, trim(fld%name), dataset_exists, ierr)
    if (.not. dataset_exists) then
       call neko_error("Checkpoint field '" // trim(fld%name) // &
            "' is missing")
    end if

    call h5dopen_f(group_id, trim(fld%name), dset_id, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)

    if (source_Xh%lxyz .eq. fld%Xh%lxyz) then
       dcount(1) = int(fld%dof%size(), hsize_t)
       doffset(1) = int(fld%msh%offset_el, hsize_t) * &
            int(fld%Xh%lxyz, hsize_t)
       call h5screate_simple_f(1, dcount, memspace, ierr)
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dread_f(dset_id, h5_neko_real, fld%x(1,1,1,1), &
            dcount, ierr, file_space_id = filespace, &
            mem_space_id = memspace, xfer_prp = plist_id)
    else
       ! Interpolation across spaces required!
       dcount(1) = int(fld%msh%nelv, hsize_t) * &
            int(source_Xh%lxyz, hsize_t)
       doffset(1) = int(fld%msh%offset_el, hsize_t) * &
            int(source_Xh%lxyz, hsize_t)
       ddim = dcount
       allocate(checkpoint_data(int(dcount(1))))
       call h5screate_simple_f(1, dcount, memspace, ierr)
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
       call h5dread_f(dset_id, h5_neko_real, checkpoint_data, ddim, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = plist_id)

       call space_interp%init(fld%Xh, source_Xh)
       call space_interp%map_host(fld%x, checkpoint_data, fld%msh%nelv, &
            fld%Xh)
       call space_interp%free()
       deallocate(checkpoint_data)
    end if

    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)

  end subroutine hdf5_checkpoint_read_field

  !> Write one nodal mesh array using its mesh distribution.
  !! @param group_id Open HDF5 payload group.
  !! @param plist_id Collective dataset transfer property list.
  !! @param h5_neko_real HDF5 datatype corresponding to `rp`.
  !! @param array Mesh-array descriptor to write.
  subroutine hdf5_checkpoint_write_mesh_array(group_id, plist_id, &
       h5_neko_real, array)
    integer(hid_t), intent(in) :: group_id, plist_id, h5_neko_real
    type(checkpoint_mesh_array_t), intent(in) :: array
    integer(hid_t) :: dset_id, filespace, memspace, attr_id, attr_space
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset, attr_dims
    integer :: nodal_shape(3)
    integer :: ierr

    dcount(1) = int(size(array%x), hsize_t)
    doffset(1) = int(array%msh%offset_el, hsize_t) * &
         int(array%Xh%lxyz, hsize_t)
    ddim(1) = int(array%msh%glb_nelv, hsize_t) * &
         int(array%Xh%lxyz, hsize_t)

    call h5screate_simple_f(1, ddim, filespace, ierr)
    call h5screate_simple_f(1, dcount, memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset, dcount, ierr)
    call h5dcreate_f(group_id, trim(array%name), h5_neko_real, &
         filespace, dset_id, ierr)

    nodal_shape = [array%Xh%lx, array%Xh%ly, array%Xh%lz]
    attr_dims(1) = size(nodal_shape)
    call h5screate_simple_f(1, attr_dims, attr_space, ierr)
    call h5acreate_f(dset_id, "NodalShape", H5T_NATIVE_INTEGER, &
         attr_space, attr_id, ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, nodal_shape, &
         attr_dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(attr_space, ierr)

    call h5dwrite_f(dset_id, h5_neko_real, array%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)

  end subroutine hdf5_checkpoint_write_mesh_array

  !> Read and optionally interpolate one nodal mesh array.
  !! @param group_id Open HDF5 payload group.
  !! @param plist_id Collective dataset transfer property list.
  !! @param h5_neko_real HDF5 datatype corresponding to `rp`.
  !! @param array Mesh-array descriptor to populate.
  !! @param source_Xh Default function space stored in the HDF5 file.
  subroutine hdf5_checkpoint_read_mesh_array(group_id, plist_id, &
       h5_neko_real, array, source_Xh)
    integer(hid_t), intent(in) :: group_id, plist_id, h5_neko_real
    type(checkpoint_mesh_array_t), intent(inout) :: array
    type(neko_space_t), target, intent(inout) :: source_Xh
    integer(hid_t) :: dset_id, filespace, memspace, attr_id, attr_space
    integer(hsize_t), dimension(1) :: dcount, doffset
    integer(hsize_t), dimension(1) :: dataset_dims, dataset_maxdims
    integer(hsize_t), dimension(1) :: attr_dims, attr_maxdims
    integer :: stored_shape(3)
    real(kind=rp), allocatable :: stored_data(:)
    type(neko_space_t), target :: stored_Xh
    type(interpolator_t) :: space_interp
    integer :: ierr
    logical :: dataset_exists, shape_exists

    call h5lexists_f(group_id, trim(array%name), dataset_exists, ierr)
    if (.not. dataset_exists) then
       call neko_error("Checkpoint mesh array '" // trim(array%name) // &
            "' is missing")
    end if

    call h5dopen_f(group_id, trim(array%name), dset_id, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)
    call h5sget_simple_extent_dims_f(filespace, dataset_dims, &
         dataset_maxdims, ierr)

    stored_shape = [source_Xh%lx, source_Xh%ly, source_Xh%lz]
    call h5aexists_f(dset_id, "NodalShape", shape_exists, ierr)
    if (shape_exists) then
       call h5aopen_f(dset_id, "NodalShape", attr_id, ierr)
    else
       ! Compatibility with development checkpoints using generic blocks.
       call h5aexists_f(dset_id, "BlockShape", shape_exists, ierr)
       if (shape_exists) then
          call h5aopen_f(dset_id, "BlockShape", attr_id, ierr)
       end if
    end if
    if (shape_exists) then
       call h5aget_space_f(attr_id, attr_space, ierr)
       call h5sget_simple_extent_dims_f(attr_space, attr_dims, &
            attr_maxdims, ierr)
       if (attr_dims(1) .ne. 3) then
          call neko_error("Checkpoint mesh-array nodal shape must have " // &
               "three dimensions")
       end if
       call h5aread_f(attr_id, H5T_NATIVE_INTEGER, stored_shape, &
            attr_dims, ierr)
       call h5sclose_f(attr_space, ierr)
       call h5aclose_f(attr_id, ierr)
    end if

    if (any(stored_shape .le. 0)) then
       call neko_error("Checkpoint mesh-array nodal shape must be positive")
    end if
    if (array%msh%gdim .eq. 3) then
       call stored_Xh%init(GLL, stored_shape(1), stored_shape(2), &
            stored_shape(3))
    else
       if (stored_shape(3) .ne. 1) then
          call neko_error("Two-dimensional checkpoint mesh arrays must " // &
               "have one nodal plane")
       end if
       call stored_Xh%init(GLL, stored_shape(1), stored_shape(2))
    end if

    if (int(dataset_dims(1), i8) .ne. &
         int(array%msh%glb_nelv, i8) * int(stored_Xh%lxyz, i8)) then
       call neko_error("Checkpoint mesh array '" // trim(array%name) // &
            "' does not match the current mesh")
    end if

    dcount(1) = int(array%msh%nelv, hsize_t) * &
         int(stored_Xh%lxyz, hsize_t)
    doffset(1) = int(array%msh%offset_el, hsize_t) * &
         int(stored_Xh%lxyz, hsize_t)
    call h5screate_simple_f(1, dcount, memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset, dcount, ierr)

    if (stored_Xh%lxyz .eq. array%Xh%lxyz) then
       call h5dread_f(dset_id, h5_neko_real, array%x, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = plist_id)
    else
       allocate(stored_data(int(dcount(1))))
       call h5dread_f(dset_id, h5_neko_real, stored_data, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = plist_id)
       call space_interp%init(array%Xh, stored_Xh)
       call space_interp%map_host(array%x, stored_data, array%msh%nelv, &
            array%Xh)
       call space_interp%free()
       deallocate(stored_data)
    end if

    call stored_Xh%free()
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)

  end subroutine hdf5_checkpoint_read_mesh_array

  !> Write one generic real array using its caller-provided distribution.
  !! @param group_id Open HDF5 payload group.
  !! @param plist_id Collective dataset transfer property list.
  !! @param h5_neko_real HDF5 datatype corresponding to `rp`.
  !! @param array Array descriptor to write.
  subroutine hdf5_checkpoint_write_array(group_id, plist_id, &
       h5_neko_real, array)
    integer(hid_t), intent(in) :: group_id, plist_id, h5_neko_real
    type(checkpoint_array_t), intent(in) :: array
    integer(hid_t) :: dset_id, filespace, memspace
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: ierr

    ddim(1) = int(array%global_count, hsize_t)
    dcount(1) = int(size(array%x), hsize_t)
    doffset(1) = int(array%offset, hsize_t)

    call h5screate_simple_f(1, ddim, filespace, ierr)
    call h5screate_simple_f(1, dcount, memspace, ierr)
    call h5dcreate_f(group_id, trim(array%name), h5_neko_real, &
         filespace, dset_id, ierr)

    if (array%replicated .and. pe_rank .ne. 0) then
       call h5sselect_none_f(filespace, ierr)
       call h5sselect_none_f(memspace, ierr)
    else
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
    end if

    call h5dwrite_f(dset_id, h5_neko_real, array%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)

  end subroutine hdf5_checkpoint_write_array

  !> Read one generic real array using its caller-provided distribution.
  !! @param group_id Open HDF5 payload group.
  !! @param plist_id Collective dataset transfer property list.
  !! @param h5_neko_real HDF5 datatype corresponding to `rp`.
  !! @param array Array descriptor to populate.
  subroutine hdf5_checkpoint_read_array(group_id, plist_id, &
       h5_neko_real, array)
    integer(hid_t), intent(in) :: group_id, plist_id, h5_neko_real
    type(checkpoint_array_t), intent(inout) :: array
    integer(hid_t) :: dset_id, filespace, memspace
    integer(hsize_t), dimension(1) :: dcount, doffset
    integer(hsize_t), dimension(1) :: dataset_dims, dataset_maxdims
    integer :: ierr
    logical :: dataset_exists

    call h5lexists_f(group_id, trim(array%name), dataset_exists, ierr)
    if (.not. dataset_exists) then
       call neko_error("Checkpoint array '" // trim(array%name) // &
            "' is missing")
    end if

    call h5dopen_f(group_id, trim(array%name), dset_id, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)
    call h5sget_simple_extent_dims_f(filespace, dataset_dims, &
         dataset_maxdims, ierr)

    if (int(dataset_dims(1), i8) .ne. array%global_count) then
       call neko_error("Checkpoint array '" // trim(array%name) // &
            "' has an incompatible global extent")
    end if

    dcount(1) = int(size(array%x), hsize_t)
    doffset(1) = int(array%offset, hsize_t)
    call h5screate_simple_f(1, dcount, memspace, ierr)

    if (array%replicated .and. pe_rank .ne. 0) then
       call h5sselect_none_f(filespace, ierr)
       call h5sselect_none_f(memspace, ierr)
    else
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
            doffset, dcount, ierr)
    end if
    call h5dread_f(dset_id, h5_neko_real, array%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)

    if (array%replicated) then
       call MPI_Bcast(array%x, size(array%x), MPI_REAL_PRECISION, 0, &
            NEKO_COMM, ierr)
    end if

  end subroutine hdf5_checkpoint_read_array

  !> Determine mesh, fields, and histories represented by a data object.
  !! @param data Data object to inspect.
  !! @param msh Associated mesh, if available.
  !! @param dof Associated degree-of-freedom map, if available.
  !! @param fp Individual fields represented by `data`.
  !! @param fsp Field series represented by `data`.
  !! @param dtlag Previous time-step sizes, if available.
  !! @param tlag Previous simulation times, if available.
  subroutine hdf5_file_determine_data(data, msh, dof, fp, fsp, dtlag, tlag)
    class(*), target, intent(in) :: data
    type(mesh_t), pointer, intent(inout) :: msh
    type(dofmap_t), pointer, intent(inout) :: dof
    type(field_ptr_t), allocatable, intent(inout) :: fp(:)
    type(field_series_ptr_t), allocatable, intent(inout) :: fsp(:)
    real(kind=rp), pointer, intent(inout) :: dtlag(:)
    real(kind=rp), pointer, intent(inout) :: tlag(:)
    integer :: i, j, fp_size, fp_cur, fsp_size, fsp_cur

    select type (data)
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

    type is (chkp_t)
       block
         type(checkpoint_payload_t), pointer :: fluid
         type(field_t), pointer :: u

         fluid => data%get_payload("fluid")
         u => fluid%find_field("u")
         if (.not. associated(u) .or. &
              .not. associated(fluid%find_field("v")) .or. &
              .not. associated(fluid%find_field("w")) .or. &
              .not. associated(fluid%find_field("p"))) then
            call neko_error("Checkpoint not initialized")
         end if
         dof => u%dof
         msh => u%msh

         fp_size = 0
         fsp_size = 0
         do i = 1, data%payload_count()
            fp_size = fp_size + data%payloads(i)%ptr%field_count()
            fsp_size = fsp_size + data%payloads(i)%ptr%series_count()
         end do

         if (fp_size .gt. 0) allocate(fp(fp_size))
         if (fsp_size .gt. 0) allocate(fsp(fsp_size))

         fp_cur = 1
         fsp_cur = 1
         do i = 1, data%payload_count()
            do j = 1, data%payloads(i)%ptr%field_count()
               fp(fp_cur)%ptr => data%payloads(i)%ptr%fields(j)%ptr
               fp_cur = fp_cur + 1
            end do
            do j = 1, data%payloads(i)%ptr%series_count()
               fsp(fsp_cur)%ptr => data%payloads(i)%ptr%series(j)%ptr
               fsp_cur = fsp_cur + 1
            end do
         end do

         call data%get_time_history(tlag, dtlag)
       end block

    class default
       call neko_log%error('Invalid data')
    end select

  end subroutine hdf5_file_determine_data

  !> Determine hdf5 real type corresponding to NEKO_REAL
  !! @note This must be called after h5open_f, otherwise
  !! the H5T_NATIVE_XYZ types has a value of 0
  subroutine hdf5_file_determine_real(H5T_NEKO_REAL)
    integer(hid_t), intent(inout) :: H5T_NEKO_REAL
    select case (rp)
    case (dp)
       H5T_NEKO_REAL = H5T_NATIVE_DOUBLE
    case (sp)
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
    character(len=LOG_SIZE) :: log_buf

    ! Set the mode for the file
    this%mode = mode

    ! Ensure precision is set and are valid.
    if (this%precision .gt. rp) then
       this%precision = rp
       call neko_warning('Requested precision is higher than working precision')
    else if (this%precision .eq. -1) then
       this%precision = rp
    end if

    fname = trim(file_get_fname(this))
    counter = this%get_counter() - this%get_start_counter()

    ! Set the configuration for MPI IO
    call h5open_f(ierr)

    mpi_info = MPI_INFO_NULL%mpi_val
    mpi_comm = NEKO_COMM%mpi_val
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

    write (log_buf, *) "Opened HDF5 file: ", trim(fname), " with counter: ", &
         counter
    call neko_log%message(log_buf, lvl = NEKO_LOG_DEBUG)

  end subroutine hdf5_file_open

  !> Close the file
  subroutine hdf5_file_close(this)
    class(hdf5_file_t), intent(inout) :: this
    integer :: ierr

    if (this%active_group_id .ne. -1_hid_t .and. &
         this%active_group_id .ne. this%file_id) then
       call h5gclose_f(this%active_group_id, ierr)
    end if
    this%active_group_id = -1_hid_t

    call h5pclose_f(this%plist_id, ierr)
    this%plist_id = -1_hid_t
    call h5fclose_f(this%file_id, ierr)
    this%file_id = -1_hid_t
    call h5close_f(ierr)

    call neko_log%message("Closed HDF5 file: " // trim(this%get_fname()), &
         lvl = NEKO_LOG_DEBUG)

  end subroutine hdf5_file_close


  !> Set the active group for HDF5 files from an input string
  !! @param this The HDF5 file object
  !! @param A string indicating the path to the group to create or open.
  subroutine hdf5_file_set_group(this, group_name_path)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in), optional :: group_name_path
    character(len=1000), allocatable :: group_name(:)

    integer(hid_t) :: current_id, group_id
    integer :: ierr, i, j, num_groups, name_len, group_loc
    logical :: group_exists


    ! Close previous active group if one is open
    if (this%active_group_id .ne. -1_hid_t .and. this%active_group_id .ne. &
         this%file_id) then
       call h5gclose_f(this%active_group_id, ierr)
    end if
    this%active_group_id = -1_hid_t

    ! Start from root location = file
    current_id = this%file_id
    ! Return the root directory if no group name is given
    if (.not. present(group_name_path)) then
       this%active_group_id = current_id
       return
    end if

    ! Split the input string into group names using "/" as a delimiter
    name_len = len(trim(group_name_path))
    ! Count how many groups
    num_groups = 1 ! There is at least one group if this was passed
    do i = 1, name_len
       if (group_name_path .eq. "/") then
          num_groups = num_groups + 1
       end if
    end do

    ! Allocate the group array and populate it
    allocate(group_name(num_groups))
    j = 1
    group_loc = 1
    do i = 1, name_len
       if (group_name_path .eq. "/") then
          group_name(group_loc) = group_name_path(j:i-1)
          group_loc = group_loc + 1
          j = i + 1
       end if
    end do
    if (j .ne. name_len) then
       group_name(group_loc) = group_name_path(j:name_len)
    end if

    ! Iterate over the groups in the path
    do i = 1, num_groups
       call h5lexists_f(current_id, trim(group_name(i)), group_exists, ierr)

       ! Only create groups if they dont exist and we are in write mode "w"
       if (group_exists) then
          call h5gopen_f(current_id, trim(group_name(i)), group_id, ierr)
       else
          if (this%mode == "r") then
             call neko_error("Group " // trim(group_name(i)) // &
                  " does not exist in file " // trim(file_get_fname(this)))
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

  subroutine hdf5_file_read_dataset(this, data_name, data, strategy)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: data_name
    class(*), intent(inout) :: data
    character(len=*), intent(in), optional :: strategy

    select type (d => data)
    type is (vector_t)
       call this%read_vector(data_name, d, strategy)
    type is (matrix_t)
       call this%read_matrix(data_name, d, strategy)
    type is (field_t)
       call neko_error("Reading a field_t is not supported yet")
    class default
       call neko_error("read_dataset not implemented for this data type")
    end select
  end subroutine hdf5_file_read_dataset

  subroutine hdf5_file_write_attribute(this, data_name, data)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: data_name
    class(*), intent(inout) :: data

    select type (d => data)
    type is (integer)
       call this%write_int_attribute(data_name, d)
    type is (real(kind=rp))
       call this%write_rp_attribute(data_name, d)
    class default
       call neko_error("write_attribute not implemented for this data type")
    end select
  end subroutine hdf5_file_write_attribute

  subroutine hdf5_file_read_attribute(this, data_name, data, exist)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: data_name
    class(*), intent(inout) :: data
    logical, intent(inout) :: exist

    select type (d => data)
    type is (integer)
       call this%read_int_attribute(data_name, d, exist)
    type is (real(kind=rp))
       call this%read_rp_attribute(data_name, d, exist)
    class default
       call neko_error("read_attribute not implemented for this data type")
    end select
  end subroutine hdf5_file_read_attribute


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

    ! Enable chunking to be able to append
    chunkdims = [max(int(max_count, hsize_t), 1_hsize_t)]
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
          call h5sget_simple_extent_dims_f(filespace, tempddims, tempmaxddims, &
               ierr)
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

    ! offset for this rank in the global vector
    doffset = [int(offset, hsize_t) + append_offset]
    ! Get the total file space (shape) of the data set
    call h5dget_space_f(dset_id, filespace, ierr)
    ! Get only the slice where my rank writes
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, &
         ierr)
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
    ! global size of the matrix
    ddims = [int(strides, hsize_t), int(total_count, hsize_t)]
    chunkdims = [int(strides, hsize_t), max(int(max_count, hsize_t), 1_hsize_t)]
    ddims_max = [int(strides, hsize_t), H5S_UNLIMITED_F]
    call h5lexists_f(this%active_group_id, trim(mat%name), dset_exists, ierr)
    if (dset_exists) then
       if (this%overwrite) then

          if (pe_rank .eq. 0) then
             call neko_warning("Dataset " // trim(mat%name) // &
                  " already exists and wil be overwritten")
          end if
          ! retrieve the dset id for the existing data set
          call h5dopen_f(this%active_group_id, trim(mat%name), dset_id, ierr)
       else
          call h5dopen_f(this%active_group_id, trim(mat%name), dset_id, ierr)
          call h5dget_space_f(dset_id, filespace, ierr)
          call h5sget_simple_extent_dims_f(filespace, tempddims, tempmaxddims, &
               ierr)
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
    ! local size of the matrix
    dcount = [int(strides, hsize_t), int(counts, hsize_t)]
    ! offset for this rank in the global matrix
    doffset = [0_hsize_t, int(offset, hsize_t) + append_offset]
    ! Get the total file space (shape) of the data set
    call h5dget_space_f(dset_id, filespace, ierr)
    ! Get only the slice where my rank writes
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, &
         ierr)
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
             call neko_warning("Overwriting dataset: " // trim(field%name))
          end if
          call h5dopen_f(this%active_group_id, trim(field%name), dset_id, ierr)
       else
          call h5dopen_f(this%active_group_id, trim(field%name), dset_id, ierr)
          call h5dget_space_f(dset_id, filespace, ierr)
          call h5sget_simple_extent_dims_f(filespace, tempddims, tempmaxddims, &
               ierr)
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
    dcount = [int(stride_ax_1, hsize_t), &
         int(stride_ax_2, hsize_t), &
         int(stride_ax_3, hsize_t), &
         int(counts, hsize_t)] ! local size of the tensor
    doffset = [0_hsize_t, 0_hsize_t, 0_hsize_t, &
         int(offset, hsize_t) + append_offset] ! offset in the global tensor
    ! Get the total file space (shape) of the data set
    call h5dget_space_f(dset_id, filespace, ierr)
    ! Get only the slice where my rank writes
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, &
         ierr)
    ! Create the corresponding memory space (buffer) for my local data
    call h5screate_simple_f(dset_rank, dcount, memspace, ierr)

    ! =======================
    ! Cast and write the data
    ! =======================
    if (this%precision == sp) then
       allocate(write_buffer_sp(field%Xh%lx, field%Xh%ly, field%Xh%lz, &
            field%msh%nelv))
       if (field%msh%nelv > 0) write_buffer_sp = real(field%x, kind=sp)
       ! Write the data
       call h5dwrite_f(dset_id, precision_hdf, write_buffer_sp, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = xf_id)
       deallocate(write_buffer_sp)
    else if (this%precision == dp) then
       allocate(write_buffer_dp(field%Xh%lx, field%Xh%ly, field%Xh%lz, &
            field%msh%nelv))
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
  subroutine hdf5_file_read_vector(this, data_name, vec, strategy)
    class(hdf5_file_t) :: this
    character(len=*), intent(in) :: data_name
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
    call h5lexists_f(this%active_group_id, trim(data_name), dset_exists, ierr)
    if (dset_exists) then
       ! Open the data set
       call h5dopen_f(this%active_group_id, trim(data_name), dset_id, ierr)
       ! Get the current rank of the dataset
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sget_simple_extent_ndims_f(filespace, temprank, ierr)
       if (temprank .ne. 1) then
          call neko_error("Dataset " // trim(data_name) // &
               " is not a rank 1 vector in file " // trim(file_get_fname(this)))
       end if
       ! Get the current shape and close the filespace
       call h5sget_simple_extent_dims_f(filespace, tempddims, tempmaxddims, &
            ierr)
       call h5sclose_f(filespace, ierr)
    else
       call neko_error("Dataset " // trim(data_name) // &
            " does not exist in current group " // trim(file_get_fname(this)))
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
    call MPI_Exscan(counts, offset, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)

    ! ===========================
    ! Set up reading the data set
    ! ===========================
    dset_rank = 1 ! rank 1 array, i.e. a vector
    dcount = [int(counts, hsize_t)] ! local size of the vector
    doffset = [int(offset, hsize_t)] ! offset for this rank in the global vector
    ! Get the total file space (shape) of the data set
    call h5dget_space_f(dset_id, filespace, ierr)
    ! Get only the slice where my rank reads
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, &
         ierr)
    ! Create the corresponding memory space (buffer) for my local data
    call h5screate_simple_f(dset_rank, dcount, memspace, ierr)

    ! =============================
    ! Allocate data. HDF5 will cast
    ! =============================
    call vec%init(counts, trim(data_name)) ! this is rp
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
  subroutine hdf5_file_read_matrix(this, data_name, mat, strategy)
    class(hdf5_file_t) :: this
    character(len=*), intent(in) :: data_name
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
    call h5lexists_f(this%active_group_id, trim(data_name), dset_exists, ierr)
    if (dset_exists) then
       ! Openr the data set
       call h5dopen_f(this%active_group_id, trim(data_name), dset_id, ierr)
       ! Get the current rank of the of the dataset
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sget_simple_extent_ndims_f(filespace, temprank, ierr)
       if (temprank .ne. 2) then
          call neko_error("Dataset " // trim(data_name) // &
               " is not a rank 2 matrix in file " // trim(file_get_fname(this)))
       end if
       ! Get the current shape and close the filespace
       call h5sget_simple_extent_dims_f(filespace, tempddims, tempmaxddims, &
            ierr)
       call h5sclose_f(filespace, ierr)
    else
       call neko_error("Dataset " // trim(data_name) &
            // " does not exist in current group " // &
            trim(file_get_fname(this)))
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
    dcount = [int(tempddims(1), hsize_t), int(counts, hsize_t)] ! local size
    doffset = [0_hsize_t, int(offset, hsize_t)] ! offset in the global matrix
    ! Get the total file space (shape) of the data set
    call h5dget_space_f(dset_id, filespace, ierr)
    ! Get only the slice where my rank reads
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, &
         ierr)
    ! Create the corresponding memory space (buffer) for my local data
    call h5screate_simple_f(dset_rank, dcount, memspace, ierr)

    ! =============================
    ! Allocate data. HDF5 will cast
    ! =============================
    call mat%init(int(tempddims(1)), counts, trim(data_name)) ! this is rp
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
  subroutine hdf5_file_write_int_attribute(this, attr_name, attr)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: attr_name
    integer, intent(in) :: attr
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
       call h5acreate_f(this%active_group_id, trim(attr_name), &
            H5T_NATIVE_INTEGER, &
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
  subroutine hdf5_file_write_rp_attribute(this, attr_name, attr)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: attr_name
    real(kind=rp), intent(in) :: attr
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

  !> Read an integer attribute
  subroutine hdf5_file_read_int_attribute(this, attr_name, attr, attr_exists)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: attr_name
    integer, intent(inout) :: attr
    logical, intent(inout) :: attr_exists
    integer :: ierr
    integer(hid_t) :: filespace, attr_id
    integer(hsize_t), dimension(1) :: dcount

    ! ====================
    ! Create the attribute
    ! ====================
    dcount = [int(1, hsize_t)]
    call h5aexists_f(this%active_group_id, trim(attr_name), attr_exists, ierr)
    if (attr_exists) then
       ! retrieve the attr id for the existing attribute
       call h5aopen_f(this%active_group_id, trim(attr_name), attr_id, ierr)
    else
       return
    end if

    ! ===========================
    ! Set up writing the data set
    ! ===========================
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, attr, dcount, ierr)

    ! =======================
    ! Clean up
    ! =======================
    call h5aclose_f(attr_id, ierr)

  end subroutine hdf5_file_read_int_attribute

  !> Read a real (kind=rp) attribute
  subroutine hdf5_file_read_rp_attribute(this, attr_name, attr, attr_exists)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: attr_name
    real(kind=rp), intent(inout) :: attr
    logical, intent(inout) :: attr_exists
    integer :: ierr
    integer(hid_t) :: precision_hdf
    integer(hid_t) :: filespace, attr_id
    integer(hsize_t), dimension(1) :: dcount

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
       return
    end if

    ! ===========================
    ! Set up writing the data set
    ! ===========================
    call h5aread_f(attr_id, precision_hdf, attr, dcount, ierr)

    ! =======================
    ! Clean up
    ! =======================
    call h5aclose_f(attr_id, ierr)

  end subroutine hdf5_file_read_rp_attribute

#else

  !> Open a HDF5 file in a given mode
  subroutine hdf5_file_open(this, mode)
    class(hdf5_file_t), intent(inout) :: this
    character(len=1), intent(in) :: mode
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_open

  !> Close the file
  subroutine hdf5_file_close(this)
    class(hdf5_file_t), intent(inout) :: this
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_close

  !> Set the active group for HDF5 files
  subroutine hdf5_file_set_group(this, group_name_path)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in), optional :: group_name_path
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_set_group

  !> Write data in HDF5 format.
  !! @param data Data object to write.
  !! @param[optional] t Simulation time.
  subroutine hdf5_file_write(this, data, t)
    class(hdf5_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_write

  !> Read data in HDF5 format.
  !! @param data Data object to populate.
  subroutine hdf5_file_read(this, data)
    class(hdf5_file_t) :: this
    class(*), target, intent(inout) :: data
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_read

  subroutine hdf5_file_write_dataset(this, data)
    class(hdf5_file_t), intent(inout) :: this
    class(*), intent(inout) :: data
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_write_dataset

  subroutine hdf5_file_read_dataset(this, data_name, data, strategy)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: data_name
    class(*), intent(inout) :: data
    character(len=*), intent(in), optional :: strategy
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_read_dataset

  subroutine hdf5_file_write_attribute(this, data_name, data)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: data_name
    class(*), intent(inout) :: data
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_write_attribute

  subroutine hdf5_file_read_attribute(this, data_name, data, exist)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: data_name
    class(*), intent(inout) :: data
    logical, intent(inout) :: exist
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_read_attribute

  subroutine hdf5_file_write_vector(this, vec)
    class(hdf5_file_t), intent(inout) :: this
    type(vector_t), intent(inout) :: vec
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_write_vector

  subroutine hdf5_file_write_matrix(this, mat)
    class(hdf5_file_t), intent(inout) :: this
    type(matrix_t), intent(inout) :: mat
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_write_matrix

  subroutine hdf5_file_write_field(this, fld)
    class(hdf5_file_t), intent(inout) :: this
    type(field_t), intent(inout) :: fld
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_write_field

  subroutine hdf5_file_read_vector(this, data_name, vec, strategy)
    class(hdf5_file_t) :: this
    character(len=*), intent(in) :: data_name
    type(vector_t), intent(inout) :: vec
    character(len=*), intent(in), optional :: strategy
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_read_vector

  subroutine hdf5_file_read_matrix(this, data_name, mat, strategy)
    class(hdf5_file_t) :: this
    character(len=*), intent(in) :: data_name
    type(matrix_t), intent(inout) :: mat
    character(len=*), intent(in), optional :: strategy
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_read_matrix

  subroutine hdf5_file_write_int_attribute(this, attr_name, attr)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: attr_name
    integer, intent(in) :: attr
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_write_int_attribute

  subroutine hdf5_file_write_rp_attribute(this, attr_name, attr)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: attr_name
    real(kind=rp), intent(in) :: attr
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_write_rp_attribute

  subroutine hdf5_file_read_int_attribute(this, attr_name, attr, attr_exists)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: attr_name
    integer, intent(inout) :: attr
    logical, intent(inout) :: attr_exists
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_read_int_attribute

  subroutine hdf5_file_read_rp_attribute(this, attr_name, attr, attr_exists)
    class(hdf5_file_t), intent(inout) :: this
    character(len=*), intent(in) :: attr_name
    real(kind=rp), intent(inout) :: attr
    logical, intent(inout) :: attr_exists
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_read_rp_attribute

#endif

end module hdf5_file
