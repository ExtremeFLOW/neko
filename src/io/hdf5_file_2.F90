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
module hdf5_file_2
  use num_types, only : rp, dp, sp
  use generic_file, only : generic_file_t
  use utils, only : neko_error, neko_warning, filename_suffix_pos, filename_split
  use logger, only : neko_log
  use vector, only : vector_t
  use matrix, only : matrix_t
  use device, only : DEVICE_TO_HOST
  use comm, only : pe_rank, NEKO_COMM
  use mpi_f08, only : MPI_INFO_NULL, MPI_Allreduce, MPI_Allgather, &
       MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, MPI_MAX, MPI_Comm_size, MPI_Exscan, &
       MPI_Barrier
#ifdef HAVE_HDF5
  use hdf5
#endif
  implicit none
  private

  !> Interface for HDF5 files
  type, public, extends(generic_file_t) :: hdf5_file_2_t
#ifdef HAVE_HDF5
    integer(hid_t) :: file_id = -1_hid_t
    integer(hid_t) :: active_group_id = -1_hid_t
    integer(hid_t) :: plist_id

#endif

    character(len=1) :: mode
    integer :: precision
    integer :: offset
    integer :: count

   contains
     procedure :: open => hdf5_file_2_open
     procedure :: close => hdf5_file_2_close
     procedure :: read => hdf5_file_2_read
     procedure :: write => hdf5_file_2_write
     procedure :: set_overwrite => hdf5_file_2_set_overwrite
     procedure :: set_active_group => hdf5_file_2_set_group 
     procedure :: get_fname => file_get_fname
     procedure :: set_precision => hdf5_file_2_set_precision
     procedure, pass(this) :: write_vector => hdf5_file_2_write_vector
     procedure, pass(this) :: write_matrix => hdf5_file_2_write_matrix
     generic :: write_dataset => write_vector, write_matrix
  end type hdf5_file_2_t

contains

  !> Return the file name with the start counter.
  function file_get_fname(this) result(base_fname)
    class(hdf5_file_2_t), intent(in) :: this
    character(len=1024) :: base_fname
    character(len=1024) :: fname
    character(len=1024) :: path, name, suffix

    fname = trim(this%get_base_fname())
    call filename_split(fname, path, name, suffix)

    write(base_fname, '(A,A,"_",I0,A)') &
         trim(path), trim(name), this%get_start_counter(), trim(suffix)

  end function file_get_fname
  
  !> Set the precision for the output (single or double)
  subroutine hdf5_file_2_set_precision(this, precision)
    class(hdf5_file_2_t), intent(inout) :: this
    integer, intent(in) :: precision
    this%precision = precision
  end subroutine hdf5_file_2_set_precision
  
  !> Set the overwrite flag for HDF5 files
  subroutine hdf5_file_2_set_overwrite(this, overwrite)
    class(hdf5_file_2_t), intent(inout) :: this
    logical, intent(in) :: overwrite
    this%overwrite = overwrite
  end subroutine hdf5_file_2_set_overwrite


#ifdef HAVE_HDF5

  !> Open a HDF5 file in a mode
  subroutine hdf5_file_2_open(this, mode)
    class(hdf5_file_2_t), intent(inout) :: this
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
    else if (this%precision .eq. 0) then
       this%precision = rp
    end if

    ! File counter management
    call this%increment_counter()
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

    write(*,*) "Opened HDF5 file: ", trim(fname), " with counter: ", counter

  end subroutine hdf5_file_2_open


  !> Close the file
  subroutine hdf5_file_2_close(this)
    class(hdf5_file_2_t), intent(inout) :: this
    integer :: ierr

    if (this%active_group_id /= -1_hid_t) then
       call h5gclose_f(this%active_group_id, ierr)
       this%active_group_id = -1_hid_t
    end if

    call h5pclose_f(this%plist_id, ierr)
    call h5fclose_f(this%file_id, ierr)
    call h5close_f(ierr)

    write(*,*) "Closed HDF5 file: ", trim(this%get_fname())

  end subroutine hdf5_file_2_close

subroutine hdf5_file_2_write_vector(this, vec)
  class(hdf5_file_2_t), intent(inout) :: this
  type(vector_t), intent(inout) :: vec 
  integer :: ierr, counts, offset, total_count, dset_rank
  integer(hid_t) :: precision_hdf
  integer(hid_t) :: xf_id, filespace, dset_id, memspace
  integer(hsize_t), dimension(1) :: dcount, doffset
  integer(hsize_t), dimension(1) :: ddims
  logical :: dset_exists
  
  ! ===============
  ! Get vector info
  ! ===============
  counts = vec%size()
  offset = 0
  total_count = 0
  call MPI_Exscan(counts, offset, 1, MPI_INTEGER, &
       MPI_SUM, NEKO_COMM, ierr)
  call MPI_Allreduce(counts, total_count, 1, MPI_INTEGER, &
       MPI_SUM, NEKO_COMM, ierr)

  ! Sync the data
  call vec%copy_from(DEVICE_TO_HOST, .true.)
  
  ! ===============
  ! Configure MPIIO
  ! ===============
  call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
  call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)
  precision_hdf = hdf5_file_2_determine_real(this%precision)

  ! =================== 
  ! Create the data set
  ! ===================
  call h5lexists_f(this%active_group_id, trim(vec%name), dset_exists, ierr)
  if (dset_exists) then
     ! retireve the dset id for the existing data set
     !call h5dopen_f(this%active_group_id, trim(vec%name), dset_id, ierr)
     call neko_error("dataset already exist in the file")
  else
    dset_rank = 1 ! rank 1 array, i.e. a vector
    ddims = [int(total_count, hsize_t)] ! global size of the vector
    ! create file space of this shape
    call h5screate_simple_f(dset_rank, ddims, filespace, ierr)
    ! create the data set with the given shape  
    call h5dcreate_f(this%active_group_id, trim(vec%name), precision_hdf, &
          filespace, dset_id, ierr)
    call h5sclose_f(filespace, ierr)
  end if

  ! =======================
  ! Write the data set
  ! =======================
  dcount = [int(counts, hsize_t)] ! local size of the vector
  doffset = [int(offset, hsize_t)] ! offset for this rank in the global vector
  ! Get the total file space (shape) of the data set
  call h5dget_space_f(dset_id, filespace, ierr)
  ! Get only the slice where my rank writes
  call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, ierr)
  ! Create the corresponding memory space (buffer) for my local data
  call h5screate_simple_f(dset_rank, dcount, memspace, ierr)
  ! Write the data
  call h5dwrite_f(dset_id, precision_hdf, vec%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = xf_id)

  ! =======================
  ! Clean up
  ! =======================
  call h5pclose_f(xf_id, ierr)
  call h5sclose_f(memspace, ierr)
  call h5sclose_f(filespace, ierr)
  call h5dclose_f(dset_id, ierr)

end subroutine hdf5_file_2_write_vector

subroutine hdf5_file_2_write_matrix(this, mat)
  class(hdf5_file_2_t), intent(inout) :: this
  type(matrix_t), intent(inout) :: mat
  integer :: ierr, counts, offset, total_count, dset_rank, strides
  integer(hid_t) :: precision_hdf
  integer(hid_t) :: xf_id, filespace, dset_id, memspace
  integer(hsize_t), dimension(2) :: dcount, doffset
  integer(hsize_t), dimension(2) :: ddims
  logical :: dset_exists
 
  ! ===============
  ! Get Matrix info
  ! ===============
  strides = mat%get_nrows()
  counts = mat%get_ncols()
  total_count = 0
  offset = 0
  call MPI_Exscan(counts, offset, 1, MPI_INTEGER, &
       MPI_SUM, NEKO_COMM, ierr)
  call MPI_Allreduce(counts, total_count, 1, MPI_INTEGER, &
       MPI_SUM, NEKO_COMM, ierr)

  ! Sync the data
  call mat%copy_from(DEVICE_TO_HOST, .true.)
  
  ! ===============
  ! Configure MPIIO
  ! ===============
  call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
  call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)
  precision_hdf = hdf5_file_2_determine_real(this%precision)

  ! =================== 
  ! Create the data set
  ! ===================
  call h5lexists_f(this%active_group_id, trim(mat%name), dset_exists, ierr)
  if (dset_exists) then
     !! retireve the dset id for the existing data set
     !call h5dopen_f(this%active_group_id, trim(mat%name), dset_id, ierr)
     call neko_error("dataset already exist in the file")
  else
    dset_rank = 2 ! rank 2 array, i.e. a matrix
    ddims = [int(strides, hsize_t), int(total_count, hsize_t)] ! global size of the matrix
    ! create file space of this shape
    call h5screate_simple_f(dset_rank, ddims, filespace, ierr)
    ! create the data set with the given shape  
    call h5dcreate_f(this%active_group_id, trim(mat%name), precision_hdf, &
          filespace, dset_id, ierr)
    call h5sclose_f(filespace, ierr)
  end if

  ! =======================
  ! Write the data set
  ! =======================
  dcount = [int(strides, hsize_t), int(counts, hsize_t)] ! local size of the matrix
  doffset = [0_hsize_t, int(offset, hsize_t)] ! offset for this rank in the global matrix
  ! Get the total file space (shape) of the data set
  call h5dget_space_f(dset_id, filespace, ierr)
  ! Get only the slice where my rank writes
  call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, ierr)
  ! Create the corresponding memory space (buffer) for my local data
  call h5screate_simple_f(dset_rank, dcount, memspace, ierr)
  ! Write the data
  call h5dwrite_f(dset_id, precision_hdf, mat%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = xf_id)

  ! =======================
  ! Clean up
  ! =======================
  call h5pclose_f(xf_id, ierr)
  call h5sclose_f(memspace, ierr)
  call h5sclose_f(filespace, ierr)
  call h5dclose_f(dset_id, ierr)

end subroutine hdf5_file_2_write_matrix


!> Write data in HDF5 format
subroutine hdf5_file_2_write(this, data, t)
  class(hdf5_file_2_t), intent(inout) :: this
  class(*), target, intent(in) :: data
  real(kind=rp), intent(in), optional :: t
  
  select type (data)
  type is (vector_t)
    call neko_error("Nothing implemented here yet")
  class default
    call neko_error("Unsupported data type for HDF5 output")
  end select

end subroutine hdf5_file_2_write

  !> Read data in HDF5 format
  subroutine hdf5_file_2_read(this, data)
    class(hdf5_file_2_t) :: this
    class(*), target, intent(inout) :: data
  
    call neko_error('file reading is not yet implemented in this type')

  end subroutine hdf5_file_2_read

  !> Set the active group for HDF5 files
  !! @param this The HDF5 file object
  !! @param An array of strings that show the path to the group to create or open.
  subroutine hdf5_file_2_set_group(this, group_name)
   class(hdf5_file_2_t), intent(inout) :: this
   character(len=*), intent(in) :: group_name(:)

   integer(hid_t) :: current_id, group_id
   integer :: ierr, i, num_groups
   logical :: group_exists

   num_groups = size(group_name)

   ! Close previous active group if one is open
   if (this%active_group_id /= -1_hid_t) then
      call h5gclose_f(this%active_group_id, ierr)
      this%active_group_id = -1_hid_t
   end if

   ! Start from root location = file
   current_id = this%file_id

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
  end subroutine hdf5_file_2_set_group

  !> Determine hdf5 real type corresponding to NEKO_REAL
  !! @note This must be called after h5open_f, otherwise
  !! the H5T_NATIVE_XYZ types has a value of 0
  function hdf5_file_2_determine_real(precision) result(H5T_NEKO_REAL)
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
  end function hdf5_file_2_determine_real

#else

  !> Write data in HDF5 format
  subroutine hdf5_file_2_write(this, data, t)
    class(hdf5_file_2_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_2_write

  !> Read data in HDF5 format
  subroutine hdf5_file_2_read(this, data)
    class(hdf5_file_2_t) :: this
    class(*), target, intent(inout) :: data
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_2_read

#endif

end module hdf5_file_2
