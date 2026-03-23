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
  use field, only : field_t
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
     integer :: precision = 0
     integer :: offset
     integer :: count

   contains
     procedure :: open => hdf5_file_2_open
     procedure :: set_active_group => hdf5_file_2_set_group
     procedure :: close => hdf5_file_2_close
     procedure :: read => hdf5_file_2_read
     procedure :: write => hdf5_file_2_write
     procedure :: set_overwrite => hdf5_file_2_set_overwrite 
     procedure :: get_fname => file_get_fname
     procedure :: set_precision => hdf5_file_2_set_precision
     procedure, pass(this) :: write_vector => hdf5_file_2_write_vector
     procedure, pass(this) :: write_matrix => hdf5_file_2_write_matrix
     procedure, pass(this) :: write_field => hdf5_file_2_write_field
     !generic :: write_dataset => write_vector, write_matrix, write_field
     procedure :: write_dataset => hdf5_file_2_write_dataset
  end type hdf5_file_2_t

contains

#ifdef HAVE_HDF5



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
