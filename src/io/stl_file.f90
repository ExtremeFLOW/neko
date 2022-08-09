! Copyright (c) 2022, The Neko Authors
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
!>  Stereolithography (STL) file 
module stl_file
  use num_types
  use generic_file
  use tri_mesh  
  use logger
  use point
  use mpi_types
  use mpi
  use comm
  use stl    
  implicit none
  private

  !>  Interface for STL files
  type, public, extends(generic_file_t) :: stl_file_t
   contains
     procedure :: read => stl_file_read
     procedure :: write => stl_file_write
  end type stl_file_t

contains

  subroutine stl_file_write(this, data, t)
    class(stl_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    call neko_log%error('Not implemented')
  end subroutine stl_file_write
  
  subroutine stl_file_read(this, data)
    class(stl_file_t) :: this
    class(*), target, intent(inout) :: data
    type(tri_mesh_t), pointer :: tri_msh => null()    
    integer :: status(MPI_STATUS_SIZE)
    integer :: fh
    type(point_t) :: p(3)
    type(stl_hdr_t) :: stl_hdr
    type(stl_triangle_t), allocatable :: stl_tri(:)
    integer :: i, p_idx, ierr

    select type(data)
    type is(tri_mesh_t)
       tri_msh => data
    class default
       call neko_log%error('Invalid data')
    end select

    call MPI_File_open(NEKO_COMM, trim(this%fname), &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_File_read_all(fh, stl_hdr, 1, MPI_STL_HEADER, status, ierr)

    if (stl_hdr%hdr(1:6) .eq. 'solid') then
       call neko_log%error('Invalid STL file (ASCII)')
    end if
    
    call tri_msh%init(stl_hdr%ntri)
    allocate(stl_tri(stl_hdr%ntri))

    call MPI_File_read_all(fh, stl_tri, stl_hdr%ntri, &
         MPI_STL_TRIANGLE, status, ierr)

    p_idx = 0
    do i = 1, stl_hdr%ntri
       p_idx = p_idx + 1
       p(1) = point_t(dble(stl_tri(i)%v1), p_idx)
       p_idx = p_idx + 1
       p(2) = point_t(dble(stl_tri(i)%v2), p_idx)
       p_idx = p_idx + 1
       p(3) = point_t(dble(stl_tri(i)%v3), p_idx)
       call tri_msh%add_element(p(1), p(2), p(3))
    end do

    deallocate(stl_tri)

    call MPI_File_close(fh, ierr)
     
  end subroutine stl_file_read
  
end module stl_file
