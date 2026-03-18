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

!> VTK Module containing utilities for VTK file handling
!!
!! References:
!! - https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
!! - https://scicomp.stackexchange.com/questions/42092/vtk-arbitrary-order-
!!           lagrange-elements-node-positions-ordering-on-reference-tri
!! - https://www.kitware.com/modeling-arbitrary-order-lagrange-finite-elements-
!!           in-the-visualization-toolkit/#:%7E:text=The%20new%20cells%20in%20VTK,may%20vary%20in%20Lagrange%20cells.
module vtk
  use utils, only: linear_index, neko_error
  implicit none
  private

  public :: vtk_ordering

contains

  !> Get the VTK node ordering for a given cell type.
  !! For Lagrange cells, returns an array mapping VTK node position to the
  !! 0-based tensor-product index (i + lx*j + lx*ly*k).
  !! @param cell_type VTK cell type (e.g. 72 for hexahedron)
  !! @param lx Number of points per edge in x-direction (polynomial order + 1)
  !! @param ly Number of points per edge in y-direction (polynomial order + 1)
  !! @param lz Number of points per edge in z-direction (polynomial order + 1)
  !! @return Array of 0-based tensor-product indices in VTK order
  function vtk_ordering(cell_type, lx, ly, lz) result(ordering)
    integer(kind=1), intent(in) :: cell_type
    integer, intent(in), optional :: lx, ly, lz
    integer, allocatable :: ordering(:)

    if (allocated(ordering)) deallocate(ordering)

    select case (cell_type)
    case (70) ! VTK_LAGRANGE_QUADRILATERAL
       if (present(lx) .and. present(ly)) then
          ordering = vtk_lagrange_quad_ordering(lx, ly)
       else
          call neko_error('lx and ly must be provided for arbitrary lagrange ' &
               // 'quadrilateral cells')
       end if
    case (72) ! VTK_LAGRANGE_HEXAHEDRON
       if (present(lx) .and. present(ly) .and. present(lz)) then
          ordering = vtk_lagrange_hex_ordering(lx, ly, lz)
       else
          call neko_error('lx, ly, and lz must be provided for arbitrary ' &
               // 'lagrange hexahedron cells')
       end if
    case default
       call neko_error('Unsupported VTK cell type in vtk_ordering')
    end select

  end function vtk_ordering

  !> Build the VTK Lagrange hexahedron node ordering for a given lx.
  !! Returns an array of size lx*ly*lz mapping VTK node position to the
  !! 0-based tensor-product index (i + lx*j + lx*ly*k).
  !! Implements VTK's PointIndexFromIJK for Lagrange hexahedra.
  !! Node ordering: 8 corners, 4 * (lx - 2 + ly - 2 + lz - 2) edge interiors,
  !! 6*(lx - 2)*(ly - 2) face interiors, (lx - 2)*(ly - 2)*(lz - 2) body
  !! interior.
  !! @param lx Number of points per edge in x-direction (polynomial order + 1)
  !! @param ly Number of points per edge in y-direction (polynomial order + 1)
  !! @param lz Number of points per edge in z-direction (polynomial order + 1)
  !! @return Array of 0-based tensor-product indices in VTK order
  pure function vtk_lagrange_hex_ordering(lx, ly, lz) result(ordering)
    integer, intent(in) :: lx, ly, lz
    integer :: ordering(lx * ly * lz)
    integer :: i, j, k, vtk_idx
    integer :: ibdy, jbdy, kbdy, nbdy, offset
    integer :: n_corners, n_edges, n_faces

    n_corners = 8
    n_edges = 4 * ((lx - 2) + (ly - 2) + (lz - 2))
    n_faces = 2 * ((lx - 2) * (ly - 2) + (lx - 2) * (lz - 2) &
         + (ly - 2) * (lz - 2))

    do concurrent (i = 1:lx, j = 1:ly, k = 1:lz)
       ibdy = merge(1, 0, i .eq. 1 .or. i .eq. lx)
       jbdy = merge(1, 0, j .eq. 1 .or. j .eq. ly)
       kbdy = merge(1, 0, k .eq. 1 .or. k .eq. lz)
       nbdy = ibdy + jbdy + kbdy

       if (nbdy .eq. 3) then
          ! Corner node
          vtk_idx = merge( &
               merge(2, 1, j .ne. 1), &
               merge(3, 0, j .ne. 1), i .ne. 1) &
               + merge(4, 0, k .ne. 1)

       else if (nbdy .eq. 2) then
          ! Edge interior node
          offset = n_corners
          if (ibdy .eq. 0) then
             vtk_idx = (i - 2) &
                  + merge((lx - 2) + (ly - 2), 0, j .ne. 1) &
                  + merge(2 * ((lx - 2) + (ly - 2)), 0, k .ne. 1) &
                  + offset
          else if (jbdy .eq. 0) then
             vtk_idx = (j - 2) &
                  + merge(lx - 2, 2 * (lx - 2) + (ly - 2), i .ne. 1) &
                  + merge(2 * ((lx - 2) + (ly - 2)), 0, k .ne. 1) &
                  + offset
          else
             vtk_idx = (k - 2) + (lz - 2) &
                  * merge(merge(2, 1, j .ne. 1), &
                  merge(3, 0, j .ne. 1), i .ne. 1) &
                  + offset + 4 * ((lx - 2) + (ly - 2))
          end if

       else if (nbdy .eq. 1) then
          ! Face interior node
          offset = n_corners + n_edges
          if (ibdy .eq. 1) then
             vtk_idx = (j - 2) + (ly - 2) * (k - 2) &
                  + merge((ly - 2) * (lz - 2), 0, i .ne. 1) &
                  + offset
          else if (jbdy .eq. 1) then
             offset = offset + 2 * (ly - 2) * (lz - 2)
             vtk_idx = (i - 2) + (lx - 2) * (k - 2) &
                  + merge((lx - 2) * (lz - 2), 0, j .ne. 1) &
                  + offset
          else if (kbdy .eq. 1) then
             offset = offset + 2 * (ly - 2) * (lz - 2) + 2 * (lx - 2) * (lz - 2)
             vtk_idx = (i - 2) + (lx - 2) * (j - 2) &
                  + merge((lx - 2) * (ly - 2), 0, k .ne. 1) &
                  + offset
          end if

       else
          ! Body interior node
          offset = n_corners + n_edges + n_faces
          vtk_idx = offset &
               + (i - 2) + (lx - 2) * ((j - 2) + (ly - 2) * (k - 2))
       end if

       ! ordering(vtk_position + 1) = tensor-product index
       ordering(vtk_idx + 1) = linear_index(i, j, k, 1, lx, ly, lz) - 1
    end do

  end function vtk_lagrange_hex_ordering

  !> Build the VTK Lagrange quadrilateral node ordering for a given lx, ly.
  !! Returns an array of size lx*ly mapping VTK node position to the
  !! 0-based tensor-product index (i + lx*j).
  !! Implements VTK's PointIndexFromIJK for Lagrange quadrilaterals.
  !! Node ordering: 4 corners, 2*(lx-2) + 2*(ly-2) edge interiors,
  !! (lx-2)*(ly-2) face interior.
  !! @param lx Number of points per edge in x-direction (polynomial order + 1)
  !! @param ly Number of points per edge in y-direction (polynomial order + 1)
  !! @return Array of 0-based tensor-product indices in VTK order
  pure function vtk_lagrange_quad_ordering(lx, ly) result(ordering)
    integer, intent(in) :: lx, ly
    integer :: ordering(lx * ly)
    integer :: i, j, vtk_idx
    integer :: ibdy, jbdy, nbdy, offset
    integer :: n_corners, n_edges

    n_corners = 4
    n_edges = 2 * ((lx - 2) + (ly - 2))

    do concurrent (i = 1:lx, j = 1:ly)
       ibdy = merge(1, 0, i .eq. 1 .or. i .eq. lx)
       jbdy = merge(1, 0, j .eq. 1 .or. j .eq. ly)
       nbdy = ibdy + jbdy

       if (nbdy .eq. 2) then
          ! Corner node
          vtk_idx = merge( &
               merge(2, 1, j .ne. 1), &
               merge(3, 0, j .ne. 1), i .ne. 1)

       else if (nbdy .eq. 1) then
          ! Edge interior node
          offset = n_corners
          if (ibdy .eq. 0) then
             vtk_idx = (i - 2) &
                  + merge((lx - 2) + (ly - 2), 0, j .ne. 1) &
                  + offset
          else
             vtk_idx = (j - 2) &
                  + merge(lx - 2, 2 * (lx - 2) + (ly - 2), i .ne. 1) &
                  + offset
          end if

       else
          ! Face interior node
          offset = n_corners + n_edges
          vtk_idx = offset + (i - 2) + (lx - 2) * (j - 2)
       end if

       ordering(vtk_idx + 1) = linear_index(i, j, 1, 1, lx, ly, 1) - 1
    end do

  end function vtk_lagrange_quad_ordering
end module vtk
