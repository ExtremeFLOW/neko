! Copyright (c) 2018-2023, The Neko Authors
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
!> Connectivity cell abstract type
module cell_cnn
  use num_types, only : i4
  use polytope_cnn, only : polytope_cnn_t
  use vertex_cnn, only : vertex_cnn_t, vertex_cnn_ptr, vertex_ncnf_cnn_t
  use edge_cnn, only : edge_cnn_t, edge_cnn_ptr, edge_aligned_cnn_t
  implicit none
  private

  public :: cell_cnn_t

  ! object information
  integer(i4), public, parameter :: NEKO_CELL_DIM = 3

  !> Base type for an abstract three-dimensional polytope (polyhedron)
  !! @details There are multiple possible realisation, so it is just
  !! an abstract type providing common functionality
  type, extends(polytope_cnn_t), abstract :: cell_cnn_t
     !> Ridges are aligned
     type(edge_aligned_cnn_t), dimension(:), allocatable :: ridge
     !> Peak pointers
     type(vertex_ncnf_cnn_t), dimension(:), allocatable :: peak
   contains
     !> Initialise cell dimension
     procedure, pass(this) :: init_dim => cell_init_dim
     !> Get pointers to ridges
     procedure, pass(this) :: rdg => cell_ridge
     !> Get pointers to peaks
     procedure, pass(this) :: pek => cell_peak
     !> Return edges shared by cells
     procedure, pass(this) :: rdg_share => cell_ridge_share
     !> Return vertices shared by cells
     procedure, pass(this) :: pek_share => cell_peak_share
  end type cell_cnn_t

contains

  !> @brief Initialise cell dimension
  subroutine cell_init_dim(this)
    class(cell_cnn_t), intent(inout) :: this

    call this%set_dim(NEKO_CELL_DIM)

    return
  end subroutine cell_init_dim

  !> @brief Return pointers to cell ridges
  !! @parameter[out]  ridge   ridge pointers array
  subroutine cell_ridge(this, ridge)
    class(cell_cnn_t), intent(in) :: this
    type(edge_aligned_cnn_t), dimension(:), allocatable, intent(out) :: ridge
    integer(i4) :: il

    allocate(ridge(this%nridge))
    do il = 1, this%nridge
       ridge(il) = this%ridge(il)
    end do

    return
  end subroutine cell_ridge

  !> @brief Return pointers to cell peaks
  !! @parameter[out]  peak   peak pointers array
  subroutine cell_peak(this, peak)
    class(cell_cnn_t), intent(in) :: this
    type(vertex_cnn_ptr), dimension(:), allocatable, intent(out) :: peak
    integer(i4) :: il

    allocate(peak(this%npeak))
    do il = 1, this%npeak
       peak(il)%obj => this%peak(il)%vertex%obj
    end do

    return
  end subroutine cell_peak

  !> @brief Return positions of ridges shared by cells
  !! @note Cells can be self-periodic
  !! @parameter[in]   other   second cell
  !! @parameter[out]  ishare  number of shared edges
  !! @parameter[out]  ridgep  integer position of shared edges
  subroutine cell_ridge_share(this, other, ishare, ridgep)
    class(cell_cnn_t), intent(in) :: this, other
    integer(i4), intent(out) :: ishare
    integer(i4), dimension(:, :), allocatable, intent(out) :: ridgep
    integer(i4) :: il, jl

    allocate(ridgep(2, this%nridge * other%nridge))
    ishare = 0
    ridgep(:,:) = 0
    do il = 1, this%nridge
       do jl = 1, other%nridge
          if (this%ridge(il)%edge%obj .eq. other%ridge(jl)%edge%obj) then
             ishare = ishare + 1
             ridgep(1,ishare) = il
             ridgep(2,ishare) = jl
          end if
       end do
    end do

    return
  end subroutine cell_ridge_share

  !> @brief Return positions of peaks shared by cells
  !! @note Cells can be self-periodic
  !! @parameter[in]   other   second cell
  !! @parameter[out]  ishare  number of shared vertices
  !! @parameter[out]  peakp   integer position of shared vertices
  pure subroutine cell_peak_share(this, other, ishare, peakp)
    class(cell_cnn_t), intent(in) :: this, other
    integer(i4), intent(out) :: ishare
    integer(i4), dimension(:, :), allocatable, intent(out) :: peakp
    integer(i4) :: il, jl

    allocate(peakp(2, this%npeak * other%npeak))
    ishare = 0
    peakp(:,:) = 0
    do il = 1, this%npeak
       do jl = 1, other%npeak
          if (this%peak(il)%vertex%obj%id() == &
               & other%peak(jl)%vertex%obj%id()) then
             ishare = ishare + 1
             peakp(1,ishare) = il
             peakp(2,ishare) = jl
          end if
       end do
    end do

    return
  end subroutine cell_peak_share

end module cell_cnn
