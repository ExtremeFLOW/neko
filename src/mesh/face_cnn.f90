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
!> Connectivity face abstract type
module face_cnn
  use num_types, only : i4
  use polytope_cnn, only : polytope_cnn_t
  use vertex_cnn, only : vertex_cnn_t, vertex_ncnf_cnn_t, vertex_ncnf_cnn_ptr
  use edge_cnn, only : edge_cnn_t, edge_cnn_ptr, edge_aligned_cnn_t
  implicit none
  private

  public :: face_cnn_t

  ! object information
  integer(i4), public, parameter :: NEKO_FACE_DIM = 2

  !> Base type for an abstract two-dimensional polytope (polygon)
  !! @details There are multiple possible realisation, so it is just
  !! an abstract type providing common functionality
  type, extends(polytope_cnn_t), abstract :: face_cnn_t
     !> Facets are aligned
     type(edge_aligned_cnn_t), dimension(:), allocatable :: facet
     !> Ridges (vertices with hanging information)
     type(vertex_ncnf_cnn_t), dimension(:), allocatable :: ridge
   contains
     !> Initialise face dimension
     procedure, pass(this) :: init_dim => face_init_dim
     !> Free face data
     procedure, pass(this) :: free => face_free
     !> Is face self-periodic
     procedure, pass(this) :: selfp => face_self_periodic
     !> Get pointers to facets
     procedure, pass(this) :: fct => face_facet
     !> Get pointers to ridges
     procedure, pass(this) :: rdg => face_ridge
     !> Return edges shared by faces
     procedure, pass(this) :: fct_share => face_facet_share
     !> Return vertices shared by faces
     procedure, pass(this) :: rdg_share => face_ridge_share
  end type face_cnn_t

contains

  !> @brief Initialise face dimension
  subroutine face_init_dim(this)
    class(face_cnn_t), intent(inout) :: this

    call this%free()

    call this%set_dim(NEKO_FACE_DIM)

    return
  end subroutine face_init_dim

  !> @brief Free face data
  subroutine face_free(this)
    class(face_cnn_t), intent(inout) :: this
    !local variables
    integer(i4) :: il

    call this%set_dim(-1)
    if (allocated(this%facet)) then
       do il = 1, this%nfacet
          call this%facet(il)%free()
       end do
       deallocate(this%facet)
    end if
    if (allocated(this%ridge)) then
       do il = 1, this%nridge
          call this%ridge(il)%free()
       end do
       deallocate(this%ridge)
    end if

    return
  end subroutine face_free

  !> @brief Check if face is self-periodic
  !! @return   selfp
  function face_self_periodic(this) result(selfp)
    class(face_cnn_t), intent(in) :: this
    logical :: selfp
    integer(i4) :: il, jl, itmp

    ! count self periodic edges
    itmp = 0
    do il = 1, this%nfacet - 1
       do jl = il + 1, this%nfacet
          selfp = this%facet(il)%edge%obj.eq.this%facet(jl)%edge%obj
          if (selfp) itmp = itmp + 1
       end do
    end do
    ! count self periodic vertices
    do il = 1, this%nridge - 1
       do jl = il + 1, this%nridge
          selfp = (this%ridge(il)%vertex%obj%id() == &
               & this%ridge(jl)%vertex%obj%id())
          if (selfp) itmp = itmp + 1
       end do
    end do
    if (itmp == 0) then
       selfp = .false.
    else
       selfp = .true.
    end if

    return
  end function face_self_periodic

  !> @brief Return pointers to face facets
  !! @parameter[out]  facet   facet pointers array
  subroutine face_facet(this, facet)
    class(face_cnn_t), intent(in) :: this
    type(edge_aligned_cnn_t), dimension(:), allocatable, intent(out) :: facet
    integer(i4) :: il

    allocate(facet(this%nfacet))
    do il = 1, this%nfacet
       facet(il) = this%facet(il)
    end do

    return
  end subroutine face_facet

  !> @brief Return pointers to face ridges
  !! @parameter[out]  ridge   ridge pointer
  !! @parameter[in]   pos     ridge position
  subroutine face_ridge(this, ridge, pos)
    class(face_cnn_t), target, intent(in) :: this
    type(vertex_ncnf_cnn_ptr), intent(out) :: ridge
    integer(i4), intent(in) :: pos

    if ((pos > 0) .and. (pos <= this%nridge)) then
       ridge%obj => this%ridge(pos)
    else
       ridge%obj => null()
    end if

    return
  end subroutine face_ridge

  !> @brief Return positions of facets shared by faces
  !! @note Faces can be self-periodic
  !! @parameter[in]   other   second face
  !! @parameter[out]  ishare  number of shared edges
  !! @parameter[out]  facetp  integer position of shared edges
  subroutine face_facet_share(this, other, ishare, facetp)
    class(face_cnn_t), intent(in) :: this, other
    integer(i4), intent(out) :: ishare
    integer(i4), dimension(:, :), allocatable, intent(out) :: facetp
    integer(i4) :: il, jl

    allocate(facetp(2, this%nfacet * other%nfacet))
    ishare = 0
    facetp(:,:) = 0
    do il = 1, this%nfacet
       do jl = 1, other%nfacet
          if (this%facet(il)%edge%obj .eq. other%facet(jl)%edge%obj) then
             ishare = ishare + 1
             facetp(1,ishare) = il
             facetp(2,ishare) = jl
          end if
       end do
    end do

    return
  end subroutine face_facet_share

  !> @brief Return positions of ridges shared by faces
  !! @note Faces can be self-periodic
  !! @parameter[in]   other   second face
  !! @parameter[out]  ishare  number of shared vertices
  !! @parameter[out]  ridgep  integer position of shared vertices
  pure subroutine face_ridge_share(this, other, ishare, ridgep)
    class(face_cnn_t), intent(in) :: this, other
    integer(i4), intent(out) :: ishare
    integer(i4), dimension(:, :), allocatable, intent(out) :: ridgep
    integer(i4) :: il, jl

    allocate(ridgep(2, this%nridge * other%nridge))
    ishare = 0
    ridgep(:,:) = 0
    do il = 1, this%nridge
       do jl = 1, other%nridge
          if (this%ridge(il)%vertex%obj%id() == &
               & other%ridge(jl)%vertex%obj%id()) then
             ishare = ishare + 1
             ridgep(1,ishare) = il
             ridgep(2,ishare) = jl
          end if
       end do
    end do

    return
  end subroutine face_ridge_share

end module face_cnn
