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
!> Connectivity edge type
module edge
  use num_types, only : i4
  use utils, only : neko_error
  use polytope, only : polytope_t
  use vertex, only : vertex_t, vertex_ptr
  use alignment_edge, only : alignment_edge_t, alignment_edge_op_set_t
  implicit none
  private

  public :: edge_t, edge_ptr, edge_aligned_t, edge_aligned_ptr

  ! object information
  integer(i4), parameter :: NEKO_EDGE_DIM = 1
  integer(i4), parameter :: NEKO_EDGE_NFACET = 2
  integer(i4), parameter :: NEKO_EDGE_NRIDGE = 0
  integer(i4), parameter :: NEKO_EDGE_NPEAK = 0

  !> Edge type for global communication
  !! @details Edge as the only realisation of one-dimensional polytope (dion)
  !! and contains unique global id and two vertices. Edge is an object with
  !! alignment.
  !! @verbatim
  !! Node numbering
  !!      1+-----+2    +----> r
  !! @endverbatim
  type, extends(polytope_t) :: edge_t
     !> Facet pointers
     type(vertex_ptr), dimension(:), allocatable :: facet
   contains
     !> Initialise edge
     procedure, pass(this) :: init => edge_init
     !> Free edge data
     procedure, pass(this) :: free => edge_free
     !> Is edged self-periodic
     procedure, pass(this) :: selfp => edge_self_periodic
     !> Get pointers to facets
     procedure, pass(this) :: fct => edge_facet
     !> Return vertices shared by edges
     procedure, pass(this) :: fct_share => edge_facet_share
     !> Get alignment info
     procedure, pass(this) :: eq_algn => edge_equal_align
     !> Edge equality including vertex information
     procedure, pass(this) :: equal => edge_equal
     generic :: operator(.eq.) => equal
  end type edge_t

  !> Pointer to an edge type
  type ::  edge_ptr
     type(edge_t), pointer :: obj
  end type edge_ptr

  !> Edge with alignment information
  type :: edge_aligned_t
     !> edge pointer
     type(edge_ptr) :: edge
     !> alignment operator
     type(alignment_edge_op_set_t) :: algn_op
   contains
     !> Initialise aligned edge
     procedure, pass(this) :: init => edge_aligned_init
     !> Free aligned edge
     procedure, pass(this) :: free => edge_aligned_free
     !> Return edge pointer
     procedure, pass(this) :: edgep => edge_aligned_edgep
     !> Return edge relative alignment
     procedure, pass(this) :: algn => edge_aligned_alignment_get
     !> Test alignment
     procedure, pass(this) :: test => edge_aligned_test
  end type edge_aligned_t

  !> Pointer to an aligned edge type
  type ::  edge_aligned_ptr
     type(edge_aligned_t), pointer :: obj
  end type edge_aligned_ptr

contains

  !> @brief Initialise edge with global id and vertices
  !! @details Vertex order defines edge orientation
  !! @parameter[in]   id     unique id
  !! @parameter[in]   vrt1, vrt2  bounding vertices
  subroutine edge_init(this, id, vrt1, vrt2)
    class(edge_t), intent(inout) :: this
    integer(i4), intent(in) :: id
    type(vertex_t), intent(in), target :: vrt1, vrt2

    call this%free()

    call this%set_dim(NEKO_EDGE_DIM)
    call this%set_nelem(NEKO_EDGE_NFACET, NEKO_EDGE_NRIDGE,&
         & NEKO_EDGE_NPEAK)
    call this%set_id(id)
    allocate(this%facet(NEKO_EDGE_NFACET))
    this%facet(1)%obj => vrt1
    this%facet(2)%obj => vrt2

    return
  end subroutine edge_init

  !> @brief Free edge data
  subroutine edge_free(this)
    class(edge_t), intent(inout) :: this
    !local variables
    integer(i4) :: il

    call this%set_dim(-1)
    if (allocated(this%facet)) then
       do il = 1, this%nfacet
          this%facet(il)%obj => null()
       end do
       deallocate(this%facet)
    end if

    return
  end subroutine edge_free

  !> @brief Check if edge is self-periodic
  !! @return   selfp
  pure function edge_self_periodic(this) result(selfp)
    class(edge_t), intent(in) :: this
    logical :: selfp

    selfp = (this%facet(1)%obj%id() == this%facet(2)%obj%id())

    return
  end function edge_self_periodic

  !> @brief Return pointers to edge facets
  !! @parameter[out]  facet   facet pointers array
  subroutine edge_facet(this, facet)
    class(edge_t), intent(in) :: this
    type(vertex_ptr), dimension(:), allocatable, intent(out) :: facet
    integer(i4) :: il

    allocate(facet(this%nfacet))
    do il = 1, this%nfacet
       facet(il)%obj => this%facet(il)%obj
    end do

    return
  end subroutine edge_facet

  !> @brief Return pointers to facets shared by edges
  !! @note Edges can be self-periodic
  !! @parameter[in]   other   second edge
  !! @parameter[out]  ishare  number of shared vertices
  !! @parameter[out]  facetp  integer position of shared vertices
  pure subroutine edge_facet_share(this, other, ishare, facetp)
    class(edge_t), intent(in) :: this, other
    integer(i4), intent(out) :: ishare
    integer(i4), dimension(:, :), allocatable, intent(out) :: facetp
    integer(i4) :: il, jl

    allocate(facetp(2, this%nfacet * other%nfacet))
    ishare = 0
    facetp(:,:) = 0
    do il = 1, this%nfacet
       do jl = 1, other%nfacet
          if (this%facet(il)%obj%id() == other%facet(jl)%obj%id()) then
             ishare = ishare + 1
             facetp(1,ishare) = il
             facetp(2,ishare) = jl
          end if
       end do
    end do

    return
  end subroutine edge_facet_share

  !> @brief Check if two edges are the same
  !! @note Alignment for self-periodic edges will not be correct
  !! @parameter[in]  other    second edge
  !! @return   equal
  subroutine edge_equal_align(this, other, equal, algn)
    class(edge_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn
    type(alignment_edge_t) :: algn_op
    integer(i4), dimension(2) :: trans = [1, 2]

    algn = -1
    ! check polygon information
    equal = this%equal_poly(other)
    if (equal) then
       ! check global id
       equal = (this%id() == other%id())
       ! check vertices getting possible orientation (doesn't work for
       ! self-periodic)
       if (equal) then
          call algn_op%init()
          select type(other)
          type is (edge_t)
             ! check all the alignment options
             do algn = 0, algn_op%nop()
                call algn_op%trns_f_i4(algn)%obj(2, trans)
                equal = (this%facet(1)%obj.eq.other%facet(trans(1))%obj).and.&
                     &(this%facet(2)%obj.eq.other%facet(trans(2))%obj)
                if (equal) return
                call algn_op%trns_inv_f_i4(algn)%obj(2, trans)
             end do
          class default
             equal = .false.
          end select
          if (.not.equal) then
             ! Something wrong; edge with the same global id should have
             ! the same type and the same facets
             call neko_error('Mismatch in edge and vertex global id')
          end if
       end if
    end if

    return
  end subroutine edge_equal_align

  !> @brief Check if two edges are the same
  !! @note No special treatment of self-periodic edges
  !! @parameter[in]  other    second edge
  !! @return   equal
  function edge_equal(this, other) result(equal)
    class(edge_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: equal
    integer(i4) :: algn

    call this%eq_algn(other, equal, algn)

    return
  end function edge_equal

  !> @brief Initialise edge with alignment information
  !! @parameter[in]   edge   edge
  !! @parameter[in]   algn   alignment
  subroutine edge_aligned_init(this, edge, algn)
    class(edge_aligned_t), intent(inout) :: this
    type(edge_t), intent(in), target :: edge
    integer(i4), intent(in) :: algn

    call this%free()
    ! set global edge pointer
    this%edge%obj => edge
    ! set relative alignment transformation
    call this%algn_op%init(algn)

    return
  end subroutine edge_aligned_init

  !> @brief Free edge with alignment information
  subroutine edge_aligned_free(this)
    class(edge_aligned_t), intent(inout) :: this

    this%edge%obj => null()
    call this%algn_op%free()

    return
  end subroutine edge_aligned_free

  !> @brief Return pointers to the edge
  !! @parameter[out]  edge   edge pointer
  subroutine edge_aligned_edgep(this, edge)
    class(edge_aligned_t), intent(in) :: this
    type(edge_ptr), intent(out) :: edge
    edge%obj => this%edge%obj
    return
  end subroutine edge_aligned_edgep

  !> @brief Get relative edge alignment
  !! @return   alignment
  pure function edge_aligned_alignment_get(this) result(alignment)
    class(edge_aligned_t), intent(in) :: this
    integer(i4) :: alignment
    alignment = this%algn_op%alignment
  end function edge_aligned_alignment_get

  !> @brief Check if two edges are properly aligned
  !! @note No special treatment of self-periodic edges, so this is not checked
  !! @parameter[in]  other    second edge
  !! @return   aligned
  function edge_aligned_test(this, other) result(aligned)
    class(edge_aligned_t), intent(in) :: this
    class(edge_t), intent(in) :: other
    logical :: aligned
    integer(i4), dimension(2) :: vrt, vrto

    ! only equal edges can be checked
    if (this%edge%obj.eq.other) then
       vrt(1) = this%edge%obj%facet(1)%obj%id()
       vrt(2) = this%edge%obj%facet(2)%obj%id()
       vrto(1) = other%facet(1)%obj%id()
       vrto(2) = other%facet(2)%obj%id()
       call this%algn_op%trns_f_i4%obj( 2, vrto)
       aligned = (vrt(1) == vrto(1)).and.(vrt(2) == vrto(2))
    else
       call neko_error('Edges not equal')
    end if

    return
  end function edge_aligned_test

end module edge
