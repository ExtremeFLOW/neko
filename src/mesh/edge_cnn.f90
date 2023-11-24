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
!> Connectivity edge types
module edge_cnn
  use num_types, only : i4
  use utils, only : neko_error
  use polytope_cnn, only : polytope_cnn_t
  use vertex_cnn, only : vertex_cab_t, vertex_cab_ptr
  use alignment_edge, only : alignment_edge_t, alignment_edge_op_set_t
  use ncnf_interpolation_edge, only : ncnf_intp_edge_op_set_t
  implicit none
  private

  public :: edge_cab_t, edge_cab_ptr, edge_aligned_cab_t, edge_3d_ncnf_cac_t,&
       & edge_3d_ncnf_cac_ptr, edge_2d_ncnf_cac_t, edge_2d_ncnf_cac_ptr

  ! object information
  integer(i4), public, parameter :: NEKO_EDGE_DIM = 1
  integer(i4), public, parameter :: NEKO_EDGE_NFACET = 2
  integer(i4), public, parameter :: NEKO_EDGE_NRIDGE = 0
  integer(i4), public, parameter :: NEKO_EDGE_NPEAK = 0

  !> Type for abstract edge object
  !! @details Edge is the only realisation of one-dimensional polytope (dion)
  !! and contains unique global id and two abstract vertices. Edge is an object
  !! with alignment.
  !! @verbatim
  !! Node numbering
  !!      f_1-----f_2    +----> r
  !! @endverbatim
  !! Its only actualisation are components of higher-dimension objects.
  type, extends(polytope_cnn_t) :: edge_cab_t
     !> Facets (abstract vertices)
     type(vertex_cab_ptr), dimension(:), allocatable :: facet
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
     !> Check equality and get relative alignment info
     procedure, pass(this) :: eq_algn => edge_equal_align
     !> Edge equality including vertex information
     procedure, pass(this) :: equal => edge_equal
     generic :: operator(.eq.) => equal
  end type edge_cab_t

  !> Pointer to an abstract edge object
  type ::  edge_cab_ptr
     type(edge_cab_t), pointer :: ptr
  end type edge_cab_ptr

  !> Abstract edge object with alignment information
  type :: edge_aligned_cab_t
     !> edge pointer
     type(edge_cab_ptr) :: edge
     !> alignment operator
     type(alignment_edge_op_set_t) :: algn_op
   contains
     !> Initialise aligned edge
     procedure, pass(this) :: init_algn => edge_aligned_init
     !> Free aligned edge
     procedure, pass(this) :: free_algn => edge_aligned_free
     !> Return edge pointer
     procedure, pass(this) :: edgep => edge_aligned_edgep
     !> Return edge relative alignment
     procedure, pass(this) :: algn => edge_aligned_alignment_get
     !> Test alignment
     procedure, pass(this) :: test => edge_aligned_test
  end type edge_aligned_cab_t

  !> Actualisation of the abstract edge for 3D nonconforming meshes
  !! @details Edge can be either independent (part of conforming interface or
  !! parent; marked 0) or hanging. The meaning of hanging depends on the mesh
  !! dimension. In 3D case hanging edges can be either ridge hanging (edge is
  !! a component of the full ridge, but neither of the two facets touching the
  !! ridge is hanging), or facet hanging. Ridge hanging edges can be either
  !! the first (marked 1) or the second (marked 2) half of the full edge.
  !! The facet hanging edges can be either located in the middle of the facet
  !! (marked 5) or at the facet border. In the last case they can correspond
  !! to the the first (marked 3) or the second (marked 4) half of the full
  !! edge. See the diagram below for clarification.
  !! @verbatim
  !! Example of edge marking for hex meshes.
  !! Hanging edge marking for three-dimensional nonconforming interface
  !!
  !!                     o                  +-------+
  !!                     :                  |\       \
  !!                     :                  2 \       \
  !!                     :                  |  +-------+
  !!                     +-------+          +  |       |
  !!                     |\       \         :\ |       |
  !!                     1 \       \        : \|       |
  !!                     |  +-------+       :  +-------+
  !!                     +  |       |       o
  !!                      \ |       |
  !!                       \|       |
  !!                        +-------+
  !!
  !!  o...............o o...............o +---3---+.......o o.......+---4---+
  !!  :               : :               : |       |       : :       |       |
  !!  :               : :               : 4       5       : :       5       4
  !!  :               : :               : |       |       : :       |       |
  !!  +---5---+       : :       +---5---+ +---5---+       : :       +---5---+
  !!  |       |       : :       |       | :               : :               :
  !!  3       5       : :       5       3 :               : :               :
  !!  |       |       : :       |       | :               : :               :
  !!  +---3---+.......o o.......+---4---+ o...............o o...............o
  !! @endverbatim
  !! There are two different 1D interpolation operators corresponding the first
  !! and the second half of the full edge (marked 1 and 2). The rest of the
  !! options (marked 3, 4, and 5) result from the 2D face interpolations.
  !! Edge is always a component part of the higher-dimension object, so
  !! position gives it's location in the object. There is no boundary
  !! information in this case.
  type, extends(edge_aligned_cab_t) :: edge_3d_ncnf_cac_t
     !> interpolation operator
     type(ncnf_intp_edge_op_set_t) :: intp_op
     !> position in the object
     integer(i4) :: position = -1
   contains
     !> Initialise aligned edge pointer and position
     procedure, pass(this) :: init_3d => edge_3d_hanging_init
     !> Free edge data
     procedure, pass(this) :: free => edge_3d_hanging_free
     !> Set hanging information
     procedure, pass(this) :: set_hng => edge_3d_hanging_set
     !> Get hanging information
     procedure, pass(this) :: hng => edge_3d_hanging_get
     !> Set position information
     procedure, pass(this) :: set_pos => edge_3d_position_set
     !> Get position information
     procedure, pass(this) :: pos => edge_3d_position_get
  end type edge_3d_ncnf_cac_t

  !> Pointer to a nonconforming edge actualisation
  type ::  edge_3d_ncnf_cac_ptr
     type(edge_3d_ncnf_cac_t), pointer :: ptr
  end type edge_3d_ncnf_cac_ptr

  !> Actualisation of the abstract edge for 2D nonconforming meshes
  !! @details Edge can be either independent (part of conforming interface or
  !! parent; marked 0) or hanging. The meaning of hanging depends on the mesh
  !! dimension. In 2D case edges are face facets and can be either the first
  !! (marked 1) or the second (marked 2) half of a full edge. See the diagram
  !! below for clarification.
  !! @verbatim
  !! Example of edge marking for quad meshes.
  !! Hanging edge marking for two-dimensional nonconforming interface
  !!  o...............o +...............o
  !!  :               : |               :
  !!  :               : 2               :
  !!  :               : |               :
  !!  +               : +               :
  !!  |               : :               :
  !!  1               : :               :
  !!  |               : :               :
  !!  +...............o o...............o
  !! @endverbatim
  !! There are two different 1D interpolation operators corresponding the first
  !! and the second half of the full edge (marked 1 and 2).
  !! For 3D meshes edge is a facet of a face, so @a position gives it's location
  !! in the cell. In addition @a boundary gives an information regarding
  !! internal/external boundary condition (0 -internal, 1 - periodic, ....).
  type, extends(edge_3d_ncnf_cac_t) :: edge_2d_ncnf_cac_t
     !> Internal/external boundary condition flag
     integer(i4) :: boundary = -1
   contains
     !> Initialise 2D edge actualisation
     procedure, pass(this) :: init_2d => edge_2d_hanging_init
     !> Free edge data
     procedure, pass(this) :: free => edge_2d_hanging_free
     !> Set boundary information
     procedure, pass(this) :: set_bnd => edge_2d_boundary_set
     !> Get boundary information
     procedure, pass(this) :: bnd => edge_2d_boundary_get
  end type edge_2d_ncnf_cac_t

  !> Pointer to a nonconforming 2D edge actualisation
  type ::  edge_2d_ncnf_cac_ptr
     type(edge_2d_ncnf_cac_t), pointer :: ptr
  end type edge_2d_ncnf_cac_ptr

contains

  !> @brief Initialise edge with global id and vertices
  !! @details Vertex order defines edge orientation
  !! @parameter[in]   id          unique id
  !! @parameter[in]   vrt1, vrt2  bounding vertices
  subroutine edge_init(this, id, vrt1, vrt2)
    class(edge_cab_t), intent(inout) :: this
    integer(i4), intent(in) :: id
    type(vertex_cab_t), intent(in), target :: vrt1, vrt2

    call this%free()

    call this%set_dim(NEKO_EDGE_DIM)
    call this%set_nelem(NEKO_EDGE_NFACET, NEKO_EDGE_NRIDGE,&
         & NEKO_EDGE_NPEAK)
    call this%set_id(id)
    ! get facet pointers
    allocate(this%facet(NEKO_EDGE_NFACET))
    this%facet(1)%ptr => vrt1
    this%facet(2)%ptr => vrt2
  end subroutine edge_init

  !> @brief Free edge data
  subroutine edge_free(this)
    class(edge_cab_t), intent(inout) :: this
    !local variables
    integer(i4) :: il

    call this%set_dim(-1)
    if (allocated(this%facet)) then
       do il = 1, this%nfacet
          this%facet(il)%ptr => null()
       end do
       deallocate(this%facet)
    end if
  end subroutine edge_free

  !> @brief Check if edge is self-periodic
  !! @return   selfp
  pure function edge_self_periodic(this) result(selfp)
    class(edge_cab_t), intent(in) :: this
    logical :: selfp

    selfp = (this%facet(1)%ptr%id() == this%facet(2)%ptr%id())
  end function edge_self_periodic

  !> @brief Return edge facet at given position
  !! @parameter[out]  facet   facet pointer
  !! @parameter[in]   pos     facet position
  subroutine edge_facet(this, facet, pos)
    class(edge_cab_t), target, intent(in) :: this
    type(vertex_cab_ptr), intent(out) :: facet
    integer(i4), intent(in) :: pos

    if ((pos > 0) .and. (pos <= NEKO_EDGE_NFACET)) then
       facet%ptr => this%facet(pos)%ptr
    else
       facet%ptr => null()
    end if
  end subroutine edge_facet

  !> @brief Return positions of facets shared by edges
  !! @note Edges can be self-periodic
  !! @parameter[in]   other   second edge
  !! @parameter[out]  ishare  number of shared vertices
  !! @parameter[out]  facetp  integer position of shared vertices
  pure subroutine edge_facet_share(this, other, ishare, facetp)
    class(edge_cab_t), intent(in) :: this, other
    integer(i4), intent(out) :: ishare
    integer(i4), dimension(:, :), allocatable, intent(out) :: facetp
    integer(i4) :: il, jl

    allocate(facetp(2, this%nfacet * other%nfacet))
    ishare = 0
    facetp(:, :) = 0
    do il = 1, this%nfacet
       do jl = 1, other%nfacet
          if (this%facet(il)%ptr%id() == other%facet(jl)%ptr%id()) then
             ishare = ishare + 1
             facetp(1, ishare) = il
             facetp(2, ishare) = jl
          end if
       end do
    end do
  end subroutine edge_facet_share

  !> @brief Check if two edges are the same
  !! @note Alignment for self-periodic edges will not be correct
  !! @parameter[in]  other    second edge
  !! @return   equal
  subroutine edge_equal_align(this, other, equal, algn)
    class(edge_cab_t), intent(in) :: this
    class(polytope_cnn_t), intent(in) :: other
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn
    type(alignment_edge_t) :: algn_op
    integer(i4), dimension(NEKO_EDGE_NFACET) :: trans

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
          type is (edge_cab_t)
             ! check all the alignment options
             trans(1) = 1
             trans(2) = 2
             do algn = 0, algn_op%nop()
                call algn_op%trns_inv_f_i4(algn)%ptr(NEKO_EDGE_NFACET, trans)
                equal = (this%facet(1)%ptr .eq. &
                     & other%facet(trans(1))%ptr) .and. &
                     & (this%facet(2)%ptr .eq. other%facet(trans(2))%ptr)
                if (equal) return
                call algn_op%trns_f_i4(algn)%ptr(NEKO_EDGE_NFACET, trans)
             end do
          class default
             equal = .false.
          end select
          if (.not. equal) then
             ! Something wrong; edge with the same global id should have
             ! the same type and the same facets
             call neko_error('Mismatch in edge and vertex global id')
          end if
       end if
    end if
  end subroutine edge_equal_align

  !> @brief Check if two edges are the same
  !! @note No special treatment of self-periodic edges
  !! @parameter[in]  other    second edge
  !! @return   equal
  function edge_equal(this, other) result(equal)
    class(edge_cab_t), intent(in) :: this
    class(polytope_cnn_t), intent(in) :: other
    logical :: equal
    integer(i4) :: algn

    call this%eq_algn(other, equal, algn)
  end function edge_equal

  !> @brief Initialise edge with alignment information
  !! @parameter[in]   edge   edge
  !! @parameter[in]   algn   alignment
  subroutine edge_aligned_init(this, edge, algn)
    class(edge_aligned_cab_t), intent(inout) :: this
    type(edge_cab_t), intent(in), target :: edge
    integer(i4), intent(in) :: algn

    call this%free_algn()
    ! set global edge pointer
    this%edge%ptr => edge
    ! set relative alignment transformation
    call this%algn_op%init(algn)
  end subroutine edge_aligned_init

  !> @brief Free edge with alignment information
  subroutine edge_aligned_free(this)
    class(edge_aligned_cab_t), intent(inout) :: this

    this%edge%ptr => null()
    call this%algn_op%free()
  end subroutine edge_aligned_free

  !> @brief Return pointer to the edge
  !! @parameter[out]  edge   edge pointer
  subroutine edge_aligned_edgep(this, edge)
    class(edge_aligned_cab_t), intent(in) :: this
    type(edge_cab_ptr), intent(out) :: edge
    edge%ptr => this%edge%ptr
  end subroutine edge_aligned_edgep

  !> @brief Get relative edge alignment
  !! @return   alignment
  pure function edge_aligned_alignment_get(this) result(alignment)
    class(edge_aligned_cab_t), intent(in) :: this
    integer(i4) :: alignment
    alignment = this%algn_op%alignment
  end function edge_aligned_alignment_get

  !> @brief Check if two edges are properly aligned
  !! @note No special treatment of self-periodic edges, so this is not checked
  !! @parameter[in]  other    second edge
  !! @return   aligned
  function edge_aligned_test(this, other) result(aligned)
    class(edge_aligned_cab_t), intent(in) :: this
    class(edge_cab_t), intent(in) :: other
    logical :: aligned
    integer(i4), dimension(NEKO_EDGE_NFACET) :: vrt, vrto

    ! only equal edges can be checked
    if (this%edge%ptr .eq. other) then
       vrt(1) = this%edge%ptr%facet(1)%ptr%id()
       vrt(2) = this%edge%ptr%facet(2)%ptr%id()
       vrto(1) = other%facet(1)%ptr%id()
       vrto(2) = other%facet(2)%ptr%id()
       call this%algn_op%trns_inv_f_i4%ptr(NEKO_EDGE_NFACET, vrto)
       aligned = (vrt(1) == vrto(1)) .and. (vrt(2) == vrto(2))
    else
       call neko_error('Edges not equal')
    end if
  end function edge_aligned_test

  !> @brief Initialise 3D edge actualisation
  !! @parameter[in]   edge    edge
  !! @parameter[in]   algn    alignment
  !! @parameter[in]   pos     position in the object
  subroutine edge_3d_hanging_init(this, edge, algn, pos)
    class(edge_3d_ncnf_cac_t), intent(inout) :: this
    type(edge_cab_t), target, intent(in) :: edge
    integer(i4), intent(in) :: algn, pos

    call this%free()
    call this%init_algn(edge, algn)
    this%position = pos
    ! Assume not hanging edge
    call this%intp_op%init(0)
  end subroutine edge_3d_hanging_init

  !> @brief Free 3D edge actualisation
  subroutine edge_3d_hanging_free(this)
    class(edge_3d_ncnf_cac_t), intent(inout) :: this
    call this%free_algn()
    call this%intp_op%free()
    this%position = -1
  end subroutine edge_3d_hanging_free

  !> @brief Set hanging information
  !! @parameter[in]   hng     hanging information
  subroutine edge_3d_hanging_set(this, hng)
    class(edge_3d_ncnf_cac_t), intent(inout) :: this
    integer(i4), intent(in) :: hng
    call this%intp_op%free()
    call this%intp_op%init(hng)
  end subroutine edge_3d_hanging_set

  !> @brief Get hanging information
  !! @return   hng
  pure function edge_3d_hanging_get(this) result(hng)
    class(edge_3d_ncnf_cac_t), intent(in) :: this
    integer(i4) :: hng
    hng = this%intp_op%hanging
  end function edge_3d_hanging_get

  !> @brief Set position information
  !! @parameter[in]   pos     position information
  pure subroutine edge_3d_position_set(this, pos)
    class(edge_3d_ncnf_cac_t), intent(inout) :: this
    integer(i4), intent(in) :: pos
    this%position = pos
  end subroutine edge_3d_position_set

  !> @brief Get position information
  !! @return   pos
  pure function edge_3d_position_get(this) result(pos)
    class(edge_3d_ncnf_cac_t), intent(in) :: this
    integer(i4) :: pos
    pos = this%position
  end function edge_3d_position_get

  !> @brief Initialise 2D edge actualisation
  !! @parameter[in]   edge    edge
  !! @parameter[in]   algn    alignment
  !! @parameter[in]   pos     position in the object
  !! @parameter[in]   bnd     boundary information
  subroutine edge_2d_hanging_init(this, edge, algn, pos, bnd)
    class(edge_2d_ncnf_cac_t), intent(inout) :: this
    type(edge_cab_t), target, intent(in) :: edge
    integer(i4), intent(in) :: algn, pos, bnd

    call this%free()
    call this%init_3d(edge, algn, pos)
    this%boundary = bnd
  end subroutine edge_2d_hanging_init

  !> @brief Free 2D edge actualisation
  subroutine edge_2d_hanging_free(this)
    class(edge_2d_ncnf_cac_t), intent(inout) :: this

    call this%edge_3d_ncnf_cac_t%free()
    this%boundary = -1
  end subroutine edge_2d_hanging_free

  !> @brief Set boundary information
  !! @parameter[in]   bnd     boundary information
  pure subroutine edge_2d_boundary_set(this, bnd)
    class(edge_2d_ncnf_cac_t), intent(inout) :: this
    integer(i4), intent(in) :: bnd
    this%boundary = bnd
  end subroutine edge_2d_boundary_set

  !> @brief Get boundary information
  !! @return   bnd
  pure function edge_2d_boundary_get(this) result(bnd)
    class(edge_2d_ncnf_cac_t), intent(in) :: this
    integer(i4) :: bnd
    bnd = this%boundary
  end function edge_2d_boundary_get

end module edge_cnn
