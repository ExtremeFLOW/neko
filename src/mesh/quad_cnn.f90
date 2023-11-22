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
!> Connectivity quadrilateral types
module quad_cnn
  use num_types, only : i4
  use utils, only : neko_error
  use polytope_cnn, only : polytope_cnn_t
  use vertex_cnn, only : vertex_cab_ptr
  use edge_cnn, only : edge_aligned_cab_t, NEKO_EDGE_NFACET
  use face_cnn, only : face_cab_t
  use alignment_quad, only : alignment_quad_t, alignment_quad_op_set_t
  use ncnf_interpolation_quad, only : ncnf_intp_quad_op_set_t
  implicit none
  private

  public :: quad_cab_t, quad_cab_ptr, quad_aligned_cab_t, quad_ncnf_crl_t,&
       & quad_ncnf_crl_ptr

  ! object information
  integer(i4), public, parameter :: NEKO_QUAD_NFACET = 4
  integer(i4), public, parameter :: NEKO_QUAD_NRIDGE = 4
  integer(i4), public, parameter :: NEKO_QUAD_NPEAK = 0

  !> Type for abstract quad object
  !! @details Quad is one of the realisations of the face containing 4 facets
  !! (edges) and 4 ridges (vertices). Quad is and object with alignment.
  !! Facets are oriented according to the numbers of ridges from the one with
  !! the smaller number towards the one with the bigger number.This corresponds
  !! to the data alignment in the data array.
  !! @verbatim
  !! Facet and ridge numbering (symmetric notation)
  !!
  !!         f_4            ^ s
  !!   r_3+-------+r_4      |
  !!      |       |         |
  !! f_1  |       |  f_2    |
  !!      |       |         |
  !!   r_1+-------+r_2      +-------> r
  !!         f_3
  !!
  !! @endverbatim
  !! Abstract quad realisations are facets of cells.
  type, extends(face_cab_t) :: quad_cab_t
   contains!> Initialise quad
     procedure, pass(this) :: init => quad_init
     !> Check equality and get relative alignment info
     procedure, pass(this) :: eq_algn => quad_equal_align
     !> Quad equality including edge and vertex information
     procedure, pass(this) :: equal => quad_equal
     generic :: operator(.eq.) => equal
  end type quad_cab_t

  !> Pointer to an abstract quad object
  type ::  quad_cab_ptr
     type(quad_cab_t), pointer :: ptr
  end type quad_cab_ptr

  !> Abstract quad object with alignment information
  type :: quad_aligned_cab_t
     !> edge pointer
     type(quad_cab_ptr) :: face
     !> alignment operator
     type(alignment_quad_op_set_t) :: algn_op
   contains
     !> Initialise aligned quad
     procedure, pass(this) :: init_algn => quad_aligned_init
     !> Free aligned quad
     procedure, pass(this) :: free_algn => quad_aligned_free
     !> Return quad pointer
     procedure, pass(this) :: facep => quad_aligned_quadp
     !> Return quad relative alignment
     procedure, pass(this) :: algn => quad_aligned_alignment_get
     !> Test alignment
     procedure, pass(this) :: test => quad_aligned_test
  end type quad_aligned_cab_t

  !> Realisation of the abstract quad for nonconforming meshes
  !! @details Quad can be either independent (part of conforming interface or
  !! parent; marked 0) or hanging. The hanging marking corresponds to the parent
  !! quad corner number adjacent to the child. See the diagram below for
  !! clarification.
  !! @verbatim
  !! Example of quad marking for hex meshes.
  !! Hanging edge marking for three-dimensional nonconforming interface
  !!
  !!  o...............o o...............o +-------+.......o o.......+-------+
  !!  :               : :               : |       |       : :       |       |
  !!  :               : :               : |   3   |       : :       |   4   |
  !!  :               : :               : |       |       : :       |       |
  !!  +-------+       : :       +-------+ +-------+       : :       +-------+
  !!  |       |       : :       |       | :               : :               :
  !!  |   1   |       : :       |   2   | :               : :               :
  !!  |       |       : :       |       | :               : :               :
  !!  +-------+.......o o.......+-------+ o...............o o...............o
  !! @endverbatim
  !! There are four different 2D interpolation operators corresponding the
  !! quarter of a parent quad occupied by the child.
  !! For 3D meshes quad is a facet of a cell, so position gives it's location
  !! in the cell.
  type, extends(quad_aligned_cab_t) :: quad_ncnf_crl_t
     !> interpolation operator
     type(ncnf_intp_quad_op_set_t) :: intp_op
     !> position in the object
     integer(i4) :: position = 0
   contains
     !> Initialise aligned quad pointer and position
     procedure, pass(this) :: init => quad_hanging_init
     !> Free quad data
     procedure, pass(this) :: free => quad_hanging_free
     !> Set hanging information
     procedure, pass(this) :: set_hng => quad_hanging_set
     !> Get hanging information
     procedure, pass(this) :: hng => quad_hanging_get
     !> Set position information
     procedure, pass(this) :: set_pos => quad_position_set
     !> Get position information
     procedure, pass(this) :: pos => quad_position_get
  end type quad_ncnf_crl_t

  !> Pointer to a nonconforming quad realisation
  type ::  quad_ncnf_crl_ptr
     type(quad_ncnf_crl_t), pointer :: ptr
  end type quad_ncnf_crl_ptr

  ! Lookup tables
  !> Facet corner to ridge
  integer, parameter, dimension(2, 4) :: fct_to_rdg = reshape((/1,3, &
       & 2,4, 1,2, 3,4 /), shape(fct_to_rdg))
  !> Facets connected to the ridge
  integer, parameter, dimension(2, 4) :: rdg_to_fct = reshape((/1,3, &
       & 2,3, 1,4, 2,4 /), shape(rdg_to_fct))
  !> Ridge to the facet corners (-1 means ridge is not part of the facet)
  integer, parameter, dimension(4, 4) :: rdg_to_fctc = reshape((/1,-1,1,-1, &
       & -1,1,2,-1, 2,-1,-1,1, -1,2,-1,2 /), shape(rdg_to_fctc))

  !> Transformation of the edge alignment (0,1) with respect to the edge
  !! position on the quad and the quad alignment
  integer, public, parameter, dimension(0:1, 4, 0:7) :: quad_to_edg_algn =&
       & reshape((/&
       & 0,1 , 0,1 , 0,1 , 0,1 , & ! I
       & 0,1 , 0,1 , 0,1 , 0,1 , & ! T
       & 0,1 , 0,1 , 1,0 , 1,0 , & ! PX
       & 0,1 , 0,1 , 1,0 , 1,0 , & ! PXT
       & 1,0 , 1,0 , 0,1 , 0,1 , & ! PYT
       & 1,0 , 1,0 , 0,1 , 0,1 , & ! PY
       & 1,0 , 1,0 , 1,0 , 1,0 , & ! PXPYT
       & 1,0 , 1,0 , 1,0 , 1,0   & ! PXPY
       &/), shape(quad_to_edg_algn))
  integer, public, parameter, dimension(0:1, 4, 0:7) :: quad_to_edg_algn_inv =&
       & reshape((/&
       & 0,1 , 0,1 , 0,1 , 0,1 , & ! I
       & 0,1 , 0,1 , 0,1 , 0,1 , & ! T
       & 0,1 , 0,1 , 1,0 , 1,0 , & ! PX
       & 1,0 , 1,0 , 0,1 , 0,1 , & ! PYT
       & 0,1 , 0,1 , 1,0 , 1,0 , & ! PXT
       & 1,0 , 1,0 , 0,1 , 0,1 , & ! PY
       & 1,0 , 1,0 , 1,0 , 1,0 , & ! PXPYT
       & 1,0 , 1,0 , 1,0 , 1,0   & ! PXPY
       &/), shape(quad_to_edg_algn))

contains
  !> @brief Initialise quad with global id and four edges
  !! @details Edges order defines face orientation
  !! @parameter[in]   id                        unique id
  !! @parameter[in]   edg1, edg2, edg3, edg4    bounding edges
  subroutine quad_init(this, id, edg1, edg2, edg3, edg4)
    class(quad_cab_t), intent(inout) :: this
    integer(i4), intent(in) :: id
    type(edge_aligned_cab_t), intent(in), target :: edg1, edg2, edg3, edg4
    integer(i4) :: il, jl, ifct, icrn
    integer(i4), dimension(NEKO_EDGE_NFACET) :: rdg
    type(vertex_cab_ptr), dimension(2) :: vrt
    logical :: equal

    call this%free()

    ! init_dim calls free
    call this%init_dim()

    call this%set_nelem(NEKO_QUAD_NFACET, NEKO_QUAD_NRIDGE,&
         & NEKO_QUAD_NPEAK)
    call this%set_id(id)
    ! get facet pointers
    allocate(this%facet(NEKO_QUAD_NFACET))
    this%facet(1) = edg1
    this%facet(2) = edg2
    this%facet(3) = edg3
    this%facet(4) = edg4

    ! THIS SHOULD BE IN DIFFERENT PLACE; It checks internal consistency of edges
    ! Check if edges are consistent. Self-periodicity is allowed, but vertices
    ! should not be messed up
    do il = 1, NEKO_QUAD_NFACET - 1
       do jl = il + 1, NEKO_QUAD_NFACET
          equal = this%facet(il)%edge%ptr .eq. this%facet(jl)%edge%ptr
       end do
    end do

    ! Get ridge pointers checking quad structure and edge orientation
    ! no special treatment of self-periodic edges
    allocate(this%ridge(NEKO_QUAD_NRIDGE))
    do il = 1, NEKO_QUAD_NRIDGE
       ! find proper vertices
       do jl = 1, 2
          ifct = rdg_to_fct(jl, il)
          icrn = rdg_to_fctc(ifct, il)
          ! mark vertices
          rdg(1) = 1
          rdg(2) = 2
          ! transformation
          call this%facet(ifct)%algn_op%trns_inv_f_i4%ptr(NEKO_EDGE_NFACET, rdg)
          ! extract vertex
          vrt(jl)%ptr => this%facet(ifct)%edge%ptr%facet(rdg(icrn))%ptr
       end do
       if (vrt(1)%ptr%id() == vrt(2)%ptr%id()) then
          this%ridge(il)%ptr => vrt(1)%ptr
       else
          call neko_error('Inconsistent edge vertices in the quad.')
       end if
    end do
  end subroutine quad_init

  !> @brief Check if two quads are the same
  !! @note Alignment for self-periodic quads will not be correct
  !! @parameter[in]  other    second quad
  !! @return   equal
  subroutine quad_equal_align(this, other, equal, algn)
    class(quad_cab_t), intent(in) :: this
    class(polytope_cnn_t), intent(in) :: other
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn
    type(alignment_quad_t) :: algn_op
    integer(i4), parameter :: sz = 3
    integer(i4), dimension(sz, sz) :: trans
    integer(i4), dimension(sz) :: work
    integer(i4), dimension(NEKO_QUAD_NFACET) :: mapf
    integer(i4), dimension(NEKO_QUAD_NRIDGE) :: mapr
    integer(i4) :: il, algnl
    logical :: equall

    algn = -1
    ! check polygon information
    equal = this%equal_poly(other)
    if (equal) then
       ! check global id
       equal = (this%id() == other%id())
       ! check edges/vertices getting possible orientation (doesn't work for
       ! self-periodic)
       if (equal) then
          call algn_op%init()
          select type(other)
          type is (quad_cab_t)
             ! check all the alignment options
             ! mark faces
             trans(1, 2) = 1
             trans(3, 2) = 2
             trans(2, 1) = 3
             trans(2, 3) = 4
             ! mark vertices
             trans(1, 1) = 1
             trans(3, 1) = 2
             trans(1, 3) = 3
             trans(3, 3) = 4
             do algn = 0, algn_op%nop()
                call algn_op%trns_inv_f_i4(algn)%ptr(sz, trans, work)
                ! get mappings for
                ! facet
                mapf(1) = trans(1, 2)
                mapf(2) = trans(3, 2)
                mapf(3) = trans(2, 1)
                mapf(4) = trans(2, 3)
                ! ridge
                mapr(1) = trans(1, 1)
                mapr(2) = trans(3, 1)
                mapr(3) = trans(1, 3)
                mapr(4) = trans(3, 3)
                equal = .true.
                ! check ridges
                do il = 1, this%nridge
                   equal = equal .and. (this%ridge(il)%ptr%id() == &
                        & other%ridge(mapr(il))%ptr%id())
                end do
                if (equal) then
                   ! check facets discarding orientation (covered by vertices)
                   do il = 1, this%nfacet
                      equal = equal .and. (this%facet(il)%edge%ptr .eq. &
                           & other%facet(mapf(il))%edge%ptr)
                   end do
                   if (equal) return
                end if
                call algn_op%trns_f_i4(algn)%ptr(sz, trans, work)
             end do
          class default
             equal = .false.
          end select
          if (.not.equal) then
             ! Something wrong; quads with the same global id should have
             ! the same type and the same facets/ridges
             call neko_error('Mismatch in quad, edge and vertex global id')
          end if
       end if
    end if
  end subroutine quad_equal_align

  !> @brief Check if two quads are the same
  !! @note No special treatment of self-periodic quads
  !! @parameter[in]  other    second quad
  !! @return   equal
  function quad_equal(this, other) result(equal)
    class(quad_cab_t), intent(in) :: this
    class(polytope_cnn_t), intent(in) :: other
    logical :: equal
    integer(i4) :: algn

    call this%eq_algn(other, equal, algn)
  end function quad_equal

  !> @brief Initialise quad with alignment information
  !! @parameter[in]   quad   quad
  !! @parameter[in]   algn   alignment
  subroutine quad_aligned_init(this, quad, algn)
    class(quad_aligned_cab_t), intent(inout) :: this
    type(quad_cab_t), intent(in), target :: quad
    integer(i4), intent(in) :: algn

    call this%free_algn()
    ! set global edge pointer
    this%face%ptr => quad
    ! set relative alignment transformation
    call this%algn_op%init(algn)
  end subroutine quad_aligned_init

  !> @brief Free quad with alignment information
  subroutine quad_aligned_free(this)
    class(quad_aligned_cab_t), intent(inout) :: this

    this%face%ptr => null()
    call this%algn_op%free()
  end subroutine quad_aligned_free

  !> @brief Return pointer to the quad
  !! @parameter[out]  quad   quad pointer
  subroutine quad_aligned_quadp(this, quad)
    class(quad_aligned_cab_t), intent(in) :: this
    type(quad_cab_ptr), intent(out) :: quad
    quad%ptr => this%face%ptr
  end subroutine quad_aligned_quadp

  !> @brief Get relative quad alignment
  !! @return   alignment
  pure function quad_aligned_alignment_get(this) result(alignment)
    class(quad_aligned_cab_t), intent(in) :: this
    integer(i4) :: alignment
    alignment = this%algn_op%alignment
  end function quad_aligned_alignment_get

  !> @brief Check if two quads are properly aligned
  !! @note No special treatment of self-periodic quads, so this is not checked
  !! @parameter[in]  other    second quad
  !! @return   aligned
  function quad_aligned_test(this, other) result(aligned)
    class(quad_aligned_cab_t), intent(in) :: this
    class(quad_cab_t), intent(in) :: other
    logical :: aligned
    integer(i4), parameter :: sz = 3
    integer(i4), dimension(sz, sz) :: elm, elmo
    integer(i4), dimension(sz) :: work

    ! only equal quads can be checked
    if (this%face%ptr.eq.other) then
       ! edges
       elm(1, 2) = this%face%ptr%facet(1)%edge%ptr%id()
       elm(3, 2) = this%face%ptr%facet(2)%edge%ptr%id()
       elm(2, 1) = this%face%ptr%facet(3)%edge%ptr%id()
       elm(2, 3) = this%face%ptr%facet(4)%edge%ptr%id()
       ! vertices
       elm(1, 1) = this%face%ptr%ridge(1)%ptr%id()
       elm(3, 1) = this%face%ptr%ridge(2)%ptr%id()
       elm(1, 3) = this%face%ptr%ridge(3)%ptr%id()
       elm(3, 3) = this%face%ptr%ridge(4)%ptr%id()
       ! edges
       elmo(1, 2) = other%facet(1)%edge%ptr%id()
       elmo(3, 2) = other%facet(2)%edge%ptr%id()
       elmo(2, 1) = other%facet(3)%edge%ptr%id()
       elmo(2, 3) = other%facet(4)%edge%ptr%id()
       ! vertices
       elmo(1, 1) = other%ridge(1)%ptr%id()
       elmo(3, 1) = other%ridge(2)%ptr%id()
       elmo(1, 3) = other%ridge(3)%ptr%id()
       elmo(3, 3) = other%ridge(4)%ptr%id()
       call this%algn_op%trns_inv_f_i4%ptr( sz, elmo, work)
       aligned = (elm(1, 1) == elmo(1, 1)) .and. (elm(1, 2) == elmo(1, 2)).and.&
            &(elm(1, 3) == elmo(1, 3)) .and. (elm(2, 1) == elmo(2, 1)) .and. &
            &(elm(2, 3) == elmo(2, 3)) .and. (elm(3, 1) == elmo(3, 1)) .and. &
            &(elm(3, 2) == elmo(3, 2)) .and. (elm(3, 3) == elmo(3, 3))
    else
       call neko_error('Quads not equal')
    end if
  end function quad_aligned_test

  !> @brief Initialise aligned quad pointer and position
  !! @parameter[in]   quad    quad
  !! @parameter[in]   algn    alignment
  !! @parameter[in]   pos     position in the object
  subroutine quad_hanging_init(this, quad, algn, pos)
    class(quad_ncnf_crl_t), intent(inout) :: this
    type(quad_cab_t), target, intent(in) :: quad
    integer(i4), intent(in) :: algn, pos

    call this%free()
    call this%init_algn(quad, algn)
    this%position = pos
    ! Assume not hanging edge
    call this%intp_op%init(0)
  end subroutine quad_hanging_init

   !> @brief Free aligned quad and hanging information
  subroutine quad_hanging_free(this)
    class(quad_ncnf_crl_t), intent(inout) :: this
    call this%free_algn()
    call this%intp_op%free()
    this%position = -1
  end subroutine quad_hanging_free

  !> @brief Set hanging information
  !! @parameter[in]   hng     hanging information
  subroutine quad_hanging_set(this, hng)
    class(quad_ncnf_crl_t), intent(inout) :: this
    integer(i4), intent(in) :: hng
    call this%intp_op%free()
    call this%intp_op%init(hng)
  end subroutine quad_hanging_set

  !> @brief Get hanging information
  !! @return   hng
  pure function quad_hanging_get(this) result(hng)
    class(quad_ncnf_crl_t), intent(in) :: this
    integer(i4) :: hng
    hng = this%intp_op%hanging
  end function quad_hanging_get

  !> @brief Set position information
  !! @parameter[in]   pos     position information
  pure subroutine quad_position_set(this, pos)
    class(quad_ncnf_crl_t), intent(inout) :: this
    integer(i4), intent(in) :: pos
    this%position = pos
  end subroutine quad_position_set

  !> @brief Get position information
  !! @return   pos
  pure function quad_position_get(this) result(pos)
    class(quad_ncnf_crl_t), intent(in) :: this
    integer(i4) :: pos
    pos = this%position
  end function quad_position_get

end module quad_cnn
