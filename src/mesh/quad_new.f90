! Copyright (c) 2018-2024, The Neko Authors
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
!> Various quad types
module quad_new
  use num_types, only : i4, dp
  use utils, only : neko_error
  use polytope, only : polytope_t, polytope_ptr
  use topology, only : topology_t, topology_component_t
  use polytope_actualisation, only : polytope_actualisation_t
  use element_new, only : element_new_t, element_component_t
  use alignment_quad, only : alignment_quad_init, alignment_quad_find
  use vertex, only : vertex_ornt_t, vertex_act_t
  use edge, only : NEKO_EDGE_TDIM, NEKO_EDGE_NFACET
  use point, only : point_t, point_ptr
  implicit none
  private

  public :: quad_tpl_t, quad_act_t, quad_elm_t

  ! object information
  integer(i4), public, parameter :: NEKO_QUAD_TDIM = 2
  integer(i4), public, parameter :: NEKO_QUAD_NFACET = 4
  integer(i4), public, parameter :: NEKO_QUAD_NRIDGE = 4
  integer(i4), public, parameter :: NEKO_QUAD_NPEAK = 0

  !> Type for quad as topology object
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
  !! Abstract quad actualisation are facets of cells. For three-dimensional
  !! meshes quads contain internal/external boundary information.
  type, extends(topology_t) :: quad_tpl_t
   contains
     !> Initialise a topology polytope
     procedure, pass(this) :: init => quad_tpl_init
     !> Test equality
     procedure, pass(this) :: equal => quad_tpl_equal
  end type quad_tpl_t

  !> Actualisation of the topology quad for 3D nonconforming meshes
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
  !! For 3D meshes quad is a facet of a cell, so @a position gives it's location
  !! in the cell.
  type, extends(polytope_actualisation_t) :: quad_act_t
   contains
     !> Initialise a polytope actualisation
     procedure, pass(this) :: init => quad_act_init
     !> Test equality and find alignment
     procedure, pass(this) :: equal_algn => quad_act_equal
     !> Test alignment
     procedure, pass(this) :: test => quad_act_test
  end type quad_act_t

  !> Type for quad as a mesh element
  !! @details For two-dimensional meshes quads are elements. For facet and
  !! ridge numbering see @ref quad_tpl_t. Local numbering of the geometrical
  !! points corresponds to the ridge numbering.
  type, extends(element_new_t) :: quad_elm_t
   contains
     !> Initialise a topology polytope
     procedure, pass(this)  :: init => quad_elm_init
     !> Test equality
     procedure, pass(this) :: equal => quad_elm_equal
     !> Return element diameter
     procedure, pass(this) :: diameter => quad_elm_diameter
     !> Return element centroid
     procedure, pass(this) :: centroid => quad_elm_centroid
     !> Return facet @a r and @s local directions with respect to the element
     procedure, pass(this) :: fct_dir => quad_elm_fct_dir
     !> Return ridge @a r local direction with respect to the element
     procedure, pass(this) :: rdg_dir => quad_elm_rdg_dir
  end type quad_elm_t

  ! Lookup tables
  !> Facet corner to ridge
  integer, public, parameter, dimension(2, 4) :: fct_to_rdg = reshape((/1, 3, &
       & 2, 4,   1, 2,   3, 4 /), shape(fct_to_rdg))
  !> Facets connected to the ridge
  integer, public, parameter, dimension(2, 4) :: rdg_to_fct = reshape((/1, 3, &
       & 2, 3,   1, 4,   2, 4 /), shape(rdg_to_fct))
  !> Ridge to the facet corners (-1 means ridge is not part of the facet)
  integer, public, parameter, dimension(4, 4) :: rdg_to_fctc = reshape((/ &
       & 1, -1, 1, -1,   -1, 1, 2, -1,   2, -1, -1, 1,   -1, 2, -1, 2 /), &
       & shape(rdg_to_fctc))

  !> Transformation of the edge alignment (0,1) with respect to the edge
  !! position on the quad and the quad alignment
  integer, public, parameter, dimension(0: 1, 4, 0: 7) :: quad_to_edg_algn =&
       & reshape((/&
       & 0, 1  ,  0, 1  , 0, 1  , 0, 1  , & ! I
       & 0, 1  ,  0, 1  , 0, 1  , 0, 1  , & ! T
       & 0, 1  ,  0, 1  , 1, 0  , 1, 0  , & ! PX
       & 0, 1  ,  0, 1  , 1, 0  , 1, 0  , & ! PXT
       & 1, 0  ,  1, 0  , 0, 1  , 0, 1  , & ! PYT
       & 1, 0  ,  1, 0  , 0, 1  , 0, 1  , & ! PY
       & 1, 0  ,  1, 0  , 1, 0  , 1, 0  , & ! PXPYT
       & 1, 0  ,  1, 0  , 1, 0  , 1, 0    & ! PXPY
       &/), shape(quad_to_edg_algn))
  integer, public, parameter, dimension(0:1, 4, 0:7) :: quad_to_edg_algn_inv =&
       & reshape((/&
       & 0, 1  ,  0, 1  ,  0, 1  , 0, 1  , & ! I
       & 0, 1  ,  0, 1  ,  0, 1  , 0, 1  , & ! T
       & 0, 1  ,  0, 1  ,  1, 0  , 1, 0  , & ! PX
       & 1, 0  ,  1, 0  ,  0, 1  , 0, 1  , & ! PYT
       & 0, 1  ,  0, 1  ,  1, 0  , 1, 0  , & ! PXT
       & 1, 0  ,  1, 0  ,  0, 1  , 0, 1  , & ! PY
       & 1, 0  ,  1, 0  ,  1, 0  , 1, 0  , & ! PXPYT
       & 1, 0  ,  1, 0  ,  1, 0  , 1, 0    & ! PXPY
       &/), shape(quad_to_edg_algn))

  !> Facet to local direction mapping
  integer, public, parameter, dimension(4) :: fct_to_dir = (/ 2, 2, 1, 1 /)

contains

  !> Initialise a polytope with boundary information
  !! @parameter[in]      id     polytope id
  !! @parameter[in]      nfct   number of facets
  !! @parameter[inout]   fct    polytope facets
  !! @parameter[in]      bnd    external boundary information
  subroutine quad_tpl_init(this, id, nfct, fct, bnd)
    class(quad_tpl_t), intent(inout) :: this
    integer(i4), intent(in) :: id, nfct, bnd
    type(topology_component_t), dimension(nfct), intent(inout) :: fct
    integer(i4) :: il, jl, ifct, icrn
    integer(i4), dimension(NEKO_EDGE_NFACET) :: rdg
    type(polytope_ptr), dimension(2) :: vrt

    call this%free()

    call this%set_tdim(NEKO_QUAD_TDIM)
    call this%set_ncomp(NEKO_QUAD_NFACET, NEKO_QUAD_NRIDGE,&
         & NEKO_QUAD_NPEAK)
    call this%set_id(id)
    ! quad can have boundary information for 3D meshes
    call this%set_bnd(bnd)
    ! get facets
    if (nfct == NEKO_QUAD_NFACET) then
       allocate (this%facet(NEKO_QUAD_NFACET))
       do il = 1, nfct
          ! There is just a single realisation of dion, so just check dimension
          if (NEKO_EDGE_TDIM == fct(il)%obj%polytope%tdim()) then
             call move_alloc(fct(il)%obj, this%facet(il)%obj)
          else
             call neko_error('Quad topology; wrong facet dimension')
          end if
       end do
    else
       call neko_error('Quad topology; Inconsistent number of facets.')
    end if

    ! Get ridge pointers checking quad structure and edge orientation
    ! no special treatment of self-periodic edges
    allocate(this%ridge(NEKO_QUAD_NRIDGE))
    do il = 1, NEKO_QUAD_NRIDGE
       allocate(vertex_ornt_t :: this%ridge(il)%obj)
       ! find proper vertices
       ! Ridge is shared by 2 facets
       do jl = 1, 2
          ifct = rdg_to_fct(jl, il)
          icrn = rdg_to_fctc(ifct, il)
          ! check face alignment
          if (this%facet(ifct)%obj%ifalgn()) then
             ! mark vertices
             rdg(1) = 1
             rdg(2) = 2
             ! transformation
             call this%facet(ifct)%obj%algn_op%trns_inv_i4(rdg, &
                  & NEKO_EDGE_NFACET, 1, 1)
             ! extract vertex
             vrt(jl)%ptr => this%facet(ifct)%obj%polytope%fct(rdg(icrn))
          else
             vrt(jl)%ptr => this%facet(ifct)%obj%polytope%fct(icrn)
          end if
       end do
       if (vrt(1)%ptr%id() == vrt(2)%ptr%id()) then
          call this%ridge(il)%obj%init(vrt(1)%ptr, -1)
       else
          call neko_error('Inconsistent edge vertices in the quad.')
       end if
    end do
  end subroutine quad_tpl_init

  !> Test equality
  !! @parameter[in]   other   polytope
  !! @return equal
  function quad_tpl_equal(this, other) result(equal)
    class(quad_tpl_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: equal
    integer(i4), parameter :: sz = 3
    integer(i4), dimension(sz, sz) :: qad1, qad2
    integer(i4) :: algn

    ! check polygon information
    equal = this%equal_poly(other)
    if (equal) then
       ! check global id
       equal = (this%id() == other%id())
       ! check vertices getting possible orientation (doesn't work for
       ! self-periodic)
       if (equal) then
          select type(other)
          type is (quad_tpl_t)
             ! other
             qad1(2, 2) = 0
             ! mark faces
             qad1(1, 2) = other%facet(1)%obj%polytope%id()
             qad1(3, 2) = other%facet(2)%obj%polytope%id()
             qad1(2, 1) = other%facet(3)%obj%polytope%id()
             qad1(2, 3) = other%facet(4)%obj%polytope%id()
             ! mark vertices
             qad1(1, 1) = other%ridge(1)%obj%polytope%id()
             qad1(3, 1) = other%ridge(2)%obj%polytope%id()
             qad1(1, 3) = other%ridge(3)%obj%polytope%id()
             qad1(3, 3) = other%ridge(4)%obj%polytope%id()
             ! this
             qad2(2, 2) = 0
             ! mark faces
             qad2(1, 2) = this%facet(1)%obj%polytope%id()
             qad2(3, 2) = this%facet(2)%obj%polytope%id()
             qad2(2, 1) = this%facet(3)%obj%polytope%id()
             qad2(2, 3) = this%facet(4)%obj%polytope%id()
             ! mark vertices
             qad2(1, 1) = this%ridge(1)%obj%polytope%id()
             qad2(3, 1) = this%ridge(2)%obj%polytope%id()
             qad2(1, 3) = this%ridge(3)%obj%polytope%id()
             qad2(3, 3) = this%ridge(4)%obj%polytope%id()
             call alignment_quad_find(equal, algn, qad1, qad2, sz, sz)
          class default
             equal = .false.
          end select
          if (.not. equal) then
             ! Something wrong; quads with the same global id should have
             ! the same type and the same facets/ridges
             call neko_error('Mismatch in quad, edge or vertex global id')
          end if
       end if
    end if
  end function quad_tpl_equal

  !> Initialise a polytope with hanging information
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  !! @parameter[in]   ifint  interpolation flag
  !! @parameter[in]   hng    hanging information
  !! @parameter[in]   pos    position in the higher order element
  subroutine quad_act_init(this, pltp, algn, ifint, hng, pos)
    class(quad_act_t), intent(inout) :: this
    class(polytope_t), target, intent(in) :: pltp
    integer(i4), intent(in) :: algn, hng, pos
    logical, intent(in) :: ifint
    logical :: ifalgn

    call this%free()

    ! There are multiple realisations of polygon, so check everything
    if (pltp%tdim() == NEKO_QUAD_TDIM .and. pltp%check_comp(NEKO_QUAD_NFACET, &
         & NEKO_QUAD_NRIDGE, NEKO_QUAD_NPEAK)) then
       ! set alignment operator
       call alignment_quad_init(algn, this%algn_op)
       ! mark non identity alignment
       ifalgn =  .not. this%algn_op%ifid()
       if (hng >= 0 .and. hng <= 4) then
          call this%init_act(pltp, ifalgn, ifint, hng, pos)
       else
          call neko_error('Inconsistent quad hanging information.')
       end if
    else
       call neko_error('Quad actualisation; wrong pointer dimensions.')
    end if
  end subroutine quad_act_init

  !> Check polytope equality and return alignment
  !! @note This subroutine does not take into account the alignment defined in
  !! @ref polytope_aligned_t, but refers directly to the topology object. This
  !! means the result is not the relative orientation between @a this and
  !! @a other, but the absolute on between @a this%polytope and @a other.
  !! @parameter[in]    other  polytope
  !! @parameter[out]   equal  polytope equality
  !! @parameter[out]   algn   alignment information
  subroutine quad_act_equal(this, other, equal, algn)
    class(quad_act_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn
    integer(i4), parameter :: sz = 3
    integer(i4), dimension(sz, sz) :: qad1, qad2
    class(polytope_t), pointer :: ptr

    algn = -1
    ! check polygon information
    equal = this%polytope%equal_poly(other)
    if (equal) then
       ! check global id
       equal = (this%polytope%id() == other%id())
       ! check vertices getting possible orientation (doesn't work for
       ! self-periodic)
       if (equal) then
          select type(other)
          type is (quad_tpl_t)
             ! other
             qad1(2, 2) = 0
             ! mark faces
             qad1(1, 2) = other%facet(1)%obj%polytope%id()
             qad1(3, 2) = other%facet(2)%obj%polytope%id()
             qad1(2, 1) = other%facet(3)%obj%polytope%id()
             qad1(2, 3) = other%facet(4)%obj%polytope%id()
             ! mark vertices
             qad1(1, 1) = other%ridge(1)%obj%polytope%id()
             qad1(3, 1) = other%ridge(2)%obj%polytope%id()
             qad1(1, 3) = other%ridge(3)%obj%polytope%id()
             qad1(3, 3) = other%ridge(4)%obj%polytope%id()
             ! this
             qad2(2, 2) = 0
             ! mark faces
             ptr => this%polytope%fct(1)
             qad2(1, 2) = ptr%id()
             ptr => this%polytope%fct(2)
             qad2(3, 2) = ptr%id()
             ptr => this%polytope%fct(3)
             qad2(2, 1) = ptr%id()
             ptr => this%polytope%fct(4)
             qad2(2, 3) = ptr%id()
             ! mark vertices
             ptr => this%polytope%rdg(1)
             qad2(1, 1) = ptr%id()
             ptr => this%polytope%rdg(2)
             qad2(3, 1) = ptr%id()
             ptr => this%polytope%rdg(3)
             qad2(1, 3) = ptr%id()
             ptr => this%polytope%rdg(4)
             qad2(3, 3) = ptr%id()
             call alignment_quad_find(equal, algn, qad1, qad2, sz, sz)
          class default
             equal = .false.
          end select
          if (.not. equal) then
             ! Something wrong; quads with the same global id should have
             ! the same type and the same facets/ridges
             call neko_error('Mismatch in quad, edge or vertex global id')
          end if
       end if
    end if
  end subroutine quad_act_equal

  !> Test alignment
  !! @parameter[in]   other   polytope actualisation
  !! @return ifalgn
  function quad_act_test(this, other) result(ifalgn)
    class(quad_act_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: ifalgn
    integer(i4), parameter :: sz = 3
    integer(i4), dimension(sz, sz) :: trans
    integer(i4), dimension(NEKO_QUAD_NFACET) :: mapf
    integer(i4), dimension(NEKO_QUAD_NRIDGE) :: mapr
    integer(i4) :: il
    class(polytope_t), pointer :: tptr, optr

    ! check polygon information; multiple possible realisations
    ifalgn = this%polytope%equal_poly(other)
    if (ifalgn) then
       ! check global id
       ifalgn = (this%polytope%id() == other%id())
       ! check vertices getting possible orientation (doesn't work for
       ! self-periodic)
       if (ifalgn) then
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
          call this%algn_op%trns_inv_i4(trans, sz, sz, 1)
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

          ifalgn = .true.
          ! check ridges
          do il = 1, this%polytope%nridge
             tptr => this%polytope%rdg(il)
             optr => other%rdg(mapr(il))
             ifalgn = ifalgn .and. (tptr%id() == optr%id())
          end do
          if (ifalgn) then
             ! check facets discarding orientation (covered by vertices)
             do il = 1, this%polytope%nfacet
                tptr => this%polytope%fct(il)
                optr => other%fct(mapf(il))
                ifalgn = ifalgn .and. (tptr%id() == optr%id())
             end do
          end if
       end if
    end if
  end function quad_act_test

  !> Initialise a polytope with geometry information
  !! @parameter[in]   id       polytope id
  !! @parameter[in]   nfct     number of facets
  !! @parameter[in]   fct      polytope facets
  !! @parameter[in]   npts     number of points
  !! @parameter[in]   pts      points
  !! @parameter[in]   gdim     geometrical dimension
  !! @parameter[in]   nrdg     number of hanging ridges
  !! @parameter[in]   rdg_hng  ridge hanging flag
  subroutine quad_elm_init(this, id, nfct, fct, npts, pts, gdim, nrdg, &
          & rdg_hng)
    class(quad_elm_t), intent(inout) :: this
    integer(i4), intent(in) :: id, nfct, npts, gdim, nrdg
    type(element_component_t), dimension(nfct), intent(inout) :: fct
    type(point_ptr), dimension(npts), intent(in) :: pts
    integer(i4), dimension(2, 3), intent(in) :: rdg_hng
    integer(i4) :: il, jl, ifct, icrn, itmp
    integer(i4), dimension(NEKO_EDGE_NFACET) :: rdg
    integer(i4), dimension(3, 2) :: hng_fct
    integer(i4), dimension(2) :: hng
    type(polytope_ptr), dimension(2) :: vrt

    call this%free()

    call this%set_tdim(NEKO_QUAD_TDIM)
    call this%set_ncomp(NEKO_QUAD_NFACET, NEKO_QUAD_NRIDGE, &
         & NEKO_QUAD_NPEAK)
    call this%set_id(id)
    call this%init_base(gdim, npts)
    ! get facets
    if (nfct == NEKO_QUAD_NFACET) then
       allocate (this%facet(NEKO_QUAD_NFACET))
       do il = 1, nfct
          ! There is just a single realisation of dion, so just check dimension
          if (NEKO_EDGE_TDIM == fct(il)%obj%polytope%tdim()) then
             ! facet position
             jl = fct(il)%obj%pos()
             call move_alloc(fct(il)%obj, this%facet(jl)%obj)
          else
             call neko_error('Quad mesh; wrong facet dimension')
          end if
       end do
    else
       call neko_error('Quad mesh; inconsistent number of facets.')
    end if

    ! Get ridge pointers checking quad structure and edge orientation
    ! no special treatment of self-periodic edges
    allocate(this%ridge(NEKO_QUAD_NRIDGE))
    do il = 1, NEKO_QUAD_NRIDGE
       allocate(vertex_act_t :: this%ridge(il)%obj)
       ! find proper vertices and collect facet hanging information
       ! Ridge is shared by 2 facets
       do jl = 1, 2
          ifct = rdg_to_fct(jl, il)
          icrn = rdg_to_fctc(ifct, il)
          ! check face alignment
          if (this%facet(ifct)%obj%ifalgn()) then
             ! mark vertices
             rdg(1) = 1
             rdg(2) = 2
             ! transformation
             call this%facet(ifct)%obj%algn_op%trns_inv_i4(rdg, &
                  & NEKO_EDGE_NFACET, 1, 1)
             ! extract vertex
             vrt(jl)%ptr => this%facet(ifct)%obj%polytope%fct(rdg(icrn))
          else
             vrt(jl)%ptr => this%facet(ifct)%obj%polytope%fct(icrn)
          end if
          ! collect facet hanging information
          hng_fct(1, jl) = this%facet(ifct)%obj%hng()
          hng_fct(2, jl) = ifct
          hng_fct(3, jl) = icrn
       end do
       ! is it a proper vertex
       if (vrt(1)%ptr%id() == vrt(2)%ptr%id()) then
          ! Vertex hanging flag depends on the edge hanging flag and the
          ! vertex position in the edge (corner).
          ! Ridge is shared by 2 facets
          do jl = 1, 2
             select case (hng_fct(1, jl))
             case (0) ! independent edge
                hng(jl) = 0 ! independent vertex
             case (1) ! face (2D) hanging; lower part
                select case(hng_fct(3, jl)) ! vertex position in the edge
                case (1)
                   hng(jl) = 0 ! independent vertex
                case (2)
                   hng(jl) = 1 ! face hanging vertex
                case default
                   call neko_error('Quad mesh; wrong corner position')
                end select
             case (2) ! face (2D) hanging; upper part
                select case(hng_fct(3, jl)) ! vertex position in the edge
                case (1)
                   hng(jl) = 1 ! face hanging vertex
                case (2)
                   hng(jl) = 0 ! independent vertex
                case default
                   call neko_error('Quad mesh; wrong corner position')
                end select
             case default
                call neko_error('Quad mesh; facet hanging flag')
             end select
          end do
          ! A single vertex can be at the same time marked as independent and
          ! face hanging. Take the max.
          itmp = maxval(hng)
          ! No need to check rdg_hng for 2D meshes.
          ! create vertex actualisation
          ! Vertex has neither alignment, nor interpolation
          call this%ridge(il)%obj%init(vrt(1)%ptr, -1, .false., itmp, il)
       else
          call neko_error('Inconsistent edge vertices in the mesh quad.')
       end if
    end do

    ! Add geometrical points
    ! Sanity check
    if (npts /= NEKO_QUAD_NRIDGE) &
         &call neko_error('Inconsistent point number in the mesh quad.')
    allocate(this%pts(NEKO_QUAD_NRIDGE))
    do il = 1, npts
       this%pts(il)%p => pts(il)%p
    end do
  end subroutine quad_elm_init

  !> Test equality
  !! @parameter[in]   other   polytope
  !! @return equal
  function quad_elm_equal(this, other) result(equal)
    class(quad_elm_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: equal

    ! check polygon information
    equal = this%equal_poly(other)
    if (equal) then
       ! check global id
       equal = (this%id() == other%id())
       ! check edges and vertices getting discarding orientation
       ! (may not work for self-periodic)
       if (equal) then
          select type(other)
          type is (quad_elm_t)
             ! geometrical dimension
             equal = (this%gdim() == other%gdim())
             if (equal) then

                call neko_error('Not finished; missing sorting routines.')

             end if
          class default
             equal = .false.
          end select
          if (.not. equal) then
             ! Something wrong; edge with the same global id should have
             ! the same type and the same facets/ridges
             call neko_error('Mismatch in class or element global id')
          end if
       end if
    end if

  end function quad_elm_equal

  !> Get element diameter
  !! @return res
  function quad_elm_diameter(this) result(res)
    class(quad_elm_t), intent(in) :: this
    real(dp) :: res
    real(dp) :: d1, d2
    integer(i4) :: il

    d1 = 0d0
    d2 = 0d0
    do il = 1, this%gdim()
       d1 = d1 + (this%pts(4)%p%x(il) - this%pts(1)%p%x(il))**2
       d2 = d2 + (this%pts(3)%p%x(il) - this%pts(2)%p%x(il))**2
    end do

    res = sqrt(max(d1, d2))
  end function quad_elm_diameter

  !> Get element centroid
  !! @return res
  function quad_elm_centroid(this) result(res)
    class(quad_elm_t), intent(in) :: this
    type(point_t) :: res
    integer(i4) :: il

    res%x(:) = 0d0
    do il = 1, this%gdim()
       res%x(il) = 0.25d0 * (this%pts(1)%p%x(il) + this%pts(2)%p%x(il) + &
            & this%pts(3)%p%x(il) + this%pts(4)%p%x(il))
    end do
  end function quad_elm_centroid

  !> Get @a r and @a s facet local directions
  !! @parameter[in]   pos          facet position
  !! @parameter[out]  dirr, dirs   local directions
  subroutine quad_elm_fct_dir(this, pos, dirr, dirs)
    class(quad_elm_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4), intent(out) :: dirr, dirs
    dirr = fct_to_dir(pos)
    dirs = -1
  end subroutine quad_elm_fct_dir

  !> Get @a r ridge local direction
  !! @parameter[in]   pos          ridge position
  !! @parameter[out]  dirr         local direction
  subroutine quad_elm_rdg_dir(this, pos, dirr)
    class(quad_elm_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4), intent(out) :: dirr
    dirr = -1
  end subroutine quad_elm_rdg_dir

end module quad_new
