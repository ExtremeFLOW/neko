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
!> Various edge types
module edge
  use num_types, only : i4
  use utils, only : neko_error
  use polytope, only : polytope_t
  use polytope_oriented, only : polytope_oriented_t
  use polytope_topology, only : polytope_topology_t, topology_object_t
  use polytope_actualisation, only : polytope_actualisation_t
  use alignment_edge, only : alignment_edge_init, alignment_edge_find
  use vertex, only : NEKO_VERTEX_TDIM
  implicit none
  private

  public :: edge_ornt_t, edge_tpl_t, edge_act_t

  ! object information
  integer(i4), public, parameter :: NEKO_EDGE_TDIM = 1
  integer(i4), public, parameter :: NEKO_EDGE_NFACET = 2
  integer(i4), public, parameter :: NEKO_EDGE_NRIDGE = 0
  integer(i4), public, parameter :: NEKO_EDGE_NPEAK = 0

  !> Type for edge as topology object
  !! @details Edge is the only realisation of one-dimensional polytope (dion)
  !! and contains unique global id and two vertices. Edge is an object
  !! with alignment.
  !! @verbatim
  !! Node numbering
  !!      f_1-----f_2    +----> r
  !! @endverbatim
  !! Its only actualisation are components of higher-dimension objects.
  !! Depending on mesh dimension edges can contain internal/external boundary
  !! information (in case of 2D meshes they are face facets). To simplify type
  !! structure I add this information here. This information is stored in
  !! @a boundary field (0 -internal, 1 - periodic, ....). In case of 3D mesh
  !! its value should be set to -1 and not used.
  type, extends(polytope_topology_t) :: edge_tpl_t
   contains
     !> Initialise a topology polytope
     procedure, pass(this) :: init => edge_tpl_init
     !> Test equality
     procedure, pass(this) :: equal => edge_tpl_equal
  end type edge_tpl_t

  !> Type for edge as aligned object
  !! @details There are two alignment operations for edges: identity and
  !! permutation
  type, extends(polytope_oriented_t) :: edge_ornt_t
   contains
     !> Initialise an aligned polytope
     procedure, pass(this) :: init => edge_ornt_init
     !> Test equality and find alignment
     procedure, pass(this) :: equal_algn => edge_ornt_equal
     !> Test alignment
     procedure, pass(this) :: test  => edge_ornt_test
  end type edge_ornt_t

  !> Actualisation of the topology edge for nonconforming meshes
  !! @details Edge can be either independent (part of conforming interface or
  !! parent; marked 0) or hanging. The meaning of hanging depends on the mesh
  !! dimension. In 2D case edges are face facets and can be either the first
  !! (marked 1) or the second (marked 2) half of a full edge. In 3D case hanging
  !! edges can be either ridge hanging (edge is a component of the full ridge,
  !! but neither of the two facets touching the ridge is hanging), or facet
  !! hanging. Ridge hanging edges can be either the first (marked 1) or the
  !! second (marked 2) half of the full edge. The facet hanging edges can be
  !! either located in the middle of the facet (marked 5) or at the facet
  !! border. In the last case they can correspond to the the first (marked 3)
  !! or the second (marked 4) half of the full edge. See the diagram below
  !! for clarification.
  !! @verbatim
  !! Example of edge marking for two-dimensional nonconforming quad meshes.
  !!  o...............o +...............o
  !!  :               : |               :
  !!  :               : 2               :
  !!  :               : |               :
  !!  +               : +               :
  !!  |               : :               :
  !!  1               : :               :
  !!  |               : :               :
  !!  +...............o o...............o
  !!
  !! Example of edge marking for three-dimensional nonconforming hex meshes.
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
  !! Edges have alignment. There are two different 1D interpolation operators
  !! corresponding the first and the second half of the full edge (marked 1
  !! and 2). The rest of the options (marked 3, 4, and 5) result from the 2D
  !! face interpolations.
  !! Edge is always a component part of the higher-dimension object, so
  !! position gives it's location in the object. It has a single dimension @a r.
  type, extends(polytope_actualisation_t) :: edge_act_t
   contains
     !> Initialise a polytope actualisation
     procedure, pass(this) :: init => edge_act_init
     !> Test equality and find alignment
     procedure, pass(this) :: equal_algn => edge_act_equal
     !> Test alignment
     procedure, pass(this) :: test => edge_act_test
  end type edge_act_t

contains

  !> Initialise a polytope with boundary information
  !! @parameter[in]      id     polytope id
  !! @parameter[in]      nfct   number of facets
  !! @parameter[inout]   fct    polytope facets
  !! @parameter[in]      bnd    external boundary information
  subroutine edge_tpl_init(this, id, nfct, fct, bnd)
    class(edge_tpl_t), intent(inout) :: this
    integer(i4), intent(in) :: id, nfct, bnd
    type(topology_object_t), dimension(nfct), intent(inout) :: fct
    integer(i4) :: il

    call this%free()

    call this%set_tdim(NEKO_EDGE_TDIM)
    call this%set_nelem(NEKO_EDGE_NFACET, NEKO_EDGE_NRIDGE,&
         & NEKO_EDGE_NPEAK)
    call this%set_id(id)
    ! edge can have boundary information for 2D meshes
    call this%set_bnd(bnd)
    ! get facets
    if (nfct == NEKO_EDGE_NFACET) then
       allocate (this%facet(NEKO_EDGE_NFACET))
       do il = 1, nfct
          ! There is just a single realisation of monon, so just check dimension
          if (NEKO_VERTEX_TDIM == fct(il)%obj%polytope%tdim()) then
             call move_alloc(fct(il)%obj, this%facet(il)%obj)
          else
             call neko_error('Edge; wrong facet dimension')
          end if
       end do
    else
       call neko_error('Edge; Inconsistent number of facets.')
    end if
  end subroutine edge_tpl_init

  !> Test equality
  !! @parameter[in]   other   polytope
  !! @return equal
  function edge_tpl_equal(this, other) result(equal)
    class(edge_tpl_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: equal
    integer(i4), dimension(NEKO_EDGE_NFACET, 1) :: edg1, edg2
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
          type is (edge_tpl_t)
             edg1(1, 1) = other%facet(1)%obj%polytope%id()
             edg1(2, 1) = other%facet(2)%obj%polytope%id()
             edg2(1, 1) = this%facet(1)%obj%polytope%id()
             edg2(2, 1) = this%facet(2)%obj%polytope%id()
             call alignment_edge_find(equal, algn, edg1, edg2, &
                  & NEKO_EDGE_NFACET, 1)
          class default
             equal = .false.
          end select
          if (.not. equal) then
             ! Something wrong; edge with the same global id should have
             ! the same type and the same facets
             call neko_error('Mismatch in class or edge and vertex global id')
          end if
       end if
    end if
  end function edge_tpl_equal

  !> Initialise a polytope with alignment information
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  subroutine edge_ornt_init(this, pltp, algn)
    class(edge_ornt_t), intent(inout) :: this
    class(polytope_t), target, intent(in) :: pltp
    integer(i4), intent(in) :: algn
    logical :: ifalgn

    call this%free()

    ! There is just a single realisation of dion, so just check dimension
    if (pltp%tdim() == NEKO_EDGE_TDIM) then
       ! set alignment operator
       call alignment_edge_init(algn, this%algn_op)
       ! mark non identity alignment
       ifalgn = .not. this%algn_op%ifid()
       call this%init_data(pltp, ifalgn)
    else
       call neko_error('Edge aligned; wrong pointer dimension.')
    end if
  end subroutine edge_ornt_init

  !> Check polytope equality and return alignment
  !! @note This subroutine does not take into account the alignment defined in
  !! @ref polytope_aligned_t, but refers directly to the topology object. This
  !! means the result is not the relative orientation between @a this and
  !! @a other, but the absolute on between @a this%polytope and @a other.
  !! @parameter[in]    other  polytope
  !! @parameter[out]   equal  polytope equality
  !! @parameter[out]   algn   alignment information
  subroutine edge_ornt_equal(this, other, equal, algn)
    class(edge_ornt_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn
    integer(i4), dimension(NEKO_EDGE_NFACET, 1) :: edg1, edg2
    class(polytope_t), pointer :: fct

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
          type is (edge_tpl_t)
             edg1(1, 1) = other%facet(1)%obj%polytope%id()
             edg1(2, 1) = other%facet(2)%obj%polytope%id()
             fct => this%polytope%fct(1)
             edg2(1, 1) = fct%id()
             fct => this%polytope%fct(2)
             edg2(2, 1) = fct%id()
             call alignment_edge_find(equal, algn, edg1, edg2, &
                  & NEKO_EDGE_NFACET, 1)
          class default
             equal = .false.
          end select
          if (.not. equal) then
             ! Something wrong; edge with the same global id should have
             ! the same type and the same facets
             call neko_error('Mismatch in class or edge and vertex global id')
          end if
       end if
    end if
  end subroutine edge_ornt_equal

  !> Test alignment
  !! @parameter[in]   other   polytope
  !! @return ifalgn
  function edge_ornt_test(this, other) result(ifalgn)
    class(edge_ornt_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: ifalgn
    integer(i4), dimension(NEKO_EDGE_NFACET, 1) :: trans
    class(polytope_t), pointer :: tfct1, tfct2, ofct1, ofct2

    ! check polygon information there is just one realisation
    ifalgn = (this%polytope%tdim() == other%tdim()) .and. &
         & (this%polytope%id() == other%id())
    ! check vertices getting possible orientation (doesn't work for
    ! self-periodic)
    if (ifalgn) then
       trans(1, 1) = 1
       trans(2, 1) = 2
       call this%algn_op%trns_inv_i4(trans, NEKO_EDGE_NFACET, 1, 1)
       tfct1 => this%polytope%fct(1)
       tfct2 => this%polytope%fct(2)
       ofct1 => other%fct(trans(1, 1))
       ofct2 => other%fct(trans(2, 1))
       ifalgn = (tfct1%id() == ofct1%id()) .and. (tfct2%id() == ofct2%id())
    end if
  end function edge_ornt_test

  !> Initialise a polytope with hanging information
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  !! @parameter[in]   ifint  interpolation flag
  !! @parameter[in]   hng    hanging information
  !! @parameter[in]   pos    position in the higher order element
  subroutine edge_act_init(this, pltp, algn, ifint, hng, pos)
    class(edge_act_t), intent(inout) :: this
    class(polytope_t), target, intent(in) :: pltp
    integer(i4), intent(in) :: algn, hng, pos
    logical, intent(in) :: ifint
    logical :: ifalgn

    call this%free()

    ! There is just a single realisation of dion, so just check dimension
    if (pltp%tdim() == NEKO_EDGE_TDIM) then
       ! set alignment operator
       call alignment_edge_init(algn, this%algn_op)
       ! mark non identity alignment
       ifalgn = .not. this%algn_op%ifid()
       if (hng >= 0 .and. hng <= 5) then
          call this%init_dat(pltp, ifalgn, ifint, hng, pos)
       else
          call neko_error('Inconsistent edge hanging information.')
       end if
    else
       call neko_error('Edge actualisation; wrong pointer dimension.')
    end if
  end subroutine edge_act_init

  !> Check polytope equality and return alignment
  !! @note This subroutine does not take into account the alignment defined in
  !! @ref polytope_aligned_t, but refers directly to the topology object. This
  !! means the result is not the relative orientation between @a this and
  !! @a other, but the absolute on between @a this%polytope and @a other.
  !! @parameter[in]    other  polytope
  !! @parameter[out]   equal  polytope equality
  !! @parameter[out]   algn   alignment information
  subroutine edge_act_equal(this, other, equal, algn)
    class(edge_act_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn
    integer(i4), dimension(NEKO_EDGE_NFACET, 1) :: edg1, edg2
    class(polytope_t), pointer :: fct

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
          type is (edge_tpl_t)
             edg1(1, 1) = other%facet(1)%obj%polytope%id()
             edg1(2, 1) = other%facet(2)%obj%polytope%id()
             fct => this%polytope%fct(1)
             edg2(1, 1) = fct%id()
             fct => this%polytope%fct(2)
             edg2(2, 1) = fct%id()
             call alignment_edge_find(equal, algn, edg1, edg2, &
                  & NEKO_EDGE_NFACET, 1)
          class default
             equal = .false.
          end select
          if (.not. equal) then
             ! Something wrong; edge with the same global id should have
             ! the same type and the same facets
             call neko_error('Mismatch in class or edge and vertex global id')
          end if
       end if
    end if
  end subroutine edge_act_equal

  !> Test alignment
  !! @parameter[in]   other   polytope
  !! @return ifalgn
  function edge_act_test(this, other) result(ifalgn)
    class(edge_act_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: ifalgn
    integer(i4), dimension(NEKO_EDGE_NFACET, 1) :: trans
    class(polytope_t), pointer :: tfct1, tfct2, ofct1, ofct2

    ! check polygon information; there is just one realisation
    ifalgn = (this%polytope%tdim() == other%tdim()) .and. &
         & (this%polytope%id() == other%id())
    ! check vertices getting possible orientation (doesn't work for
    ! self-periodic)
    if (ifalgn) then
       trans(1, 1) = 1
       trans(2, 1) = 2
       call this%algn_op%trns_inv_i4(trans, NEKO_EDGE_NFACET, 1, 1)
       tfct1 => this%polytope%fct(1)
       tfct2 => this%polytope%fct(2)
       ofct1 => other%fct(trans(1, 1))
       ofct2 => other%fct(trans(2, 1))
       ifalgn = (tfct1%id() == ofct1%id()) .and. (tfct2%id() == ofct2%id())
    end if
  end function edge_act_test

end module edge
