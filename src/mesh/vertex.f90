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
!> Various vertex types
module vertex
  use num_types, only : i4
  use utils, only : neko_error, neko_warning
  use polytope, only : polytope_t
  use polytope_aligned, only : polytope_aligned_t
  use polytope_topology, only : polytope_topology_t, topology_object_t
  use polytope_actualisation, only : polytope_actualisation_t
  implicit none
  private

  public :: vertex_tpl_t, vertex_algn_t, vertex_act_t

  ! object information
  integer(i4), public, parameter :: NEKO_VERTEX_DIM = 0
  integer(i4), public, parameter :: NEKO_VERTEX_NFACET = 0
  integer(i4), public, parameter :: NEKO_VERTEX_NRIDGE = 0
  integer(i4), public, parameter :: NEKO_VERTEX_NPEAK = 0

  !> Type for vertex as topology object
  !! @details Vertex is the only realisation of zero-dimensional polytope
  !! (monon) and contains a unique global id only. Vertex has no alignment.
  !! Its only actualisation are components of higher-dimension objects.
    type, extends(polytope_topology_t) :: vertex_tpl_t
   contains
     !> Initialise a topology polytope
     procedure, pass(this) :: init => vertex_tpl_init
     !> Test equality
     procedure, pass(this) :: equal => vertex_tpl_equal
  end type vertex_tpl_t

  !> Type for vertex as aligned object
  !! @details Vertex hes no alignment, but the type is required for
  !! completeness.
  type, extends(polytope_aligned_t) :: vertex_algn_t
   contains
     !> Initialise an aligned polytope
     procedure, pass(this) :: init => vertex_algn_init
     !> Test equality and find alignment
     procedure, pass(this) :: equal_algn => vertex_algn_equal
     !> Test alignment
     procedure, pass(this) :: test  => vertex_algn_test
  end type vertex_algn_t

  !> Actualisation of the topology vertex for nonconforming meshes
  !! @details Vertex actualisation can be either independent (located at edge,
  !! face, cell corner of all neighbours; marked 0), facet hanging (located at
  !! facet centre of parent neighbours in case of faces or cells; marked 1) or
  !! ridge hanging (located at ridge centre of parent neighbours in case of
  !! cells; marked 2). See the diagram below for clarification.
  !! @verbatim
  !! Example of vertex marking for quad/hex meshes.
  !! Hanging vertex marking for two-dimensional nonconforming interface
  !!  o...............o 0...............o
  !!  :               : |               :
  !!  :               : |               :
  !!  :               : |               :
  !!  1               : 1               :
  !!  |               : :               :
  !!  |               : :               :
  !!  |               : :               :
  !!  0...............o o...............o
  !! Hanging vertex marking for three-dimensional nonconforming interface
  !!  o...............o o...............o 0-------2.......o o.......2-------0
  !!  :               : :               : |       |       : :       |       |
  !!  :               : :               : |       |       : :       |       |
  !!  :               : :               : |       |       : :       |       |
  !!  2-------1       : :       1-------2 2-------1       : :       1-------2
  !!  |       |       : :       |       | :               : :               :
  !!  |       |       : :       |       | :               : :               :
  !!  |       |       : :       |       | :               : :               :
  !!  0-------2.......o o.......2-------0 o...............o o...............o
  !! @endverbatim
  !! Vertices have no alignment. There are no interpolation operations on
  !! vertices and vertex has no dimension. Vertex is always a component part
  !! of the higher-dimension object, so position gives its location in the
  !! object.
  type, extends(polytope_actualisation_t) :: vertex_act_t
   contains
     !> Initialise a polytope actualisation
     procedure, pass(this) :: init => vertex_act_init
     !> Test equality and find alignment
     procedure, pass(this) :: equal_algn => vertex_act_equal
     !> Test alignment
     procedure, pass(this) :: test => vertex_act_test
  end type vertex_act_t

contains

  !> Initialise a topology polytope with boundary information
  !! @parameter[in]      id     polytope id
  !! @parameter[in]      nfct   number of facets
  !! @parameter[inout]   fct    polytope facets
  !! @parameter[in]      bnd    external boundary information
  subroutine vertex_tpl_init(this, id, nfct, fct, bnd)
    class(vertex_tpl_t), intent(inout) :: this
    integer(i4), intent(in) :: id, nfct, bnd
    type(topology_object_t), dimension(nfct), intent(inout) :: fct

    call this%free()

    call this%set_tdim(NEKO_VERTEX_DIM)
    call this%set_nelem(NEKO_VERTEX_NFACET, NEKO_VERTEX_NRIDGE,&
         & NEKO_VERTEX_NPEAK)
    call this%set_id(id)
    call this%set_bnd(bnd)
    ! no facets or ridges, so nothing to do
    if (nfct /= NEKO_VERTEX_NFACET) call neko_warning('Vertex has no facets.')
  end subroutine vertex_tpl_init

  !> Test equality
  !! @parameter[in]   other   polytope
  !! @return equal
  function vertex_tpl_equal(this, other) result(equal)
    class(vertex_tpl_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: equal

    ! check polygon information
    equal = this%equal_poly(other)
    if (equal) then
       ! check global id
       equal = (this%id() == other%id())
    end if
  end function vertex_tpl_equal

  !> Initialise a polytope with alignment information
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  subroutine vertex_algn_init(this, pltp, algn)
    class(vertex_algn_t), intent(inout) :: this
    class(polytope_t), target, intent(in) :: pltp
    integer(i4), intent(in) :: algn

    call this%free()

    ! There is just a single realisation of monon, so just check dimension
    if (pltp%tdim() == NEKO_VERTEX_DIM) then
       ! vertex has no alignment
       call this%init_dat(pltp, .false.)
       ! nothing more to do here
    else
       call neko_error('Vertex aligned; wrong pointer dimension.')
    end if
  end subroutine vertex_algn_init

  !> Check polytope equality and return alignment
  !! @parameter[in]    other  polytope
  !! @parameter[out]   equal  polytope equality
  !! @parameter[out]   algn   alignment information
  subroutine vertex_algn_equal(this, other, equal, algn)
    class(vertex_algn_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn

    ! check polytope information
    equal = this%polytope%equal_poly(other)
    if (equal) then
       ! check global id
       equal = (this%polytope%id() == other%id())
    end if
    ! there is no alignment
    algn = -1
  end subroutine vertex_algn_equal

  !> Test alignment
  !! @parameter[in]   other   polytope
  !! @return ifalgn
  function vertex_algn_test(this, other) result(ifalgn)
    class(vertex_algn_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: ifalgn
    ! there is no alignment
    if (this%equal(other)) then
       ifalgn = .true.
    else
       call neko_error('Alignment test: vertices not equal')
    end if
  end function vertex_algn_test

  !> Initialise a polytope with hanging information
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  !! @parameter[in]   ifint  interpolation flag
  !! @parameter[in]   hng    hanging information
  !! @parameter[in]   pos    position in the higher order element
  subroutine vertex_act_init(this, pltp, algn, ifint, hng, pos)
    class(vertex_act_t), intent(inout) :: this
    class(polytope_topology_t), target, intent(in) :: pltp
    integer(i4), intent(in) :: algn, hng, pos
    logical, intent(in) :: ifint

    call this%free()

    ! There is just a single realisation of monon, so just check dimension
    if (pltp%tdim() == NEKO_VERTEX_DIM) then
       if (hng >= 0 .and. hng <= 2) then
          ! vertex has no alignment and cannot be interpolated
          call this%init_dat(pltp, .false., .false., hng, pos)
       else
          call neko_error('Inconsistent vertex hanging information.')
       end if
       ! nothing more to do here
    else
       call neko_error('Vertex actualisation; wrong pointer dimension.')
    end if
  end subroutine vertex_act_init

  !> Check polytope equality and return alignment
  !! @parameter[in]    other  polytope
  !! @parameter[out]   equal  polytope equality
  !! @parameter[out]   algn   alignment information
  subroutine vertex_act_equal(this, other, equal, algn)
    class(vertex_act_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn

    ! check polytope information
    equal = this%polytope%equal_poly(other)
    if (equal) then
       ! check global id
       equal = (this%polytope%id() == other%id())
    end if
    ! there is no alignment
    algn = -1
  end subroutine vertex_act_equal

  !> Test alignment
  !! @parameter[in]   other   polytope
  !! @return ifalgn
  function vertex_act_test(this, other) result(ifalgn)
    class(vertex_act_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: ifalgn
    ! there is no alignment
    if (this%equal(other)) then
       ifalgn = .true.
    else
       call neko_error('Alignmentx test: vertices not equal')
    end if
  end function vertex_act_test

end module vertex
