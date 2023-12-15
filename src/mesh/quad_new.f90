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
!> Connectivity quad types
module quad_new
  use num_types, only : i4, dp
  use utils, only : neko_error
  use polytope, only : polytope_t
  use polytope_topology, only : polytope_topology_t, topology_object_t
  use polytope_actualisation, only : polytope_actualisation_t
  use polytope_mesh, only : polytope_mesh_t, mesh_object_t
  use alignment_quad, only : alignment_quad_init, alignment_quad_set_t
  use point, only : point_t, point_ptr
  implicit none
  private

  public :: quad_tpl_t, quad_act_t, quad_msh_t

  ! object information
  integer(i4), public, parameter :: NEKO_QUAD_DIM = 2
  integer(i4), public, parameter :: NEKO_QUAD_NFACET = 4
  integer(i4), public, parameter :: NEKO_QUAD_NRIDGE = 4
  integer(i4), public, parameter :: NEKO_QUAD_NPEAK = 0

  type, extends(polytope_topology_t) :: quad_tpl_t
   contains
     !> Initialise a topology polytope
     procedure, pass(this) :: init => quad_tpl_init
     !> Test equality
     procedure, pass(this) :: equal => quad_tpl_equal
  end type quad_tpl_t

  type, extends(polytope_actualisation_t) :: quad_act_t
   contains
     !> Initialise a polytope actualisation
     procedure, pass(this) :: init => quad_act_init
     !> Test equality and find alignment
     procedure, pass(this) :: equal_algn => quad_act_equal
     !> Test alignment
     procedure, pass(this) :: test => quad_act_test
  end type quad_act_t

  type, extends(polytope_mesh_t) :: quad_msh_t
   contains
     !> Initialise a topology polytope
     procedure, pass(this)  :: init => quad_msh_init
     !> Test equality
     procedure, pass(this) :: equal => quad_msh_equal
     !> Return element diameter
     procedure, pass(this) :: diameter => quad_msh_diameter
     !> Return element centroid
     procedure, pass(this) :: centroid => quad_msh_centroid
     !> Return facet @a r and @s local directions with respect to the element
     procedure, pass(this) :: fct_dir => quad_fct_dir
     !> Return ridge @a r local direction with respect to the element
     procedure, pass(this) :: rdg_dir => quad_rdg_dir
  end type quad_msh_t

contains

  !> Initialise a polytope with boundary information
  !! @parameter[in]      id     polytope id
  !! @parameter[in]      nfct   number of facets
  !! @parameter[inout]   fct    polytope facets
  !! @parameter[in]      bnd    external boundary information
  subroutine quad_tpl_init(this, id, nfct, fct, bnd)
    class(quad_tpl_t), intent(inout) :: this
    integer(i4), intent(in) :: id, nfct, bnd
    type(topology_object_t), dimension(nfct), intent(inout) :: fct
  end subroutine quad_tpl_init

  !> Test equality
  !! @parameter[in]   other   polytope
  !! @return equal
  function quad_tpl_equal(this, other) result(equal)
    class(quad_tpl_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: equal
  end function quad_tpl_equal

  !> Initialise a polytope with hanging information
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  !! @parameter[in]   ifint  interpolation flag
  !! @parameter[in]   hng    hanging information
  !! @parameter[in]   pos    position in the higher order element
  subroutine quad_act_init(this, pltp, algn, ifint, hng, pos)
    class(quad_act_t), intent(inout) :: this
    class(polytope_topology_t), target, intent(in) :: pltp
    integer(i4), intent(in) :: algn, hng, pos
    logical, intent(in) :: ifint
  end subroutine quad_act_init

  !> Get higher dimension object direction along local @a r
  !! @return dir
  function quad_act_dirr(this) result(dir)
    class(quad_act_t), intent(in) :: this
    integer(i4) :: dir
  end function quad_act_dirr

  !> Get higher dimension object direction along local @a s
  !! @return dir
  function quad_act_dirs(this) result(dir)
    class(quad_act_t), intent(in) :: this
    integer(i4) :: dir
  end function quad_act_dirs

  !> Check polytope equality and return alignment
  !! @parameter[in]    other  polytope
  !! @parameter[out]   equal  polytope equality
  !! @parameter[out]   algn   alignment information
  subroutine quad_act_equal(this, other, equal, algn)
    class(quad_act_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn
  end subroutine quad_act_equal

  !> Test alignment
  !! @parameter[in]   other   polytope
  !! @return ifalgn
  function quad_act_test(this, other) result(ifalgn)
    class(quad_act_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: ifalgn
  end function quad_act_test

  !> Initialise a polytope with geometry information
  !! @parameter[in]   id     polytope id
  !! @parameter[in]   nfct   number of facets
  !! @parameter[in]   fct    polytope facets
  !! @parameter[in]   npts   number of points
  !! @parameter[in]   pts    points
  subroutine quad_msh_init(this, id, nfct, fct, npts, pts)
    class(quad_msh_t), intent(inout) :: this
    integer(i4), intent(in) :: id, nfct, npts
    type(mesh_object_t), dimension(nfct), intent(in) :: fct
    type(mesh_object_t), dimension(npts), intent(in) :: pts
  end subroutine quad_msh_init

  !> Test equality
  !! @parameter[in]   other   polytope
  !! @return equal
  function quad_msh_equal(this, other) result(equal)
    class(quad_msh_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: equal
  end function quad_msh_equal

  !> Get element diameter
  !! @return res
  function quad_msh_diameter(this) result(res)
    class(quad_msh_t), intent(in) :: this
    real(dp) :: res
  end function quad_msh_diameter

  !> Get element centroid
  !! @return res
  function quad_msh_centroid(this) result(res)
    class(quad_msh_t), intent(in) :: this
    type(point_t) :: res
  end function quad_msh_centroid

  !> Get @a r and @a s facet local directions
  !! @parameter[in]   pos          facet position
  !! @parameter[out]  dirr, dirs   local directions
  subroutine quad_fct_dir(this, pos, dirr, dirs)
    class(quad_msh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4), intent(out) :: dirr, dirs
    dirs = -1
  end subroutine quad_fct_dir

  !> Get @a r ridge local direction
  !! @parameter[in]   pos          ridge position
  !! @parameter[out]  dirr         local direction
  subroutine quad_rdg_dir(this, pos, dirr)
    class(quad_msh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4), intent(out) :: dirr
    dirr = -1
  end subroutine quad_rdg_dir

end module quad_new
