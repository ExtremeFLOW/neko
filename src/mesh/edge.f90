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
module edge
  use num_types, only : i4
  use utils, only : neko_error
  use polytope, only : polytope_t
  use polytope_aligned, only : polytope_aligned_t
  use polytope_topology, only : polytope_topology_t, topology_object_t
  use polytope_actualisation, only : polytope_actualisation_t
  use alignment_edge, only : alignment_edge_init, alignment_edge_set_t
  implicit none
  private

  public :: edge_algn_t, edge_tpl_t, edge_act_t

  ! object information
  integer(i4), public, parameter :: NEKO_EDGE_DIM = 1
  integer(i4), public, parameter :: NEKO_EDGE_NFACET = 2
  integer(i4), public, parameter :: NEKO_EDGE_NRIDGE = 0
  integer(i4), public, parameter :: NEKO_EDGE_NPEAK = 0

  type, extends(polytope_aligned_t) :: edge_algn_t
   contains
     !> Initialise an aligned polytope
     procedure, pass(this) :: init => edge_algn_init
     !> Test equality and find alignment
     procedure, pass(this) :: equal_algn => edge_algn_equal
     !> Test alignment
     procedure, pass(this) :: test  => edge_algn_test
  end type edge_algn_t

  type, extends(polytope_topology_t) :: edge_tpl_t
   contains
     !> Initialise a topology polytope
     procedure, pass(this) :: init => edge_tpl_init
     !> Test equality
     procedure, pass(this) :: equal => edge_tpl_equal
  end type edge_tpl_t

  type, extends(polytope_actualisation_t) :: edge_act_t
   contains
     !> Initialise a polytope actualisation
     procedure, pass(this) :: init => edge_act_init
     !> Return higher dimension object direction for local direction @a r
     procedure, pass(this) :: dirr => edge_act_dirr
     !> Return higher dimension object direction for local direction @a s
     procedure, pass(this) :: dirs => edge_act_dirs
     !> Test equality and find alignment
     procedure, pass(this) :: equal_algn => edge_act_equal
     !> Test alignment
     procedure, pass(this) :: test => edge_act_test
  end type edge_act_t

contains

  !> Initialise a polytope with alignment information
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  subroutine edge_algn_init(this, pltp, algn)
    class(edge_algn_t), intent(inout) :: this
    class(polytope_t), target, intent(in) :: pltp
    integer(i4), intent(in) :: algn
  end subroutine edge_algn_init

  !> Check polytope equality and return alignment
  !! @parameter[in]    pltp   polytope
  !! @parameter[out]   equal  polytope equality
  !! @parameter[out]   algn   alignment information
  subroutine edge_algn_equal(this, pltp, equal, algn)
    class(edge_algn_t), intent(in) :: this
    class(polytope_t), intent(in) :: pltp
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn
  end subroutine edge_algn_equal

  !> Test alignment
  !! @parameter[in]   pltp   polytope
  !! @return ifalgn
  function edge_algn_test(this, pltp) result(ifalgn)
    class(edge_algn_t), intent(in) :: this
    class(polytope_t), intent(in) :: pltp
    logical :: ifalgn
  end function edge_algn_test

  !> Initialise a polytope with boundary information
  !! @parameter[in]   nfct   number of facets
  !! @parameter[in]   fct    polytope facets
  !! @parameter[in]   bnd    external boundary information
  subroutine edge_tpl_init(this, nfct, fct, bnd)
    class(edge_tpl_t), intent(inout) :: this
    integer(i4), intent(in) :: nfct, bnd
    type(topology_object_t), dimension(nfct), intent(in) :: fct
  end subroutine edge_tpl_init

  !> Test equality
  !! @parameter[in]   pltp   polytope
  !! @return equal
  function edge_tpl_equal(this, pltp) result(equal)
    class(edge_tpl_t), intent(in) :: this
    class(polytope_t), intent(in) :: pltp
    logical :: equal
  end function edge_tpl_equal

  !> Initialise a polytope with hanging information
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  !! @parameter[in]   ifint  interpolation flag
  !! @parameter[in]   hng    hanging information
  !! @parameter[in]   pos    position in the higher order element
  subroutine edge_act_init(this, pltp, algn, ifint, hng, pos)
    class(edge_act_t), intent(inout) :: this
    class(polytope_topology_t), target, intent(in) :: pltp
    integer(i4), intent(in) :: algn, hng, pos
    logical :: ifint
  end subroutine edge_act_init

  !> Get higher dimension object direction along local @a r
  !! @return dir
  function edge_act_dirr(this) result(dir)
    class(edge_act_t), intent(in) :: this
    integer(i4) :: dir
  end function edge_act_dirr

  !> Get higher dimension object direction along local @a s
  !! @return dir
  function edge_act_dirs(this) result(dir)
    class(edge_act_t), intent(in) :: this
    integer(i4) :: dir
  end function edge_act_dirs

  !> Check polytope equality and return alignment
  !! @parameter[in]    pltp   polytope
  !! @parameter[out]   equal  polytope equality
  !! @parameter[out]   algn   alignment information
  subroutine edge_act_equal(this, pltp, equal, algn)
    class(edge_act_t), intent(in) :: this
    class(polytope_t), intent(in) :: pltp
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn
  end subroutine edge_act_equal

  !> Test alignment
  !! @parameter[in]   pltp   polytope
  !! @return ifalgn
  function edge_act_test(this, pltp) result(ifalgn)
    class(edge_act_t), intent(in) :: this
    class(polytope_t), intent(in) :: pltp
    logical :: ifalgn
  end function edge_act_test

end module edge
