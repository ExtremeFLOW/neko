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
!> Connectivity hex types
module hex_new
  use num_types, only : i4, dp
  use utils, only : neko_error
  use polytope, only : polytope_t
  use polytope_mesh, only : polytope_mesh_t, mesh_object_t
  use point, only : point_t, point_ptr
  implicit none
  private

  public :: hex_msh_t

  ! object information
  integer(i4), public, parameter :: NEKO_HEX_DIM = 3
  integer(i4), public, parameter :: NEKO_HEX_NFACET = 6
  integer(i4), public, parameter :: NEKO_HEX_NRIDGE = 12
  integer(i4), public, parameter :: NEKO_HEX_NPEAK = 8

  type, extends(polytope_mesh_t) :: hex_msh_t
   contains
     !> Initialise a topology polytope
     procedure, pass(this)  :: init => hex_msh_init
     !> Test equality
     procedure, pass(this) :: equal => hex_msh_equal
     !> Return element diameter
     procedure, pass(this) :: diameter => hex_msh_diameter
     !> Return element centroid
     procedure, pass(this) :: centroid => hex_msh_centroid
     !> Return facet @a r and @s local directions with respect to the element
     procedure, pass(this) :: fct_dir => hex_fct_dir
     !> Return ridge @a r local direction with respect to the element
     procedure, pass(this) :: rdg_dir => hex_rdg_dir
  end type hex_msh_t

contains

  !> Initialise a polytope with geometry information
  !! @parameter[in]   id     polytope id
  !! @parameter[in]   nfct   number of facets
  !! @parameter[in]   fct    polytope facets
  !! @parameter[in]   npts   number of points
  !! @parameter[in]   pts    points
  subroutine hex_msh_init(this, id, nfct, fct, npts, pts)
    class(hex_msh_t), intent(inout) :: this
    integer(i4), intent(in) :: id, nfct, npts
    type(mesh_object_t), dimension(nfct), intent(in) :: fct
    type(mesh_object_t), dimension(npts), intent(in) :: pts
  end subroutine hex_msh_init

  !> Test equality
  !! @parameter[in]   other   polytope
  !! @return equal
  function hex_msh_equal(this, other) result(equal)
    class(hex_msh_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: equal
  end function hex_msh_equal

  !> Get element diameter
  !! @return res
  function hex_msh_diameter(this) result(res)
    class(hex_msh_t), intent(in) :: this
    real(dp) :: res
  end function hex_msh_diameter

  !> Get element centroid
  !! @return res
  function hex_msh_centroid(this) result(res)
    class(hex_msh_t), intent(in) :: this
    type(point_t) :: res
  end function hex_msh_centroid

  !> Get @a r and @a s facet local directions
  !! @parameter[in]   pos          facet position
  !! @parameter[out]  dirr, dirs   local directions
  subroutine hex_fct_dir(this, pos, dirr, dirs)
    class(hex_msh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4), intent(out) :: dirr, dirs
  end subroutine hex_fct_dir

  !> Get @a r ridge local direction
  !! @parameter[in]   pos          ridge position
  !! @parameter[out]  dirr         local direction
  subroutine hex_rdg_dir(this, pos, dirr)
    class(hex_msh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4), intent(out) :: dirr
  end subroutine hex_rdg_dir

end module hex_new
