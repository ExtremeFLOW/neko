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
!> Connectivity vertex types
module vertex_cnn
  use num_types, only : i4
  use utils, only : neko_error
  use polytope_cnn, only : polytope_cnn_t
  implicit none
  private

  public :: vertex_cab_t, vertex_cab_ptr, vertex_ncnf_cac_t, vertex_ncnf_cac_ptr

  ! object information
  integer(i4), public, parameter :: NEKO_VERTEX_DIM = 0
  integer(i4), public, parameter :: NEKO_VERTEX_NFACET = 0
  integer(i4), public, parameter :: NEKO_VERTEX_NRIDGE = 0
  integer(i4), public, parameter :: NEKO_VERTEX_NPEAK = 0

  !> Type for abstract vertex object
  !! @details Vertex is the only realisation of zero-dimensional polytope
  !! (monon) and contains unique global id only. Vertex has no alignment.
  !! Its only actualisation are components of higher-dimension objects.
  type, extends(polytope_cnn_t) :: vertex_cab_t
   contains
     !> Initialise vertex
     procedure, pass(this) :: init => vertex_init
     !> vertex equality check
     procedure, pass(this) :: equal => vertex_equal
     generic :: operator(.eq.) => equal
  end type vertex_cab_t

  !> Pointer to an abstract vertex object
  type ::  vertex_cab_ptr
     type(vertex_cab_t), pointer :: ptr
  end type vertex_cab_ptr

  !> Actualisation of the abstract vertex for nonconforming meshes
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
  !! There are no operations on vertex, so no procedure pointers.
  !! Vertex is always a component part of the higher-dimension object, so
  !! position gives its location in the object.
  type :: vertex_ncnf_cac_t
     ! vertex pointer
     type(vertex_cab_ptr) :: vertex
     !> hanging information
     integer(i4) :: hanging = -1
     !> position in the object
     integer(i4) :: position = -1
   contains
     !> Initialise vertex pointer and position
     procedure, pass(this) :: init => vertex_hanging_init
     !> Free vertex data
     procedure, pass(this) :: free => vertex_hanging_free
     !> Set hanging information
     procedure, pass(this) :: set_hng => vertex_hanging_set
     !> Get hanging information
     procedure, pass(this) :: hng => vertex_hanging_get
     !> Set position information
     procedure, pass(this) :: set_pos => vertex_position_set
     !> Get position information
     procedure, pass(this) :: pos => vertex_position_get
  end type vertex_ncnf_cac_t

  !> Pointer to a nonconforming vertex actualisation
  type ::  vertex_ncnf_cac_ptr
     type(vertex_ncnf_cac_t), pointer :: ptr
  end type vertex_ncnf_cac_ptr

contains

  !> @brief Initialise vertex with global id
  !! @parameter[in]   id     unique id
  subroutine vertex_init(this, id)
    class(vertex_cab_t), intent(inout) :: this
    integer(i4), intent(in) :: id

    call this%set_dim(NEKO_VERTEX_DIM)
    call this%set_nelem(NEKO_VERTEX_NFACET, NEKO_VERTEX_NRIDGE,&
         & NEKO_VERTEX_NPEAK)
    call this%set_id(id)
  end subroutine vertex_init

  !> @brief Check if two vertices are the same
  !! @return   equal
  pure function vertex_equal(this, other) result(equal)
    class(vertex_cab_t), intent(in) :: this
    class(polytope_cnn_t), intent(in) :: other
    logical :: equal

    ! check polygon information
    equal = this%equal_poly(other)
    if (equal) then
       ! check global id
       equal = (this%id() == other%id())
    end if
  end function vertex_equal

  !> @brief Initialise vertex pointer and position
  !! @parameter[in]   vrt     vertex
  !! @parameter[in]   pos     position in the object
  subroutine vertex_hanging_init(this, vrt, pos)
    class(vertex_ncnf_cac_t), intent(inout) :: this
    type(vertex_cab_t), target, intent(in) :: vrt
    integer(i4), intent(in) :: pos
    call this%free()
    this%vertex%ptr => vrt
    this%position = pos
    ! Assume not hanging vertex
    this%hanging = 0
  end subroutine vertex_hanging_init

  !> @brief free vertex pointer and hanging information
  subroutine vertex_hanging_free(this)
    class(vertex_ncnf_cac_t), intent(inout) :: this
    this%vertex%ptr => null()
    this%hanging = -1
    this%position = -1
  end subroutine vertex_hanging_free

  !> @brief Set hanging information
  !! @parameter[in]   hng     hanging information
  subroutine vertex_hanging_set(this, hng)
    class(vertex_ncnf_cac_t), intent(inout) :: this
    integer(i4), intent(in) :: hng
    if (hng >= 0 .and. hng <= 2) then
       this%hanging = hng
    else
       call neko_error('Inconsistent vertex hanging information.')
    end if
  end subroutine vertex_hanging_set

  !> @brief Get hanging information
  !! @return   hng
  pure function vertex_hanging_get(this) result(hng)
    class(vertex_ncnf_cac_t), intent(in) :: this
    integer(i4) :: hng
    hng = this%hanging
  end function vertex_hanging_get

  !> @brief Set position information
  !! @parameter[in]   pos     position information
  pure subroutine vertex_position_set(this, pos)
    class(vertex_ncnf_cac_t), intent(inout) :: this
    integer(i4), intent(in) :: pos
    this%position = pos
  end subroutine vertex_position_set

  !> @brief Get position information
  !! @return   pos
  pure function vertex_position_get(this) result(pos)
    class(vertex_ncnf_cac_t), intent(in) :: this
    integer(i4) :: pos
    pos = this%position
  end function vertex_position_get

end module vertex_cnn
