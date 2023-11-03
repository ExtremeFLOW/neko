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
!> Connectivity cell hexahedron type
module hex_cnn
  use num_types, only : i4
  use utils, only : neko_error
  use polytope_cnn, only : polytope_cnn_t
  use vertex_cnn, only : vertex_cnn_t, vertex_cnn_ptr
  use edge_cnn, only : edge_cnn_t, edge_cnn_ptr, edge_aligned_cnn_t,&
       & NEKO_EDGE_NFACET
  use quad_cnn, only : quad_cnn_t, quad_cnn_ptr, quad_aligned_cnn_t,&
       & NEKO_QUAD_NFACET, NEKO_QUAD_NRIDGE, quad_to_edg_algn_inv
  use cell_cnn, only : cell_cnn_t
  implicit none
  private

  public :: hex_cnn_t, hex_cnn_ptr

  ! object information
  integer(i4), public, parameter :: NEKO_HEX_NFACET = 6
  integer(i4), public, parameter :: NEKO_HEX_NRIDGE = 12
  integer(i4), public, parameter :: NEKO_HEX_NPEAK = 8

  !> Hex type for global communication
  !! @details Hex is one of the realisations of cell containing 6 facets
  !! (quads), 12 ridges (edges) and 8 peaks (vertices).
  !! Facets are oriented according to the numbers of ridges and peaks
  !! (from the peak with the smaller number towards the one with the bigger
  !! number).This corresponds to the data alignment in the face data sub-array.
  !! @verbatim
  !! Facet, ridge and peak numbering (symmetric notation)
  !!
  !!              +--------+     ^ s
  !!             /        /|     |
  !!            /  f_4   / |     |
  !!    f_1 -> /        /  |     |
  !!          +--------+f_2+     +----> r
  !!          |        |  /     /
  !!          |   f_6  | /     /
  !!          |        |/     /
  !!          +--------+      t
  !!               ^
  !!               |
  !!              f_3
  !!
  !!              +--r_2---+     ^ s
  !!             /        /|     |
  !!         r_11     r_12 r_6   |
  !!           /        /  |     |
  !!          +---r_4--+   +     +----> r
  !!          |        |  /     /
  !!         r_7     r_8 r_10  /
  !!          |        |/     /
  !!          +---r_3--+     t
  !!
  !!           p_3+--------+p_4  ^ s
  !!             /        /|     |
  !!            /        / |     |
  !!           /        /  |     |
  !!       p_7+--------+p_8+p_2  +----> r
  !!          |        |  /     /
  !!          |        | /     /
  !!          |        |/     /
  !!       p_5+--------+p_6   t
  !!
  !! @endverbatim
  type, extends(cell_cnn_t) :: hex_cnn_t
     !> Facets are aligned
     type(quad_aligned_cnn_t), dimension(:), allocatable :: facet
   contains
     !> Initialise hex
     procedure, pass(this) :: init => hex_init
     !> Free hex data
     procedure, pass(this) :: free => hex_free
     !> Is hex self-periodic
     procedure, pass(this) :: selfp => hex_self_periodic
     !> Get pointers to facets
     procedure, pass(this) :: fct => hex_facet
     !> Return faces shared by hexes
     procedure, pass(this) :: fct_share => hex_facet_share
     !> Hex equality including face, edge and vertex information
     procedure, pass(this) :: equal => hex_equal
     generic :: operator(.eq.) => equal
  end type hex_cnn_t

  !> Pointer to a hex type
  type ::  hex_cnn_ptr
     type(hex_cnn_t), pointer :: obj
  end type hex_cnn_ptr

  ! Lookup tables
  !> Facet edge to ridge
  integer, parameter, public, dimension(4, 6) :: fct_to_rdg = reshape((/ &
       & 9,11,5,7 , 10,12,6,8 , 9,10,1,3 , 11,12,2,4 , 5,6,1,2 , 7,8,3,4 /), &
       & shape(fct_to_rdg))
  !> Facet corner to peak
  integer, parameter, public, dimension(4, 6) :: fct_to_pek = reshape((/ &
       & 1,3,5,7 , 2,4,6,8 , 1,2,5,6 , 3,4,7,8 , 1,2,3,4 , 5,6,7,8 /), &
       & shape(fct_to_pek))
  !> Ridge corner to peak
  integer, parameter, public, dimension(2, 12) :: rdg_to_pek = reshape((/&
       & 1,2 , 3,4 , 5,6 , 7,8 , 1,3 , 2,4 , 5,7 , 6,8 , 1,5 , 2,6 , 3,7 ,&
       & 4,8 /), shape(rdg_to_pek))
  !> Facets connected to the ridge
  integer, parameter, public, dimension(2, 12) :: rdg_to_fct = reshape((/&
       & 3,5 , 4,5 , 3,6 , 4,6 , 1,5 , 2,5 , 1,6 , 2,6 , 1,3 , 2,3 , 1,4 ,&
       & 2,4 /), shape(rdg_to_fct))
  !> Facets connected to the peak
  integer, parameter, public, dimension(3, 8) :: pek_to_fct = reshape((/&
       & 1,3,5 , 2,3,5 , 1,4,5 , 2,4,5 , 1,3,6 , 2,3,6 , 1,4,6 , 2,4,6 /), &
       & shape(pek_to_fct))
  !> Ridges connected to the peak
  integer, parameter, public, dimension(3, 8) :: pek_to_rdg = reshape((/&
       & 1,5,9 , 1,6,10 , 2,5,11 , 2,6,12 , 3,7,9 , 3,8,10 , 4,7,11 ,&
       & 4,8,12 /), shape(pek_to_rdg))
  !> Ridge to the facet edge (-1 means ridge is not part of the facet)
  integer, parameter, dimension(6, 12) :: rdg_to_fcte = reshape((/&
       & -1,-1, 3,-1, 3,-1 , -1,-1,-1, 3, 4,-1 , -1,-1, 4,-1,-1, 3 , &
       & -1,-1,-1, 4,-1, 4 ,  3,-1,-1,-1, 1,-1 , -1, 3,-1,-1, 2,-1 , &
       &  4,-1,-1,-1,-1, 1 , -1, 4,-1,-1,-1, 2 ,  1,-1, 1,-1,-1,-1 , &
       & -1, 1, 2,-1,-1,-1 ,  2,-1,-1, 1,-1,-1 , -1, 2,-1, 2,-1,-1 /),&
       & shape(rdg_to_fcte))
  !> Peak to the facet corners (-1 means peak is not part of the facet)
  integer, parameter, dimension(6, 8) :: pek_to_fctc = reshape((/&
       &  1,-1, 1,-1, 1,-1 , -1, 1, 2,-1, 2,-1 ,  2,-1,-1, 1, 3,-1 , &
       & -1, 2,-1, 2, 4,-1 ,  3,-1, 3,-1,-1, 1 , -1, 3, 4,-1,-1, 2 , &
       &  4,-1,-1, 3,-1, 3 , -1, 4,-1, 4,-1, 4 /), shape(pek_to_fctc))

contains
  !> @brief Initialise hex with global id and six quads
  !! @details Quad order is important
  !! @parameter[in]   id                              unique id
  !! @parameter[in]   qd1, qd2, qd3, qd4, qd5, qd6    bounding quads
  subroutine hex_init(this, id, qd1, qd2, qd3, qd4, qd5, qd6)
    class(hex_cnn_t), intent(inout) :: this
    integer(i4), intent(in) :: id
    type(quad_aligned_cnn_t), intent(in), target :: qd1, qd2, qd3, qd4, qd5, qd6
    integer(i4) :: il, jl, ifct, icrn
    integer(i4), parameter :: sz = 3
    integer(i4), dimension(sz, sz) :: trans
    integer(i4), dimension(sz) :: work
    integer(i4), dimension(NEKO_QUAD_NFACET) :: mapf
    integer(i4), dimension(NEKO_QUAD_NRIDGE) :: mapr
    type(vertex_cnn_ptr), dimension(3) :: vrtp
    type(edge_cnn_ptr), dimension(2) :: edgp
    type(edge_cnn_t) :: edg
    integer(i4), dimension(2) :: edg_algn
    logical :: equal

    call this%free()

    ! init_dim calls free
    call this%init_dim()

    call this%set_nelem(NEKO_HEX_NFACET, NEKO_HEX_NRIDGE,&
         & NEKO_HEX_NPEAK)
    call this%set_id(id)
    ! get facet pointers
    allocate(this%facet(NEKO_HEX_NFACET))
    this%facet(1) = qd1
    this%facet(2) = qd2
    this%facet(3) = qd3
    this%facet(4) = qd4
    this%facet(5) = qd5
    this%facet(6) = qd6

    ! THIS SHOULD BE IN DIFFERENT PLACE; It checks internal consistency of faces
    ! Check if faces are consistent. Self-periodicity is allowed, but edges and
    ! vertices should not be messed up
    do il = 1, NEKO_HEX_NFACET - 1
       do jl = il + 1, NEKO_HEX_NFACET
          equal = this%facet(il)%face%obj.eq.this%facet(jl)%face%obj
       end do
    end do

    ! Get peak pointers checking hex structure and face orientation
    ! no special treatment of self-periodic faces
    allocate(this%peak(NEKO_HEX_NPEAK))
    do il = 1, NEKO_HEX_NPEAK
       ! find proper vertices
       do jl = 1, 3
          ifct = pek_to_fct(jl, il)
          icrn = pek_to_fctc(ifct, il)
          ! check face alignment
          ! mark vertices
          trans(1, 1) = 1
          trans(3, 1) = 2
          trans(1, 3) = 3
          trans(3, 3) = 4
          ! transformation
          call this%facet(ifct)%algn_op%trns_inv_f_i4%obj( sz, trans, work)
          ! get vertex mapping
          mapr(1) = trans(1, 1)
          mapr(2) = trans(3, 1)
          mapr(3) = trans(1, 3)
          mapr(4) = trans(3, 3)
          ! extract vertex
          vrtp(jl)%obj => this%facet(ifct)%face%obj%ridge(mapr(icrn))%vertex%obj
       end do
       ! is it a proper vertex
       if ((vrtp(1)%obj%id() == vrtp(2)%obj%id()) .and. &
            & (vrtp(1)%obj%id() == vrtp(3)%obj%id())) then
          call this%peak(il)%init(vrtp(1)%obj)
       else
          call neko_error('Inconsistent face vertices in the hex.')
       end if

    end do

    ! Get ridges checking hex structure and face orientation
    ! At the same time extract edge orientation
    ! no special treatment of self-periodic faces
    allocate(this%ridge(NEKO_HEX_NRIDGE))
    do il = 1, NEKO_HEX_NRIDGE
       ! find proper edges
       do jl = 1, 2
          ifct = rdg_to_fct(jl, il)
          icrn = rdg_to_fcte(ifct, il)
          ! check face alignment
          ! mark faces
          trans(1, 2) = 1
          trans(3, 2) = 2
          trans(2, 1) = 3
          trans(2, 3) = 4
          ! transformation
          call this%facet(ifct)%algn_op%trns_inv_f_i4%obj( sz, trans, work)
          ! get edge mapping
          mapf(1) = trans(1, 2)
          mapf(2) = trans(3, 2)
          mapf(3) = trans(2, 1)
          mapf(4) = trans(2, 3)
          ! extract edge
          edgp(jl)%obj => this%facet(ifct)%face%obj%facet(mapf(icrn))%edge%obj
          ! extract alignment
          edg_algn(jl) = quad_to_edg_algn_inv( &
               & this%facet(ifct)%face%obj%facet(mapf(icrn))%algn_op%alignment,&
               & mapf(icrn), this%facet(ifct)%algn_op%alignment)
       end do
       ! is it a proper edge
       if ((edgp(1)%obj .eq. edgp(2)%obj) .and. &
            & (edg_algn(1) == edg_algn(2))) then
          call this%ridge(il)%init(edgp(1)%obj, edg_algn(1))
          ! compare with local edge to check alignment
          call edg%init(edgp(1)%obj%id(), &
               & this%peak(rdg_to_pek(1, il))%vertex%obj, &
               & this%peak(rdg_to_pek(2, il))%vertex%obj)
          equal = this%ridge(il)%test(edg)
          if (.not.equal) &
               & call neko_error('Inconsistent edge alignment in the hex.')
       else
          call neko_error('Inconsistent face edges in the hex.')
       end if
    end do

    return
  end subroutine hex_init

  !> @brief Free hex data
  subroutine hex_free(this)
    class(hex_cnn_t), intent(inout) :: this
    !local variables
    integer(i4) :: il

    call this%set_dim(-1)
    if (allocated(this%facet)) then
       do il = 1, this%nfacet
          call this%facet(il)%free()
       end do
       deallocate(this%facet)
    end if
    if (allocated(this%ridge)) then
       do il = 1, this%nridge
          call this%ridge(il)%free()
       end do
       deallocate(this%ridge)
    end if
    if (allocated(this%peak)) then
       do il = 1, this%npeak
          call this%peak(il)%free()
       end do
       deallocate(this%peak)
    end if

    return
  end subroutine hex_free

  !> @brief Check if hex is self-periodic
  !! @return   selfp
  function hex_self_periodic(this) result(selfp)
    class(hex_cnn_t), intent(in) :: this
    logical :: selfp
    integer(i4) :: il, jl, itmp

    ! count self periodic faces
    itmp = 0
    do il = 1, this%nfacet - 1
       do jl = il + 1, this%nfacet
          selfp = this%facet(il)%face%obj .eq. this%facet(jl)%face%obj
          if (selfp) itmp = itmp + 1
       end do
    end do
    ! count self periodic edges
    do il = 1, this%nridge - 1
       do jl = il + 1, this%nridge
          selfp = this%ridge(il)%edge%obj .eq. this%ridge(jl)%edge%obj
          if (selfp) itmp = itmp + 1
       end do
    end do
    ! count self periodic vertices
    do il = 1, this%npeak - 1
       do jl = il + 1, this%npeak
          selfp = (this%peak(il)%vertex%obj%id() == &
               & this%peak(jl)%vertex%obj%id())
          if (selfp) itmp = itmp + 1
       end do
    end do
    if (itmp == 0) then
       selfp = .false.
    else
       selfp = .true.
    end if

    return
  end function hex_self_periodic

  !> @brief Return pointers to hex facets
  !! @parameter[out]  facet   facet pointers array
  subroutine hex_facet(this, facet)
    class(hex_cnn_t), intent(in) :: this
    type(quad_aligned_cnn_t), dimension(:), allocatable, intent(out) :: facet
    integer(i4) :: il

    allocate(facet(this%nfacet))
    do il = 1, this%nfacet
       facet(il) = this%facet(il)
    end do

    return
  end subroutine hex_facet

  !> @brief Return positions of facets shared by hexes
  !! @note Hexes can be self-periodic
  !! @parameter[in]   other   second hex
  !! @parameter[out]  ishare  number of shared faces
  !! @parameter[out]  facetp  integer position of shared faces
  subroutine hex_facet_share(this, other, ishare, facetp)
    class(hex_cnn_t), intent(in) :: this, other
    integer(i4), intent(out) :: ishare
    integer(i4), dimension(:, :), allocatable, intent(out) :: facetp
    integer(i4) :: il, jl

    allocate(facetp(2, this%nfacet * other%nfacet))
    ishare = 0
    facetp(:,:) = 0
    do il = 1, this%nfacet
       do jl = 1, other%nfacet
          if (this%facet(il)%face%obj.eq.other%facet(jl)%face%obj) then
             ishare = ishare + 1
             facetp(1,ishare) = il
             facetp(2,ishare) = jl
          end if
       end do
    end do

    return
  end subroutine hex_facet_share

  !> @brief Check if two hexes are the same
  !! @note No special treatment of self-periodic hexes
  !! @parameter[in]  other    second hex
  !! @return   equal
  function hex_equal(this, other) result(equal)
    class(hex_cnn_t), intent(in) :: this
    class(polytope_cnn_t), intent(in) :: other
    logical :: equal

    equal = .false.
!!    call this%eq_algn(other, equal, algn)

    return
  end function hex_equal

end module hex_cnn
