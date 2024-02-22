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
!> Connectivity hex types
module hex_new
  use num_types, only : i4, dp
  use utils, only : neko_error
  use polytope, only : polytope_t, polytope_ptr
  use topology, only : topology_t
  use element_new, only : element_new_t, element_component_t
  use vertex, only : vertex_act_t, vertex_ornt_t
  use edge, only : edge_act_t, edge_tpl_t, NEKO_EDGE_NFACET
  use quad_new, only : NEKO_QUAD_TDIM, NEKO_QUAD_NFACET, NEKO_QUAD_NRIDGE, &
       & NEKO_QUAD_NPEAK, quad_to_edg_algn_inv, quad_act_t
  use point, only : point_t, point_ptr
  implicit none
  private

  public :: hex_elm_t, hex_elm_add

  ! object information
  integer(i4), public, parameter :: NEKO_HEX_TDIM = 3
  integer(i4), public, parameter :: NEKO_HEX_NFACET = 6
  integer(i4), public, parameter :: NEKO_HEX_NRIDGE = 12
  integer(i4), public, parameter :: NEKO_HEX_NPEAK = 8

  !> Type for hex as mesh element
  !! @details For three-dimensional meshes hexes are elements. Hex is one
  !! of the realisations of cell containing 6 facets (quads), 12 ridges (edges)
  !! and 8 peaks (vertices). Facets are oriented according to the numbers of
  !! ridges and peaks (from the peak with the smaller number towards the one
  !! with the bigger number). This corresponds to the data alignment in the face
  !! data sub-array.
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
  !! Local numbering of the geometrical points corresponds to the peak
  !! numbering.
  type, extends(element_new_t) :: hex_elm_t
   contains
     !> Initialise a topology polytope
     procedure, pass(this)  :: init => hex_elm_init
     !> Test equality
     procedure, pass(this) :: equal => hex_elm_equal
     !> Return element diameter
     procedure, pass(this) :: diameter => hex_elm_diameter
     !> Return element centroid
     procedure, pass(this) :: centroid => hex_elm_centroid
     !> Return facet @a r and @s local directions with respect to the element
     procedure, pass(this) :: fct_dir => hex_elm_fct_dir
     !> Return ridge @a r local direction with respect to the element
     procedure, pass(this) :: rdg_dir => hex_elm_rdg_dir
  end type hex_elm_t

  ! Lookup tables
  !> Facet edge to ridge
  integer, parameter, public, dimension(4, 6) :: fct_to_rdg = reshape((/ &
       & 9, 11, 5, 7  ,  10, 12, 6, 8  ,  9, 10, 1, 3  ,  11, 12, 2, 4  , &
       & 5, 6, 1, 2  ,  7, 8, 3, 4 /), shape(fct_to_rdg))
  !> Facet corner to peak
  integer, parameter, public, dimension(4, 6) :: fct_to_pek = reshape((/ &
       & 1, 3, 5, 7  ,  2, 4, 6, 8  ,  1, 2, 5, 6  ,  3, 4, 7, 8  , &
       & 1, 2, 3, 4  ,  5, 6, 7, 8 /), shape(fct_to_pek))
  !> Ridge corner to peak
  integer, parameter, public, dimension(2, 12) :: rdg_to_pek = reshape((/&
       & 1, 2  ,  3, 4  ,  5, 6  ,  7, 8  ,  1, 3  ,  2, 4  ,  5, 7  , &
       & 6, 8  ,  1, 5  ,  2, 6  ,  3, 7  ,  4, 8 /), shape(rdg_to_pek))
  !> Facets connected to the ridge
  integer, parameter, public, dimension(2, 12) :: rdg_to_fct = reshape((/&
       & 3, 5  ,  4, 5  ,  3, 6  ,  4, 6  ,  1, 5  ,  2, 5  ,  1, 6  , &
       & 2, 6  ,  1, 3  ,  2, 3  ,  1, 4  ,  2, 4 /), shape(rdg_to_fct))
  !> Facets connected to the peak
  integer, parameter, public, dimension(3, 8) :: pek_to_fct = reshape((/&
       & 1, 3, 5  ,  2, 3, 5  ,  1, 4, 5  ,  2, 4, 5  ,  1, 3, 6  , &
       & 2, 3, 6  ,  1, 4, 6  ,  2, 4, 6 /), shape(pek_to_fct))
  !> Ridges connected to the peak
  integer, parameter, public, dimension(3, 8) :: pek_to_rdg = reshape((/&
       & 1, 5, 9   ,  1, 6, 10  ,  2, 5, 11  ,  2, 6, 12  ,  3, 7, 9  , &
       & 3, 8, 10  ,  4, 7, 11  ,  4, 8, 12 /), shape(pek_to_rdg))
  !> Ridge to the facet edge (-1 means ridge is not part of the facet)
  integer, parameter, dimension(6, 12) :: rdg_to_fcte = reshape((/&
       & -1, -1,  3, -1,  3, -1   ,  -1, -1, -1,  3,  4, -1   , &
       & -1, -1,  4, -1, -1,  3   ,  -1, -1, -1,  4, -1,  4   , &
       &  3, -1, -1, -1,  1, -1   ,  -1,  3, -1, -1,  2, -1   , &
       &  4, -1, -1, -1, -1,  1   ,  -1,  4, -1, -1, -1,  2   , &
       &  1, -1,  1, -1, -1, -1   ,  -1,  1,  2, -1, -1, -1   , &
       &  2, -1, -1,  1, -1, -1   ,  -1,  2, -1,  2, -1, -1  /),&
       & shape(rdg_to_fcte))
  !> Peak to the facet corners (-1 means peak is not part of the facet)
  integer, parameter, dimension(6, 8) :: pek_to_fctc = reshape((/&
       &  1, -1,  1, -1,  1, -1   ,  -1, 1,  2, -1,  2, -1  , &
       &  2, -1, -1,  1,  3, -1   ,  -1, 2, -1,  2,  4, -1  , &
       &  3, -1,  3, -1, -1,  1   ,  -1, 3,  4, -1, -1,  2  , &
       &  4, -1, -1,  3, -1,  3   ,  -1, 4, -1,  4, -1,  4  /), &
       & shape(pek_to_fctc))
  !> Peak to the ridge corners (-1 means peak is not part of the ridge)
  integer, parameter, dimension(12, 8) :: pek_to_rdgc = reshape((/&
       &  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  &
       &  2, -1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1,  &
       & -1,  1, -1, -1,  2, -1, -1, -1, -1, -1,  1, -1,  &
       & -1,  2, -1, -1, -1,  2, -1, -1, -1, -1, -1,  1,  &
       & -1, -1,  1, -1, -1, -1,  1, -1,  2, -1, -1, -1,  &
       & -1, -1,  2, -1, -1, -1, -1,  1, -1,  2, -1, -1,  &
       & -1, -1, -1,  1, -1, -1,  2, -1, -1, -1,  2, -1,  &
       & -1, -1, -1,  2, -1, -1, -1,  2, -1, -1, -1,  2  /), &
       & shape(pek_to_rdgc))

  !> Facet to local direction mapping
  integer, public, parameter, dimension(2, 6) :: fct_to_dir = reshape((/&
       & 3, 2  ,  3, 2  ,  3, 1  ,  3, 1  ,  2, 1  ,  2, 1 &
       &/),shape(fct_to_dir))
  !> Ridge to local direction mapping
  integer, public, parameter, dimension(12) :: rdg_to_dir = (/ 1, 1, 1, 1  , &
       & 2, 2, 2, 2  ,  3, 3, 3, 3 /)

contains

  !> Create hex element (3D meshes only)
  !! @parameter[out]     elm    element hex
  !! @parameter[in]      id     hex id
  !! @parameter[in]      qd1, qd2, qd3, qd4, qd5, qd6   bounding topology quads
  !! @parameter[in]      qo1, qo2, qo3, qo4, qo5, qo6   quad orientation
  !! @parameter[in]      in1, in2, in3, in4, in5, in6   quad interpolation flag
  !! @parameter[in]      hn1, hn2, hn3, hn4, hn5, hn6   quad hanging position
  !! @parameter[in]      nrdg     number of hanging ridges
  !! @parameter[in]      rdg_hng  ridge hanging flag
  !! @parameter[in]      gdim     geometrical dimension
  !! @parameter[in]      pt1, pt2, pt3, pt4, p5, p6, pt7, pt8   vertex points
  subroutine hex_elm_add(elm, id, qd1, qd2, qd3, qd4, qd5, qd6, &
       & qo1, qo2, qo3, qo4, qo5, qo6, in1, in2, in3, in4, in5, in6, &
       & hn1, hn2, hn3, hn4, hn5, hn6, nrdg, rdg_hng, &
       & gdim, pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8)
    class(element_new_t), intent(inout), allocatable :: elm
    integer(i4), intent(in) :: id, qo1, qo2, qo3, qo4, qo5, qo6, hn1, hn2, &
         & hn3, hn4, hn5, hn6, nrdg, gdim
    logical, intent(in) :: in1, in2, in3, in4, in5, in6
    class(topology_t), intent(in) :: qd1, qd2, qd3, qd4, qd5, qd6
    integer(i4), dimension(2, 3), intent(in) :: rdg_hng
    type(point_t), intent(in), target :: pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8
    type(element_component_t), dimension(NEKO_HEX_NFACET) :: fct
    type(point_ptr), dimension(NEKO_HEX_NPEAK) :: ptsl

    if (allocated(elm)) deallocate(elm)
    allocate(hex_elm_t :: elm)
    ! facets
    allocate(quad_act_t :: fct(1)%obj)
    call fct(1)%obj%init(qd1, qo1, in1, hn1, 1)
    allocate(quad_act_t :: fct(2)%obj)
    call fct(2)%obj%init(qd2, qo2, in2, hn2, 2)
    allocate(quad_act_t :: fct(3)%obj)
    call fct(3)%obj%init(qd3, qo3, in3, hn3, 3)
    allocate(quad_act_t :: fct(4)%obj)
    call fct(4)%obj%init(qd4, qo4, in4, hn4, 4)
    allocate(quad_act_t :: fct(5)%obj)
    call fct(5)%obj%init(qd5, qo5, in5, hn5, 5)
    allocate(quad_act_t :: fct(6)%obj)
    call fct(6)%obj%init(qd6, qo6, in6, hn6, 6)
    ! points
    ptsl(1)%p => pt1
    ptsl(2)%p => pt2
    ptsl(3)%p => pt3
    ptsl(4)%p => pt4
    ptsl(5)%p => pt5
    ptsl(6)%p => pt6
    ptsl(7)%p => pt7
    ptsl(8)%p => pt8
    ! edge initialisation
    call elm%init(id, NEKO_HEX_NFACET, fct, NEKO_HEX_NPEAK, ptsl, gdim, &
         & nrdg, rdg_hng)
  end subroutine hex_elm_add

  !> Initialise a polytope with geometry information
  !! @parameter[in]   id       polytope id
  !! @parameter[in]   nfct     number of facets
  !! @parameter[in]   fct      polytope facets
  !! @parameter[in]   npts     number of points
  !! @parameter[in]   pts      points
  !! @parameter[in]   gdim     geometrical dimension
  !! @parameter[in]   nrdg     number of hanging ridges
  !! @parameter[in]   rdg_hng  ridge hanging flag
  subroutine hex_elm_init(this, id, nfct, fct, npts, pts, gdim, nrdg, &
          & rdg_hng)
    class(hex_elm_t), intent(inout) :: this
    integer(i4), intent(in) :: id, nfct, npts, gdim, nrdg
    type(element_component_t), dimension(nfct), intent(inout) :: fct
    type(point_ptr), dimension(npts), intent(in) :: pts
    integer(i4), dimension(2, 3), intent(in) :: rdg_hng
    integer(i4) :: il, jl, kl, ifct, icrn, itmp
    integer(i4), parameter :: sz = 3
    integer(i4), dimension(sz, sz) :: trans
    integer(i4), dimension(NEKO_QUAD_NFACET) :: mapf
    integer(i4), dimension(NEKO_QUAD_NRIDGE) :: mapr
    integer(i4), dimension(3, 3) :: hng_fct
    integer(i4), dimension(3) :: hng
    integer(i4), dimension(2) :: edg_algn
    logical :: ifint
    type(polytope_ptr), dimension(3) :: vrt
    integer(i4), dimension(NEKO_EDGE_NFACET, 1) :: edgt

    call this%free()

    call this%set_tdim(NEKO_HEX_TDIM)
    call this%set_ncomp(NEKO_HEX_NFACET, NEKO_HEX_NRIDGE,&
         & NEKO_HEX_NPEAK)
    call this%set_id(id)
    call this%init_base(gdim, npts)
    ! get facets
    if (nfct == NEKO_HEX_NFACET) then
       allocate (this%facet(NEKO_HEX_NFACET))
       do il = 1, nfct
          ! There are more than just a single realisation of cell, so check
          ! everything
          if (NEKO_QUAD_TDIM == fct(il)%obj%polytope%tdim() .and. &
               & fct(il)%obj%polytope%check_comp(NEKO_QUAD_NFACET, &
               & NEKO_QUAD_NRIDGE, NEKO_QUAD_NPEAK)) then
             ! facet position
             jl = fct(il)%obj%pos()
             call move_alloc(fct(il)%obj, this%facet(jl)%obj)
          else
             call neko_error('Hex mesh; wrong facet dimension')
          end if
       end do
    else
       call neko_error('Hex mesh; inconsistent number of facets.')
    end if


    ! Get peak pointers checking hex structure and face orientation
    ! no special treatment of self-periodic faces
    allocate(this%peak(NEKO_HEX_NPEAK))
    do il = 1, NEKO_HEX_NPEAK
       allocate(vertex_act_t :: this%peak(il)%obj)
       ! find proper vertices
       ! Peak is shared by 3 facets
       do jl = 1, 3
          ifct = pek_to_fct(jl, il)
          icrn = pek_to_fctc(ifct, il)
          ! check face alignment
          if (this%facet(ifct)%obj%ifalgn()) then
             ! mark vertices
             trans(1, 1) = 1
             trans(3, 1) = 2
             trans(1, 3) = 3
             trans(3, 3) = 4
             ! transformation
             call this%facet(ifct)%obj%algn_op%trns_inv_i4(trans, sz, sz, 1)
             ! get vertex mapping
             mapr(1) = trans(1, 1)
             mapr(2) = trans(3, 1)
             mapr(3) = trans(1, 3)
             mapr(4) = trans(3, 3)
             ! extract vertex
             vrt(jl)%ptr => this%facet(ifct)%obj%polytope%rdg(mapr(icrn))
          else
             vrt(jl)%ptr => this%facet(ifct)%obj%polytope%rdg(icrn)
          end if
          ! collect facet hanging information
          hng_fct(1, jl) = this%facet(ifct)%obj%hng()
          hng_fct(2, jl) = ifct
          hng_fct(3, jl) = icrn
       end do
       ! is it a proper vertex
       if ((vrt(1)%ptr%id() == vrt(2)%ptr%id()) .and. &
            & (vrt(1)%ptr%id() == vrt(3)%ptr%id())) then
          ! Vertex hanging flag depends on the quad hanging flag and the
          ! vertex position in the quad (corner).
          ! Peak is shared by 3 facets
          do jl = 1, 3
             select case (hng_fct(1, jl))
             case (0) ! independent face
                hng(jl) = 0 ! independent vertex
             case (1) ! face (3D) hanging; lower left part
                select case(hng_fct(3, jl)) ! vertex position in the quad
                case (1)
                   hng(jl) = 0 ! independent vertex
                case (2, 3)
                   hng(jl) = 2 ! ridge hanging vertex
                case (4)
                   hng(jl) = 1 ! face hanging vertex
                case default
                   call neko_error('Hex mesh; wrong corner position')
                end select
             case (2) ! face (3D) hanging; lower right part
                select case(hng_fct(3, jl)) ! vertex position in the quad
                case (2)
                   hng(jl) = 0 ! independent vertex
                case (1, 4)
                   hng(jl) = 2 ! ridge hanging vertex
                case (3)
                   hng(jl) = 1 ! face hanging vertex
                case default
                   call neko_error('Hex mesh; wrong corner position')
                end select
             case (3) ! face (3D) hanging; upper left part
                select case(hng_fct(3, jl)) ! vertex position in the quad
                case (3)
                   hng(jl) = 0 ! independent vertex
                case (1, 4)
                   hng(jl) = 2 ! ridge hanging vertex
                case (2)
                   hng(jl) = 1 ! face hanging vertex
                case default
                   call neko_error('Hex mesh; wrong corner position')
                end select
             case (4) ! face (3D) hanging; upper right part
                select case(hng_fct(3, jl)) ! vertex position in the quad
                case (4)
                   hng(jl) = 0 ! independent vertex
                case (2, 3)
                   hng(jl) = 2 ! ridge hanging vertex
                case (1)
                   hng(jl) = 1 ! face hanging vertex
                case default
                   call neko_error('Hex mesh; wrong corner position')
                end select
             case default
                call neko_error('Hex mesh; facet hanging flag')
             end select
          end do
          ! A single vertex can be at the same time marked as independent and
          ! face hanging. Take the max.
          itmp = maxval(hng)
          ! For 3D meshes edges and vertices can be hanging independently of
          ! faces. This happens for vertices marked as independent only.
          ! Check rdg_hng.
          if (itmp == 0 .and. nrdg > 0) then
             edg_loop1 : do kl = 1, nrdg
                ! Peak is shared by 3 ridges
                do jl = 1, 3
                   ifct = pek_to_rdg(jl, il)
                   if (rdg_hng(1, kl) == ifct) then
                      icrn = pek_to_rdgc(ifct, il)
                      select case (rdg_hng(2, kl))
                      case (0) ! independent edge
                         itmp = 0 ! independent vertex
                      case (1) ! edge (3D) hanging; lower part
                         select case(icrn) ! vertex position in the edge
                         case (1)
                            itmp = 0 ! independent vertex
                         case (2)
                            itmp = 2 ! edge hanging vertex
                         case default
                            call neko_error('Hex mesh; wrong corner position')
                         end select
                      case (2) ! edge (2D) hanging; upper part
                         select case(icrn) ! vertex position in the edge
                         case (1)
                            itmp = 2 ! edge hanging vertex
                         case (2)
                            itmp = 0 ! independent vertex
                         case default
                            call neko_error('Hex mesh; wrong corner position')
                         end select
                      case default
                         call neko_error('Hex mesh; edge hanging flag')
                      end select
                      exit edg_loop1
                   end if
                end do
             end do edg_loop1
          end if
          ! create vertex actualisation
          ! Vertex has neither alignment, nor interpolation
          call this%peak(il)%obj%init(vrt(1)%ptr, -1, .false., itmp, il)
       else
          call neko_error('Inconsistent face vertices in the hex.')
       end if
    end do


    ! Get ridges checking hex structure and face orientation
    ! At the same time extract edge orientation
    ! no special treatment of self-periodic faces
    allocate(this%ridge(NEKO_HEX_NRIDGE))
    do il = 1, NEKO_HEX_NRIDGE
       allocate(edge_act_t :: this%ridge(il)%obj)
       ! find proper edges and collect facet hanging information
       ! Ridge is shared by 2 facets
       do jl = 1, 2
          ifct = rdg_to_fct(jl, il)
          icrn = rdg_to_fcte(ifct, il)
          ! check face alignment
          if (this%facet(ifct)%obj%ifalgn()) then
             ! mark faces
             trans(1, 2) = 1
             trans(3, 2) = 2
             trans(2, 1) = 3
             trans(2, 3) = 4
             ! transformation
             call this%facet(ifct)%obj%algn_op%trns_inv_i4(trans, sz, sz, 1)
             ! get edge mapping
             mapf(1) = trans(1, 2)
             mapf(2) = trans(3, 2)
             mapf(3) = trans(2, 1)
             mapf(4) = trans(2, 3)
             ! extract edge
             vrt(jl)%ptr => this%facet(ifct)%obj%polytope%fct(mapf(icrn))
             ! extract alignment
             edg_algn(jl) = quad_to_edg_algn_inv( &
                  & this%facet(ifct)%obj%polytope%fct_algn(mapf(icrn)), &
                  & mapf(icrn), this%facet(ifct)%obj%algn_op%algn())
          else
             vrt(jl)%ptr => this%facet(ifct)%obj%polytope%fct(icrn)
             ! extract alignment
             edg_algn(jl) = quad_to_edg_algn_inv( &
                  & this%facet(ifct)%obj%polytope%fct_algn(icrn), &
                  & icrn, this%facet(ifct)%obj%algn_op%algn())
          end if
          ! collect facet hanging information
          hng_fct(1, jl) = this%facet(ifct)%obj%hng()
          hng_fct(2, jl) = ifct
          hng_fct(3, jl) = icrn
       end do
       ! is it a proper edge
       ! Quad structure was already checked, so it is enough to check id and
       ! alignment only
       if ((vrt(1)%ptr%id() == vrt(2)%ptr%id()) .and. &
            & (edg_algn(1) == edg_algn(2))) then
          ! Edge hanging flag depends on the quad hanging flag and the
          ! edge position in the quad.
          ! Ridge is shared by 2 facets
          do jl = 1, 2
             select case (hng_fct(1, jl))
             case (0) ! independent face
                hng(jl) = 0 ! independent edge
             case (1) ! face (3D) hanging; lower left part
                select case(hng_fct(3, jl)) ! edge position in the quad
                case (1, 3)
                   hng(jl) = 3 ! ridge hanging edge; lower part
                case (2, 4)
                   hng(jl) = 5 ! face hanging edge
                case default
                   call neko_error('Hex mesh; wrong edge position')
                end select
             case (2) ! face (3D) hanging; lower right part
                select case(hng_fct(3, jl)) ! edge position in the quad
                case (1, 4)
                   hng(jl) = 5 ! face hanging edge
                case (2)
                   hng(jl) = 3 ! ridge hanging edge; lower part
                case (3)
                   hng(jl) = 4 ! ridge hanging edge; upper part
                case default
                   call neko_error('Hex mesh; wrong edge position')
                end select
             case (3) ! face (3D) hanging; upper left part
                select case(hng_fct(3, jl)) ! edge position in the quad
                case (1)
                   hng(jl) = 4 ! ridge hanging edge; upper part
                case (2, 3)
                   hng(jl) = 5 ! face hanging edge
                case (4)
                   hng(jl) = 3 ! ridge hanging edge; lower part
                case default
                   call neko_error('Hex mesh; wrong edge position')
                end select
             case (4) ! face (3D) hanging; upper right part
                select case(hng_fct(3, jl)) ! edge position in the quad
                case (1, 3)
                   hng(jl) = 5 ! face hanging edge
                case (2, 4)
                   hng(jl) = 4 ! ridge hanging edge; upper part
                case default
                   call neko_error('Hex mesh; wrong edge position')
                end select
             case default
                call neko_error('Hex mesh; facet hanging flag')
             end select
          end do
          ! A single edge can be at the same time marked as independent and
          ! face hanging. Take the max.
          itmp = maxval(hng(1:2))
          ! For 3D meshes edges and vertices can be hanging independently of
          ! faces. This happens for edges marked as independent only.
          ! Check rdg_hng.
          if (itmp == 0 .and. nrdg > 0) then
             edg_loop2 : do kl = 1, nrdg
                ! Ridge is shared by 2 facets
                do jl = 1, 2
                   ifct = pek_to_rdg(jl, il)
                   if (rdg_hng(1, kl) == ifct) then
                      if (rdg_hng(2, kl) > 0 .and. rdg_hng(2, kl) < 3) then
                         itmp = rdg_hng(2, kl)
                      else
                         call neko_error('Hex mesh; wrong edge hanging flag')
                      end if
                      exit edg_loop2
                   end if
                end do
             end do edg_loop2
          end if
          ! hanging edges require interpolation
          ifint = itmp > 0
          ! create edge actualisation
          call this%ridge(il)%obj%init(vrt(1)%ptr, edg_algn(1), ifint, itmp, &
               & il)
       else
          call neko_error('Inconsistent face edges in the hex.')
       end if
       ! compare with local edge to check alignment consistency
       ! NOT SURE THIS IS NEEDED AT ALL
       edgt(1, 1) = this%peak(rdg_to_pek(1, il))%obj%polytope%id()
       edgt(2, 1) = this%peak(rdg_to_pek(2, il))%obj%polytope%id()
       call this%ridge(il)%obj%algn_op%trns_inv_i4(edgt, NEKO_EDGE_NFACET, 1, 1)
       vrt(1)%ptr => this%ridge(il)%obj%polytope%fct(1)
       vrt(2)%ptr => this%ridge(il)%obj%polytope%fct(2)
       ifint = (vrt(1)%ptr%id() == edgt(1, 1)) .and. &
            & (vrt(2)%ptr%id() == edgt(2, 1))
       if (.not. ifint) &
            & call neko_error('Inconsistent edge alignment in the hex.')
    end do

    ! Add geometrical points
    ! Sanity check
    if (npts /= NEKO_HEX_NPEAK) &
         &call neko_error('Inconsistent point number in the mesh hex.')
    allocate(this%pts(NEKO_HEX_NPEAK))
    do il = 1, npts
       this%pts(il)%p => pts(il)%p
    end do

  end subroutine hex_elm_init

  !> Test equality
  !! @parameter[in]   other   polytope
  !! @return equal
  function hex_elm_equal(this, other) result(equal)
    class(hex_elm_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: equal
    integer(i4) :: il
    integer(i4), dimension(NEKO_HEX_NFACET) :: fctt, fcto, indf
    integer(i4), dimension(NEKO_HEX_NRIDGE) :: rdgt, rdgo, indr
    integer(i4), dimension(NEKO_HEX_NPEAK) :: pekt, peko, indp

    ! check polygon information
    equal = this%equal_poly(other)
    if (equal) then
       ! check global id
       equal = (this%id() == other%id())
       ! check edges and vertices getting discarding orientation
       ! (may not work for self-periodic)
       if (equal) then
          select type(other)
          type is (hex_elm_t)
             ! geometrical dimension
             equal = (this%gdim() == other%gdim())
             if (equal) then
                ! This is the highest dimension object for 2D meshes, so we do
                ! not consider any transformations. Simply take facets/ridges/
                ! peaks and sort them. in a similar way we treat points.
                ! We simply assume the objects built of the same components
                ! should be the same.
                ! Check facets
                do il = 1, NEKO_HEX_NFACET
                   fctt(il) = this%facet(il)%obj%polytope%id()
                   fcto(il) = other%facet(il)%obj%polytope%id()
                end do
                call sorti4(fctt, indf, NEKO_HEX_NFACET)
                call sorti4(fcto, indf, NEKO_HEX_NFACET)
                do il = 1, NEKO_HEX_NFACET
                   equal = fctt(il) == fcto(il)
                   if (.not.equal) exit
                end do
                if (equal) then
                   ! Check ridges
                   do il = 1, NEKO_HEX_NRIDGE
                      rdgt(il) = this%ridge(il)%obj%polytope%id()
                      rdgo(il) = other%ridge(il)%obj%polytope%id()
                   end do
                   call sorti4(rdgt, indr, NEKO_HEX_NRIDGE)
                   call sorti4(rdgo, indr, NEKO_HEX_NRIDGE)
                   do il = 1, NEKO_HEX_NRIDGE
                      equal = rdgt(il) == rdgo(il)
                      if (.not.equal) exit
                   end do
                   if (equal) then
                      ! Check peaks
                      do il = 1, NEKO_HEX_NPEAK
                         pekt(il) = this%peak(il)%obj%polytope%id()
                         peko(il) = other%peak(il)%obj%polytope%id()
                      end do
                      call sorti4(pekt, indp, NEKO_HEX_NPEAK)
                      call sorti4(peko, indp, NEKO_HEX_NPEAK)
                      do il = 1, NEKO_HEX_NPEAK
                         equal = pekt(il) == peko(il)
                         if (.not.equal) exit
                      end do
                      if (equal) then
                         ! Check points (the same number as peaks)
                         do il = 1, NEKO_HEX_NPEAK
                            pekt(il) = this%pts(il)%p%id()
                            peko(il) = other%pts(il)%p%id()
                         end do
                         call sorti4(pekt, indp, NEKO_HEX_NPEAK)
                         call sorti4(peko, indp, NEKO_HEX_NPEAK)
                         do il = 1, NEKO_HEX_NPEAK
                            equal = pekt(il) == peko(il)
                            if (.not.equal) exit
                         end do
                      endif
                   end if
                end if
             end if
          class default
             equal = .false.
          end select
          if (.not. equal) then
             ! Something wrong; edge with the same global id should have
             ! the same type and the same facets/ridges
             call neko_error('Mismatch in class or element global id; hex')
          end if
       end if
    end if

  end function hex_elm_equal

  !> Get element diameter
  !! @return res
  function hex_elm_diameter(this) result(res)
    class(hex_elm_t), intent(in) :: this
    real(dp) :: res
    real(dp) :: d1, d2, d3, d4
    integer(i4) :: il

    d1 = 0d0
    d2 = 0d0
    d3 = 0d0
    d4 = 0d0
    do il = 1, this%gdim()
       d1 = d1 + (this%pts(8)%p%x(il) - this%pts(1)%p%x(il))**2
       d2 = d2 + (this%pts(7)%p%x(il) - this%pts(2)%p%x(il))**2
       d3 = d3 + (this%pts(6)%p%x(il) - this%pts(3)%p%x(il))**2
       d4 = d4 + (this%pts(5)%p%x(il) - this%pts(4)%p%x(il))**2
    end do

    res = sqrt(max(d1, d2, d3, d4))
  end function hex_elm_diameter

  !> Get element centroid
  !! @return res
  function hex_elm_centroid(this) result(res)
    class(hex_elm_t), intent(in) :: this
    type(point_t) :: res
    integer(i4) :: il

    res%x(:) = 0d0
    do il = 1, this%gdim()
       res%x(il) = 0.125d0 * (this%pts(1)%p%x(il) + this%pts(2)%p%x(il) + &
            & this%pts(3)%p%x(il) + this%pts(4)%p%x(il) + this%pts(5)%p%x(il) &
            & + this%pts(6)%p%x(il) + this%pts(7)%p%x(il) + this%pts(8)%p%x(il))
    end do
  end function hex_elm_centroid

  !> Get @a r and @a s facet local directions
  !! @parameter[in]   pos          facet position
  !! @parameter[out]  dirr, dirs   local directions
  subroutine hex_elm_fct_dir(this, pos, dirr, dirs)
    class(hex_elm_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4), intent(out) :: dirr, dirs
    dirr = fct_to_dir(1, pos)
    dirs = fct_to_dir(2, pos)
  end subroutine hex_elm_fct_dir

  !> Get @a r ridge local direction
  !! @parameter[in]   pos          ridge position
  !! @parameter[out]  dirr         local direction
  subroutine hex_elm_rdg_dir(this, pos, dirr)
    class(hex_elm_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4), intent(out) :: dirr
    dirr = rdg_to_dir(pos)
  end subroutine hex_elm_rdg_dir

  ! THIS SHOULD BE LOCATED IN OTHER PLACE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Heap Sort for single integer arrays
  !! @details Following p 231 Num. Rec., 1st Ed.
  !! @param[inout]   a     vector to be sorted
  !! @param[out]     ind   permutation array
  !! @param[in]      n     array size
  subroutine sorti4(a, ind, n)
    integer, intent(in) :: n
    integer(i4), intent(inout) :: a(n)
    integer, intent(out) :: ind(n)
    integer(i4) :: aa
    integer :: j, ir, i, ii, l

    do j = 1, n
       ind(j) = j
    end do

    if (n.le.1) return

    l = n/2+1
    ir = n
    do while (.true.)
       if (l.gt.1) then
          l=l-1
          aa  = a  (l)
          ii  = ind(l)
       else
          aa =   a(ir)
          ii = ind(ir)
          a(ir) =   a( 1)
          ind(ir) = ind( 1)
          ir=ir-1
          if (ir.eq.1) then
             a(1) = aa
             ind(1) = ii
             return
          endif
       endif
       i=l
       j=l+l
       do while (j .le. ir)
          if (j.lt.ir) then
             if ( a(j).lt.a(j+1) ) j=j+1
          endif
          if (aa.lt.a(j)) then
             a(i) = a(j)
             ind(i) = ind(j)
             i=j
             j=j+j
          else
             j=ir+1
          endif
       end do
       a(i) = aa
       ind(i) = ii
    end do
  end subroutine sorti4

end module hex_new
