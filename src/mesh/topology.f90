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
!> Abstract types for abstract polytope class for mesh topology
!! @details There are three classes introduced: @ref polytope_aligned_t,
!! @ref polytope_oriented_t and @ref topology_t. The first one combines polytope
!! (topology one) with alignment information. The second is similar to the firs
!! one and is introduced to define deferred initialisation interface. This one
!! is different from the initialisation interface defined in
!! @ref polytope_actualisation_t descending from @ref polytope_alignment_t
!! as well. The last one stores information of the abstract (an unique in the
!! mesh) topology objects like vertices, edges and faces.
module topology
  use num_types, only : i4, i8
  use utils, only : neko_error
  use polytope, only : polytope_t
  use alignment, only : alignment_t
  implicit none
  private

  public :: polytope_aligned_t, polytope_oriented_t
  public :: topology_component_t, topology_t, topology_ptr, topology_element_t

  !> Base type for an abstract aligned polytope
  !! @details This is an abstract type combining @ref polytope_t and an
  !! alignment operator. Note that vertices have no alignment. This type
  !! corresponds to building blocs of the higher dimension abstract objects
  !! the topology or mesh consists of.
  type, abstract :: polytope_aligned_t
     !> Polytope pointer
     class(topology_t), pointer :: polytope => null()
     !> Is the object aligned
     logical, private :: ifaligned_ = .false.
     !> Alignment operator
     class(alignment_t), allocatable :: algn_op
   contains
     !> Free polytope and aligned data
     procedure, pass(this) :: free_base => polytope_aligned_free_base
     !> Initialise alignment data
     procedure, pass(this) :: init_base => polytope_aligned_init_base
     !> Return a pointer to the polytope
     procedure, pass(this) :: polyp => polytope_aligned_ptr
     !> Return alignment flag
     procedure, pass(this) :: ifalgn => polytope_aligned_ifalgn_get
     !> Return alignment value
     procedure, pass(this) :: algn => polytope_aligned_algn_get
     !> Test equality
     procedure, pass(this) :: equal => polytope_aligned_equal
     generic :: operator(.eq.) => equal
     !> Free type
     procedure(polytope_aligned_free), pass(this), deferred :: free
     !> Test equality and find alignment
     procedure(polytope_aligned_equal_algn), pass(this), deferred :: equal_algn
     !> Test alignment
     procedure(polytope_aligned_test), pass(this), deferred :: test
  end type polytope_aligned_t

  !> Base type for an abstract oriented polytope
  !! @details This is an abstract type basically equal to
  !! @ref polytope_aligned_t This type corresponds to building blocs of
  !! the higher dimension abstract objects the mesh topology consists of.
  !! @note The reason we need this type is to define a deferred
  !! constructor. We cannot define it in the base type because the type
  !! @ref polytope_actualisation_t also descends from it  but requires a
  !! constructor with a different signature.
  type, extends(polytope_aligned_t), abstract :: polytope_oriented_t
   contains
     !> Free aligned polytope
     procedure, pass(this) :: free => polytope_oriented_free
     !> Initialise an aligned polytope
     procedure(polytope_oriented_init), pass(this), deferred :: init
  end type polytope_oriented_t

  !> Single topology object allocatable space
  type :: topology_component_t
     class(polytope_oriented_t), allocatable :: obj
  end type topology_component_t

  !> Base type for an abstract topology polytope
  !! @details This is an abstract type combining a set of lower dimension
  !! polytopes and boundary condition into the topology object.
  type, extends(polytope_t), abstract :: topology_t
     !> Polytope facets
     type(topology_component_t), dimension(:), allocatable :: facet
     !> Polytope ridges
     type(topology_component_t), dimension(:), allocatable :: ridge
     !> Internal/external boundary condition flag
     integer(i4), private :: boundary_ = -1
     !> Polytope global id used for numbering of dof for gather-scatter
     integer(i8), private :: gsid_ = -1
   contains
     !> Free polytope data
     procedure, pass(this) :: free => topology_free
     !> Return a pointer to the polytope facets
     procedure, pass(this) :: fct => topology_fct_ptr
     !> Return a pointer to the polytope ridges
     procedure, pass(this) :: rdg => topology_rdg_ptr
     !> Return boundary value
     procedure, pass(this) :: bnd => topology_bnd_get
     !> Set boundary value
     procedure, pass(this) :: set_bnd => topology_bnd_set
     !> Return communication global id
     procedure, pass(this) :: gsid => topology_gsid_get
     !> Set communication global id
     procedure, pass(this) :: set_gsid => topology_gsid_set
     !> Is polytope self-periodic?
     procedure, pass(this) :: self_periodic => topology_self_periodic
     !> Return facets shared by polytopes
     procedure, pass(this) :: fct_share => topology_facet_share
     !> Return ridges shared by polytopes
     procedure, pass(this) :: rdg_share => topology_ridge_share
     !> Return facet alignment
     procedure, pass(this) :: fct_algn => topology_fct_algn
     !> Initialise a topology polytope
     procedure(topology_init), pass(this), deferred :: init
  end type topology_t

  !> The general pointer type to the topology polytope class
  type :: topology_ptr
     class(topology_t), pointer :: ptr
  end type topology_ptr

  !> Single topology element allocatable space
  type :: topology_element_t
     class(topology_t), allocatable :: obj
  end type topology_element_t

  !> Free type
  abstract interface
     subroutine polytope_aligned_free(this)
       import polytope_aligned_t
       class(polytope_aligned_t), intent(inout) :: this
     end subroutine polytope_aligned_free
  end interface

  !> Check polytope equality and return alignment
  !! @parameter[in]    other  polytope
  !! @parameter[out]   equal  polytope equality
  !! @parameter[out]   algn   alignment information
  abstract interface
     subroutine polytope_aligned_equal_algn(this, other, equal, algn)
       import i4
       import polytope_t
       import polytope_aligned_t
       class(polytope_aligned_t), intent(in) :: this
       class(polytope_t), intent(in) :: other
       logical, intent(out) :: equal
       integer(i4), intent(out) :: algn
     end subroutine polytope_aligned_equal_algn
  end interface

  !> Test alignment
  !! @parameter[in]   other   polytope
  !! @return ifalgn
  abstract interface
     function polytope_aligned_test(this, other) result(ifalgn)
       import topology_t
       import polytope_aligned_t
       class(polytope_aligned_t), intent(in) :: this
       class(topology_t), intent(in) :: other
       logical :: ifalgn
     end function polytope_aligned_test
  end interface

  !> Initialise a polytope with alignment information for topology
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  abstract interface
     subroutine polytope_oriented_init(this, pltp, algn)
       import i4
       import topology_t
       import polytope_oriented_t
       class(polytope_oriented_t), intent(inout) :: this
       class(topology_t), target, intent(in) :: pltp
       integer(i4), intent(in) :: algn
     end subroutine polytope_oriented_init
  end interface

  !> Abstract interface to initialise a polytope with boundary information
  !! @parameter[in]      id     polytope id
  !! @parameter[in]      nfct   number of facets
  !! @parameter[inout]   fct    polytope facets
  !! @parameter[in]      bnd    external boundary information
  abstract interface
     subroutine topology_init(this, id, nfct, fct, bnd)
       import i4
       import topology_t
       import topology_component_t
       class(topology_t), intent(inout) :: this
       integer(i4), intent(in) :: id, nfct, bnd
       type(topology_component_t), dimension(nfct), intent(inout) :: fct
     end subroutine topology_init
  end interface

contains

  !> Free polytope and aligned data
  subroutine polytope_aligned_free_base(this)
    class(polytope_aligned_t), intent(inout) :: this
    this%polytope => null()
    this%ifaligned_ = .false.
    if (allocated(this%algn_op)) deallocate(this%algn_op)
  end subroutine polytope_aligned_free_base

  !> Initialise general data
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   ifalgn if non identity alignment
  subroutine polytope_aligned_init_base(this, pltp, ifalgn)
    class(polytope_aligned_t), intent(inout) :: this
    class(topology_t), target, intent(in) :: pltp
    logical, intent(in)  :: ifalgn

    this%polytope => pltp
    this%ifaligned_ = ifalgn
  end subroutine polytope_aligned_init_base

  !> @brief Return pointer to the polytope
  !! @return  poly
  function polytope_aligned_ptr(this) result(poly)
    class(polytope_aligned_t), intent(in) :: this
    class(topology_t), pointer :: poly
    poly => this%polytope
  end function polytope_aligned_ptr

  !> @brief Get polytope alignment flag
  !! @return   ifalgn
  pure function polytope_aligned_ifalgn_get(this) result(ifalgn)
    class(polytope_aligned_t), intent(in) :: this
    logical :: ifalgn
    ifalgn = this%ifaligned_
  end function polytope_aligned_ifalgn_get

  !> @brief Get polytope alignment
  !! @return   algn
  pure function polytope_aligned_algn_get(this) result(algn)
    class(polytope_aligned_t), intent(in) :: this
    integer(i4) :: algn
    if (allocated(this%algn_op)) then
       algn = this%algn_op%algn()
    else
       algn = -1
    end if
  end function polytope_aligned_algn_get

  !> Test equality
  !! @parameter[in]   other   polytope
  !! @return equal
  function polytope_aligned_equal(this, other) result(equal)
    class(polytope_aligned_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: equal
    integer(i4) :: itmp
    call this%equal_algn(other, equal, itmp)
  end function polytope_aligned_equal

  !> Free oriented polytope
  subroutine polytope_oriented_free(this)
    class(polytope_oriented_t), intent(inout) :: this

    call this%free_base()
  end subroutine polytope_oriented_free

  !> Free polytope data
  subroutine topology_free(this)
    class(topology_t), intent(inout) :: this
    integer(i4) :: il
    this%boundary_ = -1
    this%gsid_ = -1
    if (allocated(this%facet)) then
       do il = 1, this%nfacet
          call this%facet(il)%obj%free()
          deallocate(this%facet(il)%obj)
       end do
       deallocate(this%facet)
    end if
    if (allocated(this%ridge)) then
       do il = 1, this%nridge
          call this%ridge(il)%obj%free()
          deallocate(this%ridge(il)%obj)
       end do
       deallocate(this%ridge)
    end if
    call this%free_base()
  end subroutine topology_free

  !> @brief Return pointer to the polytope facet
  !! @parameter[in]   pos   polytope component position
  !! @return ptr
  function topology_fct_ptr(this, pos) result(ptr)
    class(topology_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    class(topology_t), pointer :: ptr
    if ((pos > 0) .and. (pos <= this%nfacet)) then
       ptr => this%facet(pos)%obj%polyp()
    else
       call neko_error('Wrong facet number for topology objects.')
    end if
  end function topology_fct_ptr

  !> @brief Return pointer to the polytope ridge
  !! @parameter[in]   pos   polytope component position
  !! @return ptr
  function topology_rdg_ptr(this, pos) result(ptr)
    class(topology_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    class(topology_t), pointer :: ptr
    if ((pos > 0) .and. (pos <= this%nridge)) then
       ptr => this%ridge(pos)%obj%polyp()
    else
       call neko_error('Wrong ridge number for topology objects.')
    end if
  end function topology_rdg_ptr

  !> @brief Get polytope boundary information
  !! @return   intp
  pure function topology_bnd_get(this) result(intp)
    class(topology_t), intent(in) :: this
    integer(i4) :: intp
    intp = this%boundary_
  end function topology_bnd_get

  !> @brief Set boundary value
  !! @parameter[in]   bnd     boundary information
  subroutine topology_bnd_set(this, bnd)
    class(topology_t), intent(inout) :: this
    integer(i4), intent(in) :: bnd
    this%boundary_ = bnd
  end subroutine topology_bnd_set

  !> @brief Get communication global id
  !! @return   intp
  pure function topology_gsid_get(this) result(intp)
    class(topology_t), intent(in) :: this
    integer(i8) :: intp
    intp = this%gsid_
  end function topology_gsid_get

  !> @brief Set communication global id
  !! @parameter[in]   gsid     polytope communication global id
  subroutine topology_gsid_set(this, gsid)
    class(topology_t), intent(inout) :: this
    integer(i8), intent(in) :: gsid
    this%gsid_ = gsid
  end subroutine topology_gsid_set

  !> @brief Check if polytope is self-periodic
  !! @return   selfp
  function topology_self_periodic(this) result(selfp)
    class(topology_t), intent(in) :: this
    logical :: selfp
    integer(i4) :: il, jl, itmp, algn

    ! count self periodic facets (edges or vertices)
    itmp = 0
    do il = 1, this%nfacet - 1
       do jl = il + 1, this%nfacet
          selfp = this%facet(il)%obj%equal(this%facet(jl)%obj%polytope)
          if (selfp) itmp = itmp + 1
       end do
    end do
    ! count self periodic ridges (vertices)
    do il = 1, this%nridge - 1
       do jl = il + 1, this%nridge
          selfp = (this%ridge(il)%obj%polytope%id() == &
               & this%ridge(jl)%obj%polytope%id())
          if (selfp) itmp = itmp + 1
       end do
    end do
    if (itmp == 0) then
       selfp = .false.
    else
       selfp = .true.
    end if
  end function topology_self_periodic

  !> @brief Return positions of facets shared by polytopes
  !! @note Polytopes can be self-periodic
  !! @parameter[in]   other   second polytope
  !! @parameter[out]  ishare  number of shared facets
  !! @parameter[out]  facetp  integer position of shared facets
  subroutine topology_facet_share(this, other, ishare, facetp)
    class(topology_t), intent(in) :: this, other
    integer(i4), intent(out) :: ishare
    integer(i4), dimension(:, :), allocatable, intent(out) :: facetp
    integer(i4) :: il, jl

    allocate(facetp(2, this%nfacet * other%nfacet))
    ishare = 0
    facetp(:, :) = 0
    do il = 1, this%nfacet
       do jl = 1, other%nfacet
          if (this%facet(il)%obj%equal(other%facet(jl)%obj%polytope)) then
             ishare = ishare + 1
             facetp(1, ishare) = il
             facetp(2, ishare) = jl
          end if
       end do
    end do
  end subroutine topology_facet_share

  !> @brief Return positions of ridges (vertices) shared by polytopes
  !! @note Plytopes can be self-periodic
  !! @parameter[in]   other   second polytope
  !! @parameter[out]  ishare  number of shared vertices
  !! @parameter[out]  ridgep  integer position of shared vertices
  pure subroutine topology_ridge_share(this, other, ishare, ridgep)
    class(topology_t), intent(in) :: this, other
    integer(i4), intent(out) :: ishare
    integer(i4), dimension(:, :), allocatable, intent(out) :: ridgep
    integer(i4) :: il, jl

    allocate(ridgep(2, this%nridge * other%nridge))
    ishare = 0
    ridgep(:, :) = 0
    do il = 1, this%nridge
       do jl = 1, other%nridge
          if (this%ridge(il)%obj%polytope%id() == &
               & other%ridge(jl)%obj%polytope%id()) then
             ishare = ishare + 1
             ridgep(1, ishare) = il
             ridgep(2, ishare) = jl
          end if
       end do
    end do
  end subroutine topology_ridge_share

  !> Return facet alignment
  !! @return algn
  function topology_fct_algn(this, pos) result(algn)
    class(topology_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4) :: algn
    if ((pos > 0) .and. (pos <= this%nfacet)) then
       algn = this%facet(pos)%obj%algn()
    else
       call neko_error('Wrong facet number for topology objects.')
    end if
  end function topology_fct_algn

end module topology
