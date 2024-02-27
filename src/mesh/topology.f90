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
!> Abstract type for abstract polytope class for mesh topology
module topology
  use num_types, only : i4, i8
  use utils, only : neko_error
  use polytope, only : polytope_t
  use polytope_oriented, only : polytope_oriented_t
  implicit none
  private

  public :: topology_component_t, topology_t, topology_element_t

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
     !> Return a pointer to the polytope peaks; not used
     procedure, pass(this) :: pek => topology_pek_ptr
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

  !> Single topology element allocatable space
  type :: topology_element_t
     class(topology_t), allocatable :: obj
  end type topology_element_t

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
    class(polytope_t), pointer :: ptr
    if ((pos > 0) .and. (pos <= this%nfacet)) then
       ptr => this%facet(pos)%obj%polytope
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
    class(polytope_t), pointer :: ptr
    if ((pos > 0) .and. (pos <= this%nridge)) then
       ptr => this%ridge(pos)%obj%polytope
    else
       call neko_error('Wrong ridge number for topology objects.')
    end if
  end function topology_rdg_ptr

  !> @brief Return pointer to the polytope peak; not used
  !! @parameter[in]   pos   polytope component position
  !! @return ptr
  function topology_pek_ptr(this, pos) result(ptr)
    class(topology_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    class(polytope_t), pointer :: ptr
    ptr => null()
    call neko_error('Topology objects have no peaks.')
  end function topology_pek_ptr

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
