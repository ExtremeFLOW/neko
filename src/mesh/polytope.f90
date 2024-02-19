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
!> Abstract type for abstract polytope class for mesh structure
module polytope
  use num_types, only : i4
  use utils, only : neko_error
  use entity, only : entity_t
  implicit none
  private

  public :: polytope_t, polytope_ptr

  !> Base type for an abstract polytope
  !! @details An abstract polytope of dimension @a tdim_ is a partially ordered
  !! set of its lower dimension components: facets (dim -1), ridges (dim -2),
  !! peaks (dim -3)... We do not consider polytopes with dim > 3, so numbers of
  !! facets (@a nfacet) ridges (@a nridge) and peaks (@a npeak) are taken into
  !! account only.
  type, extends(entity_t), abstract :: polytope_t
     !> Polytope dimension
     integer(i4), private :: tdim_ = -1
     !> Facet number
     integer(i4) :: nfacet = -1
     !> Ridge number
     integer(i4) :: nridge = -1
     !> Peak number
     integer(i4) :: npeak = -1
   contains
     !> Free polytope data
     procedure, pass(this) :: free_base => polytope_free_base
     !> Return polytope dimension
     procedure, pass(this) :: tdim => polytope_tdim_get
     !> Set polytope dimension
     procedure, pass(this) :: set_tdim => polytope_tdim_set
     !> Return polytope component number
     procedure, pass(this) :: ncomp => polytope_component_number_get
     !> Set polytope component number
     procedure, pass(this) :: set_ncomp => polytope_component_number_set
     !> Polytope equality on the level of dimension and component number
     procedure, pass(this) :: equal_poly => polytope_equal_poly
     !> Compare polytope component numbers
     procedure, pass(this) :: check_comp => polytope_check_component
     !> Free type
     procedure(polytope_free), pass(this), deferred :: free
     !> Pointer to the polytope facet
     procedure(polytope_pointer), pass(this), deferred :: fct
     !> Pointer to the polytope ridge
     procedure(polytope_pointer), pass(this), deferred :: rdg
     !> Pointer to the polytope peak
     procedure(polytope_pointer), pass(this), deferred :: pek
     !> Test equality
     procedure(polytope_equal), pass(this), deferred :: equal
     !> Test self-periodicity
     procedure(polytope_self_periodic), pass(this), deferred :: self_periodic
     !> Return facet alignment
     procedure(polytope_algn), pass(this), deferred :: falgn
     !> Return boundary information (topology object only)
     procedure(polytope_int), pass(this), deferred :: bnd
     !> Return communication id (topology object only)
     procedure(polytope_int), pass(this), deferred :: gsid
  end type polytope_t

  !> The general pointer type to the polytope class
  type :: polytope_ptr
     class(polytope_t), pointer :: ptr
  end type polytope_ptr

  !> Free type template
  abstract interface
     subroutine polytope_free(this)
       import polytope_t
       class(polytope_t), intent(inout) :: this
     end subroutine polytope_free
  end interface

  !> Returns polytope component pointer
  !! @parameter[in]   pos   polytope component position
  !! @return ptr
  abstract interface
     function polytope_pointer(this, pos) result(ptr)
       import i4
       import polytope_t
       class(polytope_t), intent(in) :: this
       integer(i4), intent(in) :: pos
       class(polytope_t), pointer :: ptr
     end function polytope_pointer
  end interface

  !> Test equality
  !! @parameter[in]   other   polytope
  !! @return equal
  abstract interface
     function polytope_equal(this, other) result(equal)
       import polytope_t
       class(polytope_t), intent(in) :: this, other
       logical :: equal
     end function polytope_equal
  end interface

  !> Test self-periodicity
  !! @parameter[in]   other   polytope
  !! @return equal
  abstract interface
     function polytope_self_periodic(this) result(selfp)
       import polytope_t
       class(polytope_t), intent(in) :: this
       logical :: selfp
     end function polytope_self_periodic
  end interface

  !> Return alignment
  !! @parameter[in]   pos   polytope component position
  !! @return algn
  abstract interface
     function polytope_algn(this, pos) result(algn)
       import i4
       import polytope_t
       class(polytope_t), intent(in) :: this
       integer(i4), intent(in) :: pos
       integer(i4) :: algn
     end function polytope_algn
  end interface

  !> Extract integer property
  !! @return intp
  abstract interface
     function polytope_int(this) result(intp)
       import i4
       import polytope_t
       class(polytope_t), intent(in) :: this
       integer(i4) :: intp
     end function polytope_int
  end interface

contains
  !> Free polytope data
  subroutine polytope_free_base(this)
    class(polytope_t), intent(inout) :: this
    this%tdim_ = -1
    this%nfacet = -1
    this%nridge = -1
    this%npeak = -1
  end subroutine polytope_free_base

  !> @brief Get polytope dimension
  !! @return   dim
  pure function polytope_tdim_get(this) result(dim)
    class(polytope_t), intent(in) :: this
    integer(i4) :: dim
    dim = this%tdim_
  end function polytope_tdim_get

  !> @brief Set polytope dimension
  !! @parameter[in]   dim     polytope dimension
  subroutine polytope_tdim_set(this, dim)
    class(polytope_t), intent(inout) :: this
    integer(i4), intent(in) :: dim
    this%tdim_ = dim
  end subroutine polytope_tdim_set

  !> @brief Get polytope components numbers
  !! @parameter[out]  nfacet   number of facets
  !! @parameter[out]  nridge  number of rdiges
  !! @parameter[out]  npeak   number of peaks
  pure subroutine polytope_component_number_get(this, nfacet, nridge, npeak)
    class(polytope_t), intent(in) :: this
    integer(i4), intent(out) :: nfacet, nridge, npeak
    nfacet = this%nfacet
    nridge = this%nridge
    npeak = this%npeak
  end subroutine polytope_component_number_get

  !> @brief Set polytope component numbers
  !! @note This subroutine assumes @ dim to be set before
  !! @parameter[in]   nfacet  number of facets
  !! @parameter[in]   nridge  number of rdiges
  !! @parameter[in]   npeak   number of peaks
  subroutine polytope_component_number_set(this, nfacet, nridge, npeak)
    class(polytope_t), intent(inout) :: this
    integer(i4), intent(in) :: nfacet, nridge, npeak
    ! sanity check
    if (this%tdim_ == -1) call neko_error('Polytope dimension not initialised')
    select case(this%tdim_)
    case(0)
       if ((nfacet /= 0) .or. (nridge /= 0) .or. (npeak /= 0)) then
          call neko_error('Vertex has no components.')
       end if
    case(1)
       if ((nfacet /= 2) .or. (nridge /= 0) .or. (npeak /= 0)) then
           call neko_error('Edge contains two facets only.')
       end if
    case(2)
       if ((nfacet <= 0) .or. (nridge <= 0) .or. (npeak /= 0)) then
           call neko_error('Face has no peaks.')
       end if
    case(3)
       if ((nfacet <= 0) .or. (nridge <= 0) .or. (npeak <= 0)) then
          call neko_error('Cell has all type of components.')
       end if
    case default
       call neko_error('Unsupported polytope dimension')
    end select
    this%nfacet = nfacet
    this%nridge = nridge
    this%npeak = npeak
  end subroutine polytope_component_number_set

  !> @brief Check if two polytopes have the same dimension and component numbers
  !! @parameter[in]    other  polytope
  !! @return   equal
  pure function polytope_equal_poly(this, other) result(equal)
    class(polytope_t), intent(in) :: this, other
    logical :: equal
    integer(i4) :: nfacet, nridge, npeak, nfaceto, nridgeo, npeako
    equal = (this%tdim() == other%tdim())
    if (equal) then
       call this%ncomp(nfacet, nridge, npeak)
       call other%ncomp(nfaceto, nridgeo, npeako)
       equal = (nfacet == nfaceto) .and. (nridge == nridgeo) .and. &
            & (npeak == npeako)
    end if
  end function polytope_equal_poly

  !> @brief Check the polytope component numbers
  !! @parameter[in]   nfacet  number of facets
  !! @parameter[in]   nridge  number of rdiges
  !! @parameter[in]   npeak   number of peaks
  !! @return   equal
  pure function polytope_check_component(this, nfacet, nridge, npeak) &
       & result(equal)
    class(polytope_t), intent(in) :: this
    integer(i4), intent(in) :: nfacet, nridge, npeak
    logical :: equal
    equal = (nfacet == this%nfacet) .and. (nridge == this%nridge) .and. &
         & (npeak == this%npeak)
  end function polytope_check_component

end module polytope
