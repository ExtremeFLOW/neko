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
!> Abstract type for abstract polytope class for mesh structure
module polytope
  use num_types, only : i4
  use utils, only : neko_error
  use entity, only : entity_t
  implicit none
  private

  public :: polytope_t

  !> Base type for an abstract polytope
  !! @details An abstract polytope of dimension @a dim_ is a partially ordered
  !! set of its lower dimension components: facets (dim -1), ridges (dim -2),
  !! peaks (dim -3)... I do not consider polytopes with dim > 3, so numbers of
  !! facets (@a nfacet) ridges (@a nridge) and peaks (@a npeak) are taken into
  !! account only.
  type, extends(entity_t), abstract :: polytope_t
     !> Polytope dimension
     integer(i4), private :: dim_ = -1
     !> Facet number
     integer(i4) :: nfacet = -1
     !> Ridge number
     integer(i4) :: nridge = -1
     !> Peak number
     integer(i4) :: npeak = -1
   contains
     !> Return polytope dimension
     procedure, pass(this) :: dim => polytope_dim_get
     !> Set polytope dimension
     procedure, pass(this) :: set_dim => polytope_dim_set
     !> Return polytope element number
     procedure, pass(this) :: nelem => polytope_element_number_get
     !> Set polytope element number
     procedure, pass(this) :: set_nelem => polytope_element_number_set
     !> Polytope equality on the level of dimension and element number
     procedure, pass(this) :: equal_poly => polytope_equal
     !> Pointer to the polytope facet
     procedure(polytope_ptr), pass(this), deferred :: fct
     !> Pointer to the polytope ridge
     procedure(polytope_ptr), pass(this), deferred :: rdg
     !> Pointer to the polytope peak
     procedure(polytope_ptr), pass(this), deferred :: pek
  end type polytope_t

  !> Returns polytope element pointer
  !! @parameter[in]   pos   polytope element position
  !! @return ptr
  abstract interface
     function polytope_ptr(this, pos) result(ptr)
       import i4
       import polytope_t
       class(polytope_t), intent(in) :: this
       integer(i4), intent(in) :: pos
       class(polytope_t), pointer :: ptr
     end function polytope_ptr
  end interface

contains
  !> @brief Get polytope dimension
  !! @return   dim
  pure function polytope_dim_get(this) result(dim)
    class(polytope_t), intent(in) :: this
    integer(i4) :: dim
    dim = this%dim_
  end function polytope_dim_get

  !> @brief Set polytope dimension
  !! @parameter[in]   dim     polytope dimension
  subroutine polytope_dim_set(this, dim)
    class(polytope_t), intent(inout) :: this
    integer(i4), intent(in) :: dim
    this%dim_ = dim
  end subroutine polytope_dim_set

  !> @brief Get polytope elements numbers
  !! @parameter[out]  nfacet   number of facets
  !! @parameter[out]  nridge  number of rdiges
  !! @parameter[out]  npeak   number of peaks
  pure subroutine polytope_element_number_get(this, nfacet, nridge, npeak)
    class(polytope_t), intent(in) :: this
    integer(i4), intent(out) :: nfacet, nridge, npeak
    nfacet = this%nfacet
    nridge = this%nridge
    npeak = this%npeak
  end subroutine polytope_element_number_get

  !> @brief Set polytope dimension
  !! @note This subroutine assumes @ dim to be set before
  !! @parameter[in]   nfacet  number of facets
  !! @parameter[in]   nridge  number of rdiges
  !! @parameter[in]   npeak   number of peaks
  subroutine polytope_element_number_set(this, nfacet, nridge, npeak)
    class(polytope_t), intent(inout) :: this
    integer(i4), intent(in) :: nfacet, nridge, npeak
    ! sanity check
    if (this%dim_ == -1) call neko_error('Polytope dimension not initialised')
    select case(this%dim_)
    case(0)
       if ((nfacet /= 0) .or. (nridge /= 0) .or. (npeak /= 0)) then
          call neko_error('Vertex has no elements.')
       else
          this%nfacet = 0
          this%nridge = 0
          this%npeak = 0
       end if
    case(1)
       if ((nfacet == 2) .and. (nridge == 0) .and. (npeak == 0)) then
          this%nfacet = 2
          this%nridge = 0
          this%npeak = 0
       else
           call neko_error('Edge contains two facets only.')
       end if
    case(2)
       if ((nfacet > 0) .and. (nridge > 0) .and. (npeak == 0)) then
          this%nfacet = nfacet
          this%nridge = nridge
          this%npeak = 0
       else
           call neko_error('Face has no peaks.')
       end if
    case(3)
       if ((nfacet > 0) .and. (nridge > 0) .and. (npeak > 0)) then
          this%nfacet = nfacet
          this%nridge = nridge
          this%npeak = npeak
       else
           call neko_error('Cell has all type of elements.')
       end if
    case default
       call neko_error('Unsupported polytope dimension')
    end select
  end subroutine polytope_element_number_set

  !> @brief Check if two polytopes have the same dimension and element numbers
  !! @return   equal
  pure function polytope_equal(this, other) result(equal)
    class(polytope_t), intent(in) :: this, other
    logical :: equal
    integer(i4) :: nfacet, nridge, npeak, nfaceto, nridgeo, npeako
    equal = (this%dim() == other%dim())
    if (equal) then
       call this%nelem(nfacet, nridge, npeak)
       call other%nelem(nfaceto, nridgeo, npeako)
       equal = (nfacet == nfaceto) .and. (nridge == nridgeo) .and. &
            & (npeak == npeako)
    end if
  end function polytope_equal

end module polytope
