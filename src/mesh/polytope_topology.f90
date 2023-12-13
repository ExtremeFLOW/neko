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
module polytope_topology
  use num_types, only : i4
  use utils, only : neko_error
  use polytope, only : polytope_t
  use polytope_aligned, only : topology_object_t
  implicit none
  private

  public :: polytope_topology_t

  !> Base type for an abstract topology polytope
  !! @details This is an abstract type combining a set of lower dimension
  !! polytopes and boundary condition into the topology object.
  type, extends(polytope_t), abstract :: polytope_topology_t
     !> Polytope facets
     type(topology_object_t), dimension(:), allocatable :: facet
     !> Polytope ridges
     type(topology_object_t), dimension(:), allocatable :: ridge
     !> Internal/external boundary condition flag
     integer(i4), private :: boundary_ = -1
   contains
     !> Free polytope data
     procedure, pass(this) :: free => polytope_free
     !> Return a pointer to the polytope facets
     procedure, pass(this) :: fct => polytope_fct_ptr
     !> Return a pointer to the polytope ridges
     procedure, pass(this) :: rdg => polytope_rdg_ptr
     !> Return a pointer to the polytope peaks; not used
     procedure, pass(this) :: pek => polytope_pek_ptr
     !> Return boundary value
     procedure, pass(this) :: bnd => polytope_bnd_get
     !> Initialise an aligned polytope
     procedure(polytope_topology_init), pass(this), deferred :: init
  end type polytope_topology_t

  !> Abstract interface to initialise a polytope with alignment information
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  !! @parameter[in]   bnd    external boundary information
  abstract interface
     subroutine polytope_topology_init(this, pltp, algn, bnd)
       import i4
       import polytope_t
       import polytope_topology_t
       class(polytope_topology_t), intent(inout) :: this
       class(polytope_t), target, intent(in) :: pltp
       integer(i4), intent(in) :: algn, bnd
     end subroutine polytope_topology_init
  end interface

contains

  !> Free polytope data
  subroutine polytope_free(this)
    class(polytope_topology_t), intent(inout) :: this
    integer(i4) :: il
    this%boundary_ = -1
    if (allocated(this%facet)) then
       do il = 1, this%nfacet
          call this%facet(il)%obj%free()
          deallocate(this%facet(il)%obj)
       end do
       deallocate(this%facet)
    end if
    if (allocated(this%facet)) then
       do il = 1, this%nfacet
          call this%ridge(il)%obj%free()
          deallocate(this%ridge(il)%obj)
       end do
       deallocate(this%ridge)
    end if
  end subroutine polytope_free

  !> @brief Return pointer to the polytope facet
  !! @parameter[in]   pos   polytope element position
  !! @return ptr
  function polytope_fct_ptr(this, pos) result(ptr)
    class(polytope_topology_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    class(polytope_t), pointer :: ptr
    if ((pos > 0) .and. (pos <= this%nfacet)) then
       ptr => this%facet(pos)%obj%polytope
    else
       call neko_error('Wrong facet number for topology objects.')
    end if
  end function polytope_fct_ptr

  !> @brief Return pointer to the polytope ridge
  !! @parameter[in]   pos   polytope element position
  !! @return ptr
  function polytope_rdg_ptr(this, pos) result(ptr)
    class(polytope_topology_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    class(polytope_t), pointer :: ptr
    if ((pos > 0) .and. (pos <= this%nridge)) then
       ptr => this%ridge(pos)%obj%polytope
    else
       call neko_error('Wrong ridge number for topology objects.')
    end if
  end function polytope_rdg_ptr

  !> @brief Return pointer to the polytope peak; not used
  !! @parameter[in]   pos   polytope element position
  !! @return ptr
  function polytope_pek_ptr(this, pos) result(ptr)
    class(polytope_topology_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    class(polytope_t), pointer :: ptr
    ptr => null()
    call neko_error('Topology objects have no peaks.')
  end function polytope_pek_ptr

  !> @brief Get polytope boundary information
  !! @return   bnd
  pure function polytope_bnd_get(this) result(bnd)
    class(polytope_topology_t), intent(in) :: this
    integer(i4) :: bnd
    bnd = this%boundary_
  end function polytope_bnd_get

end module polytope_topology
