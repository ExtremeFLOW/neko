! Copyright (c) 2020-2023, The Neko Authors
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
!> Distributed mesh data
module distdata
  use stack, only : stack_i4t2_t
  use tuple, only : tuple_i4_t
  use uset, only : uset_i4_t
  implicit none
  private

  type, public :: distdata_t
     type(stack_i4t2_t) :: shared_el_facet !< Elemenets with shared facets
     
     type(uset_i4_t) :: shared_facet    !< List of shared facets
     type(uset_i4_t) :: shared_edge     !< List of shared edges
     type(uset_i4_t) :: shared_point    !< List of shared points
     
     integer, allocatable :: local_to_global_facet(:)!< Local to global (facets)
     integer, allocatable :: local_to_global_edge(:) !< Local to global (edges)
     
  end type distdata_t

  public :: distdata_init, distdata_free, distdata_set_shared_el_facet, &
       distdata_set_shared_facet, distdata_set_shared_edge, &
       distdata_set_shared_point, distdata_set_local_to_global_facet, &
       distdata_set_local_to_global_edge

contains

  !> Initialise a distdata type
  subroutine distdata_init(ddata)
    type(distdata_t), intent(inout) :: ddata

    call distdata_free(ddata)
    
    call ddata%shared_el_facet%init()

    call ddata%shared_facet%init()
    call ddata%shared_edge%init()   
    call ddata%shared_point%init()
    
  end subroutine distdata_init

  !> Free a distdata type
  subroutine distdata_free(ddata)
    type(distdata_t), intent(inout) :: ddata

    call ddata%shared_el_facet%free()
    
    call ddata%shared_facet%free()
    call ddata%shared_edge%free()
    call ddata%shared_point%free()

    if (allocated(ddata%local_to_global_facet)) then
       deallocate(ddata%local_to_global_facet)
    end if
    
    if (allocated(ddata%local_to_global_edge)) then
       deallocate(ddata%local_to_global_edge)
    end if
    
  end subroutine distdata_free

  !> Mark an element's facet as shared
  subroutine distdata_set_shared_el_facet(ddata, element, side)
    type(distdata_t), intent(inout) :: ddata
    integer, intent(in), value :: element !< Element index (local numbering)
    integer, intent(in), value :: side    !< Facet index
    type(tuple_i4_t) :: t

    t%x = (/ element, side /)
    call ddata%shared_el_facet%push(t)
    
  end subroutine distdata_set_shared_el_facet

  !> Mark a facet as shared
  subroutine distdata_set_shared_facet(ddata, facet)
    type(distdata_t), intent(inout) :: ddata
    integer, value :: facet     !< Facet index (local numbering)

    call ddata%shared_facet%add(facet)
    
  end subroutine distdata_set_shared_facet
  
  !> Mark an element's edge as shared
  !! @attention only defined for elements where facet .ne. edges
  subroutine distdata_set_shared_edge(ddata, edge)
    type(distdata_t), intent(inout) :: ddata
    integer, value :: edge      !< Edge index (local numbering) 

    call ddata%shared_edge%add(edge)
    
  end subroutine distdata_set_shared_edge

  !> Mark a point as shared
  subroutine distdata_set_shared_point(ddata, point)
    type(distdata_t), intent(inout) :: ddata
    integer, value :: point !< Point index (local numbering)

    call ddata%shared_point%add(point)
    
  end subroutine distdata_set_shared_point

  !> Set local to global mapping (facets)
  subroutine distdata_set_local_to_global_facet(ddata, local, global)
    type(distdata_t), intent(inout) :: ddata
    integer, intent(in), value :: local  !< Local facet index
    integer, intent(in), value :: global !< Global facet index

    ddata%local_to_global_facet(local) = global
    
  end subroutine distdata_set_local_to_global_facet

  !> Set local to global mapping (edges)
  subroutine distdata_set_local_to_global_edge(ddata, local, global)
    type(distdata_t), intent(inout) :: ddata
    integer, intent(in) , value :: local  !< Local edge index
    integer, intent(in) , value :: global !< Global edge index

    ddata%local_to_global_edge(local) = global
    
  end subroutine distdata_set_local_to_global_edge
  
end module distdata
