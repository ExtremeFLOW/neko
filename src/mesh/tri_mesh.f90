! Copyright (c) 2022, The Neko Authors
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
!> Defines a triangular surface mesh
!! @details Mesh derived from a surface geometry 
module tri_mesh
  use tri
  use point, only : point_t
  implicit none
  private

  type, public :: tri_mesh_t
     type(tri_t), allocatable :: el(:)       !< Tetrahedron elements
     type(point_t), allocatable :: points(:) !< List of points
     integer :: nelv
     integer :: mpts
     integer, private :: melv
   contains
     procedure, pass(this) :: init => tri_mesh_init
     procedure, pass(this) :: free => tri_mesh_free
     procedure, pass(this) :: add_element => tri_mesh_add_element
  end type tri_mesh_t

contains

  !> Initialise a triangular surface mesh
  subroutine tri_mesh_init(this, nelv)
    class(tri_mesh_t), intent(inout) :: this
    integer, intent(in) :: nelv

    call this%free()

    this%nelv = nelv
    allocate(this%el(this%nelv))
    allocate(this%points(nelv * NEKO_TRI_NPTS))

    this%mpts = 0
    this%melv = 0
    
  end subroutine tri_mesh_init

  !> Deallocate a triangular surface mesh
  subroutine tri_mesh_free(this)
    class(tri_mesh_t), intent(inout) :: this

    if (allocated(this%el)) then
       deallocate(this%el)
    end if

    if (allocated(this%points)) then
       deallocate(this%points)
    end if
    
  end subroutine tri_mesh_free

  !> Add an element to a mesh
  subroutine tri_mesh_add_element(this, p1, p2, p3)
    class(tri_mesh_t), intent(inout) :: this
    type(point_t), intent(inout) :: p1, p2, p3

    this%points(this%mpts + 1) = p1
    this%points(this%mpts + 2) = p2
    this%points(this%mpts + 3) = p3

    this%melv = this%melv + 1
    call this%el(this%melv)%init(this%melv, &
         this%points(this%mpts + 1), &
         this%points(this%mpts + 2), &
         this%points(this%mpts + 3))
    this%mpts = this%mpts + 3

  end subroutine tri_mesh_add_element
  
end module tri_mesh
