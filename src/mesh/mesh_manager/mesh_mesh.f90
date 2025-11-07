! Copyright (c) 2019-2025, The Neko Authors
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
!> Implementation of the mesh type for mesh manager
module mesh_mesh
  use num_types, only : i4, i8, rp, dp
  use mesh_geom, only :  mesh_geom_t
  use mesh_conn, only :  mesh_conn_t

  implicit none
  private

  type, abstract, public :: mesh_mesh_t
     !> Topological mesh dimension
     integer(i4) :: tdim
     !> total number of local element (both V- and T-type)
     integer(i4) :: nelt
     !> global element number
     integer(i8) :: gnelt
     !> global element offset
     integer(i8) :: gnelto
     !> global element index
     integer(i8), allocatable, dimension(:) :: gidx
     !> Number of faces per element
     integer(i4) :: nfcs
     !> Geometrical information
     class(mesh_geom_t), allocatable :: geom
     !> Connectivity information
     class(mesh_conn_t), allocatable :: conn
   contains
     !> Constructor for the mesh data
     procedure, pass(this) :: init_data_base => mesh_mesh_init_data_base
     !> Constructor for the mesh data based on the other mesh type
     procedure, pass(this) :: init_type_base => mesh_mesh_init_type_base
     !> Free mesh data
     procedure, pass(this) :: free_data_base => mesh_mesh_free_data_base
     !> Free mesh type
     procedure, pass(this) :: free_base => mesh_mesh_free_base
     !> Allocate types
     procedure(mesh_mesh_init), pass(this), deferred :: init
     !> Initialise data from type
     procedure(mesh_mesh_init_type), pass(this), deferred :: init_type
     !> Free type_data
     procedure(mesh_mesh_init), pass(this), deferred :: free_data
     !> Free type
     procedure(mesh_mesh_init), pass(this), deferred :: free
  end type mesh_mesh_t

  abstract interface
     subroutine mesh_mesh_init(this)
       import mesh_mesh_t
       class(mesh_mesh_t), intent(inout) :: this
     end subroutine mesh_mesh_init

     subroutine mesh_mesh_init_type(this, mesh)
       import mesh_mesh_t
       class(mesh_mesh_t), intent(inout) :: this
       class(mesh_mesh_t), intent(inout) :: mesh
     end subroutine mesh_mesh_init_type
  end interface

contains

  !>  Initialise mesh data
  !! @param[in]    nelt       total number of local elements
  !! @param[in]    gnelt      global number of all the elements
  !! @param[in]    gnelto     global element offset
  !! @param[in]    tdim       topological mesh dimension
  !! @param[inout] gidx       global element number
  subroutine mesh_mesh_init_data_base(this, nelt, gnelt, gnelto, tdim, gidx)
    class(mesh_mesh_t), intent(inout) :: this
    integer(i4), intent(in) :: nelt, tdim
    integer(i8), intent(in) :: gnelt, gnelto
    integer(i8), allocatable, dimension(:), intent(inout)  :: gidx

    call this%free_data_base()

    this%tdim = tdim
    this%nelt = nelt
    this%gnelt = gnelt
    this%gnelto = gnelto
    ! we work with hex/quad only and there is no difference between
    ! topology and geometrical dimensions
    this%nfcs = 2 * tdim

    if (allocated(gidx)) call move_alloc(gidx, this%gidx)

  end subroutine mesh_mesh_init_data_base

  !>  Initialise mesh data based on another mesh type
  !! @param[inout] mesh   mesh data
  subroutine mesh_mesh_init_type_base(this, mesh)
    class(mesh_mesh_t), intent(inout) :: this
    class(mesh_mesh_t), intent(inout) :: mesh

    call this%free_data_base()

    if (allocated(this%geom) .and. allocated(mesh%geom)) &
         call this%geom%init_type(mesh%geom)
    if (allocated(this%conn) .and. allocated(mesh%conn)) &
         call this%conn%init_type(mesh%conn)

    this%tdim = mesh%tdim
    this%nelt = mesh%nelt
    this%gnelt = mesh%gnelt
    this%gnelto = mesh%gnelto
    this%nfcs = mesh%nfcs

    if (allocated(mesh%gidx)) call move_alloc(mesh%gidx, this%gidx)

  end subroutine mesh_mesh_init_type_base

  !> Destructor for the data in `mesh_manager_t` (base) type.
  subroutine mesh_mesh_free_data_base(this)
    class(mesh_mesh_t), intent(inout) :: this

    if (allocated(this%geom)) call this%geom%free_data()
    if (allocated(this%conn)) call this%conn%free_data()

    this%tdim = 0
    this%nelt = 0
    this%gnelt = 0
    this%gnelto = 0
    this%nfcs = 0

    if (allocated(this%gidx)) deallocate(this%gidx)

  end subroutine mesh_mesh_free_data_base

  !> Destructor for the data in `mesh_manager_t` (base) type.
  subroutine mesh_mesh_free_base(this)
    class(mesh_mesh_t), intent(inout) :: this

    call this%free_data_base()
    if (allocated(this%geom)) deallocate(this%geom)
    if (allocated(this%conn)) deallocate(this%conn)

  end subroutine mesh_mesh_free_base

end module mesh_mesh
