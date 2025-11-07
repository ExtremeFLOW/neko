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

  type, public :: mesh_mesh_t
     !> Topological mesh dimension
     integer(i4) :: tdim
     !> total number of local element (both V- and T-type)
     integer(i4) :: nelt
     !> number of local V-type elements
     integer(i4) :: nelv
     !> global element number
     integer(i8) :: gnelt
     !> global element offset
     integer(i8) :: gnelto
     !> max refinement level across the whole mesh
     integer(i4) :: level_max
     !> global element index
     integer(i8), allocatable, dimension(:) :: gidx
     !> element refinement level
     integer(i4), allocatable, dimension(:) :: level
     !> element group, (not used right now)
     integer(i4), allocatable, dimension(:) :: igrp
     !> Number of faces per element
     integer(i4) :: nfcs
     !> face curvature flag (not used right now)
     integer(i4), allocatable, dimension(:,:) :: crv
     !> face boundary condition: -1- periodic, 0-internal, 0< user specified
     integer(i4), allocatable, dimension(:,:) :: bc
     !> Geometrical information
     type(mesh_geom_t) :: geom
     !> Connectivity information
     type(mesh_conn_t) :: conn
   contains
     !> Constructor for the mesh data
     procedure, pass(this) :: init_data => mesh_mesh_init_data
     !> Constructor for the mesh data based on the other mesh type
     procedure, pass(this) :: init_type => mesh_mesh_init_type
     !> Destructor for the mesh data
     procedure, pass(this) :: free => mesh_mesh_free
  end type mesh_mesh_t

contains

  !>  Initialise mesh data
  !! @param[in]    nelt       total number of local elements
  !! @param[in]    nelv       number of V-type elements
  !! @param[in]    gnelt      global number of all the elements
  !! @param[in]    gnelto     global element offset
  !! @param[in]    level_max  max refinement level across the mesh
  !! @param[in]    tdim       topological mesh dimension
  !! @param[inout] gidx       global element number
  !! @param[inout] level      element refinement level
  !! @param[inout] igrp       element group
  !! @param[inout] crv        face curvature data
  !! @param[inout] bc         face boundary condition
  subroutine mesh_mesh_init_data(this, nelt, nelv, gnelt, gnelto, &
       level_max, tdim, gidx, level, igrp, crv, bc)
    class(mesh_mesh_t), intent(inout) :: this
    integer(i4), intent(in) :: nelt, nelv, level_max, tdim
    integer(i8), intent(in) :: gnelt, gnelto
    integer(i8), allocatable, dimension(:), intent(inout)  :: gidx
    integer(i4), allocatable, dimension(:), intent(inout) :: level, igrp
    integer(i4), allocatable, dimension(:,:), intent(inout) :: crv, bc

    call this%free()

    this%tdim = tdim
    this%nelt = nelt
    this%nelv = nelv
    this%gnelt = gnelt
    this%gnelto = gnelto
    this%level_max  = level_max
    ! we work with hex/quad only and there is no difference between
    ! topology and geometrical dimensions
    this%nfcs = 2 * tdim

    if (allocated(gidx)) call move_alloc(gidx, this%gidx)
    if (allocated(level)) call move_alloc(level, this%level)
    if (allocated(igrp)) call move_alloc(igrp, this%igrp)
    if (allocated(crv)) call move_alloc(crv, this%crv)
    if (allocated(bc)) call move_alloc(bc, this%bc)

  end subroutine mesh_mesh_init_data

  !>  Initialise mesh data based on another mesh type
  !! @param[inout] mesh   mesh data
  subroutine mesh_mesh_init_type(this, mesh)
    class(mesh_mesh_t), intent(inout) :: this
    type(mesh_mesh_t), intent(inout) :: mesh

    call this%free()

    call this%geom%init_type(mesh%geom)
    call this%conn%init_type(mesh%conn)

    this%tdim = mesh%tdim
    this%nelt = mesh%nelt
    this%nelv = mesh%nelv
    this%gnelt = mesh%gnelt
    this%gnelto = mesh%gnelto
    this%level_max = mesh%level_max
    this%nfcs = mesh%nfcs

    if (allocated(mesh%gidx)) call move_alloc(mesh%gidx, this%gidx)
    if (allocated(mesh%level)) call move_alloc(mesh%level, this%level)
    if (allocated(mesh%igrp)) call move_alloc(mesh%igrp, this%igrp)
    if (allocated(mesh%crv)) call move_alloc(mesh%crv, this%crv)
    if (allocated(mesh%bc)) call move_alloc(mesh%bc, this%bc)

  end subroutine mesh_mesh_init_type

  !> Destructor for the data in `mesh_manager_t` (base) type.
  subroutine mesh_mesh_free(this)
    class(mesh_mesh_t), intent(inout) :: this

    call this%geom%free()
    call this%conn%free()

    this%tdim = 0
    this%nelt = 0
    this%nelv = 0
    this%gnelt = 0
    this%gnelto = 0
    this%level_max = 0
    this%nfcs = 0

    if (allocated(this%gidx)) deallocate(this%gidx)
    if (allocated(this%level)) deallocate(this%level)
    if (allocated(this%igrp)) deallocate(this%igrp)
    if (allocated(this%crv)) deallocate(this%crv)
    if (allocated(this%bc)) deallocate(this%bc)

  end subroutine mesh_mesh_free

end module mesh_mesh
