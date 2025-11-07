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
!> Implementation of the mesh type for p4est mesh manager
module mesh_mesh_p4est
  use num_types, only : i4, i8, rp, dp
  use mesh_geom, only :  mesh_geom_t
  use mesh_conn, only :  mesh_conn_t
  use mesh_geom_p4est, only :  mesh_geom_p4est_t
  use mesh_conn_p4est, only :  mesh_conn_p4est_t
  use mesh_mesh, only : mesh_mesh_t

  implicit none
  private

  type, extends(mesh_mesh_t), public :: mesh_mesh_p4est_t
     !> number of local V-type elements
     integer(i4) :: nelv
     !> max refinement level across the whole mesh
     integer(i4) :: level_max
     !> element refinement level
     integer(i4), allocatable, dimension(:) :: level
     !> element group, (not used right now)
     integer(i4), allocatable, dimension(:) :: igrp
     !> face curvature flag (not used right now)
     integer(i4), allocatable, dimension(:,:) :: crv
     !> face boundary condition: -1- periodic, 0-internal, 0< user specified
     integer(i4), allocatable, dimension(:,:) :: bc
   contains
     !> Allocate types
     procedure, pass(this) :: init => mesh_mesh_init_p4est
     !> Constructor for the mesh data
     procedure, pass(this) :: init_data => mesh_mesh_init_data_p4est
     !> Constructor for the mesh data based on the other mesh type
     procedure, pass(this) :: init_type => mesh_mesh_init_type_p4est
     !> Free mesh data
     procedure, pass(this) :: free_data => mesh_mesh_free_data_p4est
     !> Free mesh type
     procedure, pass(this) :: free => mesh_mesh_free_p4est
  end type mesh_mesh_p4est_t

contains

  !> Allocate types
  subroutine mesh_mesh_init_p4est(this)
    class(mesh_mesh_p4est_t), intent(inout) :: this

    if (allocated(this%geom))then
       call this%geom%free()
       deallocate(this%geom)
    end if
    allocate(mesh_geom_p4est_t::this%geom)

    if (allocated(this%conn))then
       call this%conn%free()
       deallocate(this%conn)
    end if
    allocate(mesh_conn_p4est_t::this%conn)

  end subroutine mesh_mesh_init_p4est

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
  subroutine mesh_mesh_init_data_p4est(this, nelt, nelv, gnelt, gnelto, &
       level_max, tdim, gidx, level, igrp, crv, bc)
    class(mesh_mesh_p4est_t), intent(inout) :: this
    integer(i4), intent(in) :: nelt, nelv, level_max, tdim
    integer(i8), intent(in) :: gnelt, gnelto
    integer(i8), allocatable, dimension(:), intent(inout)  :: gidx
    integer(i4), allocatable, dimension(:), intent(inout) :: level, igrp
    integer(i4), allocatable, dimension(:,:), intent(inout) :: crv, bc

    call this%free_data()

    call this%init_data_base(nelt, gnelt, gnelto, tdim, gidx)

    this%nelv = nelv
    this%level_max  = level_max

    if (allocated(level)) call move_alloc(level, this%level)
    if (allocated(igrp)) call move_alloc(igrp, this%igrp)
    if (allocated(crv)) call move_alloc(crv, this%crv)
    if (allocated(bc)) call move_alloc(bc, this%bc)

  end subroutine mesh_mesh_init_data_p4est

  !>  Initialise mesh data based on another mesh type
  !! @param[inout] mesh   mesh data
  subroutine mesh_mesh_init_type_p4est(this, mesh)
    class(mesh_mesh_p4est_t), intent(inout) :: this
    class(mesh_mesh_t), intent(inout) :: mesh

    call this%free_data()

    call this%init_type_base(mesh)
    select type (mesh)
    type is(mesh_mesh_p4est_t)
       this%nelv = mesh%nelv
       this%level_max = mesh%level_max

       if (allocated(mesh%level)) call move_alloc(mesh%level, this%level)
       if (allocated(mesh%igrp)) call move_alloc(mesh%igrp, this%igrp)
       if (allocated(mesh%crv)) call move_alloc(mesh%crv, this%crv)
       if (allocated(mesh%bc)) call move_alloc(mesh%bc, this%bc)
    end select

  end subroutine mesh_mesh_init_type_p4est

  !> Destructor for the data in `mesh_manager_t` (base) type.
  subroutine mesh_mesh_free_data_p4est(this)
    class(mesh_mesh_p4est_t), intent(inout) :: this

    call this%free_data_base()

    this%nelv = 0
    this%level_max = 0

    if (allocated(this%level)) deallocate(this%level)
    if (allocated(this%igrp)) deallocate(this%igrp)
    if (allocated(this%crv)) deallocate(this%crv)
    if (allocated(this%bc)) deallocate(this%bc)

  end subroutine mesh_mesh_free_data_p4est

  !> Destructor for the data in `mesh_manager_t` (base) type.
  subroutine mesh_mesh_free_p4est(this)
    class(mesh_mesh_p4est_t), intent(inout) :: this

    call this%free_base()

    this%nelv = 0
    this%level_max = 0

    if (allocated(this%level)) deallocate(this%level)
    if (allocated(this%igrp)) deallocate(this%igrp)
    if (allocated(this%crv)) deallocate(this%crv)
    if (allocated(this%bc)) deallocate(this%bc)

  end subroutine mesh_mesh_free_p4est

end module mesh_mesh_p4est
