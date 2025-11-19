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
module manager_mesh_p4est
  use num_types, only : i4, i8, rp, dp
  use utils, only : neko_error
  use manager_geom, only :  manager_geom_t
  use manager_conn, only :  manager_conn_t
  use manager_geom_p4est, only :  manager_geom_p4est_t
  use manager_conn_p4est, only :  manager_conn_p4est_t
  use manager_mesh, only : manager_mesh_t

  implicit none
  private

  type, extends(manager_mesh_t), public :: manager_mesh_p4est_t
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
     ! Flag elements that can be coarsened and share the same parent
     ! family(1, lelt) - mark of a parent (not a real element number as parents
     ! do not exist on neko side)
     !                         0 - element that cannot be coarsened
     !                        >0 - family mark
     ! family(2, lelt) - vertex number shared by all the family members
     !> Family flag
     integer(i8), allocatable, dimension(:,:) :: family
   contains
     !> Allocate types
     procedure, pass(this) :: init => manager_mesh_init_p4est
     !> Constructor for the mesh data
     procedure, pass(this) :: init_data => manager_mesh_init_data_p4est
     !> Constructor for the mesh data based on the other mesh type
     procedure, pass(this) :: init_type => manager_mesh_init_type_p4est
     !> Free mesh data
     procedure, pass(this) :: free_data => manager_mesh_free_data_p4est
     !> Free mesh type
     procedure, pass(this) :: free => manager_mesh_free_p4est
  end type manager_mesh_p4est_t

  ! connectivity parameter arrays
  ! face vertices
  integer, public, parameter, dimension(4,6) :: p4_vface = reshape(&
       (/ 1,3,5,7 , 2,4,6,8 , 1,2,5,6 , 3,4,7,8 , 1,2,3,4 , 5,6,7,8 /),&
       shape(p4_vface))

  ! edge vertices
  integer, public, parameter, dimension(2,12) :: p4_vedge  = reshape(&
       (/ 1,2 , 3,4 , 5,6 , 7,8 , 1,3 , 2,4 , 5,7 , 6,8 , 1,5 , 2,6 , 3,7&
       , 4,8 /),shape(p4_vedge))

  ! edge related faces
  integer, public, parameter, dimension(2,12) :: p4_eface  = reshape(&
       (/ 3,5 , 4,5 , 3,6 , 4,6 , 1,5 , 2,5 , 1,6 , 2,6 , 1,3 , 2,3 , 1,4&
       , 2,4 /),shape(p4_eface))

  ! corner related faces
  integer, public, parameter, dimension(3,8) :: p4_cface = reshape(&
       (/ 1,3,5 , 2,3,5 , 1,4,5 , 2,4,5 , 1,3,6 , 2,3,6 , 1,4,6 , 2,4,6 /),&
       shape(p4_cface))

  ! corner related edges
  integer, public, parameter, dimension(3,8) :: p4_cedge = reshape(&
       (/ 1,5,9 , 1,6,10 , 2,5,11 , 2,6,12 , 3,7,9 , 3,8,10 , 4,7,11 ,&
        4,8,12 /),shape(p4_cedge))

  ! corner to face corner
  integer, public, parameter, dimension(6,8) :: p4_cfcrn = reshape(&
       (/ 1,-1, 1,-1, 1,-1 , -1, 1, 2,-1, 2,-1 ,  2,-1,-1, 1, 3,-1 &
       , -1, 2,-1, 2, 4,-1 ,  3,-1, 3,-1,-1, 1 , -1, 3, 4,-1,-1, 2 &
       ,  4,-1,-1, 3,-1, 3 , -1, 4,-1, 4,-1, 4 /),shape(p4_cfcrn))

  ! to calculate neighbour face corner
  integer, public, parameter, dimension(6,6) :: p4_rt =reshape( &
       (/ 1,2,2,1,1,2 , 3,1,1,2,2,1 , 3,1,1,2,2,1 , 1,3,3,1,1,2 &
       , 1,3,3,1,1,2 , 3,1,1,3,3,1 /),shape(p4_rt))
  integer, public, parameter, dimension(4,3) :: p4_qt = reshape(&
       (/ 2,3,6,7 , 1,4,5,8 , 1,5,4,8 /),shape(p4_qt))
  integer, public, parameter, dimension(4,8) :: p4_pt = reshape(&
       (/ 1,2,3,4 , 1,3,2,4 , 2,1,4,3 , 2,4,1,3 , 3,1,4,2 , 3,4,1,2 &
       , 4,2,3,1 , 4,3,2,1 /),shape(p4_pt))

contains

  !> Allocate types
  subroutine manager_mesh_init_p4est(this)
    class(manager_mesh_p4est_t), intent(inout) :: this

    if (allocated(this%geom))then
       call this%geom%free()
       deallocate(this%geom)
    end if
    allocate(manager_geom_p4est_t::this%geom)
    call this%geom%init()

    if (allocated(this%conn))then
       call this%conn%free()
       deallocate(this%conn)
    end if
    allocate(manager_conn_p4est_t::this%conn)
     call this%conn%init()

  end subroutine manager_mesh_init_p4est

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
  !! @param[inout] family     family flag
  !! @param[in]   ifsave      save component types
  subroutine manager_mesh_init_data_p4est(this, nelt, nelv, gnelt, gnelto, &
       level_max, tdim, gidx, level, igrp, crv, bc, family, ifsave)
    class(manager_mesh_p4est_t), intent(inout) :: this
    integer(i4), intent(in) :: nelt, nelv, level_max, tdim
    integer(i8), intent(in) :: gnelt, gnelto
    integer(i8), allocatable, dimension(:), intent(inout)  :: gidx
    integer(i4), allocatable, dimension(:), intent(inout) :: level, igrp
    integer(i4), allocatable, dimension(:,:), intent(inout) :: crv, bc
    integer(i8), allocatable, dimension(:,:), intent(inout) :: family
    logical, optional, intent(in) :: ifsave
    integer(i4) :: nvert, nface, nedge

    nvert = 2**tdim
    nface = 2 * tdim
    nedge = 12 * (tdim - 2)

    ! sanity check
    if ((nelt .ne. size(gidx)) .or. (nelt .ne. size(level)) .or. &
         (nelt .ne. size(igrp)) .or. &
         (nface .ne. size(crv, 1)) .or. (nelt .ne. size(crv, 2)) .or. &
         (nface .ne. size(bc, 1)) .or. (nelt .ne. size(bc, 2)) .or. &
         (2 .ne. size(family, 1)) .or. (nelt .ne. size(family, 2))) &
         call neko_error('Inconsistent array sizes; p4est%mesh')

    if (present(ifsave)) then
       call this%free_data(ifsave)
       call this%init_data_base(nelt, gnelt, gnelto, tdim, gidx, ifsave)
    else
       call this%free_data()
       call this%init_data_base(nelt, gnelt, gnelto, tdim, gidx)
    end if

    this%nelv = nelv
    this%level_max  = level_max

    if (allocated(level)) call move_alloc(level, this%level)
    if (allocated(igrp)) call move_alloc(igrp, this%igrp)
    if (allocated(crv)) call move_alloc(crv, this%crv)
    if (allocated(bc)) call move_alloc(bc, this%bc)
    if (allocated(family)) call move_alloc(family, this%family)

  end subroutine manager_mesh_init_data_p4est

  !>  Initialise mesh data based on another mesh type
  !! @param[inout] mesh   mesh data
  subroutine manager_mesh_init_type_p4est(this, mesh)
    class(manager_mesh_p4est_t), intent(inout) :: this
    class(manager_mesh_t), intent(inout) :: mesh

    call this%free_data()

    call this%init_type_base(mesh)

    select type (mesh)
    type is(manager_mesh_p4est_t)
       this%nelv = mesh%nelv
       this%level_max = mesh%level_max

       if (allocated(mesh%level)) call move_alloc(mesh%level, this%level)
       if (allocated(mesh%igrp)) call move_alloc(mesh%igrp, this%igrp)
       if (allocated(mesh%crv)) call move_alloc(mesh%crv, this%crv)
       if (allocated(mesh%bc)) call move_alloc(mesh%bc, this%bc)
       if (allocated(mesh%family)) call move_alloc(mesh%family, this%family)
    end select

  end subroutine manager_mesh_init_type_p4est

  !> Destructor for the data in `mesh_manager_t` (base) type.
  !! @param[in]   ifsave      save component types
  subroutine manager_mesh_free_data_p4est(this, ifsave)
    class(manager_mesh_p4est_t), intent(inout) :: this
    logical, optional, intent(in) :: ifsave

    if (present(ifsave)) then
       call this%free_data_base(ifsave)
    else
       call this%free_data_base()
    end if

    this%nelv = 0
    this%level_max = 0

    if (allocated(this%level)) deallocate(this%level)
    if (allocated(this%igrp)) deallocate(this%igrp)
    if (allocated(this%crv)) deallocate(this%crv)
    if (allocated(this%bc)) deallocate(this%bc)
    if (allocated(this%family)) deallocate(this%family)

  end subroutine manager_mesh_free_data_p4est

  !> Destructor for the data in `mesh_manager_t` (base) type.
  subroutine manager_mesh_free_p4est(this)
    class(manager_mesh_p4est_t), intent(inout) :: this

    call this%free_base()

    this%nelv = 0
    this%level_max = 0

    if (allocated(this%level)) deallocate(this%level)
    if (allocated(this%igrp)) deallocate(this%igrp)
    if (allocated(this%crv)) deallocate(this%crv)
    if (allocated(this%bc)) deallocate(this%bc)
    if (allocated(this%family)) deallocate(this%family)

  end subroutine manager_mesh_free_p4est

end module manager_mesh_p4est
