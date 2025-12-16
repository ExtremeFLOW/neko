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
!> Implementation of the `mesh_manager_t` type
module mesh_manager
  use num_types, only : i4, i8, rp, dp
  use json_module, only : json_file
  use nmsh, only: nmsh_mesh_t
  use mesh, only : mesh_t
  use manager_mesh, only : manager_mesh_t
  use mesh_manager_transfer, only : mesh_manager_transfer_t

  implicit none
  private

  !> Base abstract type for mesh manager.
  type, abstract, public :: mesh_manager_t
     !> Manager type name
     character(len=:), allocatable :: type_name
     !> 3rd-party software activation flag
     logical :: ifstarted
     !> AMR execution flag
     logical :: isamr
     !> mesh information
     class(manager_mesh_t), allocatable :: mesh
     !> raw mesh data from the nmsh file
     type(nmsh_mesh_t) :: nmsh_mesh
     !> data redistribution routines
     class(mesh_manager_transfer_t), allocatable :: transfer
   contains
     !> Constructor for the mesh_manager_t (base) type.
     procedure, pass(this) :: init_base => mesh_manager_init_base
     !> Free mesh manager data
     procedure, pass(this) :: free_data_base => mesh_manager_free_data_base
     !> Free mesh manager type.
     procedure, pass(this) :: free_base => mesh_manager_free_base
     !> Get element distribution to mesh file reader
     procedure, pass(this) :: elm_dst_copy => mesh_manager_elm_dst_copy
     !> Start 3rd-party software (if needed)
     procedure(mesh_manager_start), pass(this), deferred :: start
     !> Stop 3rd-party software (if needed)
     procedure(mesh_manager_free), pass(this), deferred :: stop
     !> The common constructor using a JSON object.
     procedure(mesh_manager_init), pass(this), deferred :: init
     !> Destructor.
     procedure(mesh_manager_free), pass(this), deferred :: free
     !> Import mesh data into current type
     procedure(mesh_manager_import), pass(this), deferred :: import
     !> Apply data read from mesh file to mesh manager structures
     procedure(mesh_manager_free), pass(this), deferred :: mesh_file_apply
     !> Perform refinement/coarsening on the mesh manager side
     procedure(mesh_manager_refine), pass(this), deferred :: refine
     !> Construct neko mesh type based on mesh manager data
     procedure(mesh_manager_mesh), pass(this), deferred :: mesh_construct
  end type mesh_manager_t

  abstract interface
     !> Start 3rd-party software (if needed)
     !! @param[out]  ierr  error flag
     subroutine mesh_manager_start(this, json, ierr)
       import mesh_manager_t, json_file
       class(mesh_manager_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       integer, intent(out) :: ierr
     end subroutine mesh_manager_start

     !> The common constructor using a JSON object.
     !! @param json       The JSON object for the mesh manager.
     !! @param type_name  Manager type name
     subroutine mesh_manager_init(this, json)
       import mesh_manager_t, json_file
       class(mesh_manager_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
     end subroutine mesh_manager_init

     !> Destructor, 3rd-party stopping
     subroutine mesh_manager_free(this)
       import mesh_manager_t
       class(mesh_manager_t), intent(inout) :: this
     end subroutine mesh_manager_free

     !> Importing mesh data
     subroutine mesh_manager_import(this, ifcomplete)
       import mesh_manager_t
       class(mesh_manager_t), intent(inout) :: this
       logical, intent(in) :: ifcomplete
     end subroutine mesh_manager_import

     !> Perform refinement/coarsening on the mesh manager side
     !! @param[in]   ref_mark     refinement flag
     !! @param[out]  ifmod        mesh modification flag
     subroutine mesh_manager_refine(this, ref_mark, ifmod)
       import mesh_manager_t, i4
       class(mesh_manager_t), intent(inout) :: this
       integer(i4), dimension(:), intent(in) :: ref_mark
       logical, intent(out) :: ifmod
     end subroutine mesh_manager_refine

     !> Construct neko mesh type based on mesh manager data
     !! @param[inout]   mesh     neko mesh type
     !! @param[in]      ifnmsh   use curvature information form nmsh file
     subroutine mesh_manager_mesh(this, mesh, ifnmsh)
       import mesh_manager_t, mesh_t
       class(mesh_manager_t), intent(in) :: this
       type(mesh_t), intent(inout) :: mesh
       logical, intent(in) :: ifnmsh
     end subroutine mesh_manager_mesh
  end interface

  interface
     !> Mesh manager factory. Both constructs and initializes the object.
     !! @param[inout]  object         The object to be initialised.
     !! @param[inout]  json           JSON object initialising the mesh manager.
     !! @param[in]     ifpartition    partitioning flag
     module subroutine mesh_manager_factory(object, json, ifpartition)
       class(mesh_manager_t), allocatable, intent(inout) :: object
       type(json_file), intent(inout) :: json
       logical, intent(in) :: ifpartition
     end subroutine mesh_manager_factory
  end interface

  interface
     !> Mesh manager allocator.
     !! @param[inout]  object      The object to be allocated.
     !! @param[in]     type_name   The name of the type to allocate.
     module subroutine mesh_manager_allocator(object, type_name)
       class(mesh_manager_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: type_name
     end subroutine mesh_manager_allocator
  end interface

  public :: mesh_manager_factory

contains

  !> Constructor for the `mesh_manager_t` (base) type.
  subroutine mesh_manager_init_base(this, type_name)
    class(mesh_manager_t), intent(inout) :: this
    character(len=*), intent(in) :: type_name

    this%type_name = trim(type_name)
    this%ifstarted = .false.
    this%isamr = .false.

  end subroutine mesh_manager_init_base

  !> Free mesh manager data
  subroutine mesh_manager_free_data_base(this)
    class(mesh_manager_t), intent(inout) :: this

    if (allocated(this%mesh)) call this%mesh%free()
    call this%nmsh_mesh%free()
    if (allocated(this%transfer)) call this%transfer%free()

    this%ifstarted = .false.

    if (allocated(this%type_name)) deallocate(this%type_name)

  end subroutine mesh_manager_free_data_base

  !> Free mesh manager type
  subroutine mesh_manager_free_base(this)
    class(mesh_manager_t), intent(inout) :: this

    call this%free_data_base()
    if (allocated(this%mesh)) deallocate(this%mesh)
    if (allocated(this%transfer)) deallocate(this%transfer)

  end subroutine mesh_manager_free_base

  !> Get element distribution to mesh file reader
  subroutine mesh_manager_elm_dst_copy(this)
    class(mesh_manager_t), intent(inout) :: this

    this%nmsh_mesh%gdim = this%mesh%tdim
    this%nmsh_mesh%nelt = this%mesh%nelt
    this%nmsh_mesh%gnelt = this%mesh%gnelt
    this%nmsh_mesh%offset_el = this%mesh%gnelto

  end subroutine mesh_manager_elm_dst_copy

end module mesh_manager
