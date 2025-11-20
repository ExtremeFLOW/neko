! Copyright (c) 2025, The Neko Authors
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
!> Defines a registry entry for storing and requesting temporary objects
!! This is used in the scratch registries to store temporary fields, vectors,
!! and matrices.
module registry_entry
  use field, only : field_t
  use vector, only : vector_t
  use matrix, only : matrix_t

  use dofmap, only : dofmap_t
  use utils, only : neko_error
  implicit none
  private

  type, public :: registry_entry_t
     character(len=:), allocatable :: name
     character(len=:), allocatable :: type
     logical :: allocated = .false.

     type(field_t), pointer :: field_ptr => null()
     type(vector_t), pointer :: vector_ptr => null()
     type(matrix_t), pointer :: matrix_ptr => null()

   contains
     procedure, pass(this) :: init_field => init_register_field
     procedure, pass(this) :: init_vector => init_register_vector
     procedure, pass(this) :: init_matrix => init_register_matrix
     procedure, pass(this) :: free => free_register
  end type registry_entry_t

contains

!> Initialize a register entry
  subroutine init_register_field(this, dof, name)
    class(registry_entry_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof
    character(len=*), intent(in) :: name

    if (this%allocated) then
       call neko_error("scratch_registry::init_register_field: "&
            // "Register entry is already allocated.")
    end if

    call this%free()

    allocate(this%field_ptr)
    call this%field_ptr%init(dof, trim(name))

    this%name = trim(name)
    this%type = 'field'
    this%allocated = .true.

  end subroutine init_register_field

  !> Initialize a register entry
  subroutine init_register_vector(this, n, name)
    class(registry_entry_t), intent(inout) :: this
    integer, intent(in) :: n
    character(len=*), optional, intent(in) :: name

    if (this%allocated) then
       call neko_error("scratch_registry::init_register_vector: "&
            // "Register entry is already allocated.")
    end if

    call this%free()

    allocate(this%vector_ptr)
    call this%vector_ptr%init(n)

    if (present(name)) this%name = trim(name)
    this%type = 'vector'
    this%allocated = .true.

  end subroutine init_register_vector

  !> Initialize a register entry
  subroutine init_register_matrix(this, nrows, ncols, name)
    class(registry_entry_t), intent(inout) :: this
    integer, intent(in) :: nrows, ncols
    character(len=*), optional, intent(in) :: name

    if (this%allocated) then
       call neko_error("scratch_registry::init_register_matrix: "&
            // "Register entry is already allocated.")
    end if

    call this%free()

    allocate(this%matrix_ptr)
    call this%matrix_ptr%init(nrows, ncols)

    if (present(name)) this%name = trim(name)
    this%type = 'matrix'
    this%allocated = .true.

  end subroutine init_register_matrix

  !> Free a register entry
  subroutine free_register(this)
    class(registry_entry_t), intent(inout) :: this

    if (associated(this%field_ptr)) then
       call this%field_ptr%free()
       deallocate(this%field_ptr)
    end if

    if (associated(this%vector_ptr)) then
       call this%vector_ptr%free()
       deallocate(this%vector_ptr)
    end if

    if (associated(this%matrix_ptr)) then
       call this%matrix_ptr%free()
       deallocate(this%matrix_ptr)
    end if

    if (allocated(this%name)) deallocate(this%name)
    if (allocated(this%type)) deallocate(this%type)
    this%allocated = .false.

  end subroutine free_register

end module registry_entry
