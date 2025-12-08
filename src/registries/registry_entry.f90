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
!! This is used in the registries to store temporary fields, vectors,
!! and matrices.
module registry_entry
  use num_types, only : rp
  use field, only : field_t
  use vector, only : vector_t
  use matrix, only : matrix_t

  use dofmap, only : dofmap_t
  use utils, only : neko_error
  implicit none
  private

  type, public :: registry_entry_t
     !> Name of the registry entry
     character(len=:), private, allocatable :: name
     !> Type of the registry entry; must be supproted.
     character(len=:), private, allocatable :: type
     !> Whether the entry is allocated
     logical, private :: allocated = .false.

     !> Storage. Only one of these will be allocated at a time.
     type(field_t), private, pointer :: field_ptr => null()
     type(vector_t), private, pointer :: vector_ptr => null()
     type(matrix_t), private, pointer :: matrix_ptr => null()
     real(kind=rp), private :: scalar = 0.0_rp

   contains
     !> Constructors
     procedure, pass(this) :: init_field => init_register_field
     procedure, pass(this) :: init_vector => init_register_vector
     procedure, pass(this) :: init_matrix => init_register_matrix
     procedure, pass(this) :: init_scalar => init_register_scalar
     !> Destructor
     procedure, pass(this) :: free => free_register

     !> Getters that return a pointer to the object in the entry.
     procedure, pass(this) :: get_name
     procedure, pass(this) :: get_type
     procedure, pass(this) :: get_field
     procedure, pass(this) :: get_vector
     procedure, pass(this) :: get_matrix
     procedure, pass(this) :: get_scalar
     procedure, pass(this) :: is_allocated
  end type registry_entry_t

contains

!> Initialize a register entry
  subroutine init_register_field(this, dof, name)
    class(registry_entry_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof
    character(len=*), intent(in) :: name

    if (this%allocated) then
       call neko_error("init_register_field: " &
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
       call neko_error("init_register_vector: " &
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
       call neko_error("init_register_matrix: " &
            // "Register entry is already allocated.")
    end if

    call this%free()

    allocate(this%matrix_ptr)
    call this%matrix_ptr%init(nrows, ncols)

    if (present(name)) this%name = trim(name)
    this%type = 'matrix'
    this%allocated = .true.

  end subroutine init_register_matrix

  !> Initialize a scalar register entry
  subroutine init_register_scalar(this, val, name)
    class(registry_entry_t), intent(inout) :: this
    real(kind=rp), intent(in) :: val
    character(len=*), optional, intent(in) :: name

    if (this%allocated) then
       call neko_error("init_register_scalar: " &
            // "Register entry is already allocated.")
    end if

    call this%free()

    this%scalar = val

    if (present(name)) this%name = trim(name)
    this%type = 'scalar'
    this%allocated = .true.

  end subroutine init_register_scalar

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

    this%scalar = 0.0_rp

    if (allocated(this%name)) deallocate(this%name)
    if (allocated(this%type)) deallocate(this%type)
    this%allocated = .false.

  end subroutine free_register

  !> Get the name of the registry entry
  pure function get_name(this) result(name)
    class(registry_entry_t), intent(in) :: this
    character(len=:), allocatable :: name
    name = trim(this%name)
  end function get_name

  !> Get the type of the registry entry
  pure function get_type(this) result(type)
    class(registry_entry_t), intent(in) :: this
    character(len=:), allocatable :: type
    type = trim(this%type)
  end function get_type

  !> Check if the registry entry is allocated
  pure function is_allocated(this) result(allocated)
    class(registry_entry_t), intent(in) :: this
    logical :: allocated
    allocated = this%allocated
  end function is_allocated

  !> Get the field pointer of the registry entry
  function get_field(this) result(field_ptr)
    class(registry_entry_t), target, intent(in) :: this
    type(field_t), pointer :: field_ptr
    if (this%get_type() .ne. 'field') then
       call neko_error("registry_entry::get_field: " &
            // "Registry entry is not of type 'field'.")
    end if
    field_ptr => this%field_ptr
  end function get_field

  !> Get the vector pointer of the registry entry
  function get_vector(this) result(vector_ptr)
    class(registry_entry_t), target, intent(in) :: this
    type(vector_t), pointer :: vector_ptr
    if (this%get_type() .ne. 'vector') then
       call neko_error("registry_entry::get_vector: " &
            // "Registry entry is not of type 'vector'.")
    end if
    vector_ptr => this%vector_ptr
  end function get_vector

  !> Get the matrix pointer of the registry entry
  function get_matrix(this) result(matrix_ptr)
    class(registry_entry_t), target, intent(in) :: this
    type(matrix_t), pointer :: matrix_ptr
    if (this%get_type() .ne. 'matrix') then
       call neko_error("registry_entry::get_field: " &
            // "Registry entry is not of type 'matrix'.")
    end if
    matrix_ptr => this%matrix_ptr
  end function get_matrix

  !> Get the scalar pointer of the registry entry
  function get_scalar(this) result(scalar_ptr)
    class(registry_entry_t), target, intent(in) :: this
    real(kind=rp), pointer :: scalar_ptr
    if (this%get_type() .ne. 'scalar') then
       call neko_error("registry_entry::get_field: " &
            // "Registry entry is not of type 'scalar'.")
    end if
    scalar_ptr => this%scalar
  end function get_scalar

end module registry_entry
