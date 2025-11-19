! Copyright (c) 2018-2023, The Neko Authors
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
!> Defines a registry for storing and requesting temporary fields
!! This can be used when you have a function that will be called
!! often and you don't want to create temporary fields (work arrays) inside
!! it on each call.
module scratch_registry
  use field, only : field_t
  use vector, only : vector_t
  use matrix, only : matrix_t

  use field_math, only : field_rzero
  use dofmap, only : dofmap_t
  use utils, only : neko_error
  implicit none
  private

  type :: register_t
     character(len=32) :: type = ''
     logical :: allocated = .false.

     type(field_t), pointer :: field_ptr => null()
     type(vector_t), pointer :: vector_ptr => null()
     type(matrix_t), pointer :: matrix_ptr => null()

   contains
     procedure, pass(this) :: init => init_register_field
     procedure, pass(this) :: free => free_register
  end type register_t

  type, public :: scratch_registry_t
     !> list of scratch fields
     type(register_t), private, allocatable :: registry(:)
     !> Tracks which fields are used
     logical, private, allocatable :: inuse(:)
     !> number of registered fields
     integer, private :: n_available = 0
     !> number of fields in use
     integer, private :: n_inuse = 0
     !> the size the fields array is increased by upon reallocation
     integer, private :: expansion_size = 10
     !> Dofmap
     type(dofmap_t), pointer :: dof => null()
   contains
     procedure, private, pass(this) :: expand
     !> Constructor
     procedure, pass(this) :: init => scratch_registry_init
     !> Destructor
     procedure, pass(this) :: free => scratch_registry_free
     !> Assign a dofmap to the scratch registry
     procedure, pass(this) :: set_dofmap => scratch_registry_set_dofmap
     !> Getter for n_available
     procedure, pass(this) :: get_n_available
     !> Getter for n_inuse
     procedure, pass(this) :: get_n_inuse
     !> Getter for expansion_size
     procedure, pass(this) :: get_expansion_size
     !> Return size of allocated fields
     procedure, pass(this) :: get_size
     !> Get value of inuse for a given index
     procedure, pass(this) :: get_inuse
     !> Get a new scratch field
     procedure, pass(this) :: request_field
     procedure, pass(this) :: relinquish_field_single
     procedure, pass(this) :: relinquish_field_multiple
     !> Free a field for later reuse
     generic :: relinquish_field => relinquish_field_single, &
          relinquish_field_multiple
  end type scratch_registry_t

  !> Global scratch registry
  type(scratch_registry_t), public, target :: neko_scratch_registry

contains

  !> Constructor, optionally taking initial registry and expansion
  !! size as argument
  subroutine scratch_registry_init(this, size, expansion_size, dof)
    class(scratch_registry_t), intent(inout) :: this
    integer, optional, intent(in) :: size
    integer, optional, intent(in) :: expansion_size
    type(dofmap_t), target, intent(in), optional :: dof
    integer :: i, s

    call this%free()

    s = 10
    if (present(size)) s = size
    if (present(dof)) this%dof => dof

    allocate(this%registry(s))
    allocate(this%inuse(s))
    this%inuse(:) = .false.

    this%expansion_size = 10
    if (present(expansion_size)) this%expansion_size = expansion_size

    this%n_available = 0
    this%n_inuse = 0
  end subroutine scratch_registry_init

  !> Destructor
  subroutine scratch_registry_free(this)
    class(scratch_registry_t), intent(inout):: this
    integer :: i

    if (allocated(this%registry)) then
       do i = 1, this%n_available
          call this%registry(i)%free()
       end do

       deallocate(this%registry)
    end if

    if (allocated(this%inuse)) deallocate(this%inuse)
    if (associated(this%dof)) nullify(this%dof)

    ! Reset to default values
    this%n_available = 0
    this%n_inuse = 0
    this%expansion_size = 10

  end subroutine scratch_registry_free

  !> Assign a dofmap to the scratch registry
  subroutine scratch_registry_set_dofmap(this, dof)
    class(scratch_registry_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof

    if (associated(this%dof)) then
       call neko_error("scratch_registry::set_dofmap: "&
            // "Dofmap is already assigned to scratch registry.")
    end if
    this%dof => dof
  end subroutine scratch_registry_set_dofmap

  !> Get the number of fields stored in the registry
  pure function get_n_available(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    n = this%n_available
  end function get_n_available

  pure function get_n_inuse(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n, i

    n = 0
    do i = 1 ,this%get_size()
       if (this%inuse(i)) n = n + 1
    end do
  end function get_n_inuse

  !> Get the size of the fields array
  pure function get_size(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    n = size(this%registry)
  end function get_size

  !> Get the expansion size
  pure function get_expansion_size(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    n = this%expansion_size
  end function get_expansion_size

  subroutine expand(this)
    class(scratch_registry_t), intent(inout) :: this
    type(register_t), allocatable :: temp(:)
    logical, allocatable :: temp2(:)
    integer :: i

    allocate(temp(this%get_size() + this%expansion_size))
    temp(1:this%n_available) = this%registry(1:this%n_available)

    call move_alloc(temp, this%registry)

    allocate(temp2(this%get_size() + this%expansion_size))
    temp2(1:this%n_available) = this%inuse(1:this%n_available)
    temp2(this%n_available+1:) = .false.
    this%inuse = temp2
  end subroutine expand

  !> Get a field from the registry by assigning it to a pointer
  subroutine request_field(this, f, index)
    class(scratch_registry_t), target, intent(inout) :: this
    type(field_t), pointer, intent(inout) :: f
    integer, intent(inout) :: index !< The index of the field in the inuse array
    character(len=10) :: name

    if (.not. associated(this%dof)) then
       call neko_error("scratch_registry::request_field: "&
            // "No dofmap assigned to scratch registry.")
    end if

    associate(n_available => this%n_available, n_inuse => this%n_inuse)

      do index = 1, this%get_size()
         if (.not. this%inuse(index)) then

            if (.not. this%registry(index)%allocated) then
               write(name, "(A3,I0.3)") "wrk", index
               call this%registry(index)%init('field', this%dof, trim(name))
               n_available = n_available + 1
            else if (this%registry(index)%type .ne. 'field') then
               continue
            end if

            f => this%registry(index)%field_ptr
            call field_rzero(f)
            this%inuse(index) = .true.
            this%n_inuse = this%n_inuse + 1
            return
         end if
      end do

      ! all existing fields in use, we need to expand to add a new one
      index = n_available + 1
      call this%expand()
      n_available = n_available + 1
      n_inuse = n_inuse + 1
      this%inuse(n_available) = .true.
      write (name, "(A3,I0.3)") "wrk", index
      call this%registry(n_available)%init('field', this%dof, trim(name))
      f => this%registry(n_available)%field_ptr

    end associate
  end subroutine request_field

  !> Relinquish the use of a field in the registry
  subroutine relinquish_field_single(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index !< The index of the field to free

    if (trim(this%registry(index)%type) .ne. 'field') then
       call neko_error("scratch_registry::relinquish_field_single: " &
            // "Register entry is not a field.")
    end if

    this%inuse(index) = .false.
    this%n_inuse = this%n_inuse - 1
  end subroutine relinquish_field_single

  subroutine relinquish_field_multiple(this, indices)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: indices(:) !< The indices of the field to free
    integer :: i

    do i = 1, size(indices)
       if (trim(this%registry(indices(i))%type) .ne. 'field') then
          call neko_error("scratch_registry::relinquish_field_single: " &
               // "Register entry is not a field.")
       end if

       this%inuse(indices(i)) = .false.
    end do
    this%n_inuse = this%n_inuse - size(indices)
  end subroutine relinquish_field_multiple

  logical function get_inuse(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index !< The index of the field to check

    get_inuse = this%inuse(index)
  end function get_inuse

  ! ========================================================================== !
  ! Single register component routines

  !> Initialize a register entry
  subroutine init_register_field(this, type, dof, name)
    class(register_t), intent(inout) :: this
    character(len=*), intent(in) :: type
    type(dofmap_t), target, intent(in) :: dof
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: err_msg

    if (this%allocated) then
       call neko_error("scratch_registry::init_register_field: "&
            // "Register entry is already allocated.")
    else if (trim(type) .ne. 'field') then
       write(err_msg, '(A,A,A)') 'scratch_registry::init_register_field: ', &
            'Unknown register type: ', trim(type)
       call neko_error(err_msg)
    end if

    call this%free()

    allocate(this%field_ptr)
    call this%field_ptr%init(dof, trim(name))

    this%type = trim(type)
    this%allocated = .true.

  end subroutine init_register_field

  !> Initialize a register entry
  subroutine init_register_vector(this, type, n)
    class(register_t), intent(inout) :: this
    character(len=*), intent(in) :: type
    integer, intent(in) :: n
    character(len=:), allocatable :: err_msg

    if (this%allocated) then
       call neko_error("scratch_registry::init_register_vector: "&
            // "Register entry is already allocated.")
    else if (trim(type) .ne. 'vector') then
       write(err_msg, '(A,A,A)') 'scratch_registry::init_register_vector: ', &
            'Unknown register type: ', trim(type)
       call neko_error(err_msg)
    end if

    call this%free()

    allocate(this%vector_ptr)
    call this%vector_ptr%init(n)

    this%type = trim(type)
    this%allocated = .true.

  end subroutine init_register_vector

  !> Initialize a register entry
  subroutine init_register_matrix(this, type, nrows, ncols)
    class(register_t), intent(inout) :: this
    character(len=*), intent(in) :: type
    integer, intent(in) :: nrows, ncols
    character(len=:), allocatable :: err_msg

    if (this%allocated) then
       call neko_error("scratch_registry::init_register_matrix: "&
            // "Register entry is already allocated.")
    else if (trim(type) .ne. 'matrix') then
       write(err_msg, '(A,A,A)') 'scratch_registry::init_register_matrix: ', &
            'Unknown register type: ', trim(type)
       call neko_error(err_msg)
    end if

    call this%free()

    allocate(this%matrix_ptr)
    call this%matrix_ptr%init(nrows, ncols)

    this%type = trim(type)
    this%allocated = .true.

  end subroutine init_register_matrix

  !> Free a register entry
  subroutine free_register(this)
    class(register_t), intent(inout) :: this

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

    this%allocated = .false.
    this%type = ''
  end subroutine free_register

end module scratch_registry
