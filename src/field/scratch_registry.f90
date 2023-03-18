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
  use num_types
  use field
  use utils
  implicit none
  private
  
  type :: field_ptr_t
    type(field_t), pointer :: field => NULL()
  end type field_ptr_t
  
  type, public :: scratch_registry_t
     type(field_ptr_t), private, allocatable :: fields(:)      !< list of scratch fields 
     logical, private, allocatable :: inuse(:)                 !< Tracks which fields are used
     integer, private :: nfields                      !< number of registered fields
     integer, private :: nfields_inuse                !< number of fields in use
     integer, private :: expansion_size               !< the size the fields array is increased by upon reallocation
     type(dofmap_t), pointer :: dof                   !< Dofmap
   contains
     procedure, private, pass(this) :: expand
     procedure, pass(this) :: free => scratch_registry_free !< destructor
     procedure, pass(this) :: get_nfields                   !< getter for nfields  
     procedure, pass(this) :: get_nfields_inuse             !< getter for nfields_inuse 
     procedure, pass(this) :: get_expansion_size            !< getter for expansion_size
     procedure, pass(this) :: get_size                      !< return size of allocated fields
     procedure, pass(this) :: get_inuse                     !< get value of inuse for a given index
     procedure, pass(this) :: request_field                 !< get a new scratch field
     procedure, pass(this) :: relinquish_field              !< free a field for later reuse
  end type scratch_registry_t

  interface scratch_registry_t
     procedure :: init
  end interface scratch_registry_t

  !> Global scratch registry
  type(scratch_registry_t), public, target :: neko_scratch_registry

contains
  !> Constructor, optionally taking initial registry and expansion
  !> size as argument
  type(scratch_registry_t) function init(dof, size, expansion_size) result(this)
    type(dofmap_t), target, intent(in) :: dof
    integer, optional, intent(in) :: size
    integer, optional, intent(in) :: expansion_size
    integer :: i
    
    this%dof => dof

    if (present(size)) then
       allocate (this%fields(size))
       do i= 1, size
         allocate(this%fields(i)%field)
       end do
       allocate (this%inuse(size))
    else
       allocate (this%fields(10))
       allocate (this%inuse(10))
    end if

    this%inuse(:) = .false.
    if (present(expansion_size)) then
       this%expansion_size = expansion_size
    else
       this%expansion_size = 10
    end if

    this%nfields = 0
    this%nfields_inuse = 0
  end function init

  !> Destructor
  subroutine scratch_registry_free(this)
    class(scratch_registry_t), intent(inout):: this
    integer :: i

    do i=1, this%nfields
       call field_free(this%fields(i)%field)
       deallocate(this%fields(i)%field)
    end do
    deallocate(this%fields)
    deallocate(this%inuse)
  end subroutine scratch_registry_free


  !> Get the number of fields stored in the registry
  pure function get_nfields(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    n = this%nfields  
  end function get_nfields

  pure function get_nfields_inuse(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n, i

    n = 0
    do i=1,this%get_size()
      if (this%inuse(i)) n = n + 1
    end do
  end function get_nfields_inuse

  !> Get the size of the fields array
  pure function get_size(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    n = size(this%fields)
  end function get_size

  !> Get the expansion size
  pure function get_expansion_size(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    n = this%expansion_size
  end function get_expansion_size

  subroutine expand(this)
    class(scratch_registry_t), intent(inout) :: this
    type(field_ptr_t), allocatable :: temp(:)  
    logical, allocatable :: temp2(:)  
    integer :: i

    allocate(temp(this%get_size() + this%expansion_size))
    temp(1:this%nfields) = this%fields(1:this%nfields)
    
    do i=this%nfields +1, size(temp)
      allocate(temp(i)%field)
    enddo

    call move_alloc(temp, this%fields)

    allocate(temp2(this%get_size() + this%expansion_size))
    temp2(1:this%nfields) = this%inuse(1:this%nfields)
    temp2(this%nfields+1:) = .false.
    this%inuse = temp2
  end subroutine expand


  !> Get a field from the registry by assigning it to a pointer
  subroutine request_field(this, f, index)
    class(scratch_registry_t), target, intent(inout) :: this
    type(field_t), pointer, intent(inout) :: f
    integer, intent(inout) :: index !< The index of the field in the inuse array
    character(len=10) :: name

    
    associate(nfields => this%nfields, nfields_inuse => this%nfields_inuse) 
    
    do index=1,this%get_size()
       if (this%inuse(index) .eqv. .false.) then
         write (name, "(A3,I0.3)") "wrk", index
         
         if (.not. allocated(this%fields(index)%field%x)) then
           call field_init(this%fields(index)%field, this%dof, trim(name))
           nfields = nfields + 1
         end if
         f => this%fields(index)%field
         this%inuse(index) = .true.
         this%nfields_inuse = this%nfields_inuse + 1
         return
       end if
    end do
    ! all existing fields in use, we need to expand to add a new one
    index = nfields +1
    call this%expand()
    nfields = nfields + 1
    nfields_inuse = nfields_inuse + 1
    this%inuse(nfields) = .true.
    write (name, "(A3,I0.3)") "wrk", index
    call field_init(this%fields(nfields)%field, this%dof, trim(name))
    f => this%fields(nfields)%field

    end associate
  end subroutine request_field
  
  !> Relinquish the use of a field in the registry
  subroutine relinquish_field(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index !< The index of the field to free
    
    this%inuse(index) = .false.
    this%nfields_inuse = this%nfields_inuse - 1
  end subroutine relinquish_field

  logical function get_inuse(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index !< The index of the field to check 
    
    get_inuse = this%inuse(index)
  end function get_inuse

end module scratch_registry
