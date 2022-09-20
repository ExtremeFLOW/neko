! Copyright (c) 2018-2022, The Neko Authors
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
!> Defines a registry for storing solution fields
!
module field_registry
  use num_types
  use field

  implicit none

  type :: field_registry_t
     type(field_t), allocatable :: fields(:) !< list of fields stored
     integer, private :: n !< number of registered fields
   contains
     procedure, pass(this) :: add_field => add_field
     procedure, pass(this) :: n_fields => n_fields

  end type field_registry_t

  !> Initialise the registry
  interface field_registry_t
     module procedure init
  end interface field_registry_t

  !   interface add_field
  !     module procedure add_field_external_dof
  !   end interface add_field

  !> Global field registry
  type(field_registry_t), public, target :: neko_field_registry

contains
  function init() result(this)
    type(field_registry_t) :: this

    allocate (this%fields(10))
    this%n = 0

    write (*, *) "INIT FIELD REGISTRY"

  end function init

  !  subroutine add_field_external_dof(this, dof, fld_name)
  subroutine add_field(this, dof, fld_name)
    class(field_registry_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof  !< External dofmap for the field
    character(len=*), intent(in) :: fld_name     !< Name of the field
    type(field_t) :: f

    write(*,*) "ADDING FIELD", fld_name

    call field_init(f, dof, fld_name)
    this%n = this%n + 1
    this%fields(this%n) = f
  end subroutine add_field

  pure function n_fields(this) result(n)
    class(field_registry_t), intent(in) :: this
    integer :: n

    n = this%n  
  end function n_fields

end module field_registry
