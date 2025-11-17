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
!> Implementation of the adaptive mesh refinement workflow
module amr
  use num_types, only : i4, i8, rp, dp
  use mesh_manager_transfer, only : mesh_manager_transfer_t
  use mesh_manager, only : mesh_manager_t
  use amr_reconstruct, only : amr_reconstruct_t

  implicit none
  private

  !> basic type for component restart
  type, abstract, public :: amr_restart_t
     !> Restart counter
     integer :: counter
   contains
     !> Free base type
     procedure, pass(this) :: free_base => amr_restart_free_base
     !> Restart the component
     procedure(amr_restart_type), pass(this), deferred :: restart
     !> Free type
     procedure(amr_restart_free), pass(this), deferred :: free
  end type amr_restart_t

  abstract interface
     !> Restart the component
     !! @param[in]  reconstruct   data reconstruction type
     !! @param[in]  count         restart count
     subroutine amr_restart_type(this, reconstruct, count)
       import amr_restart_t, amr_reconstruct_t
       class(amr_restart_t), intent(inout) :: this
       type(amr_reconstruct_t), intent(in) :: reconstruct
       integer, intent(in) :: count
     end subroutine amr_restart_type

     !> Free type
     subroutine amr_restart_free(this)
       import amr_restart_t
       class(amr_restart_t), intent(inout) :: this
     end subroutine amr_restart_free
  end interface

  !> Component entrance
  type, public :: amr_component_t
     class(amr_restart_t), allocatable :: cmp
  end type amr_component_t

  !> Main AMR type
  type, public :: amr_t
     !> Data reconstruction type
     type(amr_reconstruct_t) :: reconstruct
     !> Number of components
     integer :: ncomponents
     !> Components list
     class(amr_component_t), allocatable, dimension(:) :: components
     !> Restart counter
     integer :: counter
   contains
     !> Initialise type
     procedure, pass(this) :: init => amr_init
     !> Free type
     procedure, pass(this) :: free => amr_free
     !> Restart components
     procedure, pass(this) :: restart => amr_restart
     !> Add restart component
     procedure, pass(this) :: component_add => amr_component_add
  end type amr_t


contains

  !> Free base restart type
  subroutine amr_restart_free_base(this)
    class(amr_restart_t), intent(inout) :: this

    this%counter = 0
  end subroutine amr_restart_free_base

  !> Initialise amr type
  !! @param[in]  transfer     mesh manager data transfer type
  subroutine amr_init(this, transfer)
    class(amr_t), intent(inout) :: this
    class(mesh_manager_transfer_t), intent(in) :: transfer

    call this%free()

    call this%reconstruct%init(transfer)

  end subroutine amr_init

  !> Free amr type
  subroutine amr_free(this)
    class(amr_t), intent(inout) :: this
    integer :: il

    call this%reconstruct%free()
    if (allocated(this%components)) then
       do il = 1, size(this%components)
          if (allocated(this%components(il)%cmp)) then
             call this%components(il)%cmp%free()
             deallocate(this%components(il)%cmp)
          end if
       end do
       deallocate(this%components)
    end if
    this%ncomponents = 0
    this%counter = 0
  end subroutine amr_free

  !> Restart components
  subroutine amr_restart(this)
    class(amr_t), intent(inout) :: this
    integer :: il

    ! update restart counter
    this%counter = this%counter + 1

    ! restart components
    do il = 1, this%ncomponents
       call this%components(il)%cmp%restart(this%reconstruct, this%counter)
    end do

  end subroutine amr_restart

  !> Add restart component
  !! @param[inout]  component     component to be restarted
  subroutine amr_component_add(this, component)
    class(amr_t), intent(inout) :: this
    class(amr_restart_t), allocatable, intent(inout) :: component
    class(amr_component_t), allocatable, dimension(:) :: tmp
    integer :: il

    if (allocated(component)) then
       allocate(tmp(this%ncomponents + 1))
       do il = 1, this%ncomponents
          call MOVE_ALLOC(this%components(il)%cmp, tmp(il)%cmp)
       end do
       this%ncomponents = this%ncomponents + 1
       call MOVE_ALLOC(component, tmp(this%ncomponents)%cmp)
       deallocate(this%components)
       call MOVE_ALLOC(tmp, this%components)
    end if

  end subroutine amr_component_add

end module amr
