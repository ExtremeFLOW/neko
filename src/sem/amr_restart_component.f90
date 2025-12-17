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
!> Abstract type for AMR restart component
module amr_restart_component
  use amr_reconstruct, only : amr_reconstruct_t

  implicit none
  private

  !> basic type for AMR restart component
  type, abstract, public :: amr_restart_component_t
     !> Is a component listed in amr type
     logical :: listed = .false.
     !> Position in the component list
     integer :: lst_pos
     !> Component name
     character(len=:), allocatable :: cmp_name
     !> Restart counter
     integer :: counter
   contains
     !> Initialise base type
     procedure, pass(this) :: init_amr_base => amr_restart_init_base
     !> Free base type
     procedure, pass(this) :: free_amr_base => amr_restart_free_base
     !> Restart the component
     procedure(amr_restart_comp), pass(this), deferred :: amr_restart
  end type amr_restart_component_t

  abstract interface
     !> Restart the component
     !! @param[in]  reconstruct   data reconstruction type
     !! @param[in]  counter       restart counter
     subroutine amr_restart_comp(this, reconstruct, counter)
       import amr_restart_component_t, amr_reconstruct_t
       class(amr_restart_component_t), intent(inout) :: this
       type(amr_reconstruct_t), intent(in) :: reconstruct
       integer, intent(in) :: counter
     end subroutine amr_restart_comp
  end interface

contains

  !> Initialise base restart component type
  !! @param[in]  lst_pos   position in the component list
  !! @param[in]  name      component name
  subroutine amr_restart_init_base(this, lst_pos, name)
    class(amr_restart_component_t), intent(inout) :: this
    integer, intent(in) :: lst_pos
    character(len=*), intent(in) :: name

    call this%free_amr_base()

    this%listed = .true.
    this%lst_pos = lst_pos
    this%cmp_name = trim(name)

  end subroutine amr_restart_init_base

  !> Free base restart component type
  subroutine amr_restart_free_base(this)
    class(amr_restart_component_t), intent(inout) :: this

    this%listed = .false.
    this%lst_pos = 0
    if (allocated(this%cmp_name)) deallocate(this%cmp_name)
    this%counter = 0
  end subroutine amr_restart_free_base

end module amr_restart_component
