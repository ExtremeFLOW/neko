! Copyright (c) 2019-2026, The Neko Authors
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
!> Mesh topology related objects for SEM
module sem
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use mesh_manager, only : mesh_manager_t
  use mesh, only : mesh_t
  use dofmap, only : dofmap_t
  use gather_scatter, only : gs_t
  use gs_ops, only : GS_OP_ADD
  use time_state, only : time_state_t
  use amr, only : amr_t
  use amr_reconstruct, only : amr_reconstruct_t
  use amr_restart_component, only : amr_restart_component_t



  implicit none
  private

  !> Base type for a boundary condition
  type, public, extends(amr_restart_component_t) :: sem_t
!     type(amr_t) :: amr
   contains
     !> Constructor
     procedure, pass(this) :: init => sem_init
     !> Destructor
     procedure, pass(this) :: free => sem_free
     !> AMR restart
     procedure, pass(this) :: amr_restart => sem_amr_restart
  end type sem_t

contains

  !> Constructor
  !! @param dof Map of degrees of freedom.
  subroutine sem_init(this, mesh_manager, mesh)
    class(sem_t), intent(inout) :: this
    class(mesh_manager_t), intent(in) :: mesh_manager
    type(mesh_t), target, intent(in) :: mesh

!    integer :: lx
!    lx = 8

    call this%free

!    call this%amr%init(mesh_manager%transfer, mesh_manager%isamr, &
!         mesh_manager%mesh%tdim, lx)

    ! SEM has to be a first module in AMR reconstruction list
!    if (this%amr%ifamr()) then
!       call this%amr%comp_add(this, 'SEM')
!    end if

  end subroutine sem_init

  !> Destructor
  subroutine sem_free(this)
    class(sem_t), intent(inout) :: this

!    call this%amr%free()

  end subroutine sem_free

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     time          time state
  subroutine sem_amr_restart(this, reconstruct, counter, time)
    class(sem_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter
    type(time_state_t), intent(in) :: time
    character(len=LOG_SIZE) :: log_buf

    ! Was this component already restarted?
    if (this%counter .eq. counter) return

    this%counter = counter

    log_buf = 'Reconstructing SEM module'
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

  end subroutine sem_amr_restart

end module sem
