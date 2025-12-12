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
  use logger, only : neko_log, NEKO_LOG_QUIET, NEKO_LOG_INFO, &
       NEKO_LOG_VERBOSE, NEKO_LOG_DEBUG, LOG_SIZE
  use utils, only : neko_error, neko_warning
  use profiler, only : profiler_start_region, profiler_end_region
  use time_state, only : time_state_t
  use user_intf, only : user_t
  use mesh, only : mesh_t
  use mesh_manager_transfer, only : mesh_manager_transfer_t
  use mesh_manager_transfer_p4est, only : mesh_manager_transfer_p4est_t
  use mesh_manager, only : mesh_manager_t
  use amr_reconstruct, only : amr_reconstruct_t
  use amr_restart_component, only : amr_restart_component_t

  implicit none
  private

  !> Component entrance
  type, public :: amr_component_pointer_t
     class(amr_restart_component_t), pointer :: cmp
  end type amr_component_pointer_t

  !> Main AMR type
  type, public :: amr_t
     !> Data reconstruction type
     type(amr_reconstruct_t) :: reconstruct
     !> Number of components
     integer :: ncomponents
     !> Components list
     class(amr_component_pointer_t), allocatable, dimension(:) :: components
     !> Restart counter
     integer :: counter
   contains
     !> Initialise type
     procedure, pass(this) :: init => amr_init
     !> Free type
     procedure, pass(this) :: free => amr_free
     !> Add restart component
     procedure, pass(this) :: comp_add => amr_component_add
     !> Remove restart component
     procedure, pass(this) :: comp_remove => amr_component_remove
     !> Restart components
     procedure, pass(this) :: restart => amr_restart
     !> Refine/coarsen
     procedure, pass(this) :: refine => amr_refine
  end type amr_t

contains

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
    if (allocated(this%components)) deallocate(this%components)
    this%ncomponents = 0
    this%counter = 0
  end subroutine amr_free

  !> Add restart component
  !! @param[inout]  component     component to be added
  subroutine amr_component_add(this, component)
    class(amr_t), intent(inout) :: this
    class(amr_restart_component_t), target, intent(inout) :: component
    class(amr_component_pointer_t), allocatable, dimension(:) :: tmp
    integer :: il, itmp

    if (.not. allocated(this%components)) then
       this%ncomponents = 1
       allocate(this%components(1))
       this%components(il)%cmp => component
       call component%init_amr_base(this%ncomponents)
    else
       ! check if there is an empty slot
       itmp = 0
       do il = 1, this%ncomponents
          if (.not. associated(this%components(il)%cmp)) then
             itmp = il
             exit
          end if
       end do
       if (itmp .ne. 0) then
          this%components(itmp)%cmp => component
          call component%init_amr_base(itmp)
       else
          allocate(tmp(this%ncomponents + 1))
          do il = 1, this%ncomponents
             tmp(il)%cmp => this%components(il)%cmp
          end do
          this%ncomponents = this%ncomponents + 1
          tmp(this%ncomponents)%cmp => component
          call component%init_amr_base(this%ncomponents)
          deallocate(this%components)
          call MOVE_ALLOC(tmp, this%components)
       end if
    end if

  end subroutine amr_component_add

  !> Remove restart component
  !! @param[inout]  component     component to be removed
  subroutine amr_component_remove(this, component)
    class(amr_t), intent(inout) :: this
    class(amr_restart_component_t), target, intent(inout) :: component
    class(amr_component_pointer_t), allocatable, dimension(:) :: tmp
    integer :: il, itmp

    if (component%listed) then
       if (allocated(this%components)) then
          if (component%lst_pos .lt. this%ncomponents) then
             this%components(component%lst_pos)%cmp => NULL()
          else if (component%lst_pos .eq. this%ncomponents) then
             this%ncomponents = this%ncomponents - 1
             if (this%ncomponents .eq. 0) then
                deallocate(this%components)
             else
                allocate(tmp(this%ncomponents))
                do il = 1, this%ncomponents
                   tmp(il)%cmp => this%components(il)%cmp
                end do
                deallocate(this%components)
                call MOVE_ALLOC(tmp, this%components)
             end if
          else
             call neko_warning('Restart component not listed')
          end if
       else
          call neko_warning('Restart component not listed')
       end if
    end if
    component%listed = .false.
    component%lst_pos = 0

  end subroutine amr_component_remove

  !> Restart components
  !! @param[in]      user          user interface
  subroutine amr_restart(this, user)
    class(amr_t), intent(inout) :: this
    type(user_t), intent(in) :: user
    integer :: il

    ! update restart counter
    this%counter = this%counter + 1

    ! restart components
    if (allocated(this%components)) then
       do il = 1, this%ncomponents
          if (associated(this%components(il)%cmp)) &
               call this%components(il)%cmp%amr_restart(this%reconstruct, &
               this%counter)
       end do
    end if

    ! let user reconstruct fields
    call user%amr_reconstruct(this%reconstruct, this%counter)

  end subroutine amr_restart

  !> Refine/coarsen mesh
  !! @param[inout]   mesh_manager  mesh manager
  !! @param[inout]   mesh          neko mesh type
  !! @param[in]      user          user interface
  !! @param[in]      time          time state
  subroutine amr_refine(this, mesh_manager, mesh, user, time)
    class(amr_t), intent(inout) :: this
    class(mesh_manager_t), intent(inout) :: mesh_manager
    type(mesh_t), intent(inout) :: mesh
    type(user_t), intent(in) :: user
    type(time_state_t), intent(in) :: time
    integer, allocatable, dimension(:) :: ref_mark
    logical :: ifrefine, ifmod
    integer :: nelt
    character(len=LOG_SIZE) :: log_buf

    select type(transfer => this%reconstruct%transfer)
    type is (mesh_manager_transfer_p4est_t)
       nelt = transfer%nelt_neko
    end select

    ! get refinement information
    allocate(ref_mark(nelt))
    call user%amr_refine_flag(time, ref_mark, ifrefine)

    if (ifrefine) then
       call neko_log%section("Mesh refinement")

       call profiler_start_region("Mesh refinement", 30)

       ! Perform p4est refinement/coarsening
       call mesh_manager%refine(ref_mark, ifmod)

       if (ifmod) then
          write(log_buf, '(a)') 'Mesh modified; restarting solver'
          call neko_log%message(log_buf, NEKO_LOG_INFO)

          ! place to reconstruct geometry and correct mesh manager vertex
          ! position
          ! one should do gs here as well

          ! Reconstruct neko mesh
          call mesh_manager%mesh_construct(mesh, .false.)

          ! restart solver
          call this%restart(user)
       else
          write(log_buf, '(a)') 'Mesh not changed'
          call neko_log%message(log_buf, NEKO_LOG_INFO)
       end if

       call profiler_end_region("Mesh refinement", 30)

       call neko_log%end_section()
    end if

    deallocate(ref_mark)

  end subroutine amr_refine

end module amr
