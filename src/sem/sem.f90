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
!> Mesh related services provided to SEM solver
module sem
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use utils, only : neko_error
  use mesh_manager, only : mesh_manager_t
  use mesh, only : mesh_t
  use space, only : space_t, GLL
  use dofmap, only : dofmap_t
  use gather_scatter, only : gs_t
  use gs_ops, only : GS_OP_ADD
  use coefs, only : coef_t
  use time_state, only : time_state_t
  use amr_reconstruct, only : amr_reconstruct_t
  use amr_restart_component, only : amr_restart_component_t

  implicit none
  private

  !> Basic objects defining SEM discretization for given polynomial order
  type, public, extends(amr_restart_component_t) :: sem_lx_t
     !> Mesh
     type(mesh_t), pointer :: msh => null()
     !> Number of grid point in 1D (polynomial order + 1)
     integer :: lx
     !> Function space \f$ X_h \f$
     type(space_t) :: Xh
     !> Dofmap associated with \f$ X_h \f$
     type(dofmap_t) :: dm_Xh
     !> Gather-scatter associated with \f$ X_h \f$
     type(gs_t) :: gs_Xh
     !> Geometrical coefficient flag
     logical :: ifcoef = .false.
     !> Do we allocate all coefficients
     logical :: ifcoef_all = .false.
     !> Coefficients associated with \f$ X_h \f$
     type(coef_t) :: c_Xh
   contains
     !> Constructor
     procedure, pass(this) :: init => sem_lx_init
     !> Destructor
     procedure, pass(this) :: free => sem_lx_free
     !> Add empty geometrical coefficient to existing type
     procedure, pass(this) :: coef_empty_add => sem_coef_empty_add
     !> Add all geometrical coefficients to existing type
     procedure, pass(this) :: coef_all_add => sem_coef_all_add
     !> AMR restart
     procedure, pass(this) :: amr_restart => sem_lx_amr_restart
  end type sem_lx_t

  !> Basic objects defining SEM discretisation
  type, public, extends(amr_restart_component_t) :: sem_t
     !> AMR execution flag
     logical :: isamr
     !> Number of grid point in 1D (polynomial order + 1)
     integer :: lx
     !> Mesh
     type(mesh_t), pointer :: msh => null()
     !> Mesh manager
     class(mesh_manager_t), pointer :: msh_mng => null()
     !> AMR data reconstruction
     type(amr_reconstruct_t) :: amr_reconstruct
     ! Minimal grid
     type(sem_lx_t) :: grid_min
   contains
     !> Constructor
     procedure, pass(this) :: init => sem_init
     !> Finalise initialisation
     procedure, pass(this) :: finalise => sem_finalise
     !> Destructor
     procedure, pass(this) :: free => sem_free
     !> AMR restart
     procedure, pass(this) :: amr_restart => sem_amr_restart
  end type sem_t

contains

  !> Constructor
  !! @param[in]  msh          mesh
  !! @param[in]  lx           polynomial order + 1
  !! @param[in]  ifcoef       Are geometrical coefficients initialised
  !! @param[in]  ifcoef_all   Are all geometrical coefficients allocated
  subroutine sem_lx_init(this, msh, lx, ifcoef, ifcoef_all)
    class(sem_lx_t), intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    logical, intent(in) :: ifcoef, ifcoef_all

    call this%free()

    this%msh => msh

    this%lx = lx
    this%ifcoef = ifcoef
    this%ifcoef_all = ifcoef_all

    if (msh%gdim .eq. 2) then
       call this%Xh%init(GLL, lx, lx)
    else
       call this%Xh%init(GLL, lx, lx, lx)
    end if

    call this%dm_Xh%init(msh, this%Xh)
    call this%gs_Xh%init(this%dm_Xh)

    if (this%ifcoef) then
       if (this%ifcoef_all) then
          call this%c_Xh%init(this%gs_Xh)
       else
          call this%c_Xh%init(this%Xh, this%msh)
       end if
    end if

  end subroutine sem_lx_init

  !> Destructor
  subroutine sem_lx_free(this)
    class(sem_lx_t), intent(inout) :: this

    nullify(this%msh)

    this%lx = 0
    this%ifcoef = .false.
    this%ifcoef_all = .false.


    call this%c_Xh%free()
    call this%gs_Xh%free()
    call this%dm_Xh%free()
    call this%Xh%free()

    call this%free_amr_base()

  end subroutine sem_lx_free

  !> Add empty geometrical coefficient to existing type
  subroutine sem_coef_empty_add(this)
    class(sem_lx_t), intent(inout) :: this

    call this%c_Xh%free()

    this%ifcoef = .true.
    this%ifcoef_all = .false.
    call this%c_Xh%init(this%Xh, this%msh)

  end subroutine sem_coef_empty_add

  !> Add all geometrical coefficients to existing type
  subroutine sem_coef_all_add(this)
    class(sem_lx_t), intent(inout) :: this

    call this%c_Xh%free()

    this%ifcoef = .true.
    this%ifcoef_all = .true.
    call this%c_Xh%init(this%gs_Xh)

  end subroutine sem_coef_all_add

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     time          time state
  subroutine sem_lx_amr_restart(this, reconstruct, counter, time)
    class(sem_lx_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter
    type(time_state_t), intent(in) :: time
    character(len=LOG_SIZE) :: log_buf

    ! Was this component already restarted?
    if (this%counter .eq. counter) return

    this%counter = counter

    if (this%Xh%lx .lt. 1e1) then
       write(log_buf, '(A,I2)') 'Reconstructing SEM_lx; lx =  ', this%Xh%lx
    else if (this%Xh%lx .lt. 1e2) then
       write(log_buf, '(A,I2)') 'Reconstructing SEM_lx; lx =  ', this%Xh%lx
    end if
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

    ! reconstruct dofmap
    call this%dm_Xh%amr_restart(reconstruct, counter, time)

    ! reconstruct gs
    call this%gs_Xh%amr_restart(reconstruct, counter, time)

    if (this%ifcoef) call this%c_Xh%amr_restart(reconstruct, counter, time)

  end subroutine sem_lx_amr_restart

  !> Constructor
  !! @param[inout]   msh             mesh
  !! @param[in]      mesh_manager    mesh manager
  !! @param[in]      lx              polynomial order + 1
  subroutine sem_init(this, msh, mesh_manager, lx)
    class(sem_t), intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    class(mesh_manager_t), optional, target, intent(in) :: mesh_manager
    integer, optional, intent(in) :: lx

    call this%free()
    this%msh => msh

    ! for AMR reconstruction tool
    if (present(mesh_manager)) then
       this%msh_mng => mesh_manager
       this%isamr = mesh_manager%isamr
       if (.not. present(lx)) then
          call neko_error('SEM initialisation; missing lx')
       else
          this%lx = lx
       end if
    end if

    ! Add minimal grid
    call this%grid_min%init(msh, 3, .false., .false.)

  end subroutine sem_init

  !> Finalise initialisation
  !! @param[inout]   msh             mesh
  !! @param[in]      mesh_manager    mesh manager
  !! @param[in]      lx              polynomial order + 1
  subroutine sem_finalise(this)
    class(sem_t), intent(inout) :: this

    ! initialise AMR reconstruction tool
    if (this%isamr) call this%amr_reconstruct%init(this%msh_mng%transfer, &
         this%msh%gdim, this%lx)

  end subroutine sem_finalise

  !> Destructor
  subroutine sem_free(this)
    class(sem_t), intent(inout) :: this

    this%isamr = .false.
    this%lx = 0

    nullify(this%msh)
    nullify(this%msh_mng)

    call this%amr_reconstruct%free()

    call this%grid_min%free()

    call this%free_amr_base()

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
    call neko_log%section(log_buf, NEKO_LOG_VERBOSE)

    ! reconstruct minimal grid
    call this%grid_min%amr_restart(reconstruct, counter, time)

    call neko_log%end_section(lvl = NEKO_LOG_VERBOSE)

  end subroutine sem_amr_restart

end module sem
