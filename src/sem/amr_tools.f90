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
!> Various tools for AMR
module amr_tools
  use num_types, only : rp, dp
  use utils, only : neko_error
  use vector_list, only : vector_list_t
  use simulation_component, only : simulation_component_t
  use spectral_error, only : spectral_error_t
  use simcomp_executor, only : neko_simcomps
  use mesh, only : mesh_t
  use space, only : space_t, GLL
  use dofmap, only : dofmap_t
  use gather_scatter, only : gs_t
  use field, only : field_t
  use fld_file, only : fld_file_t
  use fld_file_output, only : fld_file_output_t
  use amr_reconstruct, only : amr_reconstruct_t
  use amr_restart_component, only : amr_restart_component_t

  implicit none
  private

  !> Type for spectral error indicator operations
  type, public, extends(amr_restart_component_t) :: amr_spectral_error_t
     !> Pointer to the spectral error indicator simcomp
     type(spectral_error_t), pointer :: simcomp_spectral
     !> Vector list pointing to instantaneous error vectors
     type(vector_list_t) :: eind
     !> Vector list pointing to averaged error vectors
     type(vector_list_t) :: eind_av
     !> Corresponding mesh
     type(mesh_t), pointer :: msh
     !> Minimal space
     type(space_t) :: Xh
     !> Minimal dofmap
     type(dofmap_t) :: dm_Xh
     !> Minimal gather-scatter communicator
     type(gs_t) :: gs_Xh
     !> File output for error indicator
     type(fld_file_output_t) :: eind_av_fld
   contains
     !> Initialise type
     procedure, pass(this) :: init => amr_spectral_error_init
     !> Free type
     procedure, pass(this) :: free => amr_spectral_error_free
     !> AMR restart
     procedure, pass(this) :: amr_restart => amr_spectral_error_amr_restart
  end type amr_spectral_error_t

contains

  !> Initialise spectral error indicator type
  !! @param[in]     msh       Mesh
  subroutine amr_spectral_error_init(this, msh)
    class(amr_spectral_error_t), intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    integer :: il, n_simcomps, lx

    Call this%free()

    ! Identify spectral error indicator simcomp
    n_simcomps = neko_simcomps%get_n()
    do il = 1, n_simcomps
       call simcomp_get_pointer(neko_simcomps%simcomps(il)%simcomp, &
            this%simcomp_spectral)
    end do

    if (.not. associated(this%simcomp_spectral)) call neko_error('AMR tool; &
         &spectral error indicator simcomp not found')

    ! minimal mesh
    this%msh => msh
    lx = 3
    call this%Xh%init(GLL, lx, lx, lx)
    call this%dm_Xh%init(this%msh, this%Xh)
    call this%gs_Xh%init(this%dm_Xh)

  end subroutine amr_spectral_error_init

  ! Subroutine to get simcomps pointers
  subroutine simcomp_get_pointer(simcomp, simcomp_ptr)
    class(simulation_component_t), target, intent(in) :: simcomp
    type(spectral_error_t), pointer, intent(out) :: simcomp_ptr

    nullify(simcomp_ptr)

    select type (simcomp)
       type is(spectral_error_t)
          simcomp_ptr => simcomp
       end select

  end subroutine simcomp_get_pointer

  !> Initialise spectral error indicator type
  subroutine amr_spectral_error_free(this)
    class(amr_spectral_error_t), intent(inout) :: this

    nullify(this%simcomp_spectral)
    nullify(this%msh)

    call this%gs_Xh%free()
    call this%dm_Xh%free()
    call this%Xh%free()

    call this%eind%free()
    call this%eind_av%free()

  end subroutine amr_spectral_error_free

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     tstep         time step
  subroutine amr_spectral_error_amr_restart(this, reconstruct, counter, tstep)
    class(amr_spectral_error_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter, tstep

    ! Reconstruct minimal mesh and communicator
    call this%dm_Xh%amr_restart(reconstruct, counter, tstep)
    call this%gs_Xh%amr_restart(reconstruct, counter, tstep)

  end subroutine amr_spectral_error_amr_restart

end module amr_tools
