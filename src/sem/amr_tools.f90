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
  use utils, only : NEKO_VARNAME_LEN, neko_error
  use vector_list, only : vector_list_t
  use simulation_component, only : simulation_component_t
  use spectral_error, only : spectral_error_t
  use simcomp_executor, only : neko_simcomps
  use mesh, only : mesh_t
  use space, only : space_t, GLL
  use dofmap, only : dofmap_t
  use gather_scatter, only : gs_t
  use field, only : field_t
  use field_array, only : field_array_t
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
     !> Number of error indicator fields
     integer :: nfld_smp
     !> Vector length
     integer :: nelv
     !> Corresponding mesh
     type(mesh_t), pointer :: msh
     !> Minimal space
     type(space_t) :: Xh
     !> Minimal dofmap
     type(dofmap_t) :: dm_Xh
     !> Minimal gather-scatter communicator
     type(gs_t) :: gs_Xh
     !> Field array for saving instantaneous error indicators
     type(field_array_t) :: eind_in_fld
     !> File output for instantaneous error indicator
     type(fld_file_output_t) :: eind_in_file
     !> Field array for saving averaged error indicators
     type(field_array_t) :: eind_av_fld
     !> File output for averaged error indicator
     type(fld_file_output_t) :: eind_av_file
     !> Number of fields for refinement mark assessment
     integer :: nfld_ref
     !> Field names used for refinement mark assessment
     character(NEKO_VARNAME_LEN), dimension(:), allocatable :: field_names_ref
     !> Vector list pointing to averaged error vectors used to get refinement
     !! mark
     type(vector_list_t) :: eind_av_ref
     !> Operator applied to fields
     character(:), allocatable :: opr
     !> Number of fields for stability monitoring
     integer :: nfld_mntr
     !> Field names used for refinement mark assessment
     character(NEKO_VARNAME_LEN), dimension(:), allocatable :: field_names_mntr
     !> Vector list pointing to instantaneous error vectors
     type(vector_list_t) :: eind_in_mntr
     !> Vector list pointing to averaged error vectors
     type(vector_list_t) :: eind_av_mntr
   contains
     !> Initialise type
     procedure, pass(this) :: init => amr_spectral_error_init
     !> Free type
     procedure, pass(this) :: free => amr_spectral_error_free
     !> Write instantaneous and averaged fields to the disc
     procedure, pass(this) :: sample => amr_spectral_error_sample
     !> Checkpoint averaged fields to the disc
     procedure, pass(this) :: checkpoint => amr_spectral_error_checkpoint
     !> AMR restart
     procedure, pass(this) :: amr_restart => amr_spectral_error_amr_restart
  end type amr_spectral_error_t

contains

  !> Initialise spectral error indicator type
  !! @param[in]   msh            Mesh
  !! @param[in]   fld_name_ref   List of fields for refinement mark assessment
  !! @param[in]   opr            Operator applied
  !! @param[in]   fld_name_mntr  List of fields for stability monitoring
  subroutine amr_spectral_error_init(this, msh, fld_name_ref, opr, &
       fld_name_mntr)
    class(amr_spectral_error_t), intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    character(len=NEKO_VARNAME_LEN), dimension(:), optional, intent(in) :: &
         fld_name_ref, fld_name_mntr
    character(len=*), optional, intent(in) :: opr
    integer :: il, jl, n_simcomps, lx
    type(field_t) :: fld_tmp
    logical :: found

    Call this%free()

    ! Identify spectral error indicator simcomp
    n_simcomps = neko_simcomps%get_n()
    do il = 1, n_simcomps
       call simcomp_get_pointer(neko_simcomps%simcomps(il)%simcomp, &
            this%simcomp_spectral)
    end do

    if (.not. associated(this%simcomp_spectral)) call neko_error('AMR tools; &
         &spectral error indicator simcomp not found')

    ! Number of fields in simcomp
    this%nfld_smp = this%simcomp_spectral%nfld

    ! minimal mesh
    this%msh => msh
    lx = 3
    call this%Xh%init(GLL, lx, lx, lx)
    call this%dm_Xh%init(this%msh, this%Xh)
    call this%gs_Xh%init(this%dm_Xh)

    ! Get fields for visualisation of
    ! instantaneous indicators
    call this%eind_in_fld%init(this%nfld_smp)
    do il = 1, this%nfld_smp
       call fld_tmp%init(this%dm_Xh, 'eind_in_'//&
            trim(this%simcomp_spectral%field_names(il)))
       call this%eind_in_fld%items(il)%init(fld_tmp)
       call fld_tmp%free()
    end do
    ! averaged indicators
    ! nfld_smp + 1 for refinement flag
    call this%eind_av_fld%init(this%nfld_smp + 1)
    do il = 1, this%nfld_smp
       call fld_tmp%init(this%dm_Xh, 'eind_av_'//&
            trim(this%simcomp_spectral%field_names(il)))
       call this%eind_av_fld%items(il)%init(fld_tmp)
       call fld_tmp%free()
    end do
    call fld_tmp%init(this%dm_Xh, 'rflag')
    call this%eind_av_fld%items(this%nfld_smp + 1)%init(fld_tmp)
    call fld_tmp%free()

    ! Initialise file outputs for error indicator
    ! instantaneous
    call this%eind_in_file%init(rp, "error_ind_in", this%nfld_smp)
    do il = 1, this%nfld_smp
       call this%eind_in_file%fields%assign_to_ptr(il, &
            this%eind_in_fld%items(il)%field)
    end do
    select type(file => this%eind_in_file%file_%file_type)
    type is (fld_file_t)
       file%write_mesh = .true.
    end select
    ! averaged
    call this%eind_av_file%init(rp, "error_ind_av", this%nfld_smp + 1)
    do il = 1, this%nfld_smp + 1
       call this%eind_av_file%fields%assign_to_ptr(il, &
            this%eind_av_fld%items(il)%field)
    end do
    select type(file => this%eind_av_file%file_%file_type)
    type is (fld_file_t)
       file%write_mesh = .true.
    end select

    ! Averaged fields used for refinement mark assessment
    if (present(fld_name_ref)) then
       this%nfld_ref = size(fld_name_ref)
       if (this%nfld_ref .gt. 0) then
          ! save field names
          allocate(this%field_names_ref(this%nfld_ref))
          do il = 1, this%nfld_ref
             this%field_names_ref(il) = fld_name_ref(il)
          end do

          ! check if listed fields are present in simcomp list and map them
          ! to the averaged field
          call this%eind_av_ref%init(this%nfld_ref)
          do il = 1, this%nfld_ref
             found = .false.
             do jl = 1, this%nfld_smp
                if (trim(this%field_names_ref(il)) .eq. &
                     trim(this%simcomp_spectral%field_names(jl))) then
                   this%eind_av_ref%items(il)%ptr => &
                        this%simcomp_spectral%eind_av%items(jl)%ptr
                   found = .true.
                   exit
                end if
             end do
             if (.not. found) call neko_error('AMR tools; refinement field not &
                  &mapped')
          end do
       end if
       if (present(opr)) then
          this%opr = trim(opr)
       else
          this%opr = 'max'
       end if
    end if

    ! Fields used for stability monitoring
    if (present(fld_name_mntr)) then
       this%nfld_mntr = size(fld_name_mntr)
       if (this%nfld_mntr .gt. 0) then
          ! save field names
          allocate(this%field_names_mntr(this%nfld_mntr))
          do il = 1, this%nfld_mntr
             this%field_names_mntr(il) = fld_name_mntr(il)
          end do

          ! check if listed fields are present in simcomp list and map them
          ! to the averaged field
          call this%eind_in_mntr%init(this%nfld_mntr)
          call this%eind_av_mntr%init(this%nfld_mntr)
          do il = 1, this%nfld_mntr
             found = .false.
             do jl = 1, this%nfld_smp
                if (trim(this%field_names_mntr(il)) .eq. &
                     trim(this%simcomp_spectral%field_names(jl))) then
                   this%eind_in_mntr%items(il)%ptr => &
                        this%simcomp_spectral%eind%items(jl)%ptr
                   this%eind_av_mntr%items(il)%ptr => &
                        this%simcomp_spectral%eind_av%items(jl)%ptr
                   found = .true.
                   exit
                end if
             end do
             if (.not. found) call neko_error('AMR tools; monitor field not &
                  &mapped')
          end do
       end if
    end if

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

    this%nfld_smp = 0
    this%nelv = 0
    this%nfld_ref = 0
    this%nfld_mntr = 0

    if (allocated(this%field_names_ref)) deallocate(this%field_names_ref)
    if (allocated(this%opr)) deallocate(this%opr)
    if (allocated(this%field_names_mntr)) deallocate(this%field_names_mntr)

    call this%gs_Xh%free()
    call this%dm_Xh%free()
    call this%Xh%free()

    call this%eind_in_fld%free()
    call this%eind_av_fld%free()

    call this%eind_av_ref%free()
    call this%eind_in_mntr%free()
    call this%eind_av_mntr%free()

  end subroutine amr_spectral_error_free

  !> Write instantaneous and averaged fields to the disc
  !! @param[in]     ref_mark     Refinement flag
  !! @param[in]     time         Simulation time
  subroutine amr_spectral_error_sample(this, ref_mark, time)
    class(amr_spectral_error_t), intent(inout) :: this
    integer, dimension(:), intent(in) :: ref_mark
    real(rp), intent(in) :: time
    integer :: il, jl

    ! instantaneous indicator
    do il = 1, this%nfld_smp
       do jl = 1, this%simcomp_spectral%nelv
          this%eind_in_fld%items(il)%field%x(:, :, :, jl) = &
               this%simcomp_spectral%eind%items(il)%ptr%x(jl)
       end do
    end do
    call this%eind_in_file%sample(time)

    ! averaged indicator
    do il = 1, this%nfld_smp
       do jl = 1, this%simcomp_spectral%nelv
          this%eind_av_fld%items(il)%field%x(:, :, :, jl) = &
               this%simcomp_spectral%eind_av%items(il)%ptr%x(jl)
       end do
    end do
    do jl = 1, this%simcomp_spectral%nelv
       this%eind_av_fld%items(this%nfld_smp + 1)%field%x(:, :, :, jl) = &
            ref_mark(jl)
    end do
    call this%eind_av_file%sample(time)

  end subroutine amr_spectral_error_sample

  !> Checkpoint averaged fields to the disc
  subroutine amr_spectral_error_checkpoint(this)
    class(amr_spectral_error_t), intent(inout) :: this
    integer :: il, jl
    real(rp) :: time

    ! averaged indicator
    do il = 1, this%nfld_smp
       do jl = 1, this%simcomp_spectral%nelv
          this%eind_av_fld%items(il)%field%x(:, :, :, jl) = &
               this%simcomp_spectral%eind_av%items(il)%ptr%x(jl)
       end do
    end do
    this%eind_av_fld%items(this%nfld_smp + 1)%field%x(:, :, :, :) = 0.0_rp
    time = this%simcomp_spectral%time_start
    call this%eind_av_file%sample(time)

  end subroutine amr_spectral_error_checkpoint

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     tstep         time step
  subroutine amr_spectral_error_amr_restart(this, reconstruct, counter, tstep)
    class(amr_spectral_error_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter, tstep

    this%nelv = this%simcomp_spectral%nelv

    ! Reconstruct minimal mesh and communicator
    call this%dm_Xh%amr_restart(reconstruct, counter, tstep)
    call this%gs_Xh%amr_restart(reconstruct, counter, tstep)

    ! instantaneous and averaged fields for visualisation
    call this%eind_in_fld%amr_reallocate(reconstruct, counter, tstep)
    call this%eind_av_fld%amr_reallocate(reconstruct, counter, tstep)

  end subroutine amr_spectral_error_amr_restart

end module amr_tools
