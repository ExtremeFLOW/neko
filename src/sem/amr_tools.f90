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
  use time_state, only : time_state_t
  use amr_reconstruct, only : amr_reconstruct_t
  use amr_restart_component, only : amr_restart_component_t

  implicit none
  private

  integer, parameter, public :: AMR_OP_MAX = 1, AMR_OP_ELEN = 2

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
     !> Get averaging time
     procedure, pass(this) :: get_collect_time => &
          amr_spectral_error_get_collect_time
     !> Reset averaged variables
     procedure, pass(this) :: reset_average => amr_spectral_error_reset_average
     !> Get combined averaged error indicator
      procedure, pass(this) :: err_av_get => amr_spectral_error_err_av_get
     !> AMR restart
     procedure, pass(this) :: amr_restart => amr_spectral_error_amr_restart
  end type amr_spectral_error_t

contains

  !> Initialise spectral error indicator type
  !! @param[in]   msh            Mesh
  !! @param[in]   fld_name_ref   List of fields for refinement mark assessment
  !! @param[in]   fld_name_mntr  List of fields for stability monitoring
  subroutine amr_spectral_error_init(this, msh, fld_name_ref, fld_name_mntr)
    class(amr_spectral_error_t), intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    character(len=NEKO_VARNAME_LEN), dimension(:), optional, intent(in) :: &
         fld_name_ref, fld_name_mntr
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
    ! Array length
    this%nelv = this%simcomp_spectral%nelv

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

  !> Get averaging time
  pure function amr_spectral_error_get_collect_time(this) result(time)
    class(amr_spectral_error_t), intent(in) :: this
    real(dp) :: time

    time = this%simcomp_spectral%get_collect_time()

  end function amr_spectral_error_get_collect_time

  !> Reset averaged variables
  subroutine amr_spectral_error_reset_average(this)
    class(amr_spectral_error_t), intent(inout) :: this
    real(dp) :: time

    this%simcomp_spectral%time_start = this%simcomp_spectral%time_previous
    call this%simcomp_spectral%average_reset()

  end subroutine amr_spectral_error_reset_average

  !> Get combined averaged error indicator
  !! @param[in]     opr      operator
  !! @param[out]    errind   combined averaged error indicator
  !! @param[in]     nelv     number of elements
  subroutine amr_spectral_error_err_av_get(this, opr, errind, nelv)
    class(amr_spectral_error_t), intent(in) :: this
    integer, intent(in) :: opr, nelv
    real(rp), dimension(nelv), intent(out) :: errind
    integer :: il, jl
    real(rp) :: rtmp
    real(rp), parameter :: nbig = -huge(0.0_rp)

    if (nelv .ne. this%nelv) call neko_error('AMR tools; &
         &inconsistent number of elements')

    select case(opr)
    case (AMR_OP_MAX)
       ! take max from all the listed fields
       do il = 1, this%nelv
          rtmp = nbig
          do jl = 1, this%nfld_ref
             rtmp = max(rtmp, this%eind_av_ref%items(jl)%ptr%x(il))
          end do
          errind(il) = rtmp
       end do
    case (AMR_OP_ELEN)
       ! Euclidean length of vector
       do il = 1, this%nelv
          rtmp = 0.0_rp
          do jl = 1, this%nfld_ref
             rtmp = rtmp + this%eind_av_ref%items(jl)%ptr%x(il)**2
          end do
          errind(il) = sqrt(rtmp)
       end do
    case default
       call neko_error('AMR tools; wrong operator')
    end select

  end subroutine amr_spectral_error_err_av_get

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     time          time state
  subroutine amr_spectral_error_amr_restart(this, reconstruct, counter, time)
    class(amr_spectral_error_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter
    type(time_state_t), intent(in) :: time

    this%nelv = this%simcomp_spectral%nelv

    ! Reconstruct minimal mesh and communicator
    call this%dm_Xh%amr_restart(reconstruct, counter, time)
    call this%gs_Xh%amr_restart(reconstruct, counter, time)

    ! instantaneous and averaged fields for visualisation
    call this%eind_in_fld%amr_reallocate(reconstruct, counter, time)
    call this%eind_av_fld%amr_reallocate(reconstruct, counter, time)

  end subroutine amr_spectral_error_amr_restart

end module amr_tools
