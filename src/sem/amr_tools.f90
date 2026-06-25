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
  use comm, only : NEKO_COMM, pe_rank
  use utils, only : NEKO_VARNAME_LEN, neko_error
  use logger, only : neko_log, NEKO_LOG_QUIET, NEKO_LOG_INFO, &
       NEKO_LOG_VERBOSE, NEKO_LOG_DEBUG, LOG_SIZE
  use vector_list, only : vector_list_t
  use simulation_component, only : simulation_component_t
  use spectral_error, only : spectral_error_t
  use simcomp_executor, only : neko_simcomps
  use mesh, only : mesh_t
  use sem, only : sem_lx_t
  use field, only : field_t
  use field_array, only : field_array_t
  use fld_file, only : fld_file_t
  use fld_file_output, only : fld_file_output_t
  use time_state, only : time_state_t
  use gs_ops, only : GS_OP_ADD, GS_OP_MIN, GS_OP_MAX
  use amr_reconstruct, only : amr_reconstruct_t, amr_flg_none, amr_flg_h_ref, &
       amr_flg_h_crs
  use amr_restart_component, only : amr_restart_component_t
  use mpi_f08, only : MPI_INTEGER, MPI_SUM, MPI_IN_PLACE, MPI_Allreduce

  implicit none
  private

  integer, parameter, public :: AMR_OP_MAX = 1, AMR_OP_ELEN = 2

  integer, parameter :: lx = 3, dim = 3, nvrt = 8, nedg = 12, nfcs = 6

  ! vertex position in the minimal grid box
  integer, parameter, dimension(dim, nvrt) :: vert =reshape([&
       1, 1, 1,  3, 1, 1,  1, 3, 1,  3, 3, 1,  &
       1, 1, 3,  3, 1, 3,  1, 3, 3,  3, 3, 3], shape(vert))

  ! edge position in the minimal grid box
  integer, parameter, dimension(dim, nedg) :: edge =reshape([&
       2, 1, 1,  2, 3, 1,  2, 1, 3,  2, 3, 3,  &
       1, 2, 1,  3, 2, 1,  1, 2, 3,  3, 2, 3,  &
       1, 1, 2,  3, 1, 2,  1, 3, 2,  3, 3, 2], shape(edge))

  ! face position in the minimal grid box
  integer, parameter, dimension(dim, nfcs) :: face =reshape([&
       1, 2, 2,  3, 2, 2,  2, 1, 2,  2, 3, 2,  2, 2, 1,  2, 2, 3], shape(face))

  ! faces sharing vertex in the minimal grid box
  integer, parameter, dimension(dim, 3, nvrt) :: fsv =reshape([ &
       1, 2, 2,  2, 1, 2,  2, 2, 1, &
       3, 2, 2,  2, 1, 2,  2, 2, 1, &
       1, 2, 2,  2, 3, 2,  2, 2, 1, &
       3, 2, 2,  2, 3, 2,  2, 2, 1, &
       1, 2, 2,  2, 1, 2,  2, 2, 3, &
       3, 2, 2,  2, 1, 2,  2, 2, 3, &
       1, 2, 2,  2, 3, 2,  2, 2, 3, &
       3, 2, 2,  2, 3, 2,  2, 2, 3], shape(fsv))

  ! faces sharing edge in the minimal grid box
  integer, parameter, dimension(dim, 2, nedg) :: fse =reshape([ &
       2, 1, 2,  2, 2, 1,   2, 3, 2,  2, 2, 1, &
       2, 1, 2,  2, 2, 3,   2, 3, 2,  2, 2, 3, &
       1, 2, 2,  2, 2, 1,   3, 2, 2,  2, 2, 1, &
       1, 2, 2,  2, 2, 3,   3, 2, 2,  2, 2, 3, &
       1, 2, 2,  2, 1, 2,   3, 2, 2,  2, 1, 2, &
       1, 2, 2,  2, 3, 2,   3, 2, 2,  2, 3, 2], shape(fse))

  ! pairs of opposing faces in the minimal grid box
  integer, parameter, dimension(dim, 2, dim) :: fop =reshape([ &
       1, 2, 2,  3, 2, 2,   2, 1, 2,  2, 3, 2,   2, 2, 1,  2, 2, 3], shape(fop))

  ! edge vertices
  integer, public, parameter, dimension(2, 12) :: vedge  = reshape([ &
       1, 2,  3, 4,  5, 6,  7, 8,  1, 3,  2, 4, &
       5, 7,  6, 8,  1, 5,  2, 6,  3, 7,  4, 8], shape(vedge))

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
     !> Minimal grid
     type(sem_lx_t) :: grid_min
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

  ! Specific tools for managing refinement on the user side
  public :: amr_ref_mark_check, amr_nonconf_int_remove

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
    integer :: il, jl, n_simcomps
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

    ! Add minimal grid
    this%msh => msh
    call this%grid_min%init(this%msh, lx, .false., .false.)

    ! Get fields for visualisation of
    ! instantaneous indicators
    call this%eind_in_fld%init(this%nfld_smp)
    do il = 1, this%nfld_smp
       call fld_tmp%init(this%grid_min%dm_Xh, 'eind_in_'//&
            trim(this%simcomp_spectral%field_names(il)))
       call this%eind_in_fld%items(il)%init(fld_tmp)
       call fld_tmp%free()
    end do
    ! averaged indicators
    ! nfld_smp + 1 for refinement flag
    call this%eind_av_fld%init(this%nfld_smp + 1)
    do il = 1, this%nfld_smp
       call fld_tmp%init(this%grid_min%dm_Xh, 'eind_av_'//&
            trim(this%simcomp_spectral%field_names(il)))
       call this%eind_av_fld%items(il)%init(fld_tmp)
       call fld_tmp%free()
    end do
    call fld_tmp%init(this%grid_min%dm_Xh, 'rflag')
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

    call this%grid_min%free()

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

    ! Reconstruct minimal grid
    call this%grid_min%amr_restart(reconstruct, counter, time)

    ! instantaneous and averaged fields for visualisation
    call this%eind_in_fld%amr_reallocate(reconstruct, counter, time)
    call this%eind_av_fld%amr_reallocate(reconstruct, counter, time)

  end subroutine amr_spectral_error_amr_restart

  ! Set vertex in the box
  subroutine amr_vertex_set(box, ivrt, const)
    integer, dimension(lx, lx, lx), intent(inout) :: box
    integer, intent(in) :: ivrt, const

    box(vert(1, ivrt), vert(2, ivrt), vert(3, ivrt)) = const

  end subroutine amr_vertex_set

  ! Extract vertex from the box
  subroutine amr_vertex_get(box, ivrt, const)
    integer, dimension(lx, lx, lx), intent(in) :: box
    integer, intent(in) :: ivrt
    integer, intent(out) :: const

    const = box(vert(1, ivrt), vert(2, ivrt), vert(3, ivrt))

  end subroutine amr_vertex_get

  ! Set edge in the box
  subroutine amr_edge_set(box, iedg, const)
    integer, dimension(lx, lx, lx), intent(inout) :: box
    integer, intent(in) :: iedg, const

    box(edge(1, iedg), edge(2, iedg), edge(3, iedg)) = const

  end subroutine amr_edge_set

  ! Extract edge from the box
  subroutine amr_edge_get(box, iedg, const)
    integer, dimension(lx, lx, lx), intent(in) :: box
    integer, intent(in) :: iedg
    integer, intent(out) :: const

    const = box(edge(1, iedg), edge(2, iedg), edge(3, iedg))

  end subroutine amr_edge_get

  ! Set face in the box
  subroutine amr_face_set(box, ifcs, const)
    integer, dimension(lx, lx, lx), intent(inout) :: box
    integer, intent(in) :: ifcs, const

    box(face(1, ifcs), face(2, ifcs), face(3, ifcs)) = const

  end subroutine amr_face_set

  ! Extract face from the box
  subroutine amr_face_get(box, ifcs, const)
    integer, dimension(lx, lx, lx), intent(in) :: box
    integer, intent(in) :: ifcs
    integer, intent(out) :: const

    const = box(face(1, ifcs), face(2, ifcs), face(3, ifcs))

  end subroutine amr_face_get

  ! Set vertex sharing faces in the box
  subroutine amr_vertex_face_set(box, ivrt, const)
    integer, dimension(lx, lx, lx), intent(inout) :: box
    integer, intent(in) :: ivrt, const

    box(fsv(1, 1, ivrt), fsv(2, 1, ivrt), fsv(3, 1, ivrt)) = const
    box(fsv(1, 2, ivrt), fsv(2, 2, ivrt), fsv(3, 2, ivrt)) = const
    box(fsv(1, 3, ivrt), fsv(2, 3, ivrt), fsv(3, 3, ivrt)) = const

  end subroutine amr_vertex_face_set

  ! Extract vertex sharing faces in the box
  subroutine amr_vertex_face_get(box, ivrt, const)
    integer, dimension(lx, lx, lx), intent(in) :: box
    integer, intent(in) :: ivrt
    integer, dimension(3), intent(out) :: const

    const(1) = box(fsv(1, 1, ivrt), fsv(2, 1, ivrt), fsv(3, 1, ivrt))
    const(2) = box(fsv(1, 2, ivrt), fsv(2, 2, ivrt), fsv(3, 2, ivrt))
    const(3) = box(fsv(1, 3, ivrt), fsv(2, 3, ivrt), fsv(3, 3, ivrt))

  end subroutine amr_vertex_face_get

  ! Set edge sharing faces in the box
  subroutine amr_edge_face_set(box, iedg, const)
    integer, dimension(lx, lx, lx), intent(inout) :: box
    integer, intent(in) :: iedg, const

    box(fse(1, 1, iedg), fse(2, 1, iedg), fse(3, 1, iedg)) = const
    box(fse(1, 2, iedg), fse(2, 2, iedg), fse(3, 2, iedg)) = const

  end subroutine amr_edge_face_set

  ! Extract edge sharing faces in the box
  subroutine amr_edge_face_get(box, iedg, const)
    integer, dimension(lx, lx, lx), intent(in) :: box
    integer, intent(in) :: iedg
    integer, dimension(2), intent(out) :: const

    const(1) = box(fse(1, 1, iedg), fse(2, 1, iedg), fse(3, 1, iedg))
    const(2) = box(fse(1, 2, iedg), fse(2, 2, iedg), fse(3, 2, iedg))

  end subroutine amr_edge_face_get

  ! Extract opposite face from the minimal grid box
  subroutine amr_oposite_face_extract(box, face_pair)
    integer, dimension(lx, lx, lx), intent(in) :: box
    integer, dimension(2, dim), intent(out) :: face_pair
    integer :: il, jl

    do il = 1, dim
       do jl = 1, 2
          face_pair(jl, il) = box(fop(1, jl, il), fop(2, jl, il), &
               fop(3, jl, il))
       end do
    end do
  end subroutine amr_oposite_face_extract

  !> Make coarsening flag consistent within families.
  !! @param[in]     nelv        local element number
  !! @param[in]     ref_family  element family information for coarsening
  !! @param[inout]  ref_mark    refinement flag
  !! @param[inout]  grid_min    minimal grid
  subroutine amr_family_crs_check(nelv, ref_family, ref_mark, grid_min)
    integer, intent(in) :: nelv
    integer, dimension(2, nelv), intent(in) :: ref_family
    integer, dimension(nelv), intent(inout) :: ref_mark
    type(sem_lx_t), intent(inout) :: grid_min
    real(rp), dimension(lx, lx, lx, nelv) :: exchange
    integer, dimension(lx, lx, lx, nelv) :: mark_max
    integer :: ntot, itmp, il

    ntot = lx * lx* lx * nelv

    ! Make coarsening flag consistent within families. Elements in the family
    ! can be marked for coarsening only if all the members are flagged this
    ! way.
    do il = 1, nelv
       exchange(:, :, :, il) = real(ref_mark(il), rp)
    end do
    call grid_min%gs_Xh%gs_op_vector(exchange, ntot, GS_OP_MAX)
    mark_max = nint(exchange)
    do il = 1, nelv
       if (ref_family(1, il) .gt. 0 .and. ref_mark(il) .eq. amr_flg_h_crs) then
          call amr_vertex_get(mark_max(:, :, :, il), ref_family(2, il), itmp)
          if (itmp .ge. amr_flg_none) ref_mark(il) = amr_flg_none
       end if
    end do
  end subroutine amr_family_crs_check


  !> Balance refinement with 2:1 criterion. It is a crude one.
  !! @param[in]     nelv        local element number
  !! @param[in]     ref_level   element refinement level (current or predicted)
  !! @param[in]     ref_family  element family information for coarsening
  !! @param[inout]  ref_mark    refinement flag
  !! @param[inout]  grid_min    minimal grid
  !! @param[inout]  iter        max/performed number of iterations
  !! @param[out]    nmod        global number of (modified, not fixed) elements
  subroutine amr_ref_balance(nelv, ref_level, ref_family, ref_mark, grid_min, &
       iter, nmod)
    integer, intent(in) :: nelv
    integer, dimension(nelv), intent(in) :: ref_level
    integer, dimension(2, nelv), intent(in) :: ref_family
    integer, dimension(nelv), intent(inout) :: ref_mark
    type(sem_lx_t), intent(inout) :: grid_min
    integer, intent(inout) ::  iter
    integer, dimension(2), intent(out) :: nmod
    integer, dimension(2) :: nmod_iter
    integer, dimension(nelv) :: level_tmp
    real(rp), dimension(lx, lx, lx, nelv) :: exchange
    integer, dimension(lx, lx, lx, nelv) :: level_max
    integer :: ntot, iter_max, itmp, ierr, il
    character(len=LOG_SIZE) :: log_buf

    ntot = lx * lx* lx * nelv
    nmod(:) = 0
    iter_max = iter
    iter = 0

    do
       iter = iter + 1
       nmod_iter(:) = 0

       ! Make coarsening flag consistent within families.
       call amr_family_crs_check(nelv, ref_family, ref_mark, grid_min)

       ! Get predicted refinement level
       do il = 1, nelv
          select case(ref_mark(il))
          case(amr_flg_h_crs)
             level_tmp(il) = ref_level(il) - 1
          case(amr_flg_h_ref)
             level_tmp(il) = ref_level(il) + 1
          case(amr_flg_none)
             level_tmp(il) = ref_level(il)
          end select
       end do

       ! Distribute element refinement level
       do il = 1, nelv
          exchange(:, :, :, il) = real(level_tmp(il), rp)
       end do
       ! possible problem for nonconforming interfaces; vertices/edges
       call grid_min%gs_Xh%gs_op_vector(exchange, ntot, GS_OP_MAX)
       level_max = nint(exchange)

       ! Check refinement level difference across interfaces
       do il = 1, nelv
          itmp = maxval(level_max(:, :, :, il))
          ! We modify amr_flg_none and amr_flg_h_crs only
          if (itmp - level_tmp(il) .gt. 1) then
             if (ref_mark(il) .ne. amr_flg_h_ref) then
                nmod_iter(1) = nmod_iter(1) + 1
                select case(ref_mark(il))
                case(amr_flg_h_crs)
                   ref_mark(il) = amr_flg_none
                case(amr_flg_none)
                   ref_mark(il) = amr_flg_h_ref
                end select
             else
                nmod_iter(2) = nmod_iter(2) + 1
             end if
          end if
       end do

       ! global number of modified elements
       call MPI_Allreduce(MPI_IN_PLACE, nmod_iter, 2, MPI_INTEGER, MPI_SUM, &
            NEKO_COMM, ierr)

       nmod(1) = nmod(1) + nmod_iter(1)
       nmod(2) = nmod_iter(2)
       if (nmod_iter(1) .eq. 0) then
          exit
       else
          if (iter .eq. iter_max) then
             write(log_buf, '(A,I9)') 'AMR tool balance; balancing not &
                  &finalised: ', nmod_iter(1)
             call neko_log%message(log_buf, NEKO_LOG_INFO)
             exit
          end if
       end if
    end do

  end subroutine amr_ref_balance

  !> Fix edge/vertex connected regions for single iteration
  !! @param[in]     nelv        local element number
  !! @param[in]     ref_level   element refinement level (current or predicted)
  !! @param[in]     ifcurrent   is a refinement level a current one
  !! @param[in]     flag_exc    flag marking excluded elements
  !! @param[inout]  ref_mark    refinement flag
  !! @param[inout]  grid_min    minimal grid
  !! @param[out]    nmod        global number of modified elements
  subroutine amr_ref_connect_fix(nelv, ref_level, ifcurrent, flag_exc, &
       ref_mark, grid_min, nmod)
    integer, intent(in) :: nelv
    integer, dimension(nelv), intent(in) :: ref_level
    logical, intent(in) :: ifcurrent
    logical, dimension(nelv), intent(inout) :: flag_exc
    integer, dimension(nelv), intent(inout) :: ref_mark
    type(sem_lx_t), intent(inout) :: grid_min
    integer, intent(out) :: nmod
    real(rp), dimension(lx, lx, lx, nelv) :: exchange
    integer, dimension(lx, lx, lx, nelv) :: level_min, level_max, mult, &
         local_mark
    ! vertex sharing faces
    integer, dimension(3) :: fvert
    ! edge sharing faces
    integer, dimension(2) :: fedge
    integer :: ntot, mult_min_vrt, mult_min_edg, lmax, lmin, itmp, ierr, il, jl
    logical, dimension(nvrt) :: ifvert

    ! NOTICE. Following algorithm is not fully general, as it e.g., may not
    ! work with conforming edges with multiplicity higher than 6.

    ntot = lx * lx* lx * nelv
    nmod  = 0
    if (ifcurrent) then
       ! take into account nonconforming faces
       mult_min_edg = 6
       mult_min_vrt = 7
    else
       ! conforming faces only
       mult_min_edg = 2
       mult_min_vrt = 2
    end if

    ! Distribute element refinement level
    do il = 1, nelv
       exchange(:, :, :, il) = real(ref_level(il), rp)
    end do
    call grid_min%gs_Xh%gs_op_vector(exchange, ntot, GS_OP_MIN)
    level_min = nint(exchange)
    do il = 1, nelv
       exchange(:, :, :, il) = real(ref_level(il), rp)
    end do
    call grid_min%gs_Xh%gs_op_vector(exchange, ntot, GS_OP_MAX)
    level_max = nint(exchange)

    ! Get multiplicity of highest refinement level through edges and vertices
    mult(:, :, :, :) = 0
    do il = 1, nelv
       ! Skip excluded elements
       if (.not. flag_exc(il)) then
          ! Refined regions connected by edge
          do jl = 1, nedg
             call amr_edge_get(level_max(:, :, :, il), jl, lmax)
             call amr_edge_get(level_min(:, :, :, il), jl, lmin)
             if (lmax .ne. lmin .and. ref_level(il) .gt. lmin) then
                call amr_edge_set(mult(:, :, :, il), jl, 1)
             end if
          end do

          ! Refined regions connected by vertex
          do jl = 1, nvrt
             call amr_vertex_get(level_max(:, :, :, il), jl, lmax)
             call amr_vertex_get(level_min(:, :, :, il), jl, lmin)
             if (lmax .ne. lmin .and. ref_level(il) .gt. lmin) then
                call amr_vertex_set(mult(:, :, :, il), jl, 1)
             end if
          end do
       end if
    end do

    exchange = real(mult, rp)
    call grid_min%gs_Xh%gs_op_vector(exchange, ntot, GS_OP_ADD)
    mult = nint(exchange)

    ! Mark problematic elements
    local_mark(:, :, :, :) = amr_flg_none
    do il = 1, nelv
       ! Skip excluded elements
       if (.not. flag_exc(il)) then
          ifvert(:) = .true.
          ! Refined regions connected by edge
          ! Start with edges to exclude vertices
          do jl = 1, nedg
             call amr_edge_get(level_max(:, :, :, il), jl, lmax)
             call amr_edge_get(level_min(:, :, :, il), jl, lmin)
             call amr_edge_get(mult(:, :, :, il), jl, itmp)
             call amr_edge_face_get(level_min(:, :, :, il), jl, fedge)
             if (lmax .ne. lmin .and. ref_level(il) .gt. lmin .and. &
                  itmp .gt. mult_min_edg .and. &
                  maxval(fedge) .lt. ref_level(il)) then
                flag_exc(il) = .true.
                ifvert(vedge(1, jl)) = .false.
                ifvert(vedge(2, jl)) = .false.
                call amr_edge_face_set(local_mark(:, :, :, il), jl, &
                     amr_flg_h_ref)
             end if
          end do

          ! Refined regions connected by vertex
          do jl = 1, nvrt
             ! Skip marked vertices
             if (ifvert(jl)) then
                call amr_vertex_get(level_max(:, :, :, il), jl, lmax)
                call amr_vertex_get(level_min(:, :, :, il), jl, lmin)
                call amr_vertex_get(mult(:, :, :, il), jl, itmp)
                call amr_vertex_face_get(level_min(:, :, :, il), jl, fvert)
                if (lmax .ne. lmin .and. ref_level(il) .gt. lmin .and. &
                     itmp .gt. mult_min_vrt .and.&
                     maxval(fvert) .lt. ref_level(il)) then
                   flag_exc(il) = .true.
                   call amr_vertex_face_set(local_mark(:, :, :, il), jl, &
                        amr_flg_h_ref)
                end if
             end if
          end do
       end if
    end do

    exchange = real(local_mark, rp)
    call grid_min%gs_Xh%gs_op_vector(exchange, ntot, GS_OP_MAX)
    local_mark = nint(exchange)

    if (ifcurrent) then
       do il = 1, nelv
          ! Skip excluded elements
          if (.not. flag_exc(il)) then
             faces1 : do jl = 1, nfcs
                call amr_face_get(level_max(:, :, :, il), jl, lmax)
                call amr_face_get(local_mark(:, :, :, il), jl, itmp)
                if (lmax .gt. ref_level(il) .and. itmp .eq. amr_flg_h_ref &
                     .and. ref_mark(il) .ne. amr_flg_h_ref) then
                   nmod = nmod + 1
                   ref_mark(il) = amr_flg_h_ref
                   exit faces1
                end if
             end do faces1
          end if
       end do
    else
       do il = 1, nelv
          ! Skip excluded elements
          if (.not. flag_exc(il)) then
             faces2 : do jl = 1, nfcs
                call amr_face_get(level_max(:, :, :, il), jl, lmax)
                call amr_face_get(local_mark(:, :, :, il), jl, itmp)
                if (lmax .gt. ref_level(il) .and. itmp .eq. amr_flg_h_ref &
                     .and. ref_mark(il) .ne. amr_flg_h_ref) then
                   nmod = nmod + 1
                   select case(ref_mark(il))
                   case(amr_flg_h_crs)
                      ref_mark(il) = amr_flg_none
                   case(amr_flg_none)
                      ref_mark(il) = amr_flg_h_ref
                   end select
                   exit faces2
                end if
             end do faces2
          end if
       end do
    end if

    ! global number of modified elements
    call MPI_Allreduce(MPI_IN_PLACE, nmod, 1, MPI_INTEGER, MPI_SUM, &
         NEKO_COMM, ierr)

  end subroutine amr_ref_connect_fix

  !> Remove single element gaps
  !! @param[in]     nelv        local element number
  !! @param[in]     ref_level   element refinement level (current or predicted)
  !! @param[inout]  ref_mark    refinement flag
  !! @param[inout]  grid_min    minimal grid
  !! @param[out]    nmod        global number of (modified, not fixed) elements
  subroutine amr_ref_gap_remove(nelv, ref_level, ref_mark, grid_min, nmod)
    integer, intent(in) :: nelv
    integer, dimension(nelv), intent(in) :: ref_level
    integer, dimension(nelv), intent(inout) :: ref_mark
    type(sem_lx_t), intent(inout) :: grid_min
    integer, dimension(2), intent(out) :: nmod
    real(rp), dimension(lx, lx, lx, nelv) :: exchange
    integer, dimension(lx, lx, lx, nelv) :: level_max
    integer, dimension(2, 3) :: fpair_max
    integer :: ntot, il, jl, ierr
    logical :: ifgap

    ntot = lx * lx* lx * nelv
    nmod(:) = 0

    ! Distribute element refinement level
    do il = 1, nelv
       exchange(:, :, :, il) = real(ref_level(il), rp)
    end do
    ! possible problem for nonconforming interfaces; faces
    call grid_min%gs_Xh%gs_op_vector(exchange, ntot, GS_OP_MAX)
    level_max = nint(exchange)

    do il = 1, nelv
       call amr_oposite_face_extract(level_max(:, :, :, il), fpair_max)
       ifgap = .false.
       do jl = 1, 3
          if (fpair_max(1, jl) .gt. ref_level(il) .and. &
               fpair_max(2, jl) .gt. ref_level(il)) then
             ifgap = .true.
             exit
          end if
       end do
       if (ifgap) then
          if (ref_mark(il) .ne. amr_flg_h_ref) then
             nmod(1) = nmod(1) + 1
             select case(ref_mark(il))
             case(amr_flg_h_crs)
                ref_mark(il) = amr_flg_none
             case(amr_flg_none)
                ref_mark(il) = amr_flg_h_ref
             end select
          else
             nmod(2) = nmod(2) + 1
          end if
       end if
    end do

    ! global number of modified elements
    call MPI_Allreduce(MPI_IN_PLACE, nmod, 2, MPI_INTEGER, MPI_SUM, &
         NEKO_COMM, ierr)

  end subroutine amr_ref_gap_remove

  !> Check refinement flag for possible topology problems
  !! @param[in]     nelv        local element number
  !! @param[in]     ref_level   element refinement level
  !! @param[in]     ref_family  element family information for coarsening
  !! @param[inout]  ref_mark    refinement flag
  !! @param[inout]  grid_min    minimal grid
  !! @param[inout]  iter        max/performed number of iterations
  !! @param[in]     iterb       max number of iterations for balancing
  !! @param[out]    nmod        global number of modified elements
  subroutine amr_ref_mark_check(nelv, ref_level, ref_family, ref_mark, &
       grid_min, iter, iterb, nmod)
    integer, intent(in) :: nelv, iterb
    integer, dimension(nelv), intent(in) :: ref_level
    integer, dimension(2, nelv), intent(in) :: ref_family
    integer, dimension(nelv), intent(inout) :: ref_mark
    type(sem_lx_t), intent(inout) :: grid_min
    integer, intent(inout) ::  iter
    integer, intent(out) :: nmod
    integer, dimension(nelv) :: level_tmp
    logical, dimension(nelv) :: flag_exc
    integer :: ntot, iter_max, iter_balance, nmod_iter, itmp, il
    integer, dimension(2) :: nmod2
    character(len=LOG_SIZE) :: log_buf
    logical :: refine

    ntot = lx * lx* lx * nelv
    nmod  = 0
    iter_max = iter
    iter = 0

    ! Pressure preconditioner was found to be sensitive to the refinement
    ! topology, so some of the configurations should be corrected. The goal
    ! is to avoid single unrefined elements surrounded by refined ones and
    ! refined regions connected by edge or vertex only. This is just an
    ! approximate procedure, as mesh manager can have additional constraints
    ! not applied here. However, it is not an intention to replicate all
    ! possible constraints, but to fix the most obvious problems. In general
    ! order of applying operations can matter.

    flag_exc(:) = .false.

    ! Connect through face refinement regions sharing edge or vertex
    ! Start with current refinement level to remove existing problems
    ! This step can be exact
    call amr_ref_connect_fix(nelv, ref_level, .true., flag_exc, ref_mark, &
         grid_min, itmp)
    nmod = nmod + itmp

    do
       iter = iter + 1
       nmod_iter = 0
       iter_balance = iterb

       ! Balance the refinement flag
       call amr_ref_balance(nelv, ref_level, ref_family, ref_mark, grid_min, &
            iter_balance, nmod2)
       if (nmod2(2) .ne. 0) then
          write(log_buf, '(A,I9,1X,I3)') 'AMR tool mark check; tree cannot be &
               &balanced: ', nmod2(2), iter_balance
          call neko_log%message(log_buf, NEKO_LOG_INFO)
       end if

       ! Fix refined regions connected by edge/vertex only
       ! Use predicted refinement level to avoid possible future problems
       ! This is approximate only as not all additional constraints may be
       ! taken into account
       do il = 1, nelv
          select case(ref_mark(il))
          case(amr_flg_h_crs)
             level_tmp(il) = ref_level(il) - 1
          case(amr_flg_h_ref)
             level_tmp(il) = ref_level(il) + 1
          case(amr_flg_none)
             level_tmp(il) = ref_level(il)
          end select
       end do
       call amr_ref_connect_fix(nelv, level_tmp, .false., flag_exc, &
            ref_mark, grid_min, itmp)
       nmod_iter = nmod_iter + itmp

       ! Update predicted refinement level.
       do il = 1, nelv
          select case(ref_mark(il))
          case(amr_flg_h_crs)
             level_tmp(il) = ref_level(il) - 1
          case(amr_flg_h_ref)
             level_tmp(il) = ref_level(il) + 1
          case(amr_flg_none)
             level_tmp(il) = ref_level(il)
          end select
       end do

       ! Remove single unrefined element surrounded by refined ones
       call amr_ref_gap_remove(nelv, level_tmp, ref_mark, grid_min, nmod2)
       nmod_iter = nmod_iter + nmod2(1)

       nmod = nmod + nmod_iter
       if (nmod_iter .eq. 0) then
          exit
       else
          if (iter .eq. iter_max) then
             write(log_buf, '(A,I9)') 'AMR tool mark check; correcting not &
                  &finalised: ', nmod_iter
             call neko_log%message(log_buf, NEKO_LOG_INFO)
             exit
          end if
       end if
    end do

  end subroutine amr_ref_mark_check

  !> Remove nonconforming interfaces around marked elements
  !! @param[in]     nelv        local element number
  !! @param[in]     ref_level   element refinement level
  !! @param[in]     ref_flag    flag marking problematic elements
  !! @param[inout]  ref_mark    refinement flag
  !! @param[inout]  grid_min    minimal grid
  !! @param[in]     icomm       communication type
  !! @param[in]     ifrefall    refine only flag
  !! @param[out]    nmod        global number of modified elements
  subroutine amr_nonconf_int_remove(nelv, ref_level, ref_flag, ref_mark, &
       grid_min, icomm, ifrefall, nmod)
    integer, intent(in) :: nelv, icomm
    integer, dimension(nelv), intent(in) :: ref_level
    logical, dimension(nelv), intent(in) :: ref_flag
    integer, dimension(nelv), intent(inout) :: ref_mark
    type(sem_lx_t), intent(inout) :: grid_min
    logical, intent(in) :: ifrefall
    integer, intent(out) :: nmod
    real(rp), dimension(lx, lx, lx, nelv) :: exchange
    integer, dimension(lx, lx, lx, nelv) :: level
    integer :: ntot, itmp, il, jl

    ntot = lx * lx* lx * nelv
    nmod  = 0

    ! Distribute refinement level of marked elements
    level(:, :, :, :) = -1
    ! Depending on communication type fill
    select case(icomm)
    case (1)
       ! vertices only
       do il = 1, nelv
          if (ref_flag(il)) then
             do jl = 1, nvrt
                call amr_vertex_set(level(:, :, :, il), jl, ref_level(il))
             end do
          end if
       end do
    case (2)
       ! edges only
       do il = 1, nelv
          if (ref_flag(il)) then
             do jl = 1, nedg
                call amr_edge_set(level(:, :, :, il), jl, ref_level(il))
             end do
          end if
       end do
    case (3)
       ! faces only
       do il = 1, nelv
          if (ref_flag(il)) then
             do jl = 1, nfcs
                call amr_face_set(level(:, :, :, il), jl, ref_level(il))
             end do
          end if
       end do
    case (4)
       ! all
       do il = 1, nelv
          if (ref_flag(il)) level(:, :, :, il) = ref_level(il)
       end do
    case default
       call neko_error('AMR tools removing interf.; wrong communication type')
    end select
    exchange = real(level, rp)
    call grid_min%gs_Xh%gs_op_vector(exchange, ntot, GS_OP_MAX)
    level = nint(exchange)

    ! Mark neighbours of flagged elements
    if (ifrefall) then
       do il = 1, nelv
          if (.not. ref_flag(il)) then
             itmp = maxval(level(:, :, :, il))
             if (itmp .gt. ref_level(il)) then
                nmod = nmod + 1
                ref_mark(il) = amr_flg_h_ref
             end if
          end if
       end do
    else
       do il = 1, nelv
          if (.not. ref_flag(il)) then
             itmp = maxval(level(:, :, :, il))
             if (itmp .gt. ref_level(il)) then
                nmod = nmod + 1
                select case(ref_mark(il))
                case(amr_flg_h_crs)
                   ref_mark(il) = amr_flg_none
                case(amr_flg_none)
                   ref_mark(il) = amr_flg_h_ref
                end select
             end if
          end if
       end do
    end if

  end subroutine amr_nonconf_int_remove

end module amr_tools
