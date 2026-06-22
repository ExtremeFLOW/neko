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
  public :: amr_remove_gaps

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

    ! Add minimal grid
    this%msh => msh
    lx = 3
    call this%grid_min%init(this%msh, 3, .false., .false.)

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


  !> Remove refinement gaps of the size of a single element
  !! @param[in]     nelv        local element number
  !! @param[in]     ref_level   element refinement level
  !! @param[inout]  ref_mark    refinement flag
  !! @param[inout]  grid_min    minimal grid
  !! @param[inout]  iter        max/performed number of iterations
  !! @param[out]    nmod        global number of modified elements
  subroutine amr_remove_gaps(nelv, ref_level, ref_mark, grid_min, iter, nmod)
    integer, intent(in) :: nelv
    integer, dimension(nelv), intent(in) :: ref_level
    integer, dimension(nelv), intent(inout) :: ref_mark
    type(sem_lx_t), intent(inout) :: grid_min
    integer, intent(inout) ::  iter
    integer, intent(out) :: nmod
    integer, parameter :: lx = 3
    real(rp), dimension(lx, lx, lx, nelv) :: exchange
    integer, dimension(lx, lx, lx, nelv) :: level, mark_max, mark_min
    ! pairs of opposing faces
    integer, dimension(2, 3) :: fpair_lev, fpair_max, fpair_min
    integer :: ntot, il, jl, iter_max, nmod_iter, ierr
    logical :: refine

    ntot = lx * lx* lx * nelv
    nmod  = 0
    iter_max = iter

    ! Get current refinement level across interfaces; no interpolation
    do il = 1, nelv
       exchange(:, :, :, il) = ref_level(il)
    end do
    ! take care of nonconforming interfaces
    call grid_min%gs_Xh%interp%scale_children(exchange, 0.25_rp, 0.5_rp)
    call grid_min%gs_Xh%gs_op_vector(exchange, ntot, GS_OP_ADD)
    level = nint(exchange)
    ! subtract element's refinement level
    do il = 1, nelv
       level(:, :, :, il) = level(:, :, :, il) - ref_level(il)
    end do

    ! This part can be iterated
    iter = 0
    do
       iter = iter + 1
       ! Get min/max refinement mark across interfaces; no interpolation
       ! Possibly needed to take care of nonconforming interfaces
!       do il = 1, nelv
!          exchange(:, :, :, il) = ref_mark(il)
!       end do
!       call grid_min%gs_Xh%gs_op_vector(exchange, ntot, GS_OP_MIN)
!       mark_min = nint(exchange)
       do il = 1, nelv
          exchange(:, :, :, il) = ref_mark(il)
       end do
       call grid_min%gs_Xh%gs_op_vector(exchange, ntot, GS_OP_MAX)
       mark_max = nint(exchange)

       nmod_iter = 0

       do il = 1, nelv

!          call amr_oposite_face_extract(mark_min(:, :, :, il), fpair_min)

          refine = .false.

          ! gap in refinement level
          ! during first iteration only
          if (iter .eq. 1) then
             if (ref_mark(il) .eq. amr_flg_none) then
                call amr_oposite_face_extract(level(:, :, :, il), fpair_lev)
                do jl = 1, 3
                   if (fpair_lev(1, jl) .gt. ref_level(il) .and. &
                        fpair_lev(2, jl) .gt. ref_level(il) .and. &
                        fpair_max(1, jl) .eq. amr_flg_none .and. &
                        fpair_max(2, jl) .eq. amr_flg_none) then
                      refine = .true.
                   end if
                end do
             end if
          end if

          ! gap in refinement mark
          if (ref_mark(il) .eq. amr_flg_none) then
             call amr_oposite_face_extract(mark_max(:, :, :, il), fpair_max)
             do jl = 1, 3
                if (fpair_lev(1, jl) .eq. ref_level(il) .and. &
                     fpair_lev(2, jl) .eq. ref_level(il) .and. &
                     fpair_max(1, jl) .eq. amr_flg_h_ref .and. &
                     fpair_max(2, jl) .eq. amr_flg_h_ref) then
                   refine = .true.
                end if
             end do
          end if

          ! other tests?
          if (refine) then
             nmod_iter = nmod_iter + 1
             ref_mark(il) = amr_flg_h_ref
          end if
       end do

       ! global number of modified elements
       call MPI_Allreduce(MPI_IN_PLACE, nmod_iter, 1, MPI_INTEGER, MPI_SUM, &
            NEKO_COMM, ierr)

       nmod = nmod + nmod_iter
       if (nmod_iter .eq. 0 .or. iter .eq. iter_max) exit
    end do

  end subroutine amr_remove_gaps

  ! Extract opposite face from the minimal grid box
  subroutine amr_oposite_face_extract(box, face_pair)
    integer, parameter :: lx = 3, nface = 3, dim = 3
    integer, dimension(lx, lx, lx), intent(in) :: box
    integer, dimension(2, nface), intent(out) :: face_pair
    ! pairs of opposing faces in the minimal grid box
    integer, parameter, dimension(dim, 2, nface) :: fop =reshape([1, 2, 2, &
         3, 2, 2, 2, 1, 2, 2, 3, 2, 2, 2, 1, 2, 2, 3], shape(fop))
    integer :: il, jl

    do il = 1, nface
       do jl = 1, 2
          face_pair(jl, il) = box(fop(1, jl, il), fop(2, jl, il), &
               fop(3, jl, il))
       end do
    end do
  end subroutine amr_oposite_face_extract

end module amr_tools
