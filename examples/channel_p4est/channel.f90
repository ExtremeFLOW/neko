! Test for immersed boundary method
!
module user
  use neko
  use amr_reconstruct, only : amr_reconstruct_t, amr_flg_none, amr_flg_h_ref, &
       amr_flg_h_crs
  use amr_tools, only : amr_spectral_error_t, AMR_OP_MAX, AMR_OP_ELEN, &
       amr_ref_mark_check, amr_nonconf_int_remove, amr_thrsh_get
  use math, only : glmax
  use vector_list, only : vector_list_t

  use fld_file, only : fld_file_t
  use fld_file_output, only : fld_file_output_t

  use mpi_f08, only : MPI_Wtime, MPI_Bcast, MPI_Allreduce, MPI_IN_PLACE, &
       MPI_LOGICAL, MPI_LOR
  implicit none

  ! Global user variables
  ! sponge variables
  ! base flow velocities
  type(field_t) :: u_bf, v_bf, w_bf
  ! sponge profile
  type(field_t) :: spng_prf
  ! temporary field
  type(field_t) :: ftmp
  ! sponge parameters; starting point, width, strength
  real(rp), parameter :: spng_st = 4.0_rp, spng_wdth = 2.0_rp, &
       spng_str = 0.75_rp
  ! error indicator
  real(rp), dimension(:), allocatable :: errind
  type(amr_spectral_error_t) :: amr_error_ind_tool
  ! flags for refinement and error indicator collection phases
  logical :: if_refined = .false., if_eind = .false.
  ! Last refinement marking wall time
  real(dp) :: last_wall_time
  ! Max wall time period in seconds for averaging of error indicator
  real(dp), parameter :: wall_period = 60 * 60 * 0.5_dp
  ! Averaging time of error indicator in simulation units
  real(dp), parameter :: time_int = 2.0_dp
  ! Threshold for instantaneous to average indicator ratio to trigger rescue
  ! refinement
  real(rp), parameter :: pr_ratio = 5.0_rp
  ! Max accepted pressure error value for rescue refinement; HAND TUNED
  real(rp), parameter :: pr_err_cutoff = 0.05_rp
  ! min/max refinement levels
  integer, parameter :: ref_level_min = 0, ref_level_max = 2
  ! Fixed refinement thresholds for refinement and coarsening
  real(rp), parameter :: ref_thr = 1.0D-05, crs_thr = 1.0D-07
  ! Global max number of elements
  integer, parameter :: glb_max_elem = 20000

  ! for debugging
  type(fld_file_output_t) :: tst_fld
  real(rp) :: t

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initialize => user_initialize
    user%source_term => user_source_terms
    user%compute => user_compute
    user%finalize => user_finalize
    user%amr_refine_flag => amr_refine_flag
    user%amr_reconstruct => amr_reconstructing
  end subroutine user_setup

  ! User-defined initialization called just before time loop starts
  ! I use it to set up AMR refinement output
  subroutine user_initialize(time)
    type(time_state_t), intent(in) :: time
    integer :: il, n_simcomps
    real(kind=rp) :: t, rtmp
    type(field_t), pointer :: u, v, w, p
    character(len=NEKO_VARNAME_LEN), dimension(3) :: fld_name_ref
    character(len=NEKO_VARNAME_LEN), dimension(1) :: fld_name_mntr

    u => neko_registry%get_field('u')
    v => neko_registry%get_field('v')
    w => neko_registry%get_field('w')
    p => neko_registry%get_field('p')

    ! initialise sponge variables
    call u_bf%init(u%dof, 'u_bf')
    call v_bf%init(u%dof, 'v_bf')
    call w_bf%init(u%dof, 'w_bf')
    call spng_prf%init(u%dof, 'spng_prf')
    call ftmp%init(u%dof)

    do il = 1, u%dof%size()
       u_bf%x(il, 1, 1, 1) = 1.0_rp
       v_bf%x(il, 1, 1, 1) = 0.0_rp
       w_bf%x(il, 1, 1, 1) = 0.0_rp

       rtmp = (u%dof%x(il, 1, 1, 1) - spng_st) / spng_wdth
       spng_prf%x(il, 1, 1, 1) = spng_str * step(rtmp)
    end do

    ! initialise error indicator
    fld_name_ref(1) = 'u'
    fld_name_ref(2) = 'v'
    fld_name_ref(3) = 'w'
    fld_name_mntr(1) = 'p'
    call amr_error_ind_tool%init(u%msh, fld_name_ref = fld_name_ref, &
         fld_name_mntr = fld_name_mntr)

    ! allocate error indicator array
    allocate(errind(amr_error_ind_tool%nelv))

    ! initialise file output for debugging
    call tst_fld%init(rp, "testing_refine", 4)
    call tst_fld%fields%assign_to_ptr(1, p)
    call tst_fld%fields%assign_to_ptr(2, u)
    call tst_fld%fields%assign_to_ptr(3, v)
    call tst_fld%fields%assign_to_ptr(4, w)
    select type(file => tst_fld%file_%file_type)
    type is (fld_file_t)
       file%write_mesh = .true.
       file%skip_pressure = .false.
       file%skip_velocity = .false.
    end select

    ! testing refinement frequency based on the wall time
    last_wall_time = MPI_WTIME()

  end subroutine user_initialize

  ! Calculate case specific sponge
  subroutine user_source_terms(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: rhs_u, rhs_v, rhs_w, u, v, w

    return

    if (scheme_name .eq. 'fluid') then

       rhs_u => rhs%get_by_index(1)
       rhs_v => rhs%get_by_index(2)
       rhs_w => rhs%get_by_index(3)

       u => neko_registry%get_field('u')
       v => neko_registry%get_field('v')
       w => neko_registry%get_field('w')

       ! this is just a hack; should be done differently in the future
       call field_sub3(ftmp, u_bf, u)
       call field_col2(ftmp, spng_prf)
       call field_add2(rhs_u, ftmp)

       call field_sub3(ftmp, v_bf, v)
       call field_col2(ftmp, spng_prf)
       call field_add2(rhs_v, ftmp)

       call field_sub3(ftmp, w_bf, w)
       call field_col2(ftmp, spng_prf)
       call field_add2(rhs_w, ftmp)

    end if

  end subroutine user_source_terms

  ! User-defined routine called at the end of every time step
  subroutine user_compute(time)
    type(time_state_t), intent(in) :: time

  end subroutine user_compute

  ! User-defined finalization routine called at the end of the simulation
  subroutine user_finalize(time)
    type(time_state_t), intent(in) :: time

    ! Deallocate the fields
    call u_bf%free()
    call v_bf%free()
    call w_bf%free()
    call spng_prf%free()
    call ftmp%free()

    ! Checkpoint averaged indicator fields before finalising
    call amr_error_ind_tool%checkpoint()

    call amr_error_ind_tool%free()

    deallocate(errind)

  end subroutine user_finalize

  ! Set refinement flag
  subroutine amr_refine_flag(time, nelv, ref_level, family, ref_mark, ifrefine)
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: nelv
    integer, dimension(nelv), intent(in) :: ref_level
    integer, dimension(2, nelv), intent(in) :: family
    integer, dimension(nelv), intent(inout) :: ref_mark
    logical, intent(inout) :: ifrefine
    integer :: il, ierr, iter, iterb, nmod, el_ref, el_dlt
    real(rp) :: pr_in_max, pr_av_max, ref_thr_var, rtmp, rgnelv
    real(dp) :: time_av
    logical :: ifrescue, ifwall
    logical, dimension(nelv) :: refine_flag
    character(len=LOG_SIZE) :: log_buf
    ! Number of intervals and max iteration count for threshold calculation
    integer, parameter :: nint = 500, it_max = 50
    ! Ratios for global new elements and acceptable error
    real(rp), parameter :: ratio_el = 0.2, ratio_dlt = 0.01
    ! geometrical context
    real(rp) :: x_pos
    ! Define region with refinement enabled
    real(rp), parameter :: x_out = 4.0

    real(rp) :: dst

    ! geometry based refinement
    associate(dm_Xh => amr_error_ind_tool%grid_min%dm_Xh)
      if (time%tstep .eq. 5) then
         ! flag elements; just geometrical context
         do il = 1, nelv
            ! sphere
            !dst = sqrt((dm_Xh%x(2, 2, 2, il) -1.0_rp)**2 + &
            !     (dm_Xh%y(2, 2, 2, il) - 0.5_rp)**2 + &
            !     (dm_Xh%z(2, 2, 2, il) - 0.5_rp)**2)
            !if (dst .le. 0.3_rp) then
            !   ref_mark(il) = amr_flg_h_ref
            !end if
            ! box
            !if (dm_Xh%x(2, 2, 2, il) .ge. 0.75_rp .and. &
            !     dm_Xh%x(2, 2, 2, il) .le. 1.25_rp .and. &
            !     dm_Xh%y(2, 2, 2, il) .ge. 0.25_rp .and. &
            !     dm_Xh%y(2, 2, 2, il) .le. 0.75_rp .and. &
            !     dm_Xh%z(2, 2, 2, il) .ge. 0.25_rp .and. &
            !     dm_Xh%z(2, 2, 2, il) .le. 0.75_rp) then
            !   ref_mark(il) = amr_flg_h_ref
            !end if
            ! column
            if (dm_Xh%x(2, 2, 2, il) .ge. 5.75_rp .and. &
                 dm_Xh%x(2, 2, 2, il) .le. 6.25_rp .and. &
                 dm_Xh%y(2, 2, 2, il) .ge. 0.25_rp .and. &
                 dm_Xh%y(2, 2, 2, il) .le. 0.75_rp .and. &
                 dm_Xh%z(2, 2, 2, il) .ge. 0.0_rp .and. &
                 dm_Xh%z(2, 2, 2, il) .le. 0.75_rp) then
               ref_mark(il) = amr_flg_h_ref
            end if
            ! section
            !if (dm_Xh%x(2, 2, 2, il) .ge. 0.75_rp .and. &
            !     dm_Xh%x(2, 2, 2, il) .le. 1.25_rp) then
            !   ref_mark(il) = amr_flg_h_ref
            !end if
         end do
         ifrefine = .true.

      end if
    end associate

    ! for monitoring
    if (ifrefine) then
       ! save error indicator output
       t = time%t
       call amr_error_ind_tool%sample(ref_mark, t)

       ! for debugging
       call tst_fld%sample(t)
    end if

    ! set refinement flag
    if_refined = .false.

  end subroutine amr_refine_flag

  subroutine amr_reconstructing(reconstruct, counter, time)
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter
    type(time_state_t), intent(in) :: time

    call u_bf%amr_restart(reconstruct, counter, time)
    call v_bf%amr_restart(reconstruct, counter, time)
    call w_bf%amr_restart(reconstruct, counter, time)
    call spng_prf%amr_restart(reconstruct, counter, time)
    call ftmp%amr_reallocate(reconstruct, counter, time)

    call amr_error_ind_tool%amr_restart(reconstruct, counter, time)

    ! reallocate arrays
    if (reconstruct%nold .ne. reconstruct%nnew) then
       if (allocated(errind)) then
          deallocate(errind)
          allocate(errind(reconstruct%nnew))
          errind( :) = 0.0_rp
       end if
    end if

    ! mark refinement step
    if_refined = .true.
    if_eind = .false.

    ! for debugging
    call tst_fld%sample(t)

  end subroutine amr_reconstructing

  ! Smooth step function, with zero derivatives at 0 and 1
  function step(x)
    ! Taken from Simson code
    ! x<=0 : step(x) = 0
    ! x>=1 : step(x) = 1
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: step

    if (x.le.0.02_rp) then
       step = 0.0_rp
    else
       if (x.le.0.98_rp) then
          step = 1._rp/( 1._rp + exp(1._rp/(x-1._rp) + 1._rp/x) )
       else
          step = 1._rp
       end if
    end if

  end function step

end module user
