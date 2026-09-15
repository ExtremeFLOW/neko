! Two-dimensional flow past circular cylinder cavity
!
! Note that the domain is actually 3D with width one element. In order
! to prevent any instability in the z direction, the w velocity is
! set to zero at every step. This is needed for higher Reynolds numbers.
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
  real(rp), parameter :: spng_st = 24.0_rp, spng_wdth = 6.0_rp, &
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
  real(dp), parameter :: time_int = 0.01_dp
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
    character(len=NEKO_VARNAME_LEN), dimension(2) :: fld_name_ref
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

    end if

  end subroutine user_source_terms

  ! User-defined routine called at the end of every time step
  subroutine user_compute(time)
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: w

    ! set the w component to zero to avoid any 3D instability
    ! in this quasi-2D flow
    w => neko_registry%get_field("w")
    call field_rzero(w)

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

    ! Get global max of instantaneous and averages pressure error indicators
    pr_in_max = glmax(amr_error_ind_tool%eind_in_mntr%items(1)%ptr%x, nelv)
    pr_av_max = glmax(amr_error_ind_tool%eind_av_mntr%items(1)%ptr%x, nelv)
    ! check if rescue refinement is needed; HAND TUNED
    ifrescue = (pr_in_max .gt. (min(pr_av_max, pr_err_cutoff) * pr_ratio ))

    if (ifrescue) then
       ! rescue refinement is based on instantaneous presser error indicator
       errind(:) = amr_error_ind_tool%eind_in_mntr%items(1)%ptr%x(:)

       ! flag problematic elements
       ref_thr_var = min(pr_av_max, pr_err_cutoff)
       do il = 1, nelv
          if(errind(il) .gt. ref_thr_var) then
             refine_flag(il) = .true.
          else
             refine_flag(il) = .false.
          endif
       enddo

       ! remove nonconforming interfaces next to marked elements
       call amr_nonconf_int_remove(nelv, ref_level, refine_flag, ref_mark, &
            amr_error_ind_tool%grid_min, 1, .true., nmod)

       ifrefine = .true.
       write(log_buf, '(A,3E15.7)') 'Rescue refinement: ', pr_in_max, &
            pr_av_max, ref_thr_var
       call neko_log%message(log_buf)
    else
       ! Averaging time for error indicator; simulation time
       time_av = amr_error_ind_tool%get_collect_time()

       ! refinement frequency based on the wall time
       ifwall = ((MPI_WTIME() - last_wall_time) .gt. wall_period)
       call MPI_Bcast(ifwall, 1, MPI_LOGICAL, 0, NEKO_COMM, ierr)

       if (time_av .gt. time_int .or. ifwall) then
          ! refinement frequency based on the wall time
          last_wall_time = MPI_WTIME()

          if (if_eind .and. .not. if_refined) then
             ! error indicator collection phase already completed, but no
             ! refinement performed
             if_eind = .false.

             ! reset averaged error indicator
             call amr_error_ind_tool%reset_average()
          else
             ! new error indicator collection phase completed
             if_eind = .true.

             ! get combined error indicator
             call amr_error_ind_tool%err_av_get(AMR_OP_ELEN, errind, nelv)

             ! Variable threshold for refinement based on assumed number of
             ! elements to be refined
             ! Global max number of elements
             rtmp = real(glb_max_elem, rp)
             ! global free space ratio
             rgnelv = real(amr_error_ind_tool%msh%glb_nelv, rp)
             rtmp = sqrt(max(0.0_rp, rtmp - rgnelv) / rtmp)
             ! take into account number of children
             rtmp = rtmp / 7.0_rp
             ! Expected number of elements to be refined and allowed difference
             el_ref = int(rgnelv * ratio_el * rtmp)
             el_dlt = max(int(el_ref * ratio_dlt), 3)
             ! Elements at max refinement level cannot be refined, so
             ! shouldn't be counted
             do il = 1, nelv
                if(ref_level(il) .lt. ref_level_max) then
                   refine_flag(il) = .true.
                else
                   refine_flag(il) = .false.
                endif
             enddo
             call amr_thrsh_get(el_ref, el_dlt, errind, refine_flag, nelv, &
                  nint, it_max, .false., ref_thr_var, nmod, iter)

             ! Do not refine forever
             ref_thr_var = max(ref_thr_var, ref_thr)

             write(log_buf, '(A,2E15.7)') 'Ref./crs. thresholds: ', &
                  ref_thr_var, crs_thr
             call neko_log%message(log_buf)

             ! for geometry based refinement
             associate(dm_Xh => amr_error_ind_tool%grid_min%dm_Xh, &
                  msh => amr_error_ind_tool%grid_min%msh)

               ! flag elements
               ! Error indicator based part
               do il = 1, nelv
                  if (errind(il) .gt. ref_thr_var .and. &
                       ref_level(il) .lt. ref_level_max) then
                     ref_mark(il) = amr_flg_h_ref
                     ifrefine = .true.
                  else if (errind(il) .lt. crs_thr .and. &
                       ref_level(il) .gt. ref_level_min) then
                     ref_mark(il) = amr_flg_h_crs
                     ifrefine = .true.
                  else
                     ref_mark(il) = amr_flg_none
                  end if
               end do
             end associate

             ! Global test
             call MPI_Allreduce(MPI_IN_PLACE, ifrefine, 1, MPI_LOGICAL, &
                  MPI_LOR, NEKO_COMM, ierr)

             if (ifrefine) then
                ! check consistency of refinement regions
                iter = 10
                iterb = 10
                call amr_ref_mark_check(nelv, ref_level, family, ref_mark, &
                     amr_error_ind_tool%grid_min, iter, iterb, nmod)
                write(log_buf, '(A,2I9)') 'Refinement mark check: ', iter, nmod
                call neko_log%message(log_buf)
             end if
          end if
       end if
    end if

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
