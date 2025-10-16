! Two-dimensional lid-driven cavity
!
! Time-integration of the lid-driven cavity with smoothened
! belt velocity to fulfil the continuity equation in the corners.
!
! Note that the domain is actually 3D with width one element. In order
! to prevent any instability in the z direction, the w velocity is
! set to zero at every step. This is needed for higher Reynolds numbers.
!
module user
  use neko
  implicit none

  ! Global user variables
  type(field_t) :: w1
  type(field_t) :: temp1,temp2
  type(field_t) :: vort1,vort2,vort3

  type(file_t) output_file ! output file
  type(vector_t) :: vec_out ! will store our output data
  integer :: ipostproc ! frequency of the output

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%dirichlet_conditions => dirichlet_conditions
    user%compute => compute
    user%initialize =>initialize
    user%finalize => finalize
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    character(len=50) :: mess

    ! read postprocessing interval
    call json_get(params, "case.fluid.ipostproc", ipostproc)
    write(mess,*) "postprocessing steps : ",ipostproc
    call neko_log%message(mess)
  end subroutine startup

  ! user-defined Dirichlet boundary condition
  subroutine dirichlet_conditions(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: u, v, w
    integer :: i
    real(kind=rp) :: x

    real(kind=rp) lsmoothing

    if (time%tstep .ne. 1) return

    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")

    lsmoothing = 0.1_rp ! length scale of smoothing at the edges

    do i = 1, bc%msk(0)
       x = u%dof%x(bc%msk(i), 1, 1, 1)
       u%x(bc%msk(i), 1, 1, 1) = &
            step( x/lsmoothing ) * step( (1._rp - x)/lsmoothing )

       v%x(bc%msk(i), 1, 1, 1) = 0
       w%x(bc%msk(i), 1, 1, 1) = 0
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, u%size(), HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(v%x, v%x_d, v%size(), HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(w%x, w%x_d, w%size(), HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine dirichlet_conditions

  ! User-defined initialization called just before time loop starts
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: u

    ! Initialize the file object and create the output.csv file
    ! in the working directory
    call output_file%init("ekin.csv")
    call output_file%set_header("# t,Ekin,enst")
    call vec_out%init(2) ! Initialize our vector with 2 elements (Ekin, enst)

    ! initialize work arrays for postprocessing
    u => neko_field_registry%get_field("u")
    call w1%init(u%dof, 'work1')
    call temp1%init(u%dof)
    call temp2%init(u%dof)
    call vort1%init(u%dof)
    call vort2%init(u%dof)
    call vort3%init(u%dof)


    ! call usercheck also for tstep=0
    call compute(time)

  end subroutine initialize

  ! User-defined routine called at the end of every time step
  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    integer :: ntot, i
    real(kind=rp) :: ekin, enst
    type(field_t), pointer :: u, v, w
    type(coef_t), pointer :: coef


    ! do the postprocessing only every 50 time steps (including step 0)
    if (mod(time%tstep, ipostproc) .ne. 0) return

    coef => neko_user_access%case%fluid%c_Xh
    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")

    ntot = u%dof%size()

    ! compute the kinetic energy
    ! field_ routines operate on the configured backend
    call field_col3(w1, u, u, ntot)
    call field_addcol3(w1, v, v, ntot)
    call field_addcol3(w1, w, w, ntot)
    ! For glsc2 we need to call the correct backend
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ekin = 0.5_rp * device_glsc2(w1%x_d, coef%B_d,ntot) / coef%volume
    else
       ekin = 0.5_rp * glsc2(w1%x, coef%B, ntot) / coef%volume
    end if

    ! compute enstrophy
    ! (the factor of 0.5 depends on the definition of enstrophy. We
    ! follow the reference paper by the HiOCFD4 workshop, but the book
    ! by Pope for instance would not include this factor)
    call curl(vort1,vort2,vort3, u, v, w, temp1, temp2, coef)
    call field_col3(w1, vort1, vort1, ntot)
    call field_addcol3(w1, vort2, vort2, ntot)
    call field_addcol3(w1, vort3, vort3, ntot)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       enst = 0.5_rp * device_glsc2(w1%x_d, coef%B_d, ntot) / coef%volume
    else
       enst = 0.5_rp * glsc2(w1%x, coef%B, ntot) / coef%volume
    end if


    ! output all this to file
    call neko_log%message("Writing csv file")
    vec_out%x = (/ekin, enst/)
    call output_file%write(vec_out, time%t)

    ! set the w component to zero to avoid any 3D instability
    ! in this quasi-2D flow
    call field_rzero(w)

  end subroutine compute

  ! User-defined finalization routine called at the end of the simulation
  subroutine finalize(time)
    type(time_state_t), intent(in) :: time

    ! Deallocate the fields
    call w1%free()

    ! Deallocate output file and vector
    call file_free(output_file)
    call vec_out%free

  end subroutine finalize

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
