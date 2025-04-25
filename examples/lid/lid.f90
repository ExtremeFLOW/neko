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
    user%fluid_user_ic => user_ic
    user%fluid_user_if => user_bc
    user%user_check => user_calc_quantities
    user%user_init_modules => user_initialize
    user%user_finalize_modules => user_finalize
  end subroutine user_setup

  ! user-defined boundary condition
  subroutine user_bc(u, v, w, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
    real(kind=rp), intent(inout) :: u, v,w
    real(kind=rp), intent(in) :: x, y, z
    real(kind=rp), intent(in) :: nx, ny, nz
    integer, intent(in) :: ix, iy, iz, ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    real(kind=rp) lsmoothing
    lsmoothing = 0.1_rp ! length scale of smoothing at the edges

    u = step( x/lsmoothing ) * step( (1._rp-x)/lsmoothing )
    v = 0._rp
    w = 0._rp

  end subroutine user_bc

  ! User-defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u, v, w, p
    type(json_file), intent(inout) :: params

    ! set all zero velocity and pressure
    u = 0._rp
    v = 0._rp
    w = 0._rp
    p = 0._rp
  end subroutine user_ic

  ! User-defined initialization called just before time loop starts
  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u, v, w, p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    integer tstep
    character(len=50) :: mess

    ! Initialize the file object and create the output.csv file
    ! in the working directory
    output_file = file_init("ekin.csv")
    call output_file%set_header("# t,Ekin,enst")
    call vec_out%init(2) ! Initialize our vector with 2 elements (Ekin, enst)

    ! initialize work arrays for postprocessing
    call w1%init(u%dof, 'work1')
    call temp1%init(u%dof)
    call temp2%init(u%dof)
    call vort1%init(u%dof)
    call vort2%init(u%dof)
    call vort3%init(u%dof)

    ! read postprocessing interval
    call json_get(params, "case.fluid.ipostproc", ipostproc)
    write(mess,*) "postprocessing steps : ",ipostproc
    call neko_log%message(mess)

    ! call usercheck also for tstep=0
    tstep = 0
    call user_calc_quantities(t, tstep, u, v, w, p, coef, params)

  end subroutine user_initialize

  ! User-defined routine called at the end of every time step
  subroutine user_calc_quantities(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    type(field_t), intent(inout) :: u, v, w, p
    integer :: ntot, i
    real(kind=rp) :: ekin, enst

    ! do the postprocessing only every 50 time steps (including step 0)
    if (mod(tstep,ipostproc).ne.0) return

    ntot = u%dof%size()

    ! compute the kinetic energy
    ! field_ routines operate on the configured backend
    call field_col3(w1,u,u,ntot)
    call field_addcol3(w1,v,v,ntot)
    call field_addcol3(w1,w,w,ntot)
    ! For glsc2 we need to call the correct backend
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ekin = 0.5_rp * device_glsc2(w1%x_d,coef%B_d,ntot) / coef%volume
    else
       ekin = 0.5_rp * glsc2(w1%x,coef%B,ntot) / coef%volume
    end if

    ! compute enstrophy
    ! (the factor of 0.5 depends on the definition of enstrophy. We
    ! follow the reference paper by the HiOCFD4 workshop, but the book
    ! by Pope for instance would not include this factor)
    call curl(vort1,vort2,vort3, u, v, w, temp1, temp2, coef)
    call field_col3(w1,vort1,vort1,ntot)
    call field_addcol3(w1,vort2,vort2,ntot)
    call field_addcol3(w1,vort3,vort3,ntot)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       enst = 0.5_rp * device_glsc2(w1%x_d,coef%B_d,ntot) / coef%volume
    else
       enst = 0.5_rp * glsc2(w1%x,coef%B,ntot) / coef%volume
    end if


    ! output all this to file
    call neko_log%message("Writing csv file")
    vec_out%x = (/ekin, enst/)
    call output_file%write(vec_out, t)

    ! set the w component to zero to avoid any 3D instability
    ! in this quasi-2D flow
    w = 0._rp
    ! a perhaps more verbose alternative would be:
    !   use field_math, only: field_rzero
    !   call field_rzero(w, ntot)

  end subroutine user_calc_quantities

  ! User-defined finalization routine called at the end of the simulation
  subroutine user_finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    ! Deallocate the fields
    call w1%free()

    ! Deallocate output file and vector
    call file_free(output_file)
    call vec_out%free

  end subroutine user_finalize

  ! Smooth step function, with zero derivatives at 0 and 1
  function step(x)
    ! Taken from Simson code
    ! x<=0 : step(x) = 0
    ! x>=1 : step(x) = 1
    real(kind=rp) :: step, x

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
