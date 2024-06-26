! Lid-driven cavity
!
! Time-integration of the lid-driven cavity with smoothened
! belt velocity to fulfil continuity equation.
!
module user
  use neko
  implicit none

  ! Global user variables
  type(field_t) :: w1

  type(file_t) output_file ! output file
  type(vector_t) :: vec_out    ! will store our output data
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
    real(kind=rp), intent(inout) :: u
    real(kind=rp), intent(inout) :: v
    real(kind=rp), intent(inout) :: w
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    real(kind=rp) lsmoothing
    lsmoothing = 0.05_rp    ! length scale of smoothing at the edges

    u = step( x/lsmoothing ) * step( (1._rp-x)/lsmoothing )
    v = 0._rp
    w = 0._rp
    
  end subroutine user_bc
  
  ! User-defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    u = 0._rp
    v = 0._rp
    w = 0._rp
    p = 0._rp
  end subroutine user_ic
  
  ! User-defined initialization called just before time loop starts
  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    integer tstep

    ! Initialize the file object and create the output.csv file
    ! in the working directory
    output_file = file_init("output.csv")
    call output_file%set_header("t,Ekin,enst")
    call vec_out%init(2) ! Initialize our vector with 2 elements (Ekin, enst)

    ! initialize work arrays for postprocessing
    call w1%init(u%dof, 'work1')

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
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(field_t), pointer :: om1, om2, om3
    integer :: ntot, i
    real(kind=rp) :: e1, e2

    if (mod(tstep,50).ne.0) return

    om1 => neko_field_registry%get_field("omega_x")
    om2 => neko_field_registry%get_field("omega_y")
    om3 => neko_field_registry%get_field("omega_z")
    
    ntot = u%dof%size()

    call col3(w1%x,u%x,u%x,ntot)
    call addcol3(w1%x,v%x,v%x,ntot)
    call addcol3(w1%x,w%x,w%x,ntot)
    e1 = 0.5 * glsc2(w1%x,coef%B,ntot) / coef%volume

    call col3(w1%x,om1%x,om1%x,ntot)
    call addcol3(w1%x,om2%x,om2%x,ntot)
    call addcol3(w1%x,om3%x,om3%x,ntot)
    e2 = 0.5 * glsc2(w1%x,coef%B,ntot) / coef%volume

    ! Output all this to file
    call neko_log%message("Writing csv file")
    vec_out%x = (/e1, e2/)
    call output_file%write(vec_out, t)

  end subroutine user_calc_quantities

  ! Smooth step function
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

end module user
