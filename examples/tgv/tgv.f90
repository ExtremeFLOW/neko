! Taylor-Green vortex (TGV)
!
! Time-integration of the Taylor-Green vortex up to specified time, evaluating
! total kinetic energy and enstrophy. The resulting flow fields may be
! visualised using Paraview or VisIt by opening the field0.nek5000 file.
!
module user
  use neko
  use simcomp_test, only : simcomp_test_t
  implicit none

  ! Global user variables
  type(field_t) :: w1

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%fluid_user_ic => user_ic
    user%user_mesh_setup => user_mesh_scale
    user%user_check => user_calc_quantities
    user%user_init_modules => user_initialize
    user%user_finalize_modules => user_finalize
  end subroutine user_setup

  ! Rescale mesh
  subroutine user_mesh_scale(msh)
    type(mesh_t), intent(inout) :: msh
    integer :: i, p, nvert
    real(kind=rp) :: d
    d = 4._rp

    ! original mesh has size 0..8 to be mapped onto -pi..pi
    ! will be updated later to a method giving back the vertices of the mesh
    nvert = size(msh%points)
    do i = 1, nvert
       msh%points(i)%x(1) = (msh%points(i)%x(1) - d) / d * pi
       msh%points(i)%x(2) = (msh%points(i)%x(2) - d) / d * pi
       msh%points(i)%x(3) = (msh%points(i)%x(3) - d) / d * pi
    end do

  end subroutine user_mesh_scale

  ! User-defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i, ntot
    real(kind=rp) :: uvw(3)

    ! u%dof%size() gives the total number of collocation points per rank
    ntot = u%dof%size()
    do i = 1, ntot
       uvw = tgv_ic(u%dof%x(i,1,1,1),u%dof%y(i,1,1,1),u%dof%z(i,1,1,1))
       u%x(i,1,1,1) = uvw(1)
       v%x(i,1,1,1) = uvw(2)
       w%x(i,1,1,1) = uvw(3)
    end do
    p = 0._rp
  end subroutine user_ic

  function tgv_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: ux, uy, uz
    real(kind=rp) :: uvw(3)

    uvw(1) = sin(x)*cos(y)*cos(z)
    uvw(2) = -cos(x)*sin(y)*cos(z)
    uvw(3) = 0._rp
  end function tgv_ic

  ! User-defined initialization called just before time loop starts
  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    real(kind=rp) dt
    integer tstep
    integer :: n_simcomps
    type(json_file) :: comp_subdict
    type(json_core) :: core
    type(json_value), pointer :: simcomp_object
    logical :: found
    type(simcomp_test_t), allocatable :: my_simcomp


    ! ! Probably we can put this into some %add function in the simcomp_exector.
    ! ! However, it won't know the new simcomp type, so one would have to create
    ! ! it here and inside %add there would be an allocate(..., source=...)
    ! n_simcomps = neko_simcomps%n_simcomps
    ! ! User simcomps always executed last
    ! neko_simcomps%order(n_simcomps + 1) = n_simcomps + 1
    ! allocate(simcomp_test_t::neko_simcomps%simcomps(n_simcomps + 1)%simcomp)
    ! neko_simcomps%n_simcomps = n_simcomps + 1

    ! Extracting json. I put an array under "user_simcomps", here we just
    ! assume there is only 1 item there.
    call params%get_core(core)
    call params%get('case.user_simcomps', simcomp_object, found)
    ! Extract the 1st item and create a separate json from it
    call json_extract_item(core, simcomp_object, 1, comp_subdict)
    ! Need to add order for the constructor of the simcomp, but the value is
    ! irrelevant.
    ! call comp_subdict%add("order", 1)
    ! ! Have to assume the case is the same, because case is not passed to this
    ! ! routine, but this is probably very future-proof.
    ! call neko_simcomps%simcomps(n_simcomps + 1)%simcomp%init(comp_subdict, &
    !                                                          neko_simcomps%simcomps(n_simcomps)%simcomp%case)

    ! Allocate a simulation component
    allocate(my_simcomp)
    call my_simcomp%init(comp_subdict, neko_simcomps%simcomps(1)%simcomp%case)
    ! Add the simulation component to the list
    call neko_simcomps%add(my_simcomp)

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
    type(field_t), pointer :: omega_x, omega_y, omega_z
    integer :: ntot, i
    real(kind=rp) :: vv, sum_e1(1), e1, e2, sum_e2(1), oo, e3

    if (mod(tstep, 50) .ne. 0) return

    omega_x => neko_field_registry%get_field("omega_x")
    omega_y => neko_field_registry%get_field("omega_y")
    omega_z => neko_field_registry%get_field("omega_z")

    ntot = u%dof%size()

!    Option 1:
!    sum_e1 = 0._rp
!    sum_e2 = 0._rp
!    do i = 1, ntot
!       vv = u%x(i,1,1,1)**2 + v%x(i,1,1,1)**2 + w%x(i,1,1,1)**2
!       oo = om1%x(i,1,1,1)**2 + om2%x(i,1,1,1)**2 + om3%x(i,1,1,1)**2
!       sum_e1 = sum_e1 + vv*coef%B(i,1,1,1)
!       sum_e2 = sum_e2 + oo*coef%B(i,1,1,1)
!    end do
!    e1 = 0.5 * glsum(sum_e1,1) / coef%volume
!    e2 = 0.5 * glsum(sum_e2,1) / coef%volume

!    Option 2:
!    do i = 1, ntot
!       w1%x(i,1,1,1) = u%x(i,1,1,1)**2 + v%x(i,1,1,1)**2 + w%x(i,1,1,1)**2
!       w2%x(i,1,1,1) = om1%x(i,1,1,1)**2 + om2%x(i,1,1,1)**2 + om3%x(i,1,1,1)**2
!    end do
!    e1 = 0.5 * glsc2(w1%x,coef%B,ntot) / coef%volume
!    e2 = 0.5 * glsc2(w2%x,coef%B,ntot) / coef%volume

!    Option 3:
!    w1%x = u%x**2 + v%x**2 + w%x**2
!    w2%x = om1%x**2 + om2%x**2 + om3%x**2
!    e1 = 0.5 * glsc2(w1%x,coef%B,ntot) / coef%volume
!    e2 = 0.5 * glsc2(w2%x,coef%B,ntot) / coef%volume

!    Option 4:
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col3(w1%x_d, u%x_d, u%x_d, ntot)
       call device_addcol3(w1%x_d, v%x_d, v%x_d, ntot)
       call device_addcol3(w1%x_d, w%x_d, w%x_d, ntot)
       e1 = 0.5 * device_glsc2(w1%x_d, coef%B_d, ntot) / coef%volume

       call device_col3(w1%x_d, omega_x%x_d, omega_x%x_d, ntot)
       call device_addcol3(w1%x_d, omega_x%x_d, omega_y%x_d, ntot)
       call device_addcol3(w1%x_d, omega_z%x_d, omega_z%x_d, ntot)
       e2 = 0.5 * device_glsc2(w1%x_d, coef%B_d, ntot) / coef%volume
    else
       call col3(w1%x, u%x, u%x, ntot)
       call addcol3(w1%x, v%x, v%x, ntot)
       call addcol3(w1%x, w%x, w%x, ntot)
       e1 = 0.5 * glsc2(w1%x, coef%B, ntot) / coef%volume

       call col3(w1%x, omega_x%x, omega_x%x, ntot)
       call addcol3(w1%x, omega_y%x, omega_y%x, ntot)
       call addcol3(w1%x, omega_z%x, omega_z%x, ntot)
       e2 = 0.5 * glsc2(w1%x, coef%B, ntot) / coef%volume
    end if

    if (pe_rank .eq. 0) &
      & write(*,'(a,e18.9,a,e18.9,a,e18.9)') &
      & 'POST: t:', t, ' Ekin:', e1, ' enst:', e2

  end subroutine user_calc_quantities

  ! User-defined finalization routine called at the end of the simulation
  subroutine user_finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    ! Deallocate the fields
    call w1%free()

  end subroutine user_finalize

end module user
