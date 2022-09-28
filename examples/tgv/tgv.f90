! Taylor-Green vortex (TGV)
!
! Time-integration of the Taylor-Green vortex up to specified time, evaluating
! total kinetic energy and enstrophy. The resulting flow fields may be
! visualised using Paraview or VisIt by opening the field0.nek5000 file.
!
module user
  use neko
  implicit none

  ! Global user variables
  type(field_t) :: om1, om2, om3, w1, w2

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%fluid_usr_ic => user_ic
    user%usr_msh_setup => user_mesh_scale
    user%usr_chk => user_calc_quantities
    user%user_init_modules => user_initialize
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
    type(param_t), intent(inout) :: params
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

    uvw(1)   = sin(x)*cos(y)*cos(z)
    uvw(2)   = -cos(x)*sin(y)*cos(z)
    uvw(3)   = 0._rp
  end function tgv_ic

  ! User-defined initialization called just before time loop starts
  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params

    real(kind=rp) dt
    integer tstep

    ! initialize work arrays for postprocessing
    call field_init(om1, u%dof, 'omega1')
    call field_init(om2, u%dof, 'omega2')
    call field_init(om3, u%dof, 'omega3')
    call field_init(w1, u%dof, 'work1')
    call field_init(w2, u%dof, 'work1')

    ! call usercheck also for tstep=0
    dt = params%dt
    tstep = 0
    call user_calc_quantities(t, dt, tstep, u, v, w, p, coef)

  end subroutine user_initialize
 
  ! User-defined routine called at the end of every time step
  subroutine user_calc_quantities(t, dt, tstep,u, v, w, p, coef)
    real(kind=rp), intent(in) :: t, dt
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    integer :: ntot, i
    real(kind=rp) :: vv, sum_e1(1), e1, e2, sum_e2(1), oo

    if (mod(tstep,50).ne.0) return

    ntot = u%dof%size()

    call curl(om1, om2, om3, u, v, w, w1, w2, coef)

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
    call col3(w1%x,u%x,u%x,ntot)
    call addcol3(w1%x,v%x,v%x,ntot)
    call addcol3(w1%x,w%x,w%x,ntot)
    e1 = 0.5 * glsc2(w1%x,coef%B,ntot) / coef%volume

    call col3(w1%x,om1%x,om1%x,ntot)
    call addcol3(w1%x,om2%x,om2%x,ntot)
    call addcol3(w1%x,om3%x,om3%x,ntot)
    e2 = 0.5 * glsc2(w1%x,coef%B,ntot) / coef%volume
      
    if (pe_rank .eq. 0) &
         &  write(*,'(a,e18.9,a,e18.9,a,e18.9)') &
         &  'POST: t:', t, ' Ekin:', e1, ' enst:', e2
    
  end subroutine user_calc_quantities

end module user
