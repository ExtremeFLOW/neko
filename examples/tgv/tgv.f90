module user
  use neko
  implicit none
contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_usr_ic => user_ic
    u%usr_msh_setup => user_mesh_scale
    u%usr_chk => usr_calc_quantities
  end subroutine user_setup

  ! Normalize mesh
  
  subroutine user_mesh_scale(msh)
    type(mesh_t), intent(inout) :: msh
    type(point_t), pointer :: pt
    integer :: i, p, npts
    real(kind=rp) :: pi, d
    pi = 4d0*atan(1d0)
    d = 4d0

    npts = size(msh%points)
    do i = 1, npts
       msh%points(i)%x(1) = (msh%points(i)%x(1) - d) / d * pi
       msh%points(i)%x(2) = (msh%points(i)%x(2) - d) / d * pi
       msh%points(i)%x(3) = (msh%points(i)%x(3) - d) / d * pi
    end do
    
  end subroutine user_mesh_scale



  ! User defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(param_t), intent(inout) :: params
    integer :: i
    real(kind=rp) :: uvw(3)

    do i = 1, u%dof%size()
       uvw = tgv_ic(u%dof%x(i,1,1,1),u%dof%y(i,1,1,1),u%dof%z(i,1,1,1))
       u%x(i,1,1,1) = uvw(1)
       v%x(i,1,1,1) = uvw(2)
       w%x(i,1,1,1) = uvw(3)
    end do
    p = real(0d0,rp)
  end subroutine user_ic
  
  function tgv_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: ux, uy, uz
    real(kind=rp) :: uvw(3)
    integer e,eg

    uvw(1)   = sin(x)*cos(y)*cos(z)
    uvw(2)   = -cos(x)*sin(y)*cos(z)
    uvw(3)   = real(0d0,rp)
  end function tgv_ic

  subroutine usr_calc_quantities( t, dt, tstep,u, v, w, p, coef)
    real(kind=rp), intent(in) :: t, dt
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(field_t) :: omg1, omg2, omg3, w1, w2
    integer :: n, i
    real(kind=rp) :: vv, sum_e1(1), e1, e2, sum_e2(1), oo
    n = u%dof%n_dofs

    if (mod(tstep,50).ne.0) return

    sum_e1 = 0d0
    sum_e2 = 0d0
    call field_init(omg1, u%dof, 'omg1')
    call field_init(omg2, u%dof, 'omg2')
    call field_init(omg3, u%dof, 'omg3')
    call field_init(w1, u%dof, 'work1')
    call field_init(w2, u%dof, 'work1')

    call curl(omg1, omg2, omg3, u, v, w, w1, w2, coef)
    do i = 1,n
       vv = u%x(i,1,1,1)**2 + v%x(i,1,1,1)**2 + w%x(i,1,1,1)**2
       oo = omg1%x(i,1,1,1)**2 + omg2%x(i,1,1,1)**2 + omg3%x(i,1,1,1)**2 
       sum_e1 = sum_e1 + vv*coef%B(i,1,1,1) 
       sum_e2 = sum_e2 + oo*coef%B(i,1,1,1) 
    enddo
    e1 = 0.5 * glsum(sum_e1,1) / coef%volume
    e2 = 0.5 * glsum(sum_e2,1) / coef%volume
    if (pe_rank .eq. 0) write(*,*) 'Time:', t, 'Kinetic energy:', e1, 'Enstrophy:', e2

    call field_free(omg1)
    call field_free(omg2)
    call field_free(omg3)
    call field_free(w1)
    call field_free(w2)

  end subroutine usr_calc_quantities


end module user
