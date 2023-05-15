! Martin Karp 13/3-2023
module user
  use neko
  implicit none
contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_user_ic => user_ic
    u%user_mesh_setup => user_mesh_scale
  end subroutine user_setup

  ! Rescale mesh, we create a mesh with some refinement close to the wall. mesh size (4*pi,2*delta,4/3*pi) 
  ! New mesh can easily be genreated with genmeshbox
  ! OBS refinement is not smooth and the constant valuas are a bit ad hoc. 
  ! Stats converge close to reference DNS
  subroutine user_mesh_scale(msh)
    type(mesh_t), intent(inout) :: msh
    integer :: i, p, nvert
    real(kind=rp) :: d, y, viscous_layer, visc_el_h, el_h, center_el_h, dist_from_wall
    integer :: el_in_visc_lay, el_in_y

    el_in_y = 18
    el_in_visc_lay = 2
    viscous_layer = 0.0888889
    el_h = 2.0_rp/el_in_y
    visc_el_h = viscous_layer/el_in_visc_lay
    center_el_h = (1.0_rp-viscous_layer)/(el_in_y/2-el_in_visc_lay)
    

    nvert = size(msh%points)
    do i = 1, nvert
       msh%points(i)%x(1) = pi*msh%points(i)%x(1) 
       y = msh%points(i)%x(2) 
       if ((1-abs(y)) .le. (el_in_visc_lay*el_h)) then
          dist_from_wall = (1-abs(y))/el_h*visc_el_h
       else
          dist_from_wall = viscous_layer + (1-abs(y)-el_in_visc_lay*el_h)/el_h*center_el_h
       end if
       if (y .gt. 0) msh%points(i)%x(2) = 1.0_rp - dist_from_wall
       if (y .lt. 0) msh%points(i)%x(2) = -1.0_rp + dist_from_wall
       msh%points(i)%x(3) = pi*msh%points(i)%x(3) 
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
       uvw = channel_ic(u%dof%x(i,1,1,1),u%dof%y(i,1,1,1),u%dof%z(i,1,1,1))
       u%x(i,1,1,1) = uvw(1)
       v%x(i,1,1,1) = uvw(2)
       w%x(i,1,1,1) = uvw(3)
    end do
  end subroutine user_ic

  ! Kind of brute force with rather large initial disturbances 
  function channel_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: ux, uy, uz, eps, Re_tau, yp, Re_b, alpha, beta
    real(kind=rp) :: C, k, kx, kz

      Re_tau = 180
      C      = 5.17
      k      = 0.41
      Re_b   = 2860

      yp = (1-y)*Re_tau
      if (y.lt.0) yp = (1+y)*Re_tau
      
      ! Reichardt function
      ux  = 1/k*log(1.0+k*yp) + (C - (1.0/k)*log(k)) * &
            (1.0 - exp(-yp/11.0) - yp/11*exp(-yp/3.0))
      ux  = ux * Re_tau/Re_b

      eps = 3e-2
      kx  = 16
      kz  = 7

      alpha = kx * 2*PI/(pi*4)
      beta  = kz * 2*PI/(pi*4/3)

      ! add perturbation to trigger turbulence 
      uvw(1)  = ux  + eps*beta  * sin(alpha*x)*cos(beta*z) 
      uvw(2)  =       eps       * sin(alpha*x)*sin(beta*z)
      uvw(3)  =      -eps*alpha * cos(alpha*x)*sin(beta*z)

  end function channel_ic

end module user
