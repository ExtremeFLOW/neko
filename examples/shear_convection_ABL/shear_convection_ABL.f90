! Linnea Huusko 27/2-2025
module user
  use neko
  implicit none

  real(kind=rp) :: u_geo

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    u_geo = 10.0_rp

    user%initial_conditions => user_ic
  end subroutine user_setup

  ! User defined initial condition
  subroutine user_ic(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: u, v, w, s
    type(field_t), pointer :: p
    type(dofmap_t), pointer :: dof
    integer :: i
    real(kind=rp) :: x, y, z
    real(kind=rp) :: eps, kx, ky, lx, ly, alpha, beta, gamma, delta, PI
    real(kind=rp) :: z1, z2, ze
    real(kind=rp) :: theta0, theta1, alpha_theta, gamma_theta

    PI = (4.*atan(1.))

    kx = 5.0_rp
    ky = 7.0_rp
    lx = 40.0_rp
    ly = 40.0_rp
    eps = 0.4_rp

    ! parameters for temperature profile
    alpha_theta = 0.08_rp
    gamma_theta = 0.003_rp
    z1 = 500.0_rp
    z2 = 550.0_rp
    theta0 = 300.0_rp
    theta1 = theta0 + alpha_theta * (z2 - z1)

    ! parameters for TKE profile
    ze = 600.0_rp


    alpha = kx * PI / 1500.0_rp
    beta = ky * PI / 1500.0_rp
    gamma = lx * PI / 1500.0_rp
    delta = ly * PI / 1500.0_rp
    if (scheme_name .eq. 'fluid') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")
       do i = 1, u%dof%size()
          u%x(i,1,1,1) = u_geo
          v%x(i,1,1,1) = 0.0_rp
          w%x(i,1,1,1) = 0.0_rp
          x = u%dof%x(i,1,1,1)
          y = u%dof%y(i,1,1,1)
          z = u%dof%z(i,1,1,1)
          if (z .le. 50) then ! Small perturbation to help get turbulence started
             u%x(i,1,1,1) = u%x(i,1,1,1) + eps*(sin(alpha*x)*sin(beta*y)) &
                  + eps*(sin(gamma*x)*sin(delta*y))
             v%x(i,1,1,1) = v%x(i,1,1,1) + eps*(sin(alpha*x)*sin(beta*y)) &
                  + eps*(sin(gamma*x)*sin(delta*y))
             w%x(i,1,1,1) = w%x(i,1,1,1) - eps*(alpha * cos(alpha*x)*sin(beta*y)) &
                  - eps*(gamma * cos(gamma*x)*sin(delta*y))
          endif
       end do
    else !scalar
       s => fields%get(scheme_name)
       if (scheme_name .eq. 'temperature') then
          do i = 1, s%dof%size()
             z = s%dof%z(i,1,1,1)
             if (z .le. z1) then
                s%x(i,1,1,1) = theta0
             elseif (z .le. z2) then
                s%x(i,1,1,1) = theta0 + alpha_theta * (z - z1)
             else
                s%x(i,1,1,1) = theta1 + gamma_theta * (z - z2)
             endif
          end do
       elseif (scheme_name .eq. 'TKE') then
          do i = 1, s%dof%size()
             z = s%dof%z(i,1,1,1)
             if (z .le. ze) then
                s%x(i,1,1,1) = 0.4_rp * (1 - z/ze)*(1 - z/ze)*(1 - z/ze)
             else
                s%x(i,1,1,1) = 0.0_rp
             endif
          end do
       endif
    endif
  end subroutine user_ic

end module user
