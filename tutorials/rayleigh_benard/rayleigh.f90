module user
  use neko
  implicit none

  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Re = 0
  real(kind=rp) :: Pr = 0

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%scalar_user_ic => set_ic
    u%material_properties => set_material_properties
  end subroutine user_setup

  subroutine set_material_properties(t, tstep, rho, mu, cp, lambda, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params
    real(kind=rp) :: Re

    call json_get(params, "case.fluid.Ra", Ra)
    call json_get(params, "case.fluid.Pr", Pr)

    Re = sqrt(Ra / Pr)
    mu = 1.0_rp / Re
    lambda = mu / Pr
    rho = 1.0_rp
    cp = 1.0_rp
  end subroutine set_material_properties

  !> User initial condition
  subroutine set_ic(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: rand, x, y, z, pert, p1, p2, p3, p4

    do i = 1, s%dof%size()
       s%x(i,1,1,1) = 1-s%dof%z(i,1,1,1)
    end do
    
    do e = 1, s%msh%nelv
       do k = 2,s%Xh%lx-1
          do j = 2,s%Xh%lx-1
             do i = 2,s%Xh%lx-1

                rand = 1_rp
                x = s%dof%x(i,j,k,e)
                y = s%dof%y(i,j,k,e)
                z = s%dof%z(i,j,k,e)
                p1 = 1_rp*sin(3_rp*pi*x)
                p2 = 0.6_rp*sin(7_rp*pi*x)
                p3 = 0.3_rp*sin(24_rp*pi*x)
                p4 = 0.1_rp*sin(5_rp*pi*x)
                pert = (p1 + p2 + p3 + p4)
                s%x(i,j,k,e) = 1_rp-z + 0.01_rp*pert*rand

            end do
          end do
       end do
    end do

    if ((NEKO_BCKND_DEVICE .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_memcpy(s%x, s%x_d, s%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if


  end subroutine set_ic

end module user
