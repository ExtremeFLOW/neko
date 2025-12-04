! Copyright (c) 2024-2025, The Neko Authors
! All rights reserved.
!
! Zalesak disk test case - slotted disk rotation
!
module user
  use neko
  implicit none

  real(kind=rp) :: eps
  real(kind=rp) :: gamma, u_max

contains

  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%compute => compute
    user%source_term => source_term
    user%material_properties => material_properties
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  !> User startup routine
  subroutine startup(params)
    type(json_file), intent(inout) :: params

    call json_get(params, "case.scalar.epsilon", eps)
    call json_get(params, "case.scalar.gamma", gamma)
    u_max = 1.0_rp
  end subroutine startup

  !> User compute routine called at each timestep
  subroutine compute(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v
    integer :: i
    real(kind=rp) :: x, y, u_tmp(1)

    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")

    u_tmp = 0
    do i = 1, u%size()
      x = u%dof%x(i,1,1,1)
      y = u%dof%y(i,1,1,1)
      ! Solid body rotation velocity field
      u%x(i,1,1,1) = pi*(0.5_rp - y)/3.14_rp
      v%x(i,1,1,1) = pi*(x - 0.5_rp)/3.14_rp
      u_tmp(1) = max(u_tmp(1), sqrt(u%x(i,1,1,1)**2 + v%x(i,1,1,1)**2))
    end do
    u_max = glmax(u_tmp, 1)
  end subroutine compute

  !> User source term for the scalar equation
  subroutine source_term(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time

    integer :: i
    type(field_t), pointer :: rhs_s, s
    real(kind=rp) :: absgrad
    integer :: ind(4)
    type(field_t), pointer :: work1, work2, work3, work4
    type(coef_t), pointer :: coef

    if (scheme_name .ne. 's') return

    ! Get the right-hand side field for the scalar
    rhs_s => rhs%get_by_index(1)

    ! Get the scalar field and coefficients
    s => neko_field_registry%get_field('s')
    coef => neko_user_access%case%fluid%c_Xh

    ! Request scratch fields
    call neko_scratch_registry%request_field(work1, ind(1))
    call neko_scratch_registry%request_field(work2, ind(2))
    call neko_scratch_registry%request_field(work3, ind(3))
    call neko_scratch_registry%request_field(work4, ind(4))

    ! Compute gradient of scalar field
    call grad(work1%x, work2%x, work3%x, s%x, coef)

    ! Apply gather-scatter and multiplicity
    call coef%gs_h%op(work1, GS_OP_ADD)
    call coef%gs_h%op(work2, GS_OP_ADD)
    call coef%gs_h%op(work3, GS_OP_ADD)
    call col2(work1%x, coef%mult, work4%size())
    call col2(work2%x, coef%mult, work4%size())
    call col2(work3%x, coef%mult, work4%size())

    ! Compute normalized gradient and apply phase field forcing
    do i = 1, work4%size()
      absgrad = sqrt(work1%x(i,1,1,1)**2 + work2%x(i,1,1,1)**2 + &
                     work3%x(i,1,1,1)**2)
      if (absgrad == 0.0_rp) then
        print *, 'warning, absgrad==', absgrad
        absgrad = 1e21_rp
      end if

      work1%x(i,1,1,1) = -s%x(i,1,1,1)*(1.0_rp - s%x(i,1,1,1))* &
                         (work1%x(i,1,1,1)/absgrad)
      work2%x(i,1,1,1) = -s%x(i,1,1,1)*(1.0_rp - s%x(i,1,1,1))* &
                         (work2%x(i,1,1,1)/absgrad)
      work3%x(i,1,1,1) = -s%x(i,1,1,1)*(1.0_rp - s%x(i,1,1,1))* &
                         (work3%x(i,1,1,1)/absgrad)
    end do

    ! Compute divergence of the normalized gradient
    call dudxyz(work4%x, work1%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call copy(rhs_s%x, work4%x, work4%size())
    call dudxyz(work4%x, work2%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call add2(rhs_s%x, work4%x, work4%size())
    call dudxyz(work4%x, work3%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call add2(rhs_s%x, work4%x, work4%size())

    ! Scale by gamma * u_max
    absgrad = gamma * u_max
    call cmult(rhs_s%x, absgrad, work4%size())

    ! Release scratch fields
    call neko_scratch_registry%relinquish_field(ind)
  end subroutine source_term

  !> User material properties
  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time
    real(kind=rp) :: delta

    delta = u_max * gamma
    if (scheme_name .eq. "fluid") then
      call field_cfill(properties%get("fluid_rho"), 1.0_rp)
      call field_cfill(properties%get("fluid_mu"), 1.0_rp)
    else if (scheme_name .eq. "s") then
      call field_cfill(properties%get('s_cp'), 1.0_rp)
      call field_cfill(properties%get('s_lambda'), eps*delta)
    end if
  end subroutine material_properties

  !> User initial conditions for scalar field
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type(field_t), pointer :: s
    real(kind=rp) :: rad, center_x, center_y, radius, notch_width, &
                     notch_height, notch_x, notch_y
    integer :: i

    if (scheme_name .ne. 's') return

    s => fields%items(1)%ptr

    ! Zalesak disk parameters
    center_x = 0.5_rp
    center_y = 0.75_rp
    radius = 0.15_rp

    ! Notch dimensions from the problem description
    notch_width = 0.05_rp
    notch_height = 0.05_rp

    ! Set the circular disk with notch using smooth tanh transitions
    do i = 1, s%dof%size()
      ! Create disk with tanh profile
      rad = sqrt((s%dof%x(i,1,1,1) - center_x)**2 + &
                 (s%dof%y(i,1,1,1) - center_y)**2)
      s%x(i,1,1,1) = 0.5_rp*(1.0_rp + tanh((radius - rad)/(2.0_rp*eps)))

      ! Cut out notch with smooth tanh transition
      notch_x = 0.5_rp*(1.0_rp + tanh((abs(s%dof%x(i,1,1,1) - center_x) - &
                notch_width/2.0_rp)/(2.0_rp*eps)))
      notch_y = 0.5_rp*(1.0_rp + tanh((s%dof%y(i,1,1,1) - &
                (center_y + notch_height))/(2.0_rp*eps)))

      ! Combine disk with notch
      s%x(i,1,1,1) = s%x(i,1,1,1) * (notch_x + (1.0_rp - notch_x)*notch_y)
    end do

    if ((NEKO_BCKND_DEVICE .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
        .or. (NEKO_BCKND_OPENCL .eq. 1)) then
      call device_memcpy(s%x, s%x_d, s%dof%size(), &
                         HOST_TO_DEVICE, sync=.false.)
    end if
  end subroutine initial_conditions

end module user
