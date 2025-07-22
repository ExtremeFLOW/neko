module user
  use neko
  implicit none

  !> Variables to store the Rayleigh and Prandlt numbers
  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Re = 0
  real(kind=rp) :: mu = 0
  real(kind=rp) :: Pr = 0

  !> =============================================

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%scalar_user_bc => set_scalar_boundary_conditions
    user%material_properties => set_material_properties
    user%user_startup => startup
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params

    call json_get(params, "case.fluid.Ra", Ra)
    call json_get(params, "case.scalar.Pr", Pr)
    Re = sqrt(Ra / Pr)
    mu = 1.0_rp / Re
  end subroutine startup

  ! Used here for demonstration purposes. Since the properties are
  ! actually const, it is better to set them directly in the startup routine,
  ! by adding the appropriate entries to the user file.
  subroutine set_material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time

    if (scheme_name .eq. "fluid") then
       call field_cfill(properties%get_by_name("fluid_rho"), 1.0_rp)
       call field_cfill(properties%get_by_name("fluid_mu"), mu)
    else if (scheme_name .eq. "scalar") then
       call field_cfill(properties%get_by_name("scalar_cp"), 1.0_rp)
       call field_cfill(properties%get_by_name("scalar_lambda"), mu / Pr)
    end if
  end subroutine set_material_properties


  subroutine set_scalar_boundary_conditions(dirichlet_field_list, dirichlet_bc, &
       coef, t, tstep)
    type(field_list_t), intent(inout) :: dirichlet_field_list
    type(field_dirichlet_t), intent(in) :: dirichlet_bc
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! Only do this at the first time step since our BCs are constants.
    if (tstep .ne. 1) return

    dof => dirichlet_field_list%dof(1)

    ! We know that it is the scalar calling the routine
    associate(s => dirichlet_field_bc_list%items(1)%ptr)

    ! We can set the whole field, it doesn't matter.
    ! Under the hood only the values on the bc will be set after a masked copy.
    s = 1.0_rp - dof%z

    end associate

  end subroutine set_scalar_boundary_conditions

  subroutine set_initial_conditions_for_s(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params
    integer :: i, j, k, e
    real(kind=rp) :: rand, r,z

    !> Initialize with rand perturbations on temperature
    call rzero(s%x,s%dof%size())
    do i = 1, s%dof%size()
       s%x(i,1,1,1) = 1-s%dof%z(i,1,1,1)
    end do
    ! perturb not on element boundaries
    ! Maybe not necessary, but lets be safe
    do e = 1, s%msh%nelv
       do k = 2,s%Xh%lx-1
          do j = 2,s%Xh%lx-1
             do i = 2,s%Xh%lx-1

                !call random_number(rand)
                !Somewhat random
                rand = cos(real(e+s%msh%offset_el,rp)*real(i*j*k,rp))
                r = sqrt(s%dof%x(i,j,k,e)**2+s%dof%y(i,j,k,e)**2)
                z = s%dof%z(i,j,k,e)
                s%x(i,j,k,e) = 1-z + 0.0001*rand*s%dof%x(i,j,k,e)*&
                     sin(3*pi*r/0.05_rp)*sin(10*pi*z)
             end do
          end do
       end do
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x, s%x_d, s%dof%size(), &
            HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine set_initial_conditions_for_s
end module user
