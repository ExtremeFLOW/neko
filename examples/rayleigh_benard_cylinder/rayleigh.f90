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
    user%dirichlet_conditions => dirichlet_conditions
    user%material_properties => material_properties
    user%startup => startup
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
  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time

    if (scheme_name .eq. "fluid") then
       call field_cfill(properties%get("fluid_rho"), 1.0_rp)
       call field_cfill(properties%get("fluid_mu"), mu)
    else if (scheme_name .eq. "temperature") then
       call field_cfill(properties%get("temperature_cp"), 1.0_rp)
       call field_cfill(properties%get("temperature_lambda"), mu / Pr)
    end if
  end subroutine material_properties


  !> Sets the stratification of temperature along the cylinder walls.
  subroutine dirichlet_conditions(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: s
    integer :: i

    ! Only do this at the first time step since our BCs are constants.
    if (time%tstep .ne. 1) return

    ! We know that it is the scalar calling the routine
    s => fields%get("temperature")

    do i = 1, bc%msk(0)
       s%x(bc%msk(i),1,1,1) = &
            1.0_rp - s%dof%z(bc%msk(i),1,1,1)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x, s%x_d, s%size(), HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine dirichlet_conditions

  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    integer :: i, j, k, e
    real(kind=rp) :: rand, r,z
    type(field_t), pointer :: s

    ! See scalar.name in the case file, makes sure that we only
    ! run this for the scalar field.
    if (scheme_name .ne. 'temperature') return

    s => fields%items(1)%ptr

    !> Initialize with rand perturbations on temperature
    do i = 1, s%dof%size()
       s%x(i,1,1,1) = 1 - s%dof%z(i,1,1,1)
    end do

    ! perturb not on element boundaries
    ! Maybe not necessary, but lets be safe
    do e = 1, s%msh%nelv
       do k = 2, s%Xh%lx - 1
          do j = 2, s%Xh%lx - 1
             do i = 2, s%Xh%lx - 1

                !call random_number(rand)
                !Somewhat random
                rand = cos(real(e + s%msh%offset_el,rp) * real(i*j*k, rp))
                r = sqrt(s%dof%x(i,j,k,e)**2 + s%dof%y(i,j,k,e)**2)
                z = s%dof%z(i,j,k,e)
                s%x(i,j,k,e) = 1.0_rp - z + 0.0001_rp*rand*s%dof%x(i,j,k,e)* &
                     sin(3.0_rp*pi*r/0.05_rp) * sin(10.0_rp*pi*z)
             end do
          end do
       end do
    end do

  end subroutine initial_conditions
end module user
