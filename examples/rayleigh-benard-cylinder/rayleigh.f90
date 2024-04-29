module user
  use neko
  implicit none

  !> Variables to store the Rayleigh and Prandlt numbers
  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Re = 0
  real(kind=rp) :: Pr = 0

  !> =============================================

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%scalar_user_ic => set_initial_conditions_for_s
    u%user_dirichlet_update => set_scalar_boundary_conditions
    u%material_properties => set_material_properties
  end subroutine user_setup

  subroutine set_material_properties(t, tstep, rho, mu, cp, lambda, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params
    real(kind=rp) :: Re

    call json_get(params, "case.fluid.Ra", Ra)
    call json_get(params, "case.scalar.Pr", Pr)


    Re = sqrt(Ra / Pr)
    mu = 1.0_rp / Re
    lambda = mu / Pr
    rho = 1.0_rp
    cp = 1.0_rp
  end subroutine set_material_properties

  ! Boudnary condition for the scalar
  subroutine set_scalar_boundary_conditions(field_bc_list, bc_bc_list, coef, t, tstep, which_solver)
    type(field_list_t), intent(inout) :: field_bc_list
    type(bc_list_t), intent(inout) :: bc_bc_list
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    character(len=*), intent(in) :: which_solver

    integer :: i
    real(kind=rp) :: z

    ! Only do this at the first time step since our BCs are constants.
    if (tstep .ne. 1) return

    ! Check that we are being called by `fluid`
    if (trim(which_solver) .eq. "fluid") then
    else if (trim(which_solver) .eq. "scalar") then

       associate( s => field_bc_list%items(1)%ptr, s_bc => bc_bc_list%bc(1)%bcp)

         !
         ! Perform operations on the scalar field here
         ! Note that we are checking if the field is allocated, in
         ! case the boundary is empty.
         !
         if (allocated(s%x)) then

            do i = 1, s_bc%msk(0)
               z = s_bc%dof%z(s_bc%msk(i), 1, 1, 1)
               s%x(s_bc%msk(i), 1, 1, 1) = 1.0_rp - z
            end do

         end if
       end associate

    end if

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
