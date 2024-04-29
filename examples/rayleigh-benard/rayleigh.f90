module user
  use neko
  implicit none

  real(kind=rp) :: Ra
  real(kind=rp) :: Pr
  real(kind=rp) :: ta2 = 0

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%scalar_user_ic => set_ic
    u%fluid_user_f_vector => forcing
    u%material_properties => set_material_properties
    u%user_dirichlet_update => dirichlet_update
  end subroutine user_setup

  subroutine set_material_properties(t, tstep, rho, mu, cp, lambda, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params
    real(kind=rp) :: Re

    call json_get(params, "case.fluid.Ra", Ra)
    call json_get(params, "case.scalar.Pr", Pr)

    Re = 1.0_rp / Pr

    mu = 1.0_rp / Re
    lambda = mu / Pr
    rho = 1.0_rp
    cp = 1.0_rp
  end subroutine set_material_properties

  ! Boundary condition for the scalar
  subroutine dirichlet_update(field_bc_list, bc_bc_list, coef, t, tstep, which_solver)
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

    if (trim(which_solver) .eq. "scalar") then
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

  end subroutine dirichlet_update

  !> User initial condition
  subroutine set_ic(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: rand, z

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
                z = s%dof%z(i,j,k,e)
                s%x(i,j,k,e) = 1-z + 0.0001* rand*&
                                     sin(4*pi/4.5*s%dof%x(i,j,k,e)) &
                * sin(4*pi/4.5*s%dof%y(i,j,k,e))

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

  !> Forcing
  subroutine forcing(f, t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: i
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: rapr, ta2pr
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    s => neko_field_registry%get_field('s')
    rapr = Ra*Pr
    ta2pr = ta2*Pr

    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_cmult2(f%u_d,v%x_d,Ta2Pr,f%dm%size())
       call device_cmult2(f%v_d,u%x_d,Ta2Pr,f%dm%size())
       call device_cmult2(f%w_d,s%x_d,rapr,f%dm%size())
    else
       call cmult2(f%u,v%x,Ta2Pr,f%dm%size())
       call cmult2(f%v,u%x,Ta2Pr,f%dm%size())
       call cmult2(f%w,s%x,rapr,f%dm%size())
    end if
  end subroutine forcing
end module user
