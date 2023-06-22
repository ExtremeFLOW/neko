module user
  use neko
  use json_module, only : json_file
  implicit none

  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Pr = 0
  real(kind=rp) :: ta2 = 0

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules => set_Pr
    u%fluid_user_ic => set_ic
    u%fluid_user_f_vector => forcing
    u%scalar_user_bc => scalar_bc
  end subroutine user_setup

  subroutine scalar_bc(s, x, y, z, nx, ny, nz, ix, iy, iz, ie)
    real(kind=rp), intent(inout) :: s
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    ! If we set scalar_bcs(*) = 'user' instead 
    ! this will be used instead on that zone
    s = 1.0_rp-z
  end subroutine scalar_bc
 
  !> User initial condition
  subroutine set_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    type(field_t), pointer :: s
    integer :: i, e, k, j
    real(kind=rp) :: rand, z
    s => neko_field_registry%get_field('s')

    call rzero(u%x,u%dof%size())
    call rzero(v%x,v%dof%size())
    call rzero(w%x,w%dof%size())
    
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

    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_memcpy(s%x,s%x_d,s%dof%size(),HOST_TO_DEVICE)
    end if


  end subroutine set_ic

  subroutine set_Pr(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    logical :: found
    real(kind=rp) :: Re
    ! Reset the relevant nondimensional parameters
    ! Pr = input Pr
    ! Ra = input Re
    ! Re = 1/Pr

    call params%get('case.fluid.Re', Ra, found)
    call params%get('case.scalar.Pr', Pr, found)
    Re = 1._rp / Pr
    call params%update('case.fluid.Re', Re, found)
    return
  end subroutine set_Pr



  !> Forcing
  subroutine forcing(f, t)
    class(source_t), intent(inout) :: f
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
