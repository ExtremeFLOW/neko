module user
  use neko
  implicit none

  !> Variables to store the Rayleigh and Prandlt numbers
  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Pr = 0

  !> Arrays asociated with Method#1 for nusselt calculation
  type(field_t) :: work_field ! Field to perform operations
  type(field_t) :: uzt ! u_z * T
  real(kind=rp) :: bar_uzt ! Volume integral

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules => user_initialize
    u%user_finalize_modules => user_finalize
    u%fluid_user_ic => set_initial_conditions_for_u_and_s
    u%scalar_user_bc => set_scalar_boundary_conditions
    u%fluid_user_f_vector => set_bousinesq_forcing_term
    u%user_check => calculate_nusselt
  end subroutine user_setup

  subroutine set_scalar_boundary_conditions(s, x, y, z, nx, ny, nz, ix, iy, iz, ie)
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
    ! This will be used on all zones without labels
    ! e.g. the ones hardcoded to 'v', 'w', etcetc
    s = 1.0_rp - z
  end subroutine set_scalar_boundary_conditions

  subroutine set_initial_conditions_for_u_and_s(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(param_t), intent(inout) :: params
    type(field_t), pointer :: s
    integer :: i, j, k, e
    real(kind=rp) :: rand, r,z
    s => neko_field_registry%get_field('s')

    call rzero(u%x,u%dof%size())
    call rzero(v%x,v%dof%size())
    call rzero(w%x,w%dof%size())
    call rzero(s%x,w%dof%size())
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
    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_memcpy(s%x,s%x_d,s%dof%size(),HOST_TO_DEVICE)
    end if

  end subroutine set_initial_conditions_for_u_and_s

  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params

    ! Reset the relevant nondimensional parameters
    Pr = params%Pr
    Ra = params%Re
    params%Re = sqrt(Ra / Pr)

    ! Initialize variables related to nusselt calculation
    call field_init(work_field, u%dof, 'work_field')
    call field_init(uzt, u%dof, 'uzt')


  end subroutine user_initialize


  subroutine user_finalize(t, param)
    real(kind=rp) :: t
    type(param_t), intent(inout) :: param

    ! Finalize variables related to nusselt calculation
    call field_free(work_field)
    call field_free(uzt)

  
end subroutine user_finalize

  !> Forcing
  subroutine set_bousinesq_forcing_term(f, t)
    class(source_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: i
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: rapr, ta2pr
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    s => neko_field_registry%get_field('s')

    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_rzero(f%u_d,f%dm%size())
       call device_rzero(f%v_d,f%dm%size())
       call device_copy(f%w_d,s%x_d,f%dm%size())
    else
       call rzero(f%u,f%dm%size())
       call rzero(f%v,f%dm%size())
       call copy(f%w,s%x,f%dm%size())
    end if
  end subroutine set_bousinesq_forcing_term

  subroutine calculate_nusselt(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(field_t), pointer :: s
    integer :: n,ntot

    s => neko_field_registry%get_field('s')
    n = size(coef%B)
    ntot = coef%dof%size()

    !> ------ Method #1 for nusselt calculation -----
    !> Nu_v = 1 + sqrt(Ra) * <u_z * T>_{v,t}
    !> Steps:
       !1.    Calculate the convective current T*u_z
       !2.    Get <u_z * T>_{v}(t) (Average in volume)
       !!2.1. Multiply field with mass matrix
       !!2.2. Perform a global sum to get the integral
       !!2.3. Normalize with total volume to get the average
       !3.    Get <u_z * T>_{v,t} by averaging time signal
       !4.    Multiply average by sqrt(Ra) and sum 1
    !> Steps 1. and 2. are done here. Do 3. and 4. in post  
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_col3(uzt%x_d,w%x_d,s%x_d,n)             !1.  
       call device_col3(work_field%x_d,uzt%x_d,coef%B_d,n) !2.1.   
       bar_uzt = device_glsum(work_field%x_d,n)            !2.2.
       bar_uzt = bar_uzt / coef%volume                     !2.3. 
    else
       call col3(uzt%x,w%x,s%x,n)                          !1.
       call col3(work_field%x,uzt%x,coef%B,n)              !2.1.
       bar_uzt = glsum(work_field%x,n)                     !2.2.
       bar_uzt = bar_uzt / coef%volume                     !2.3.
    end if

    if (pe_rank .eq. 0) &
    &  write(*,*) &
    &  'Space average convective current = ', bar_uzt

    if (pe_rank .eq. 0) &
    &  write(*,*) &
    &  'n = ', n

    if (pe_rank .eq. 0) &
    &  write(*,*) &
    &  'ntot = ', ntot
  end subroutine calculate_nusselt


end module user
