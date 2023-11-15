module user
  use neko
  use time_based_controller, only : time_based_controller_t
  use mean_field, only : mean_field_t
  use rbc, only : rbc_t
  implicit none

  !> Variables to store the Rayleigh and Prandlt numbers
  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Re = 0
  real(kind=rp) :: Pr = 0

  !> Type for run time
  type(rbc_t) :: runtime_rbc
     
  !> Mesh deformation
  real(kind=rp) :: delta_mesh = 0.0_rp !How much to deform to the top and bottom plate. =0 means no deformation

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_mesh_setup => deform_mesh
    u%user_init_modules => user_initialize
    u%user_finalize_modules => user_finalize
    u%fluid_user_ic => set_initial_conditions_for_u_and_s
    u%scalar_user_bc => set_scalar_boundary_conditions
    u%fluid_user_f_vector => set_bousinesq_forcing_term
    u%user_check => check
    u%material_properties => set_material_properties
  end subroutine user_setup

  subroutine set_material_properties(t, tstep, rho, mu, cp, lambda, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params

    call json_get(params, "case.fluid.Ra", Ra)
    call json_get(params, "case.scalar.Pr", Pr)

    Re = sqrt(Ra / Pr)
    mu = 1.0_rp / Re
    lambda = mu / Pr
    rho = 1.0_rp
    cp = 1.0_rp

  end subroutine set_material_properties


  subroutine deform_mesh(msh)
    type(mesh_t), intent(inout) :: msh
    msh%apply_deform => redistribute_elements
  end subroutine deform_mesh
  
  subroutine redistribute_elements(msh, x, y, z, lx, ly, lz)
    class(mesh_t) :: msh
    integer, intent(in) :: lx, ly, lz
    real(kind=rp), intent(inout) :: x(lx, lx, lx, msh%nelv)
    real(kind=rp), intent(inout) :: y(lx, lx, lx, msh%nelv)
    real(kind=rp), intent(inout) :: z(lx, lx, lx, msh%nelv)
    type(tuple_i4_t) :: el_and_facet
    real(kind=rp) :: th
    integer :: e, i, j ,k, l,  facet, nxyze
    real(kind=rp) :: Betaz, x_pt, y_pt, z_pt, z_pt_def
    
    nxyze = lx*ly*lz*msh%nelv
    if (delta_mesh .gt. 1e-8_rp) then
       Betaz = delta_mesh 
       do i = 1, nxyze
          z_pt = z(i,1,1,1)
          z_pt_def = 0.5*(tanh(Betaz*(2*z_pt-1.0))/tanh(Betaz) + 1.0)
          z(i,1,1,1) = z_pt_def
       end do
    end if
  end subroutine redistribute_elements
  
  subroutine set_scalar_boundary_conditions(s, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
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
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    !> Variables for bias
    real(kind=rp) :: arg, bias
    
    arg  = -tstep*0.0005
    bias = x*0.2*exp(arg)
    bias = 0

    ! This will be used on all zones without labels
    ! e.g. the ones hardcoded to 'v', 'w', etcetc
    s = 1.0_rp - z + bias


  end subroutine set_scalar_boundary_conditions

  subroutine set_initial_conditions_for_u_and_s(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    type(field_t), pointer :: s
    integer :: i, j, k, e
    real(kind=rp) :: rand, r,z
    s => neko_field_registry%get_field('s')

    !> Initialize with zero velocity
    call rzero(u%x,u%dof%size())
    call rzero(v%x,v%dof%size())
    call rzero(w%x,w%dof%size())

    !> Initialize with rand perturbations on temperature
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
    if (NEKO_BCKND_DEVICE .eq. 1) then
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
    type(json_file), intent(inout) :: params


    call runtime_rbc%init(t, u, v, w, p, coef, params)
    

  end subroutine user_initialize


  subroutine user_finalize(t, param)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: param

    call runtime_rbc%free()
  
  end subroutine user_finalize

  subroutine set_bousinesq_forcing_term(f, t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: i
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: rapr, ta2pr
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    s => neko_field_registry%get_field('s')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(f%u_d,f%dm%size())
       call device_rzero(f%v_d,f%dm%size())
       call device_copy(f%w_d,s%x_d,f%dm%size())
    else
       call rzero(f%u,f%dm%size())
       call rzero(f%v,f%dm%size())
       call copy(f%w,s%x,f%dm%size())
    end if
  end subroutine set_bousinesq_forcing_term

  subroutine check(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    logical :: get_spec_err_ind = .true.


    !> Determine if spectral error indicator needs to be calculated
    call json_get_or_default(params, 'case.rbc_sampler.get_spec_err_ind',&
                                          get_spec_err_ind, .true.)


    if (runtime_rbc%sample_control%check(t, tstep, .false.)) then
       
       !> Calculate at sampling time
       call runtime_rbc%calculate(t, tstep, coef, params, Ra, Pr, get_spec_err_ind)
       call runtime_rbc%update_stats(t)
       call runtime_rbc%get_integral_quantities(t, tstep, coef)

       !> Register the execution of the controller
       call runtime_rbc%sample_control%register_execution()

    end if

    if (runtime_rbc%field_write_control%check(t, tstep, .false.)) then
    
       !> Write fields
       call runtime_rbc%sync()
       call runtime_rbc%write_fields_and_stats(t)
       
       !> Register the execution of the controller
       call runtime_rbc%field_write_control%register_execution()
    
    end if

    !> Detmine if data should be streamed
    if (runtime_rbc%stream_data) then
       !> Stream data if controller works
       if (runtime_rbc%data_stream_control%check(t, tstep, .false.)) then

          !> stream data
          call runtime_rbc%dstream%stream(u,v,w,p,coef)
          
          !> Register the execution of the controller
          call runtime_rbc%data_stream_control%register_execution()
    
       end if
    end if
    
  end subroutine check

end module user
