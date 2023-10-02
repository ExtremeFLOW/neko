module user
  use neko
  use time_based_controller, only : time_based_controller_t
  implicit none

  !> Variables to store the Rayleigh and Prandlt numbers
  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Re = 0
  real(kind=rp) :: Pr = 0

  !> Fields and variables for calculations
  type(field_t) , target:: work_field ! Field to perform operations
  type(field_t) , target:: uzt ! u_z * T
  type(field_t) , target:: dtdx ! Derivative of scalar wrt x
  type(field_t) , target:: dtdy ! Detivative of scalar wrt y
  type(field_t) , target:: dtdz ! Derivative of scalar wrt z
  type(field_t) , target:: dtdn ! Derivative of scalar wrt normal
  type(field_t) , target:: div_dtdX ! divergence of heat flux

  type(field_t) , target:: dudx ! Derivative wrt x
  type(field_t) , target:: dudy ! Detivative wrt y
  type(field_t) , target:: dudz ! Derivative wrt z
  type(field_t) , target:: dvdx ! Derivative wrt x
  type(field_t) , target:: dvdy ! Detivative wrt y
  type(field_t) , target:: dvdz ! Derivative wrt z
  type(field_t) , target:: dwdx ! Derivative wrt x
  type(field_t) , target:: dwdy ! Detivative wrt y
  type(field_t) , target:: dwdz ! Derivative wrt z

  type(field_t) , target:: eps_k ! Detivative wrt y
  type(field_t) , target:: eps_t ! Derivative wrt z
  
  type(field_t) , target:: mass_area_top ! mass matrix for area on top wall
  type(field_t) , target:: mass_area_bot ! mass matrix for area on bottom wall
  type(field_t) , target:: mass_area_side ! mass matrix for area on top wall
  type(field_t) , target:: bdry_mask ! mass matrix for area on top wall
  type(field_t) , target:: bm1 ! mass matrix for area on bottom wall
  type(field_t) , target:: normal_x ! normal vector in x (non zero at the boundaries)
  type(field_t) , target:: normal_y ! normal vector in y (non zero at the boundaries)
  type(field_t) , target:: normal_z ! normal vector in z (non zero at the boundaries)
  real(kind=rp) :: bar_uzt ! Volume integral
  real(kind=rp) :: bar_div_dtdX ! Volume integral
  real(kind=rp) :: top_wall_bar_dtdz ! area integral
  real(kind=rp) :: bot_wall_bar_dtdz ! area integral
  real(kind=rp) :: top_wall_area ! area of top and bottom wall
  real(kind=rp) :: bot_wall_area ! area of top and bottom wall
  real(kind=rp) :: side_wall_area ! area of top and bottom wall

  real(kind=rp) :: heat_bot_wall
  real(kind=rp) :: heat_top_wall
  real(kind=rp) :: heat_side_wall
  real(kind=rp) :: heat_balance
  
  !> field list for outputs
  type(field_list_t) :: area_l
  type(field_list_t) :: dtdX_l
  type(field_list_t) :: eps_l

  !> Boundary conditions  
  integer :: istep = 1

  !> Variables to write extra files
  character(len=NEKO_FNAME_LEN) :: fname_dtdX
  character(len=NEKO_FNAME_LEN) :: fname_dtdn
  character(len=NEKO_FNAME_LEN) :: fname_area
  character(len=NEKO_FNAME_LEN) :: fname_bm1
  character(len=NEKO_FNAME_LEN) :: fname_eps
  type(file_t) :: mf_dtdn
  type(file_t) :: mf_dtdX
  type(file_t) :: mf_eps
  type(file_t) :: mf_area
  type(file_t) :: mf_bm1

  !> List of elements and facets in uper and lower boundary
  type(stack_i4t4_t) :: top_wall_facet
  type(stack_i4t4_t) :: bot_wall_facet
  type(stack_i4t4_t) :: side_wall_facet
  type(stack_i4t4_t) :: wall_facet
  type(tuple4_i4_t)  :: facet_el_type_0

  !> Spectral error indicator
  type(spectral_error_indicator_t) :: spectral_error_indicator

  !> Mesh deformation
  real(kind=rp) :: delta_mesh = 0.0_rp !How much to deform to the top and bottom plate. =0 means no deformation
    
  !> Calculation controllers
  type(time_based_controller_t), allocatable :: controllers(:)
    
    

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
    u%user_check => calculate_nusselt
  end subroutine user_setup

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
    
    arg  = -istep*0.0005
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
    type(json_file), intent(inout) :: params

    !> Parameters to populate the list of elements and facets
    real(kind=rp) :: stack_size(1) !Do it like this so it is cast
    real(kind=rp) :: stack_size_global
    integer :: e, facet, typ, e_counter, facet_counter, i,j,k,n
    real(kind=rp) :: normal(3)
    real(kind=rp) :: area_mass, area_rank(1),area_global(1)
    type(tuple4_i4_t), pointer :: wall_facet_list(:)
    type(tuple4_i4_t), pointer :: sidewall_facet_list(:)
    real(kind=rp) :: voll(1), voll_temp(1)
    integer :: lx, ly, lz, nones
    real(kind=rp) :: ones(u%Xh%lx,u%Xh%ly,6,u%msh%nelv)

    !> Controller data
    real(kind=rp)     :: monitor_write_par
    character(len=:), allocatable  :: monitor_write_control
    real(kind=rp)     :: T_end


    !> Recalculate the non dimensional parameters
    call json_get(params, 'case.scalar.Pr', Pr)
    call json_get(params, 'case.fluid.Re', Re)
    Ra = (Re**2)*Pr
    if (pe_rank.eq.0) write(*,*) 'Rayleigh Number is Ra=', Ra
    
    !> Initialize variables related to nusselt calculation
    call work_field%init( u%dof, 'work_field')
    call uzt%init( u%dof, 'uzt')
    call dtdx%init( u%dof, 'dtdx')
    call dtdy%init( u%dof, 'dtdy')
    call dtdz%init( u%dof, 'dtdz')
    call dtdn%init( u%dof, 'dtdn')
    
    call dudx%init( u%dof, 'dudx')
    call dudy%init( u%dof, 'dudy')
    call dudz%init( u%dof, 'dudz')
    call dvdx%init( u%dof, 'dvdx')
    call dvdy%init( u%dof, 'dvdy')
    call dvdz%init( u%dof, 'dvdz')
    call dwdx%init( u%dof, 'dwdx')
    call dwdy%init( u%dof, 'dwdy')
    call dwdz%init( u%dof, 'dwdz')


    call eps_k%init( u%dof, 'eps_k')
    call eps_t%init( u%dof, 'eps_t')

    call div_dtdX%init( u%dof, 'div_dtdX')
    call mass_area_top%init( u%dof, 'mat')
    call mass_area_bot%init( u%dof, 'mab')
    call mass_area_side%init( u%dof, 'masd')
    call bdry_mask%init( u%dof, 'bdry_mask')
    call bm1%init( u%dof, 'mass_mat')
    call normal_x%init( u%dof, 'normal_x')
    call normal_y%init( u%dof, 'normal_y')
    call normal_z%init( u%dof, 'normal_z')

    !> Initialize the file
    fname_dtdX = 'dtdX.fld'
    fname_eps = 'eps.fld'
    fname_dtdn = 'dtdn.fld'
    fname_area = 'area.fld'
    fname_bm1 = 'bm1.fld'
    mf_dtdX =  file_t(fname_dtdX)
    mf_eps =  file_t(fname_eps)
    mf_dtdn =  file_t(fname_dtdn)
    mf_area =  file_t(fname_area)
    mf_bm1 =  file_t(fname_bm1)

    !> Initialize the list of pointers to write out files
    call list_init3(area_l,mass_area_top,mass_area_bot, &
                         mass_area_side)
    call list_init3(dtdX_l,dtdx,dtdy, &
                         dtdz)
    call list_init3(eps_l,eps_t,eps_k, &
                         dtdz)

    !> Initialize list to identify relevant facets in boundaries
    call wall_facet%init()
    call top_wall_facet%init()
    call bot_wall_facet%init()
    call side_wall_facet%init()

    !> Populate the list with upper and lower and wall facets 
    do e = 1, u%msh%nelv !Go over all elements
      do facet = 1, 6 ! Go over all facets of hex element
        normal = coef%get_normal(1,1,1,e,facet) ! Get facet normal
        typ =  u%msh%facet_type(facet,e)             ! Get facet type
        if (typ.ne.0) then !then it is a boundary facet
          if ((normal(3)).ge.0.999) then !then it is on the top plate
            facet_el_type_0%x = (/facet, e, typ, 0/)
            call wall_facet%push(facet_el_type_0)
            call top_wall_facet%push(facet_el_type_0)
          else if ((normal(3)).le.(-0.999)) then !then it is on the bottom
            facet_el_type_0%x = (/facet, e, typ, 0/)
            call wall_facet%push(facet_el_type_0)
            call bot_wall_facet%push(facet_el_type_0)
          else !then it is on the sidewall
            facet_el_type_0%x = (/facet, e, typ, 0/)
            call wall_facet%push(facet_el_type_0)
            call side_wall_facet%push(facet_el_type_0)
          end if
        end if
      end do
    end do
    !! Determine the number of facets in the walls
    stack_size = real(wall_facet%top_,kind=rp)
    stack_size_global = glsum(stack_size,1)
    if (pe_rank .eq. 0) then
       write(*,*) 'Facets at the wall in this rank: ', stack_size
       write(*,*) 'Facets at the wall global: ', stack_size_global
    end if

    ! Map from facet data to field data to have an easy way to perform operations with col3
    !! Area mass matrices
    call map_from_facets_to_field(mass_area_top, coef%area, top_wall_facet)
    call map_from_facets_to_field(mass_area_bot, coef%area, bot_wall_facet)
    call map_from_facets_to_field(mass_area_side, coef%area, side_wall_facet)
    !! Normal vectors in the required facets
    call map_from_facets_to_field(normal_x, coef%nx, wall_facet)
    call map_from_facets_to_field(normal_y, coef%ny, wall_facet)
    call map_from_facets_to_field(normal_z, coef%nz, wall_facet)
    
    !! Create a mask with the boundary facet data
    n = size(coef%B)
    nones = u%Xh%lx*u%Xh%ly*6*u%msh%nelv
    call rone(ones, nones)
    call map_from_facets_to_field(bdry_mask, ones, wall_facet)

    top_wall_area = glsum(mass_area_top%x,n)
    bot_wall_area = glsum(mass_area_bot%x,n)
    side_wall_area = glsum(mass_area_side%x,n)

    !! Verify validity of averaging routines
    bar_uzt = 0_rp
    top_wall_bar_dtdz = 0_rp
    bot_wall_bar_dtdz = 0_rp
    call average_from_weights(bar_uzt, bdry_mask, &
                        work_field, mass_area_side%x, side_wall_area)
    call average_from_weights(top_wall_bar_dtdz, bdry_mask, &
                        work_field, mass_area_top%x, top_wall_area)
    call average_from_weights(bot_wall_bar_dtdz, bdry_mask, &
                        work_field, mass_area_bot%x, bot_wall_area)
    !! Write it out
    if (pe_rank .eq. 0) then
       write(*,*) 'Area of top wall * normal= ', top_wall_area
       write(*,*) 'Area of bottom wall * normal= ', bot_wall_area
       write(*,*) 'Area of side wall= ', side_wall_area
       write(*,*) 'Validity check of averaging #1 (should be 1): ', bar_uzt 
       write(*,*) 'Validity check of averaging #2 (should be 1): ', top_wall_bar_dtdz 
       write(*,*) 'Validity check of averaging #3 (should be 1): ', bot_wall_bar_dtdz 
    end if 

    !> Store volume mass matrix for post processing
    n = size(coef%B) 
    call rzero(bm1%x,u%dof%size())
    call copy(bm1%x,coef%B,n)
    voll = glsum(bm1%x,n)
    if (pe_rank .eq. 0) then
       write(*,*) 'volume=', voll
    end if 
    !! rescale mass if volume is too small
    do i = 1,n
       bm1%x(i,1,1,1)=bm1%x(i,1,1,1) * 1e0
    end do
    voll = glsum(bm1%x,n)
    if (pe_rank .eq. 0) then
       write(*,*) 'rescaled volume=', voll
    end if 

    !> Perform IO
    call mf_area%write(area_l,t)
    call mf_bm1%write(bm1,t)      

    !> Initialize the calculation controller
    !!> Take data from case file
    call json_get(params, 'case.end_time', T_end)
    call json_get(params, 'case.monitor.output_control', &
                                     monitor_write_control)
    call json_get(params, 'case.monitor.calc_frequency', &
                                     monitor_write_par)    
    !!> Calculate relevant parameters and restart                     
    allocate(controllers(1))
    call controllers(1)%init(T_end, monitor_write_control, &
                             monitor_write_par)
    if (controllers(1)%nsteps .eq. 0) then
       controllers(1)%nexecutions = &
               int(t / controllers(1)%time_interval) + 1
    end if
    
    !> Initialize the spectral error indicator
    call spectral_error_indicator%init(u,v,w,coef)


  end subroutine user_initialize



  subroutine user_finalize(t, param)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: param

    ! Finalize variables related to nusselt calculation
    call work_field%free()
    call uzt%free()
    call dtdx%free()
    call dtdy%free()
    call dtdz%free()
    call dtdn%free()
    call div_dtdX%free()

    call dudx%free()
    call dudy%free()
    call dudz%free()
    call dvdx%free()
    call dvdy%free()
    call dvdz%free()
    call dwdx%free()
    call dwdy%free()
    call dwdz%free()

    call eps_k%free()
    call eps_t%free()
    
    call mass_area_top%free()
    call mass_area_bot%free()
    call mass_area_side%free()
    call bdry_mask%free()
    call bm1%free()
    call normal_x%free()
    call normal_y%free()
    call normal_z%free()
    
    ! Finilize list that contains uper and lower wall facets
    call wall_facet%free()
    call top_wall_facet%free()
    call bot_wall_facet%free()
    call side_wall_facet%free()
 
    ! Finilize the field lists
    call list_final3(area_l)
    call list_final3(dtdX_l)
    call list_final3(eps_l)

    ! deallocate the calc controllers
    if (allocated(controllers)) then
       deallocate(controllers)       
    end if
    
    !> Finalize the spectral error indicator
    call spectral_error_indicator%free()
  
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
    type(json_file), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(field_t), pointer :: s
    integer :: n,ntot, i,j,k, facet, e, lx,ly,lz
    integer :: index(4)
    type (tuple_i4_t) :: facet_el
    real(kind=rp) :: normal(3)
    real(kind=rp) :: voll(1), voll_temp(1)
    ! output control
    logical :: verify_bc = .false. ! write boundary conditions
    real(kind=rp) :: dt
    
    !> This value is used for breaking symtetries in bc
    istep = istep + 1

    !> Get the control parameters
    call json_get(params, 'case.timestep', dt)
    call json_get(params, 'case.monitor.verify_bc', verify_bc)

  if (controllers(1)%check(t, tstep, .false.)) then

    s => neko_field_registry%get_field('s')
    n = size(coef%B)
    ntot = coef%dof%size()
    lx = u%Xh%lx
    ly = u%Xh%ly
    lz = u%Xh%lz

    !> ------ for method #1 for nusselt calculation -----
    !> Nu_v = 1 + sqrt(Ra) * <u_z * T>_{v,t}
    !> Steps:
       !1.    Calculate the convective current T*u_z
       !2.    Get <u_z * T>_{v}(t) (Average in volume)
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_col3(uzt%x_d,w%x_d,s%x_d,n)               
    else
       call col3(uzt%x,w%x,s%x,n)                          
    end if
    !> ------ For method #2 for nusselt calculation -----
    !> Steps:
       !1.    Calculate the heat flux in z dtdz
       !2.    Get <dtdz>_{A}(t) (Average in area of plates)
    !!> Calculate derivatives. Automatically in device with opgrad
    call dudxyz (dtdx%x, s%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (dtdy%x, s%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (dtdz%x, s%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    
    call dudxyz (dudx%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (dudy%x, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (dudz%x, u%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    
    call dudxyz (dvdx%x, v%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (dvdy%x, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (dvdz%x, v%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    
    call dudxyz (dwdx%x, w%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (dwdy%x, w%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (dwdz%x, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    !> --------------------------------------------------

    !> Get the required averages minding the appropiate weights
    !> This function deals with device pointers internally
    bar_uzt = 0_rp
    top_wall_bar_dtdz = 0_rp
    bot_wall_bar_dtdz = 0_rp
    call average_from_weights(bar_uzt, uzt, &
                        work_field, coef%B, coef%volume)
    call average_from_weights(top_wall_bar_dtdz, dtdz, &
                        work_field, mass_area_top%x, top_wall_area)
    call average_from_weights(bot_wall_bar_dtdz, dtdz, &
                        work_field, mass_area_bot%x, bot_wall_area)
 
    !> write variables to monitor
    !! Integral quantities
    if (pe_rank .eq. 0) then
       open(10,file="nusselt.txt",position="append")
       write(10,*) t,'', bar_uzt, '', top_wall_bar_dtdz, '', &
                   bot_wall_bar_dtdz
       close(10)
    end if

    !> Write the derivatives
    call device_memcpy(dtdx%x,dtdx%x_d, n,DEVICE_TO_HOST)
    call device_memcpy(dtdy%x,dtdy%x_d, n,DEVICE_TO_HOST)
    call device_memcpy(dtdz%x,dtdz%x_d, n,DEVICE_TO_HOST)
    call mf_dtdX%write(dtdX_l,t)


    !> Calculate dissipations
    call calculate_thermal_dissipation(eps_t, dtdx,dtdy,dtdz, &
                                       work_field,Ra,Pr)
    call calculate_kinetic_dissipation(eps_k, dudx,dudy,dudz, &
                                       dvdx,dvdy,dvdz, &
                                       dwdx,dwdy,dwdz, &
                                       work_field,Ra,Pr)

    !> Write the derivatives
    call device_memcpy(eps_t%x,eps_t%x_d, n,DEVICE_TO_HOST)
    call device_memcpy(eps_k%x,eps_k%x_d, n,DEVICE_TO_HOST)
    call mf_eps%write(eps_l,t)



    ! Calculate some extra parameters to verify the boundary conditions
    if (verify_bc.eqv..true.) then
       
       !find the normal derivatives to verify boundaries
       call project_to_vector(dtdn, dtdx, dtdy, dtdz, &
                              normal_x%x, normal_y%x, normal_z%x)

       !>Calculate the divergence of the heat flux
       call divergence_of_field(div_dtdX, dtdx, dtdy, dtdz, &
                                work_field, coef)
       !!> volume average of the divergence
       call average_from_weights(bar_div_dtdX, div_dtdX, &
                        work_field, coef%B, coef%volume)

       !> perform IO               
       call device_memcpy(dtdn%x,dtdn%x_d, n,DEVICE_TO_HOST)
       call mf_dtdn%write(dtdn,t)

       !> get the heat flux averages this way. It should match
       call average_from_weights(heat_bot_wall, dtdn, &
                        work_field, mass_area_bot%x, bot_wall_area)
       call average_from_weights(heat_top_wall, dtdn, &
                        work_field, mass_area_top%x, top_wall_area)

       call average_from_weights(heat_side_wall, dtdn, &
                        work_field, mass_area_side%x, side_wall_area)

       heat_balance = heat_bot_wall + heat_top_wall + heat_side_wall

       if (pe_rank .eq. 0) then
          open(20,file="bc_heat_balance.txt",position="append")
          write(20,*) t,'', heat_top_wall, '', heat_bot_wall, '', &
                      heat_side_wall, '', heat_balance, '', &
                      bar_div_dtdX
          close(20)
       end if

    end if

    !> Get the spectral error indicators for the mesh
    call spectral_error_indicator%get_indicators(coef)
    call spectral_error_indicator%write_as_field(t)

  call controllers(1)%register_execution()
  end if
  end subroutine calculate_nusselt


  !==================================================
  !> Supporting subroutines
  !==================================================


  subroutine list_init3(list,uu,vv,ww)
    type(field_list_t), intent(inout) :: list
    type(field_t) , target:: uu
    type(field_t) , target:: vv
    type(field_t) , target:: ww
    !> Initialize field lists
    allocate(list%fields(3))
    list%fields(1)%f => uu
    list%fields(2)%f => vv
    list%fields(3)%f => ww
  end subroutine list_init3


  subroutine list_final3(list)
    type(field_list_t), intent(inout) :: list
    !> Deallocate field lists
    deallocate(list%fields)
  end subroutine list_final3

  
  subroutine map_from_facets_to_field(u, area, wall_facet)
    !> this is needed because the facet data is arrayed differently than filed data
    !> u(lx,ly,lz,e) is the field where you want to put the facet data
    !> area(lx,ly,facet,e) is the array in facet form
    !> wall_facet is the list of facets that should be mapped
    type(field_t), intent(inout) :: u
    real(kind=rp) :: area(u%Xh%lx,u%Xh%ly,6,u%msh%nelv)
    type(stack_i4t4_t) :: wall_facet

    !> Parameters to populate the list of elements and facets
    integer :: e, facet, e_counter, facet_counter, i,j,k,n
    type(tuple4_i4_t), pointer :: wall_facet_list(:)
    integer :: lx, ly, lz

    wall_facet_list => wall_facet%array()
    !! Initialize as 0. This serves as a mask for nodes that are not integrated
    call rzero(u%x,u%dof%size())
    !! Fill up the top and bottom mass matrices
    lx = u%Xh%lx
    ly = u%Xh%ly
    lz = u%Xh%lz
    do e_counter = 1, wall_facet%top_
      e = wall_facet_list(e_counter)%x(2)
      facet = wall_facet_list(e_counter)%x(1)
      select case(facet)
         case(1)
            do j = 1 , ly
               do k = 1 , lz
                  u%x(1,j,k,e) = area(j,k,facet,e) 
               end do
            end do
         case(2)
            do j = 1 , ly
               do k = 1 , lz
                  u%x(lx,j,k,e) = area(j,k,facet,e) 
               end do
            end do
         case(3)
            do i = 1 , lx
               do k = 1 , lz
                  u%x(i,1,k,e) = area(i,k,facet,e)    
               end do
            end do
         case(4)
            do i = 1 , lx
               do k = 1 , lz
                  u%x(i,ly,k,e) = area(i,k,facet,e)       
               end do
            end do
         case(5)
            do i = 1 , lx
               do j = 1 , ly
                  u%x(i,j,1,e) = area(i,j,facet,e)
               end do
            end do
         case(6)
            do i = 1 , lx
               do j = 1 , ly
                  u%x(i,j,lz,e) = area(i,j,facet,e)
               end do
            end do
      end select
    end do


    ! Put it also on device
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_memcpy(u%x,u%x_d, &
                          u%dof%size(),HOST_TO_DEVICE)
    end if

  end subroutine map_from_facets_to_field


  subroutine average_from_weights(avrg, field, work_field, weights, sum_weights)
    real(kind=rp), intent(inout) :: avrg
    type(field_t), intent(inout) :: field
    type(field_t), intent(inout) :: work_field
    real(kind=rp), dimension(field%Xh%lxyz,field%msh%nelv), intent(inout) :: weights
    real(kind=rp), intent(inout) :: sum_weights
    integer :: n
    type(c_ptr) :: weights_d
    
    n = field%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then 
       weights_d = device_get_ptr(weights)
       call device_col3(work_field%x_d,field%x_d,weights_d,n)      
       avrg = device_glsum(work_field%x_d,n)                  
       avrg = avrg / abs(sum_weights)            
    else
       call col3(work_field%x,field%x,weights,n)      
       avrg = glsum(work_field%x,n)                  
       avrg = avrg / abs(sum_weights)            
    end if

  end subroutine average_from_weights


  subroutine project_to_vector(proj, uu, vv, ww, nx, ny, nz)
    type(field_t), intent(inout) :: proj
    type(field_t), intent(inout) :: uu
    type(field_t), intent(inout) :: vv
    type(field_t), intent(inout) :: ww
    real(kind=rp), dimension(proj%Xh%lxyz,proj%msh%nelv), intent(inout) :: nx
    real(kind=rp), dimension(proj%Xh%lxyz,proj%msh%nelv), intent(inout) :: ny
    real(kind=rp), dimension(proj%Xh%lxyz,proj%msh%nelv), intent(inout) :: nz
    integer :: n
    type(c_ptr) :: nx_d
    type(c_ptr) :: ny_d
    type(c_ptr) :: nz_d
    
    n = proj%dof%size()

    ! Find the normal temperature derivative at the boundaries
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       nx_d = device_get_ptr(nx)
       ny_d = device_get_ptr(ny)
       nz_d = device_get_ptr(nz)
       call device_rzero(proj%x_d,n)
       call device_vdot3(proj%x_d, uu%x_d,vv%x_d,ww%x_d, &
                         nx_d,ny_d,nz_d, n) 
    else
       call rzero(proj%x,n)
       call vdot3(proj%x, uu%x,vv%x,ww%x, &
                         nx,ny,nz, n) 
    end if

  end subroutine project_to_vector


  subroutine divergence_of_field(div_u, u, v, w, work_field, coef)
    type(field_t), intent(inout) :: div_u
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: work_field
    type(coef_t), intent(inout) :: coef
    integer :: n
    
    n = div_u%dof%size()
    
    !dudx in div_u
    call dudxyz (div_u%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    !dvdy in work_field
    call dudxyz (work_field%x, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    ! Add dudx + dvdy    
    if (NEKO_BCKND_DEVICE .eq. 1) then 
          call device_add2(div_u%x_d,work_field%x_d,n)
    else
          call add2(div_u%x,work_field%x,n)
    end if
    ! dwdz in work_dielf
    call dudxyz (work_field%x, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    ! Add (dudx + dvdy) + dwdz 
    if (NEKO_BCKND_DEVICE .eq. 1) then 
          call device_add2(div_u%x_d,work_field%x_d,n)
    else
          call add2(div_u%x,work_field%x,n)
    end if
     
  end subroutine divergence_of_field


  subroutine calculate_thermal_dissipation(eps_t, dtdx,dtdy,dtdz, work_field,Ra,Pr)
    type(field_t), intent(inout) :: eps_t
    type(field_t), intent(inout) :: dtdx
    type(field_t), intent(inout) :: dtdy
    type(field_t), intent(inout) :: dtdz
    type(field_t), intent(inout) :: work_field
    real(kind=rp), intent(in) :: Ra
    real(kind=rp), intent(in) :: Pr
     
    real(kind=rp) :: sqrtRaPr_inv
    integer :: n

    n = eps_t%dof%size()
    sqrtRaPr_inv = 1_rp/sqrt(Ra*Pr)

    !Thermal dissipation epsT=(grad T)**2/sqrt(RaPr)
    !jorg did it as epsT=((dtdx)**2+(dtdy)**2+(dtdz)**2)/sqrt(RaPr)
    if (NEKO_BCKND_DEVICE .eq. 1) then 
        call device_rzero(eps_t%x_d,n)
        call device_addsqr2s2(eps_t%x_d, dtdx%x_d, sqrtRaPr_inv, n)
        call device_addsqr2s2(eps_t%x_d, dtdy%x_d, sqrtRaPr_inv, n)
        call device_addsqr2s2(eps_t%x_d, dtdz%x_d, sqrtRaPr_inv, n)
    else
        call rzero(eps_t%x,n)
        call addsqr2s2(eps_t%x, dtdx%x, sqrtRaPr_inv, n)
        call addsqr2s2(eps_t%x, dtdy%x, sqrtRaPr_inv, n)
        call addsqr2s2(eps_t%x, dtdz%x, sqrtRaPr_inv, n)
  
    end if

  end subroutine calculate_thermal_dissipation
  
  
  
  subroutine calculate_kinetic_dissipation(eps_k, dudx,dudy,dudz, &
                                                  dvdx,dvdy,dvdz, &
                                                  dwdx,dwdy,dwdz, &
                                                  work_field,Ra,Pr)
    type(field_t), intent(inout) :: eps_k
    type(field_t), intent(inout) :: dudx
    type(field_t), intent(inout) :: dudy
    type(field_t), intent(inout) :: dudz
    type(field_t), intent(inout) :: dvdx
    type(field_t), intent(inout) :: dvdy
    type(field_t), intent(inout) :: dvdz
    type(field_t), intent(inout) :: dwdx
    type(field_t), intent(inout) :: dwdy
    type(field_t), intent(inout) :: dwdz
    type(field_t), intent(inout) :: work_field
    real(kind=rp), intent(in) :: Ra
    real(kind=rp), intent(in) :: Pr
     
    real(kind=rp) :: sqrtPr_Ra, c1 = 1, c2 = 1
    integer :: n

    n = eps_t%dof%size()
    sqrtPr_Ra = sqrt(Pr/Ra)*0.5_rp

    !Energy dissipation epsv=0.5*(du_i/dx_j+du_j/dx_i)**2*sqrt(Pr/Ra)
    if (NEKO_BCKND_DEVICE .eq. 1) then 
        call device_rzero(eps_k%x_d,n)
        
        call device_add3s2(work_field%x_d,dudx%x_d,dudx%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, sqrtPr_Ra, n)
        call device_add3s2(work_field%x_d,dudy%x_d,dvdx%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, sqrtPr_Ra, n)
        call device_add3s2(work_field%x_d,dudz%x_d,dwdx%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, sqrtPr_Ra, n)

        call device_add3s2(work_field%x_d,dvdx%x_d,dudy%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, sqrtPr_Ra, n)
        call device_add3s2(work_field%x_d,dvdy%x_d,dvdy%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, sqrtPr_Ra, n)
        call device_add3s2(work_field%x_d,dvdz%x_d,dwdy%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, sqrtPr_Ra, n)

        call device_add3s2(work_field%x_d,dwdx%x_d,dudz%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, sqrtPr_Ra, n)
        call device_add3s2(work_field%x_d,dwdy%x_d,dvdz%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, sqrtPr_Ra, n)
        call device_add3s2(work_field%x_d,dwdz%x_d,dwdz%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, sqrtPr_Ra, n)
    else
        call rzero(eps_k%x,n)
        
        call add3s2(work_field%x,dudx%x,dudx%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, sqrtPr_Ra, n)
        call add3s2(work_field%x,dudy%x,dvdx%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, sqrtPr_Ra, n)
        call add3s2(work_field%x,dudz%x,dwdx%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, sqrtPr_Ra, n)

        call add3s2(work_field%x,dvdx%x,dudy%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, sqrtPr_Ra, n)
        call add3s2(work_field%x,dvdy%x,dvdy%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, sqrtPr_Ra, n)
        call add3s2(work_field%x,dvdz%x,dwdy%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, sqrtPr_Ra, n)

        call add3s2(work_field%x,dwdx%x,dudz%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, sqrtPr_Ra, n)
        call add3s2(work_field%x,dwdy%x,dvdz%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, sqrtPr_Ra, n)
        call add3s2(work_field%x,dwdz%x,dwdz%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, sqrtPr_Ra, n)

    end if

  end subroutine calculate_kinetic_dissipation



end module user
