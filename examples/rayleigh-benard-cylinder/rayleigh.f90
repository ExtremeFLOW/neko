module user
  use neko
  implicit none

  !> Variables to store the Rayleigh and Prandlt numbers
  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Pr = 0

  !> Arrays asociated with Method#1 for nusselt calculation
  integer :: calc_frequency = 500 ! How frequently should we calculate Nu
  type(field_t) :: work_field ! Field to perform operations
  type(field_t) :: uzt ! u_z * T
  type(field_t) :: dtdx ! Derivative of scalar wrt x
  type(field_t) :: dtdy ! Detivative of scalar wrt y
  type(field_t) :: dtdz ! Derivative of scalar wrt z
  type(field_t) :: mass_area_top ! mass matrix for area on top wall
  type(field_t) :: mass_area_bot ! mass matrix for area on bottom wall
  real(kind=rp) :: bar_uzt ! Volume integral
  real(kind=rp) :: top_wall_bar_dtdz ! area integral
  real(kind=rp) :: bot_wall_bar_dtdz ! area integral
  real(kind=rp) :: top_wall_area ! area of top and bottom wall
  real(kind=rp) :: bot_wall_area ! area of top and bottom wall


  !> List of elements and facets in uper and lower boundary
  type(stack_i4t4_t) :: wall_facet
  type(tuple4_i4_t)  :: facet_el_type_0


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

    !Debugging info
    !if (pe_rank .eq. 0) then
    !  write(*,*) 'Element:', ie, 'in boundary'
    !end if 
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
    type(param_t), intent(inout) :: params

    !> Parameters to populate the list of elements and facets
    real(kind=rp) :: stack_size(1) !Do it like this so it is cast
    real(kind=rp) :: stack_size_global
    integer :: e, facet, typ, e_counter, facet_counter, i,j,k,n
    real(kind=rp) :: normal(3)
    real(kind=rp) :: area_mass, area_rank(1),area_global(1)
    type(tuple4_i4_t), pointer :: wall_facet_list(:)
    
    !> Reset the relevant nondimensional parameters
    Pr = params%Pr
    Ra = params%Re
    params%Re = sqrt(Ra / Pr)

    !> Initialize variables related to nusselt calculation
    call field_init(work_field, u%dof, 'work_field')
    call field_init(uzt, u%dof, 'uzt')
    call field_init(dtdx, u%dof, 'dtdx')
    call field_init(dtdy, u%dof, 'dtdy')
    call field_init(dtdz, u%dof, 'dtdz')
    call field_init(mass_area_top, u%dof, 'mat')
    call field_init(mass_area_bot, u%dof, 'mab')

    !> Initialize list with upper and lower wall facets
    call wall_facet%init()

    !> Populate the list with upper and lower wall facets 
    do e = 1, u%msh%nelv !Go over all elements
      do facet = 1, 6 ! Go over all facets of hex element
        normal = coef_get_normal(coef,1,1,1,e,facet) ! Get facet normal
        typ =  u%msh%facet_type(facet,e)             ! Get facet type
        if (typ.ne.0.and.abs(normal(3)).ge.1-1e-5) then
          write(*,*) 'In boundary: facet=', facet, 'and element=', e
          write(*,*) 'BC facet type ', typ 
          write(*,*) 'normal to this facet is: ', normal   
          facet_el_type_0%x = (/facet, e, typ, 0/)
          call wall_facet%push(facet_el_type_0)
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

    !> Fill up arrays to serve as mass matrix for integration
    !! we want <dt/dz>_A,t at z=0,1.
    !! Average in Area is: (dt/dz * normal * Area_mass) / Area_total
    wall_facet_list => wall_facet%array()
    !! Initialize the mass as 0. This serves as a mask for nodes that are not integrated
    call rzero(mass_area_top%x,u%dof%size())
    call rzero(mass_area_bot%x,u%dof%size())
    !! Fill up the top and bottom mass matrices
    do e_counter=1,wall_facet%top_
      do i = 1, u%Xh%lx
        do j = 1 , u%Xh%lx
           facet = wall_facet_list(e_counter)%x(1)
           e = wall_facet_list(e_counter)%x(2)
           normal = coef_get_normal(coef,1,1,1,e,facet) ! Get facet normal
           area_mass = get_area_mass(coef,i,j,1,e,facet) ! We are interested only in facet 5 and 6 for our case
           select case(facet) !Based on function: index_is_on_facet
              case(5)
                 mass_area_bot%x(i,j,1,e) = area_mass * normal(3)
              case(6)
                 mass_area_top%x(i,j,u%Xh%lz,e) = area_mass * normal(3)
           end select                
         end do
      end do 
    end do
    n = size(coef%B)
    top_wall_area = glsum(mass_area_top%x,n)
    bot_wall_area = glsum(mass_area_bot%x,n)
    !! Write it out
    if (pe_rank .eq. 0) then
       write(*,*) 'Area of top wall * normal= ', top_wall_area
       write(*,*) 'Area of bottom wall * normal= ', bot_wall_area
    end if 

    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_memcpy(mass_area_top%x,mass_area_top%x_d, &
                          mass_area_top%dof%size(),HOST_TO_DEVICE)
       call device_memcpy(mass_area_bot%x,mass_area_bot%x_d, &
                          mass_area_bot%dof%size(),HOST_TO_DEVICE)
    end if

  end subroutine user_initialize


  pure function get_area_mass(coef, i, j, k, e, facet) result(area_mass)
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: i, j, k, e, facet
    real(kind=rp) :: area_mass
      
    select case (facet)               
      case(1,2)
        area_mass = coef%area(j, k, facet, e)
      case(3,4)
        area_mass = coef%area(i, k, facet, e)
      case(5,6)
        area_mass = coef%area(i, j, facet, e)
      end select
  end function get_area_mass


  subroutine user_finalize(t, param)
    real(kind=rp) :: t
    type(param_t), intent(inout) :: param

    ! Finalize variables related to nusselt calculation
    call field_free(work_field)
    call field_free(uzt)
    call field_free(dtdx)
    call field_free(dtdy)
    call field_free(dtdz)
    call field_free(mass_area_top)
    call field_free(mass_area_bot)
    
    ! Finilize list that contains uper and lower wall facets
    call wall_facet%free()

  end subroutine user_finalize

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
    integer :: n,ntot, i,j,k, facet, e
    integer :: index(4)
    type (tuple_i4_t) :: facet_el
    real(kind=rp) :: normal(3)

    if (mod(tstep,calc_frequency).ne.0) return

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



    !> ------ Method #2 for nusselt calculation -----
    ! Calculate derivatives. Automatically in device with opgrad
    call opgrad(dtdx%x,dtdy%x,dtdz%x,s%x,coef) 

    if (NEKO_BCKND_DEVICE .eq. 1) then 
       ! Calculate for top wall
       call device_col3(work_field%x_d,dtdz%x_d,mass_area_top%x_d,n)      
       top_wall_bar_dtdz = device_glsum(work_field%x_d,n)                  
       top_wall_bar_dtdz = top_wall_bar_dtdz / abs(top_wall_area)            
       ! Calculate for bot wall
       call device_col3(work_field%x_d,dtdz%x_d,mass_area_bot%x_d,n)      
       bot_wall_bar_dtdz = device_glsum(work_field%x_d,n)                  
       bot_wall_bar_dtdz = bot_wall_bar_dtdz / abs(bot_wall_area)            
    else
       ! Calculate for top wall
       call col3(work_field%x,dtdz%x,mass_area_top%x,n)      
       top_wall_bar_dtdz = glsum(work_field%x,n)                  
       top_wall_bar_dtdz = top_wall_bar_dtdz / abs(top_wall_area)            
       ! Calculate for bot wall
       call col3(work_field%x,dtdz%x,mass_area_bot%x,n)      
       bot_wall_bar_dtdz = glsum(work_field%x,n)                  
       bot_wall_bar_dtdz = bot_wall_bar_dtdz / abs(bot_wall_area)            
    end if
    
    !> Include variables to monitor
    if (pe_rank .eq. 0) then
       open(10,file="nusselt.txt",position="append")
       write(10,*) t,'', bar_uzt, '', top_wall_bar_dtdz, '', &
                   bot_wall_bar_dtdz
       close(10)
    end if
  

  end subroutine calculate_nusselt

end module user
