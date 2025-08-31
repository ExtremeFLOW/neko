! Two-dimensional lid-driven cavity
!
! Time-integration of the lid-driven cavity with smoothened
! belt velocity to fulfil the continuity equation in the corners.
!
! Note that the domain is actually 3D with width one element. In order
! to prevent any instability in the z direction, the w velocity is
! set to zero at every step. This is needed for higher Reynolds numbers.
!
module user
  use neko
  implicit none

  ! Global user variables
  type(field_t) :: w1
  type(field_t) :: temp1,temp2
  type(field_t) :: vort1,vort2,vort3

  type(file_t) output_file ! output file
  type(vector_t) :: vec_out ! will store our output data
  integer :: ipostproc ! frequency of the output

  real(kind=rp) :: eps
  real(kind=rp) :: gamma, u_max

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%dirichlet_conditions => dirichlet_conditions
    user%compute => compute
    user%initialize =>initialize
    user%finalize => finalize
    user%source_term => source_term
    user%material_properties => material_properties
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    character(len=50) :: mess

    ! read postprocessing interval
    call json_get(params, "case.fluid.ipostproc", ipostproc)
    write(mess,*) "postprocessing steps : ",ipostproc
    call neko_log%message(mess)

    call json_get(params, "case.scalar.epsilon",eps)
    call json_get(params, "case.scalar.gamma", gamma)
    u_max = 1.0_rp
  end subroutine startup

  ! user-defined Dirichlet boundary condition
  subroutine dirichlet_conditions(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: u, v, w
    integer :: i
    real(kind=rp) :: x

    real(kind=rp) lsmoothing

    if (time%tstep .ne. 1) return

    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")

    lsmoothing = 0.1_rp ! length scale of smoothing at the edges

    do i = 1, bc%msk(0)
       x = u%dof%x(bc%msk(i), 1, 1, 1)
       u%x(bc%msk(i), 1, 1, 1) = &
            step( x/lsmoothing ) * step( (1._rp - x)/lsmoothing )

       v%x(bc%msk(i), 1, 1, 1) = 0
       w%x(bc%msk(i), 1, 1, 1) = 0
    end do

  end subroutine dirichlet_conditions

  ! User-defined initialization called just before time loop starts
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: u

    ! Initialize the file object and create the output.csv file
    ! in the working directory (overwrite any existing file)
    call output_file%init("ekin.csv", overwrite=.true.)
    call output_file%set_header("# t,Ekin,enst")
    call vec_out%init(2) ! Initialize our vector with 2 elements (Ekin, enst)

    ! initialize work arrays for postprocessing
    u => neko_field_registry%get_field("u")
    call w1%init(u%dof, 'work1')
    call temp1%init(u%dof)
    call temp2%init(u%dof)
    call vort1%init(u%dof)
    call vort2%init(u%dof)
    call vort3%init(u%dof)


    ! call usercheck also for tstep=0
    call compute(time)

  end subroutine initialize

  ! User-defined routine called at the end of every time step
  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    integer :: ntot, i
    real(kind=rp) :: ekin, enst
    type(field_t), pointer :: u, v, w
    type(coef_t), pointer :: coef


    ! do the postprocessing only every 50 time steps (including step 0)
    if (mod(time%tstep, ipostproc) .ne. 0) return

    coef => neko_user_access%case%fluid%c_Xh
    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")

    ntot = u%dof%size()

    ! compute the velocity magnitude in w1 (u^2 + v^2 + w^2)
    call field_col3(w1, u, u, ntot)
    call field_addcol3(w1, v, v, ntot)
    call field_addcol3(w1, w, w, ntot)
    
    ! compute maximum velocity magnitude over the domain
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(w1%x, w1%x_d, w1%size(), &
            DEVICE_TO_HOST, sync=.true.)
    end if
    u_max = sqrt(glmax(w1%x, ntot))

    ! compute the kinetic energy (reuse w1 which already has velocity magnitude squared)
    ! For glsc2 we need to call the correct backend
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ekin = 0.5_rp * device_glsc2(w1%x_d, coef%B_d,ntot) / coef%volume
    else
       ekin = 0.5_rp * glsc2(w1%x, coef%B, ntot) / coef%volume
    end if

    ! compute enstrophy
    ! (the factor of 0.5 depends on the definition of enstrophy. We
    ! follow the reference paper by the HiOCFD4 workshop, but the book
    ! by Pope for instance would not include this factor)
    call curl(vort1,vort2,vort3, u, v, w, temp1, temp2, coef)
    call field_col3(w1, vort1, vort1, ntot)
    call field_addcol3(w1, vort2, vort2, ntot)
    call field_addcol3(w1, vort3, vort3, ntot)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       enst = 0.5_rp * device_glsc2(w1%x_d, coef%B_d, ntot) / coef%volume
    else
       enst = 0.5_rp * glsc2(w1%x, coef%B, ntot) / coef%volume
    end if


    ! output all this to file
    call neko_log%message("Writing csv file")
    vec_out%x = (/ekin, enst/)
    call output_file%write(vec_out, time%t)

    ! set the w component to zero to avoid any 3D instability
    ! in this quasi-2D flow
    call field_rzero(w)

  end subroutine compute

  ! User-defined finalization routine called at the end of the simulation
  subroutine finalize(time)
    type(time_state_t), intent(in) :: time

    ! Deallocate the fields
    call w1%free()

    ! Deallocate output file and vector
    call file_free(output_file)
    call vec_out%free

  end subroutine finalize

  ! User-defined source term (forcing function)
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

    ! Only apply forcing to the scalar equation (phase)
    if (scheme_name .ne. 'phase') return

    ! Get the right-hand side field for the scalar
    rhs_s => rhs%get_by_index(1)  ! scalar field
    
    ! Get the scalar field and coefficients
    s => neko_field_registry%get_field('phase')
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
       absgrad = sqrt(work1%x(i,1,1,1)**2+work2%x(i,1,1,1)**2+work3%x(i,1,1,1)**2)
       if (absgrad == 0.0_rp) then 
          print *, 'warning, absgrad==', absgrad
          absgrad = 1e21_rp
       end if
          
       work1%x(i,1,1,1) = - s%x(i,1,1,1)*(1.0_rp-s%x(i,1,1,1))*(work1%x(i,1,1,1)/absgrad)
       work2%x(i,1,1,1) = - s%x(i,1,1,1)*(1.0_rp-s%x(i,1,1,1))*(work2%x(i,1,1,1)/absgrad)
       work3%x(i,1,1,1) = - s%x(i,1,1,1)*(1.0_rp-s%x(i,1,1,1))*(work3%x(i,1,1,1)/absgrad)
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

  ! User-defined material properties
  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time
    real(kind=rp) :: delta, lambda_val, mu_val

    if (scheme_name .eq. "fluid") then
      call field_cfill(properties%get("fluid_rho"), 1.0_rp)
      mu_val = 1.0_rp / 3000.0_rp  ! 1/Re
      call field_cfill(properties%get("fluid_mu"), mu_val)
    else if (scheme_name .eq. "phase") then
      delta = u_max * gamma
      call field_cfill(properties%get('phase_cp'), 1.0_rp)
      call field_cfill(properties%get('phase_lambda'), eps*delta)
    end if

  end subroutine material_properties

  ! User-defined initial conditions for phase field
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    
    type(field_t), pointer :: phase
    integer :: i, j, k, e
    real(kind=rp) :: x, y, d_x
    real(kind=rp) :: x_interface
    
    ! Only apply to phase field
    if (scheme_name .ne. 'phase') return
    
    phase => fields%items(1)%ptr
    
    ! Interface location (horizontal interface at y = 0.5)
    x_interface = 0.5_rp
    
    ! Initialize phase field using signed distance function
    ! φ(x,0) = 1/2 [1 - tanh(d(x)/(2ε))]
    ! where d(x) is the signed distance to the interface
    do e = 1, phase%msh%nelv
      do k = 1, phase%Xh%lx
        do j = 1, phase%Xh%lx
          do i = 1, phase%Xh%lx
            x = phase%dof%x(i,j,k,e)
            y = phase%dof%y(i,j,k,e)
            
            ! Signed distance to horizontal interface at y = x_interface
            ! Positive above interface (fluid 1), negative below (fluid 0)
            d_x = y - x_interface
            
            ! Phase field initial condition: φ(x,0) = 1/2 [1 - tanh(d(x)/(2ε))]
            ! ε is the nominal half-thickness of the diffuse interface
            phase%x(i,j,k,e) = 0.5_rp * (1.0_rp - tanh(d_x / (2.0_rp * eps)))
            
          end do
        end do
      end do
    end do
    
  end subroutine initial_conditions

  ! Smooth step function, with zero derivatives at 0 and 1
  function step(x)
    ! Taken from Simson code
    ! x<=0 : step(x) = 0
    ! x>=1 : step(x) = 1
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: step

    if (x.le.0.02_rp) then
       step = 0.0_rp
    else
       if (x.le.0.98_rp) then
          step = 1._rp/( 1._rp + exp(1._rp/(x-1._rp) + 1._rp/x) )
       else
          step = 1._rp
       end if
    end if

  end function step

end module user
