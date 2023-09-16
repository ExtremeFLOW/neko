module user
  use neko
  implicit none

  !> Variables to store the Rayleigh and Prandlt numbers
  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Re = 0
  real(kind=rp) :: Pr = 0


  !> ========== Needed for Probes =================
  
  !> Probe type
  type(probes_t) :: pb

  !> Output variables
  type(file_t) :: fout
  type(matrix_t) :: mat_out

  !> Case IO parameters  
  integer            :: n_fields
  character(len=:), allocatable  :: output_file

  !> Output control
  logical :: write_output = .false.

  !> =============================================

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules => user_initialize
    u%user_finalize_modules => user_finalize
    u%fluid_user_ic => set_initial_conditions_for_u_and_s
    u%scalar_user_bc => set_scalar_boundary_conditions
    u%fluid_user_f_vector => set_bousinesq_forcing_term
    u%user_check => check
  end subroutine user_setup
 
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
    
    ! This will be used on all zones without labels
    ! e.g. the ones hardcoded to 'v', 'w', etcetc
    s = 1.0_rp - z

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

    !> Log variable
    character(len=LOG_SIZE) :: log_buf ! For logging status

    !> Support variables for probes 
    integer :: i
    type(matrix_t) :: mat_coords

    !> Recalculate the non dimensional parameters
    call json_get(params, 'case.scalar.Pr', Pr)
    call json_get(params, 'case.fluid.Re', Re)
    Ra = (Re**2)*Pr
    write(log_buf,*) 'Rayleigh Number is Ra=', Ra
    call neko_log%message(log_buf)
    
!    if (pe_rank.eq.0) write(*,*) 

    !> ========== Needed for Probes =================
    
    !> Read the output information
    call json_get(params, 'case.probes.output_file', output_file) 

    !> Probe set up
    !! Read probe info and initialize the controller, arrays, etc.
    call pb%init(t, params)
    !! Perform the set up of gslib_findpts
    call pb%setup(coef)
    !! Find the probes in the mesh. Map from xyz -> rst
    call pb%map(coef)
    !> Write a summary of the probe info
    call pb%show()

    !> Initialize the output
    fout = file_t(trim(output_file))
    call mat_out%init(pb%n_probes, pb%n_fields)

    !> Write coordinates in output file (imitate nek5000)
    !! Initialize the arrays
    call mat_coords%init(pb%n_probes,3)
    !! Array them as rows
    call transpose(mat_coords%x, pb%n_probes, pb%xyz, 3)
    !! Write the data to the file
    call fout%write(mat_coords)
    !! Free the memory 
    call mat_coords%free

    !> ==============================================

  end subroutine user_initialize

  subroutine user_finalize(t, param)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: param

    !> ========== Needed for Probes =================
    
    call pb%free
    call mat_out%free
    call file_free(fout)
    
    !> ==============================================
     
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

    !> ========== Needed for Probes =================

    !> Interpolate the desired fields
    call pb%interpolate(t,tstep, write_output)
    !! Write if the interpolate function returs write_output=.true.
    if (write_output) then
       mat_out%x = pb%out_fields
       call fout%write(mat_out, t)
       write_output = .false.
    end if

    !> ==============================================

  end subroutine check

end module user
