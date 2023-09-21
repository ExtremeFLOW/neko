module data_streamer_2

  use num_types, only: rp
  use field, only: field_t
  use coefs, only: coef_t
  use device 
  use comm
  use mpi_types
  use neko_config

  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  type, public :: data_streamer_2_t
       
     !> Define if the execution is asyncrhonous
     integer :: if_asynch

     !> global element numbers
     integer, allocatable :: lglel(:)

     contains
       !> Constructor
       procedure, pass(this) :: init => data_streamer_2_init
       !> Destructor
       procedure, pass(this) :: free => data_streamer_2_free
       !> Stream data
       procedure, pass(this) :: stream => data_streamer_2_stream

  end type data_streamer_2_t

contains

  !> Constructor
  subroutine data_streamer_2_init(this, coef, if_asynch)
    class(data_streamer_2_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: if_asynch
    integer :: nelb, nelb2, nelv, nelgv,npts,e

    !Allocate and initialize the global element number
    allocate(this%lglel(coef%msh%nelv))
    do e = 1, coef%msh%nelv
       this%lglel(e) = e + coef%msh%offset_el
    end do

    !Assign if the streaming is asynchronous
    this%if_asynch = if_asynch

    !Assign the set up parameters
    nelv  = coef%msh%nelv
    npts  = coef%Xh%lx*coef%Xh%ly*coef%Xh%lz
    nelgv = coef%msh%glb_nelv
    nelb  = coef%msh%offset_el
    ! Alternative way to get nelb:
    !nelb = elem_running_sum(nelv)
    !nelb = nelb - nelv

!#ifdef HAVE_ADIOS2
    call fortran_adios2_setup(npts,nelv,nelb,nelgv, &
                nelgv,coef%dof%x,coef%dof%y,  &
                coef%dof%z,if_asynch,NEKO_COMM)
!#else
!    call neko_error('NEKO needs to be built with ADIOS2 support')
!#endif


  end subroutine data_streamer_2_init

  !> Destructor
  subroutine data_streamer_2_free(this)
    class(data_streamer_2_t), intent(inout) :: this

    if (allocated(this%lglel))        deallocate(this%lglel)

!#ifdef HAVE_ADIOS2
    call fortran_adios2_finalize()
!#else
!    call neko_error('NEKO needs to be built with ADIOS2 support')
!#endif
    
  end subroutine data_streamer_2_free
  
  !> streamer
  subroutine data_streamer_2_stream(this,u,v,w,p,coef)
    class(data_streamer_2_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    integer :: nelv, npts
    
    nelv  = coef%msh%nelv
    npts  = coef%Xh%lx*coef%Xh%ly*coef%Xh%lz
    
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       ! Move the data to the CPU to be able to write it
       call device_memcpy(u%x,  u%x_d, nelv*npts,DEVICE_TO_HOST)
       call device_memcpy(v%x,  v%x_d, nelv*npts,DEVICE_TO_HOST)
       call device_memcpy(w%x,  w%x_d, nelv*npts,DEVICE_TO_HOST)
       call device_memcpy(p%x,  p%x_d, nelv*npts,DEVICE_TO_HOST)
    end if

!#ifdef HAVE_ADIOS2
    call fortran_adios2_stream(this%lglel,p%x, u%x, v%x, w%x, coef%B, u%x)
!#else
!    call neko_error('NEKO needs to be built with ADIOS2 support')
!#endif

  end subroutine data_streamer_2_stream
  
  !> Supporting function to calculate the element number offset.  
  function elem_running_sum(nelv) result(rbuff)

    integer, intent(in) :: nelv
    integer ::  ierr,xbuff,wbuff,rbuff

    xbuff = nelv  ! running sum
    wbuff = nelv  ! working buff
    rbuff = 0   ! recv buff
     
    call mpi_scan(xbuff,rbuff,1,mpi_integer,mpi_sum,NEKO_COMM,ierr)

  end function elem_running_sum

  !> Interfaces to C++ code:
  subroutine fortran_adios2_setup(npts,nelv,nelb,nelgv, nelgt,x,y,z,asynch,comm)  
     use, intrinsic :: ISO_C_BINDING
     implicit none
     real(kind=rp), dimension(:,:,:,:), intent(inout) :: x
     real(kind=rp), dimension(:,:,:,:), intent(inout) :: y
     real(kind=rp), dimension(:,:,:,:), intent(inout) :: z
     integer, intent(in) :: npts, nelv, nelb,nelgv, nelgt, asynch
     type(MPI_COMM) :: comm
     interface
        ! C-definition is: int func(double *data, const int nx, const int ny)
        subroutine c_adios2_setup(npts,nelv,nelb,nelgv,nelgt,x,y,z,asynch,comm) bind(C,name="adios2_setup_")
           use, intrinsic :: ISO_C_BINDING
           implicit none
           integer(kind=C_INT) :: npts
           integer(kind=C_INT) :: nelv
           integer(kind=C_INT) :: nelb
           integer(kind=C_INT) :: nelgv
           integer(kind=C_INT) :: nelgt
           real(kind=C_DOUBLE), intent(INOUT) :: x(*)
           real(kind=C_DOUBLE), intent(INOUT)  :: y(*)
           real(kind=C_DOUBLE), intent(INOUT)  :: z(*)
           integer(kind=C_INT) :: asynch
           type(*) :: comm
        end subroutine c_adios2_setup
     end interface
     call c_adios2_setup(npts,nelv,nelb,nelgv, nelgt,x,y,z,asynch,comm)  
  end subroutine fortran_adios2_setup

  subroutine fortran_adios2_finalize()  
     use, intrinsic :: ISO_C_BINDING
     implicit none
     interface
        ! C-definition is: int func(double *data, const int nx, const int ny)
        subroutine c_adios2_finalize() bind(C,name="adios2_finalize_")
           use, intrinsic :: ISO_C_BINDING
           implicit none
        end subroutine c_adios2_finalize
     end interface
     call c_adios2_finalize()  
  end subroutine fortran_adios2_finalize
  
  subroutine fortran_adios2_stream(lglel,p,u,v,w,bm1,t)  
     use, intrinsic :: ISO_C_BINDING
     implicit none
     integer, dimension(:), intent(inout) :: lglel
     real(kind=rp), dimension(:,:,:,:), intent(inout) :: p
     real(kind=rp), dimension(:,:,:,:), intent(inout) :: u
     real(kind=rp), dimension(:,:,:,:), intent(inout) :: v
     real(kind=rp), dimension(:,:,:,:), intent(inout) :: w
     real(kind=rp), dimension(:,:,:,:), intent(inout) :: bm1
     real(kind=rp), dimension(:,:,:,:), intent(inout) :: t

     interface
        ! C-definition is: int func(double *data, const int nx, const int ny)
        subroutine c_adios2_stream(lglel,p,u,v,w,bm1,t) bind(C,name="adios2_stream_")
           use, intrinsic :: ISO_C_BINDING
           implicit none
           integer(kind=C_INT), intent(INOUT) :: lglel(*)
           real(kind=C_DOUBLE), intent(INOUT) :: p(*)
           real(kind=C_DOUBLE), intent(INOUT)  :: u(*)
           real(kind=C_DOUBLE), intent(INOUT)  :: v(*)
           real(kind=C_DOUBLE), intent(INOUT)  :: w(*)
           real(kind=C_DOUBLE), intent(INOUT)  :: bm1(*)
           real(kind=C_DOUBLE), intent(INOUT)  :: t(*)
        end subroutine c_adios2_stream
     end interface
     call c_adios2_stream(lglel,p,u,v,w,bm1,t)  
  end subroutine fortran_adios2_stream

end module data_streamer_2


module user
  use neko
  use data_streamer_2

  implicit none

  !> Variables to store the Rayleigh and Prandlt numbers
  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Re = 0
  real(kind=rp) :: Pr = 0

  type(data_streamer_2_t) :: dstream

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
    integer :: nelb, nelb2, nelv, nelgv,npts, if_asynch

    !> Recalculate the non dimensional parameters
    call json_get(params, 'case.scalar.Pr', Pr)
    call json_get(params, 'case.fluid.Re', Re)
    Ra = (Re**2)*Pr
    write(log_buf,*) 'Rayleigh Number is Ra=', Ra
    call neko_log%message(log_buf)
    
    call dstream%init(coef, 1)


  end subroutine user_initialize

  subroutine user_finalize(t, param)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: param
         
    call dstream%free()
     
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
    
    
    integer ::  lglel(u%msh%nelv)
    integer ::  n, npts, nelv

    call dstream%stream(u,v,w,p,coef)

  end subroutine check

end module user
