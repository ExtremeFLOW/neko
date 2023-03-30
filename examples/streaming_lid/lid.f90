! Lid-driven cavity
!
! Time-integration of the lid-driven cavity with smoothened
! belt velocity to fulfil continuity equation.
!
module user
  use neko
  use comm
  implicit none

  ! Global user variables
  type(field_t) :: om1, om2, om3, w1, w2

 contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%fluid_user_ic => user_ic
    user%fluid_user_if => user_bc
    user%user_check => user_calc_quantities
    user%user_init_modules => user_initialize
  end subroutine user_setup

  ! user-defined boundary condition
  subroutine user_bc(u, v, w, x, y, z, nx, ny, nz, ix, iy, iz, ie)
    real(kind=rp), intent(inout) :: u
    real(kind=rp), intent(inout) :: v
    real(kind=rp), intent(inout) :: w
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

    real(kind=rp) lsmoothing
    lsmoothing = 0.05_rp    ! length scale of smoothing at the edges

    u = step( x/lsmoothing ) * step( (1._rp-x)/lsmoothing )
    v = 0._rp
    w = 0._rp
    
  end subroutine user_bc
  
  ! User-defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(param_t), intent(inout) :: params
    u = 0._rp
    v = 0._rp
    w = 0._rp
    p = 0._rp
  end subroutine user_ic
  
  ! User-defined initialization called just before time loop starts
  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params

    integer tstep

    ! initialize work arrays for postprocessing
    call field_init(om1, u%dof, 'omega1')
    call field_init(om2, u%dof, 'omega2')
    call field_init(om3, u%dof, 'omega3')
    call field_init(w1, u%dof, 'work1')
    call field_init(w2, u%dof, 'work1')

    ! call usercheck also for tstep=0
    tstep = 0
    call user_calc_quantities(t, tstep, u, v, w, p, coef, params)

  end subroutine user_initialize

  ! User-defined routine called at the end of every time step
  subroutine user_calc_quantities(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    integer :: ntot, i
    real(kind=rp) :: e1, e2

    ! For streaming
    integer ::  nelb, nelv, nelgv, npts
    integer ::  lglel(u%msh%nelv)
    integer :: freq, endstream
    integer :: stream, finalized
    character(len=NEKO_FNAME_LEN) :: fname
    character(len=LOG_SIZE) :: log_buf 
    
    stream = 1
    freq = 100
    endstream = 1000

    !!!!!!!!!!!!!! Streaming section !!!!!!

    if (stream .eq. 1) then

      nelv  = u%msh%nelv
      npts  = u%Xh%lx**3
      nelgv = u%msh%glb_nelv


      if (tstep.eq.1) then

        write(log_buf, '(A)') &
              'Initialize ADIOS2' 

        write(log_buf, '(A,I5.2)') &
              'The communicator is:', NEKO_COMM 
      
        call neko_log%message(log_buf)

        nelb = elem_running_sum(nelv)
        
        nelb = nelb - nelv
    
        write(*,*) 'my sum is', nelb
        write(*,*) 'my number of points is', npts
        write(*,*) 'my number of elements is', nelv
        write(*,*) 'total number of elements is', nelgv

        call adios2_setup(npts,nelv,nelb,nelgv, &
                nelgv,u%dof%x,u%dof%y,  &
                u%dof%z,1,NEKO_COMM)


      else if (mod(tstep+1,freq).eq.0 .and. (tstep+1) .le.endstream) then


        ! Move the data to the CPU to be able to write it
        call device_memcpy(u%x,  u%x_d, &
                       nelv*npts,                         &
                       DEVICE_TO_HOST)

        !Write the fields into a file
        call adios2_stream(lglel,u%x, u%x, & 
                        u%x, u%x, coef%B, u%x)
      
      end if

      if ((tstep+1) .eq. (endstream)) then 
        
        call adios2_finalize()

      end if

    end if


    if (mod(tstep,50).ne.0) return

    ntot = u%dof%size()

    call curl(om1, om2, om3, u, v, w, w1, w2, coef)

    call col3(w1%x,u%x,u%x,ntot)
    call addcol3(w1%x,v%x,v%x,ntot)
    call addcol3(w1%x,w%x,w%x,ntot)
    e1 = 0.5 * glsc2(w1%x,coef%B,ntot) / coef%volume

    call col3(w1%x,om1%x,om1%x,ntot)
    call addcol3(w1%x,om2%x,om2%x,ntot)
    call addcol3(w1%x,om3%x,om3%x,ntot)
    e2 = 0.5 * glsc2(w1%x,coef%B,ntot) / coef%volume
      
    if (pe_rank .eq. 0) &
         &  write(*,'(a,e18.9,a,e18.9,a,e18.9)') &
         &  'POST: t:', t, ' Ekin:', e1, ' enst:', e2

  end subroutine user_calc_quantities

  ! Smooth step function
  function step(x)
  ! Taken from Simson code
  ! x<=0 : step(x) = 0
  ! x>=1 : step(x) = 1
    real(kind=rp) :: step, x

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

  function elem_running_sum(nelv) result(rbuff)

    integer, intent(in) :: nelv
    integer ::  ierr,xbuff,wbuff,rbuff

    xbuff = nelv  ! running sum
    wbuff = nelv  ! working buff
    rbuff = 0   ! recv buff

    call mpi_scan(xbuff,rbuff,1,mpi_integer,mpi_sum,NEKO_COMM,ierr)

  end function elem_running_sum


end module user
