module user
  use neko
  use comm
  implicit none
contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_user_ic => user_ic
    u%user_mesh_setup => user_mesh_scale
    u%user_check => usr_calc_quantities
  end subroutine user_setup

  ! Normalize mesh
  
  subroutine user_mesh_scale(msh)
    type(mesh_t), intent(inout) :: msh
    type(point_t), pointer :: pt
    integer :: i, p, npts
    real(kind=rp) :: pi, d
    pi = 4d0*atan(1d0)
    d = 4d0

    npts = size(msh%points)
    do i = 1, npts
       msh%points(i)%x(1) = (msh%points(i)%x(1) - d) / d * pi
       msh%points(i)%x(2) = (msh%points(i)%x(2) - d) / d * pi
       msh%points(i)%x(3) = (msh%points(i)%x(3) - d) / d * pi
    end do
    
  end subroutine user_mesh_scale



  ! User defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(param_t), intent(inout) :: params
    integer :: i
    real(kind=rp) :: uvw(3)

    do i = 1, u%dof%size()
       uvw = tgv_ic(u%dof%x(i,1,1,1),u%dof%y(i,1,1,1),u%dof%z(i,1,1,1))
       u%x(i,1,1,1) = uvw(1)
       v%x(i,1,1,1) = uvw(2)
       w%x(i,1,1,1) = uvw(3)
    end do
    p = real(0d0,rp)
  end subroutine user_ic
  
  function tgv_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: ux, uy, uz
    real(kind=rp) :: uvw(3)
    real(kind=rp), parameter :: zero = 0d0
    integer e,eg

    uvw(1)   = sin(x)*cos(y)*cos(z)
    uvw(2)   = -cos(x)*sin(y)*cos(z)
    uvw(3)   = zero
  end function tgv_ic

  subroutine usr_calc_quantities( t, tstep,u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    real(kind=rp) :: restemp_a, res_a
    integer, intent(in) :: tstep    
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(file_t) :: mf
    type(cpr_t) :: cpr_u
    integer :: i, freq, maxstream
    character(len=NEKO_FNAME_LEN) :: fname
    character(len=LOG_SIZE) :: log_buf 
    ! Definitions for the l2norm in device
    real(kind=rp) :: ev(u%Xh%lx, u%Xh%lx, u%Xh%lx,u%msh%nelv) 
    real(kind=rp) :: tmp(u%Xh%lx, u%Xh%lx, u%Xh%lx,u%msh%nelv) 

    !---- For running sum - compression test
    integer ::  nelb, nelv, nelgv, npts
    integer ::  lglel(u%msh%nelv)
    !----

    nelv  = u%msh%nelv
    npts  = u%Xh%lx**3
    nelgv = u%msh%glb_nelv


    call neko_log%section('Compression')       
   
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

    end if

    freq = 50
    maxstream = 2000

    if (tstep.eq.maxstream) then
      call adios2_finalize()
    end if

    call neko_log%end_section()       


    if (mod(tstep+1,freq).ne.0.or.tstep.ge.maxstream) return

    call neko_log%section('Stream data')       


    ! synchronize CPU
    ! Move the data to the CPU to be able to write it
    call device_memcpy(u%x,  u%x_d, &
                       nelv*npts,                         &
                       DEVICE_TO_HOST)

    !Write the fields into a file
    call adios2_stream(lglel,u%x, u%x, & 
                      u%x, u%x, coef%B, u%x)

    call neko_log%end_section()       

  end subroutine usr_calc_quantities

     
  function elem_running_sum(nelv) result(rbuff)

    integer, intent(in) :: nelv
    integer ::  ierr,xbuff,wbuff,rbuff

    xbuff = nelv  ! running sum
    wbuff = nelv  ! working buff
    rbuff = 0   ! recv buff

    call mpi_scan(xbuff,rbuff,1,mpi_integer,mpi_sum,NEKO_COMM,ierr)

  end function elem_running_sum

end module user
