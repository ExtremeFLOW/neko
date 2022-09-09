module neko_intf
  use neko
  use, intrinsic :: iso_c_binding
  implicit none
  private

contains

  subroutine neko_intf_init() bind(c, name='intf_init')
    character(len=LOG_SIZE) :: log_buf
    character(10) :: time
    character(8) :: date

    call neko_init()

    if (pe_rank .eq. 0) then
       write(*,*) ''
       write(*,*) '   _  __  ____  __ __  ____ '
       write(*,*) '  / |/ / / __/ / //_/ / __ \'
       write(*,*) ' /    / / _/  / ,<   / /_/ /'
       write(*,*) '/_/|_/ /___/ /_/|_|  \____/ '
       write(*,*) ''
       write(*,*) '(version: ', trim(NEKO_VERSION),')'
       write(*,*) trim(NEKO_BUILD_INFO)
       write(*,*) ''
    end if

    call date_and_time(time=time, date=date)           
    call neko_log%section("Session Information")       
    write(log_buf, '(A,A,A,A,1x,A,1x,A,A,A,A,A)') 'Start time: ',&
         time(1:2),':',time(3:4), '/', date(1:4),'-', date(5:6),'-',date(7:8)
    call neko_log%message(log_buf)
    write(log_buf, '(a)') 'Running on: '
    if (pe_size .lt. 1e1)  then
       write(log_buf(13:), '(i1,a)') pe_size, ' MPI '
       if (pe_size .eq. 1) then
          write(log_buf(19:), '(a)') 'rank'
       else
          write(log_buf(19:), '(a)') 'ranks'
       end if
    else if (pe_size .lt. 1e2) then
       write(log_buf(13:), '(i2,a)') pe_size, ' MPI ranks'
    else if (pe_size .lt. 1e3) then
       write(log_buf(13:), '(i3,a)') pe_size, ' MPI ranks'
    else if (pe_size .lt. 1e4) then
       write(log_buf(13:), '(i4,a)') pe_size, ' MPI ranks'
    else if (pe_size .lt. 1e5) then
       write(log_buf(13:), '(i5,a)') pe_size, ' MPI ranks'
    else
       write(log_buf(13:), '(i6,a)') pe_size, ' MPI ranks'
    end if
    call neko_log%message(log_buf)

    write(log_buf, '(a)') 'Bcknd type: '
    if (NEKO_BCKND_SX .eq. 1) then
       write(log_buf(13:), '(a)') 'SX-Aurora'
    else if (NEKO_BCKND_XSMM .eq. 1) then
       write(log_buf(13:), '(a)') 'CPU (libxsmm)'
    else if (NEKO_BCKND_CUDA .eq. 1) then
       write(log_buf(13:), '(a)') 'Accelerator (CUDA)'
    else if (NEKO_BCKND_HIP .eq. 1) then
       write(log_buf(13:), '(a)') 'Accelerator (HIP)'
    else if (NEKO_BCKND_OPENCL .eq. 1) then
       write(log_buf(13:), '(a)') 'Accelerator (OpenCL)'
    else
       write(log_buf(13:), '(a)') 'CPU'
    end if
    call neko_log%message(log_buf)

    if (NEKO_BCKND_HIP .eq. 1 .or. NEKO_BCKND_CUDA .eq. 1 .or. &
         NEKO_BCKND_OPENCL .eq. 1) then
       write(log_buf, '(a)') 'Dev. name : '
       call device_name(log_buf(13:))
       call neko_log%message(log_buf)
    end if

    write(log_buf, '(a)') 'Real type : '
    select case (rp)
    case (real32)
       write(log_buf(13:), '(a)') 'single precision'
    case (real64)
       write(log_buf(13:), '(a)') 'double precision'
    case (real128)
       write(log_buf(13:), '(a)') 'quad precision'
    end select
    call neko_log%message(log_buf)       
    call neko_log%end()
    call neko_log%newline

  end subroutine neko_intf_init

end module neko_intf
