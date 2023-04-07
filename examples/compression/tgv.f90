module user
  use neko
  use comm
  implicit none
  
  character(len=NEKO_FNAME_LEN) :: fname_u
  character(len=NEKO_FNAME_LEN) :: fname_v
  character(len=NEKO_FNAME_LEN) :: fname_w
  character(len=NEKO_FNAME_LEN) :: fname_p
  type(file_t) :: mf_u
  type(file_t) :: mf_v
  type(file_t) :: mf_w
  type(file_t) :: mf_p

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
    integer, intent(in) :: tstep    
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
  
    ! This step is here only temporarly, to check decompressed data 
    if (tstep.eq.1) then 
      fname_u = 'rct_u.fld'
      fname_v = 'rct_v.fld'
      fname_w = 'rct_w.fld'
      fname_p = 'rct_p.fld'
      mf_u =  file_t(fname_u)
      mf_v =  file_t(fname_v)
      mf_w =  file_t(fname_w)
      mf_p =  file_t(fname_p)
    end if

    call in_situ_initialize( t, tstep,u, v, w, p, coef, params)

    call in_situ_compress( t, tstep,u, v, w, p, coef, params)


  end subroutine usr_calc_quantities


  subroutine in_situ_initialize( t, tstep,u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep    
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    !
    character(len=LOG_SIZE) :: log_buf 
    integer ::  nelb, nelv, nelgv, npts

    if (tstep.ne.1) return
    
    ! Calculate case parameters
    nelv  = u%msh%nelv
    npts  = u%Xh%lx**3
    nelgv = u%msh%glb_nelv

    write(log_buf, '(A)') &
            'Initialize ADIOS2' 

    write(log_buf, '(A,I5.2)') &
           'The communicator is:', NEKO_COMM 
      
    call neko_log%message(log_buf)

    nelb = in_situ_elem_running_sum(nelv)
        
    nelb = nelb - nelv
    
      !write(*,*) 'my sum is', nelb
      !write(*,*) 'my number of points is', npts
      !write(*,*) 'my number of elements is', nelv
      !write(*,*) 'total number of elements is', nelgv

    call adios2_setup(npts,nelv,nelb,nelgv, &
              nelgv,u%dof%x,u%dof%y,  &
              u%dof%z,NEKO_COMM)

  end subroutine


  subroutine in_situ_compress( t, tstep,u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep    
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    !
    type(cpr_t) :: cpr_u
    type(cpr_t) :: cpr_v
    type(cpr_t) :: cpr_w
    type(cpr_t) :: cpr_p
    integer i, nelv
    integer ::  lglel(u%msh%nelv)
    !
    integer compression_frequency
    integer reconstruct_and_write

    compression_frequency = 10 
    reconstruct_and_write = 1

    if (mod(tstep,compression_frequency).ne.0) return

    call neko_log%section('Compression step')       

    ! Initialize the fields and transform matrices
    call cpr_init(cpr_u,u)
    call cpr_init(cpr_v,v)
    call cpr_init(cpr_w,w)
    call cpr_init(cpr_p,p)

    ! truncate the spectral coefficients
    call cpr_truncate_wn(cpr_u,coef)
    call cpr_truncate_wn(cpr_v,coef)
    call cpr_truncate_wn(cpr_w,coef)
    call cpr_truncate_wn(cpr_p,coef)

    ! Write global element number ==== This is temporal
    do i = 1, u%msh%nelv
      lglel(i) = i
    end do
   
    ! Write with adios2
    call in_situ_adios2_write(lglel, cpr_u,cpr_v,cpr_w,cpr_p)

    ! Check the compression ratio
    call in_situ_compression_check_cr(lglel, cpr_u,cpr_v, &
                                          cpr_w,cpr_p)

    if (reconstruct_and_write .eq. 1) then
      call cpr_goto_space(cpr_u,'phys') !< 'spec' / 'phys'
      call cpr_goto_space(cpr_v,'phys') !< 'spec' / 'phys'
      call cpr_goto_space(cpr_w,'phys') !< 'spec' / 'phys'
      call cpr_goto_space(cpr_p,'phys') !< 'spec' / 'phys'

      call in_situ_mf_write(t,tstep,lglel, cpr_u,cpr_v,cpr_w,cpr_p)
  
      call in_situ_compression_check_error(lglel, cpr_u,cpr_v, &
                                       cpr_w,cpr_p,            &
                                       u, v, w, p, coef, params)
    end if

    ! Free the memory allocated for the fields
    call cpr_free(cpr_u)
    call cpr_free(cpr_v)
    call cpr_free(cpr_w)
    call cpr_free(cpr_p)

    call neko_log%end_section()       

  end subroutine

  
  subroutine in_situ_compression_check_error(lglel, cpr_u,cpr_v, &
                                      cpr_w,cpr_p, &
                                      u, v, w, p, coef, params)
    type(coef_t), intent(in) :: coef
    type(param_t), intent(in) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(cpr_t), intent(inout) :: cpr_u
    type(cpr_t), intent(inout) :: cpr_v
    type(cpr_t), intent(inout) :: cpr_w
    type(cpr_t), intent(inout) :: cpr_p
    integer, intent(inout) ::  lglel(cpr_u%fld%msh%nelv)
    real(kind=rp) :: restemp_a, res_a
    real(kind=rp) :: ev(u%Xh%lx, u%Xh%lx, u%Xh%lx,u%msh%nelv) 
    integer :: nelv, nelgv, npts


    nelv  = u%msh%nelv
    npts  = u%Xh%lx**3
    nelgv = u%msh%glb_nelv
    
    call device_memcpy(cpr_u%fldhat%x,  cpr_u%fldhat%x_d, &
               (cpr_u%fld%msh%nelv)*(cpr_u%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)
    call device_memcpy(cpr_v%fldhat%x,  cpr_v%fldhat%x_d, &
               (cpr_v%fld%msh%nelv)*(cpr_v%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)

    call device_memcpy(cpr_w%fldhat%x,  cpr_w%fldhat%x_d, &
               (cpr_w%fld%msh%nelv)*(cpr_w%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)

    call device_memcpy(cpr_p%fldhat%x,  cpr_p%fldhat%x_d, &
               (cpr_p%fld%msh%nelv)*(cpr_p%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)
    
    call device_memcpy(u%x,  u%x_d, &
               (cpr_u%fld%msh%nelv)*(cpr_u%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)
    call device_memcpy(v%x,  v%x_d, &
               (cpr_v%fld%msh%nelv)*(cpr_v%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)

    call device_memcpy(w%x,  w%x_d, &
               (cpr_w%fld%msh%nelv)*(cpr_w%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)

    call device_memcpy(p%x,  p%x_d, &
               (cpr_p%fld%msh%nelv)*(cpr_p%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)
   
    ! Get error vector 
    call sub3(ev,cpr_u%fldhat%x,u%x,nelv*npts)
    restemp_a = glsc3(ev,coef%B,ev,nelv*npts)
    res_a = sqrt(restemp_a)/sqrt(coef%volume)
    if (pe_rank .eq. 0) then
      write(*,*) 'Reconstruction error - u = ', res_a
    end if
    
    call sub3(ev,cpr_v%fldhat%x,v%x,nelv*npts)
    restemp_a = glsc3(ev,coef%B,ev,nelv*npts)
    res_a = sqrt(restemp_a)/sqrt(coef%volume)
    if (pe_rank .eq. 0) then
      write(*,*) 'Reconstruction error - v = ', res_a
    end if
    
    call sub3(ev,cpr_w%fldhat%x,w%x,nelv*npts)
    restemp_a = glsc3(ev,coef%B,ev,nelv*npts)
    res_a = sqrt(restemp_a)/sqrt(coef%volume)
    if (pe_rank .eq. 0) then
      write(*,*) 'Reconstruction error - w = ', res_a
    end if
    
    call sub3(ev,cpr_p%fldhat%x,p%x,nelv*npts)
    restemp_a = glsc3(ev,coef%B,ev,nelv*npts)
    res_a = sqrt(restemp_a)/sqrt(coef%volume)
    if (pe_rank .eq. 0) then
      write(*,*) 'Reconstruction error - p = ', res_a
    end if

  end subroutine
  
  
  subroutine in_situ_compression_check_cr(lglel, cpr_u,cpr_v, &
                                          cpr_w,cpr_p)
    type(cpr_t), intent(inout) :: cpr_u
    type(cpr_t), intent(inout) :: cpr_v
    type(cpr_t), intent(inout) :: cpr_w
    type(cpr_t), intent(inout) :: cpr_p
    integer, intent(inout) ::  lglel(cpr_u%fld%msh%nelv)
    real(kind=rp) :: restemp_a, res_a
    real(kind=rp) :: ev(cpr_u%fld%Xh%lx, cpr_u%fld%Xh%lx, &
                        cpr_u%fld%Xh%lx,cpr_u%fld%msh%nelv) 
    integer :: i, nelv, nelgv, npts


    nelv  = cpr_u%fld%msh%nelv
    npts  = cpr_u%fld%Xh%lx**3
    nelgv = cpr_u%fld%msh%glb_nelv
    
    call device_memcpy(cpr_u%fldhat%x,  cpr_u%fldhat%x_d, &
               (cpr_u%fld%msh%nelv)*(cpr_u%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)
    call device_memcpy(cpr_v%fldhat%x,  cpr_v%fldhat%x_d, &
               (cpr_v%fld%msh%nelv)*(cpr_v%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)

    call device_memcpy(cpr_w%fldhat%x,  cpr_w%fldhat%x_d, &
               (cpr_w%fld%msh%nelv)*(cpr_w%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)

    call device_memcpy(cpr_p%fldhat%x,  cpr_p%fldhat%x_d, &
               (cpr_p%fld%msh%nelv)*(cpr_p%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)
     
    ! Count how many coefficients have been set to zero
    do i=1,nelv*npts
      if (cpr_u%fldhat%x(i,1,1,1) .le. 1e-15) then
        ev(i,1,1,1) = 1
      else
        ev(i,1,1,1) = 0
      end if
    end do
    res_a = glsum(ev,nelv*npts)
    if (pe_rank .eq. 0) then
    write(*,*) '% Data discarded - u = ', res_a/(nelgv*npts)*100
    end if
    
    ! Count how many coefficients have been set to zero
    do i=1,nelv*npts
      if (cpr_v%fldhat%x(i,1,1,1) .le. 1e-15) then
        ev(i,1,1,1) = 1
      else
        ev(i,1,1,1) = 0
      end if
    end do
    res_a = glsum(ev,nelv*npts)
    if (pe_rank .eq. 0) then
      write(*,*) '% Data discarded - v = ', res_a/(nelgv*npts)*100
    end if

    ! Count how many coefficients have been set to zero
    do i=1,nelv*npts
      if (cpr_w%fldhat%x(i,1,1,1) .le. 1e-15) then
        ev(i,1,1,1) = 1
      else
        ev(i,1,1,1) = 0
      end if
    end do
    res_a = glsum(ev,nelv*npts)
    if (pe_rank .eq. 0) then
      write(*,*) '% Data discarded - w = ', res_a/(nelgv*npts)*100
    end if
    
    ! Count how many coefficients have been set to zero
    do i=1,nelv*npts
      if (cpr_p%fldhat%x(i,1,1,1) .le. 1e-15) then
        ev(i,1,1,1) = 1
      else
        ev(i,1,1,1) = 0
      end if
    end do
    res_a = glsum(ev,nelv*npts)
    if (pe_rank .eq. 0) then
      write(*,*) '% Data discarded - p = ', res_a/(nelgv*npts)*100
    end if
  end subroutine

  subroutine in_situ_adios2_write(lglel, cpr_u,cpr_v,cpr_w,cpr_p)
    type(cpr_t), intent(inout) :: cpr_u
    type(cpr_t), intent(inout) :: cpr_v
    type(cpr_t), intent(inout) :: cpr_w
    type(cpr_t), intent(inout) :: cpr_p
    integer, intent(inout) ::  lglel(cpr_u%fld%msh%nelv)


    ! Put data in cpu to let adios2 write it
    call device_memcpy(cpr_u%fldhat%x,  cpr_u%fldhat%x_d, &
               (cpr_u%fld%msh%nelv)*(cpr_u%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)
    call device_memcpy(cpr_v%fldhat%x,  cpr_v%fldhat%x_d, &
               (cpr_v%fld%msh%nelv)*(cpr_v%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)

    call device_memcpy(cpr_w%fldhat%x,  cpr_w%fldhat%x_d, &
               (cpr_w%fld%msh%nelv)*(cpr_w%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)

    call device_memcpy(cpr_p%fldhat%x,  cpr_p%fldhat%x_d, &
               (cpr_p%fld%msh%nelv)*(cpr_p%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)
    
    !Write data to file though adios2
    call adios2_update(lglel, cpr_p%fldhat%x, cpr_u%fldhat%x, &
            cpr_v%fldhat%x, cpr_w%fldhat%x, cpr_u%fldhat%x)

  end subroutine


  subroutine in_situ_mf_write(t,tstep,lglel, cpr_u,cpr_v,cpr_w,cpr_p)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep    
    type(cpr_t), intent(inout) :: cpr_u
    type(cpr_t), intent(inout) :: cpr_v
    type(cpr_t), intent(inout) :: cpr_w
    type(cpr_t), intent(inout) :: cpr_p
    integer, intent(inout) ::  lglel(cpr_u%fld%msh%nelv)


    ! Put data in cpu to let adios2 write it
    call device_memcpy(cpr_u%fldhat%x,  cpr_u%fldhat%x_d, &
               (cpr_u%fld%msh%nelv)*(cpr_u%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)
    call device_memcpy(cpr_v%fldhat%x,  cpr_v%fldhat%x_d, &
               (cpr_v%fld%msh%nelv)*(cpr_v%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)

    call device_memcpy(cpr_w%fldhat%x,  cpr_w%fldhat%x_d, &
               (cpr_w%fld%msh%nelv)*(cpr_w%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)

    call device_memcpy(cpr_p%fldhat%x,  cpr_p%fldhat%x_d, &
               (cpr_p%fld%msh%nelv)*(cpr_p%fld%Xh%lx**3), &
                       DEVICE_TO_HOST)
   
    call mf_u%write(cpr_u%fldhat,t)
    call mf_v%write(cpr_v%fldhat,t)
    call mf_w%write(cpr_w%fldhat,t)
    call mf_p%write(cpr_p%fldhat,t)


  end subroutine


  function in_situ_elem_running_sum(nelv) result(rbuff)

    integer, intent(in) :: nelv
    integer ::  ierr,xbuff,wbuff,rbuff

    xbuff = nelv  ! running sum
    wbuff = nelv  ! working buff
    rbuff = 0   ! recv buff

    call mpi_scan(xbuff,rbuff,1,mpi_integer,mpi_sum,NEKO_COMM,ierr)

  end function in_situ_elem_running_sum

end module user
