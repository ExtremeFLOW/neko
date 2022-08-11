module user
  use neko
  implicit none
contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_usr_ic => user_ic
    u%usr_msh_setup => user_mesh_scale
    u%usr_chk => usr_calc_quantities
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

  subroutine usr_calc_quantities( t, dt, tstep,u, v, w, p, coef)
    real(kind=rp), intent(in) :: t, dt
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(cpr_t) :: cpr_u
    integer :: i
    character(len=LOG_SIZE) :: log_buf 


    if (mod(tstep,10).ne.0) return
    
    call neko_log%section('Compression')       

    ! Initialize the fields and transform matrices
    call cpr_init(cpr_u,u)

    ! truncate the spectral coefficients
    call cpr_truncate_wn(cpr_u,coef)
  
    ! just to check, go to physical space and compare
    call cpr_goto_space(cpr_u,'phys') !< 'spec' / 'phys'
    ! chech that the copy is fine in one entry
    do i = 1, 10
      write(log_buf, '(A,E15.7,A,E15.7)') &
            'u value:', cpr_u%fld%x(i,1,1,10), &
            ' reconstructed u value:', &
            cpr_u%fldhat(i,1,1,10)
      call neko_log%message(log_buf)
    enddo

    ! Free the memory allocated for the fields
    call cpr_free(cpr_u)


  end subroutine usr_calc_quantities


end module user
