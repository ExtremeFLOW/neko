
module user
  use neko
  use fld_file_output
  use FST
  implicit none

  type(FST_t) :: FST_obj
  character(len=LOG_SIZE) :: log_buf ! For logging status
  class(point_zone_t), pointer :: fstbox

  type(fld_file_output_t) :: fout

  integer :: counter = 0

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules     => initialize ! Initialize parameters
    u%user_check            => check
    u%fluid_user_f_vector   => forcing
    u%user_finalize_modules => finalize ! Finalize
  end subroutine user_setup

  subroutine initialize(t, u, v, w, p, coef, params)
    implicit none

    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    type(box_point_zone_t), pointer :: bpz
 
    type(field_t), pointer :: fu,fv,fw

    fstbox => neko_point_zone_registry%get_point_zone("fstbox")

    select type(fstbox)
    type is (box_point_zone_t)
       bpz => fstbox
       call neko_log%message("Found box zone " // bpz%name)
       write (log_buf, *) "xmin", bpz%xmin, "xmax", bpz%xmax
       call neko_log%message(log_buf)
    class default
       call neko_error("Wrong point zone type")
    end select

    call fout%init(sp, "forcing", 3)

    call neko_field_registry%add_field(u%dof, 'fu')
    call neko_field_registry%add_field(u%dof, 'fv')
    call neko_field_registry%add_field(u%dof, 'fw')
    fu => neko_field_registry%get_field('fu')
    fv => neko_field_registry%get_field('fv')
    fw => neko_field_registry%get_field('fw')

    call rzero(fu%x, fu%dof%size())
    call rzero(fv%x, fv%dof%size())
    call rzero(fw%x, fw%dof%size())

    !
    ! ========== FST ============
    !
    ! Initialize all FST parameters
    call FST_obj%init(params, -0.05_rp, 0.05_rp, bpz%xmin, bpz%xmax,&
         100000.0_rp, 0.002_rp, 0.002_rp)

    ! Generate FST (mostly files for postprocessing)
    call FST_obj%generate(u)

    ! ===========================

  end subroutine initialize
! usrcheck, this is called at the end of every time step

  subroutine check(t, tstep,u, v, w, p, coef, param)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: param

  end subroutine check

  ! Define the forcing here
  subroutine forcing(f,t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t

    ! Access velocity fields without having them as arguments
    type(field_t), pointer :: u, v, w, fu, fv, fw
    integer :: nlindex(4), i, lx, j, k, l, e

    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    fu => neko_field_registry%get_field('fu')
    fv => neko_field_registry%get_field('fv')
    fw => neko_field_registry%get_field('fw')

    lx = u%Xh%lx

    do i = 1, f%dm%size() !fstbox%size
       nlindex = nonlinear_index(i, lx, lx, lx)
       j = nlindex(1)
       k = nlindex(2)
       l = nlindex(3)
       e = nlindex(4)
       call FST_forcing_local(FST_obj, u%dof%x(j,k,l,e), u%dof%y(j,k,l,e), &
            u%dof%z(j,k,l,e), t, u%x(j,k,l,e), v%x(j,k,l,e), w%x(j,k,l,e), &
            ! I put the chi in now
            f%u(j,k,l,e), f%v(j,k,l,e), f%w(j,k,l,e),f%chi(j,k,l,e), i)
    end do

    counter = counter + 1 
    if (mod(counter, 10) .eq. 0) then
       
       call copy(fu%x, f%u, fu%dof%size())
       call copy(fv%x, f%v, fu%dof%size())
       call copy(fw%x, f%w, fu%dof%size())

       fout%fields%fields(1)%f => fu
       fout%fields%fields(2)%f => fv
       fout%fields%fields(3)%f => fw
       call fout%sample(t)

    end if

    !write (*,*) "HHUHUH Bef", f%u(1:5,1,1,1)
    !call FST_obj%FST_forcing_zone(f,u,v,w,t,fstbox)

  end subroutine forcing

  ! Free relevant objects
  subroutine finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    call FST_obj%free()

  end subroutine finalize

end module user
