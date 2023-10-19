
module user
  use neko
  implicit none

  character(len=LOG_SIZE) :: log_buf ! For logging status

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules     => initialize ! Initialize parameters
    u%user_check            => check
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

    class(*), pointer :: pz
    type(box_point_zone_t), pointer :: bpz => null()
    type(sphere_point_zone_t), pointer :: spz => null()


    pz => neko_point_zone_registry%get_point_zone("myzone")

    write (*,*) type(pz)

    select type(pz)
    type is (box_point_zone_t)
       bpz => pz
       write (*,*) bpz%xmin, bpz%xmax, bpz%ymin, bpz%ymax, bpz%zmin, bpz%zmax
    type is (sphere_point_zone_t)
       spz => pz
       write (*,*) spz%x0, spz%y0, spz%z0, spz%radius
    class default
       call neko_error("not found")
    end select

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

  ! Free relevant objects
  subroutine finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

  end subroutine finalize


end module user
