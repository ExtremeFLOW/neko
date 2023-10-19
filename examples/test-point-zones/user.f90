
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
    u%fluid_user_ic         => set_ic
  end subroutine user_setup

  !> User initial condition
  subroutine set_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    type(field_t), pointer :: myfield
    class(point_zone_t), pointer :: pz
    type(box_point_zone_t), pointer :: bpz => null()
    type(sphere_point_zone_t), pointer :: spz => null()
    real(kind=rp) :: x,y,z
    integer :: i,j, nlindex(4)

    myfield => neko_field_registry%get_field('s')

    do i = 1, neko_point_zone_registry%n_point_zones()

       pz => neko_point_zone_registry%get_point_zone(i)

       select type(pz)
       type is (box_point_zone_t)
          bpz => pz
          write (*, '(A,A,A)') "Point zone ", trim(bpz%name), " is a box"
          write (*,'(6(F10.6," "))') bpz%xmin, bpz%xmax, bpz%ymin, bpz%ymax, bpz%zmin, bpz%zmax

       type is (sphere_point_zone_t)
          spz => pz
          write (*, '(A,A,A)') "Point zone ", trim(spz%name), " is a sphere"
          write (*,'(4(F10.6, " "))') spz%x0, spz%y0, spz%z0, spz%radius
       class default
          call neko_error("not found")
       end select

       do j = 1, pz%size
          nlindex = nonlinear_index(pz%mask(j), u%Xh%lx, u%Xh%lx, u%Xh%lx)
          myfield%x(nlindex(1), nlindex(2), nlindex(3), nlindex(4)) = 10.0_rp
       end do

    end do

  end subroutine set_ic

  subroutine initialize(t, u, v, w, p, coef, params)
    implicit none

    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

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
