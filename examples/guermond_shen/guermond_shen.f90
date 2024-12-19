module user
  use neko
  implicit none

  type(file_t) :: csv
  type(vector_t) :: linf_error

contains
  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_user_f_vector => set_source
    u%user_init_modules => user_initialize
    u%user_finalize_modules => user_finalize
    u%user_check => user_compute
  end subroutine user_setup

  !> Set source term
  subroutine set_source(f, t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    real(kind=rp) :: x, y
    integer :: i
    real(kind=rp) :: pi=4.0_rp * DATAN(1.0_rp)

    do i = 1, f%dm%size()
       x = f%dm%x(i,1,1,1)
       y = f%dm%y(i,1,1,1)

       f%u(i,1,1,1) = -pi * sin(t) * sin(pi*x) * sin(pi*y) - &
                      2 * pi**3 * cos(pi*x)**2 * sin(t) * sin(2*pi*y) + &
                      pi * cos(t) * sin(pi*x)**2 * sin(2*pi*y) + &
                      6 * pi**3 * sin(t) * sin(pi*x)**2 * sin(2*pi*y)
       f%v(i,1,1,1) = pi * cos(pi*x) * cos(pi*y) * sin(t) + &
                      2 * pi**3 * cos(pi*y)**2 * sin(t) * sin(2*pi*x) - &
                      pi * cos(t) * sin(pi*y)**2 * sin(2*pi*x) - &
                      6 * pi**3 * sin(t) * sin(pi*y)**2 * sin(2*pi*x)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(f%u, f%u_d, f%dm%size(), &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(f%v, f%v_d, f%dm%size(), &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(f%w, f%w_d, f%dm%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine set_source

  ! The exact solution based on the applied forcing.
  subroutine true_solution(u, v, p, x, y, t)
    real(kind=rp), intent(inout) :: u, v, p
    real(kind=rp), intent(in) :: x, y
    real(kind=rp), intent(in) :: t

    real(kind=rp) :: pi=4.0_rp * DATAN(1.0_rp)

    u =  pi * sin(t) * sin(2 * pi * y) * sin(pi *x) **2
    v = -pi * sin(t) * sin(2 * pi * x) * sin(pi *y) **2
    p =       sin(t) * sin(pi * y)     * cos(pi *x)

  end subroutine true_solution

  subroutine user_compute(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p

    real(kind=rp) :: true_u, true_v, true_p, x, y
    integer :: i
    type(field_t), pointer :: error_u, error_v, error_p

    error_u => neko_field_registry%get_field("error_u")
    error_v => neko_field_registry%get_field("error_v")
    error_p => neko_field_registry%get_field("error_p")

    do i = 1, u%size()
       x = u%dof%x(i,1,1,1)
       y = u%dof%y(i,1,1,1)

       call true_solution(true_u, true_v, true_p, x, y, t)

       error_u%x(i,1,1,1) = u%x(i,1,1,1) - true_u
       error_v%x(i,1,1,1) = v%x(i,1,1,1) - true_v
       error_p%x(i,1,1,1) = p%x(i,1,1,1) - true_p

    end do

    linf_error%x(1) = glmax(error_u%x, u%size())
    linf_error%x(2) = glmax(error_v%x, u%size())
    linf_error%x(3) = glmax(error_p%x, u%size())

    select type (csv_file => csv%file_type)
      type is (csv_file_t)
       if (pe_rank .eq. 0) then
          call csv_file%write(linf_error)
       end if
    end select

  end subroutine user_compute

  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    character(len=:), allocatable :: fname

    real(kind=rp) dt
    integer tstep

    ! Redundant if the field_writer is used, because it will register the fields
    call neko_field_registry%add_field(u%dof, "error_u", ignore_existing=.true.)
    call neko_field_registry%add_field(u%dof, "error_v", ignore_existing=.true.)
    call neko_field_registry%add_field(u%dof, "error_p", ignore_existing=.true.)

    call json_get_or_default(params, "case.error_file", fname, "error.csv")
    csv = file_t(trim(fname))
    call linf_error%init(3)

  end subroutine user_initialize

  subroutine user_finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    call linf_error%free()
  end subroutine user_finalize

end module user
