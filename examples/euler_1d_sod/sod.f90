! Sod's shock tube problem
!
module user
  use neko
  use fluid_scheme_compressible, only : fluid_scheme_compressible_t
  implicit none

  ! Global user variables
  type(field_t) :: w1

  type(file_t) output_file ! output file
  type(vector_t) :: vec_out ! will store our output data
  logical :: case_ready = .false.

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initialize => initialize
    user%initial_conditions => initial_conditions
    user%material_properties => material_properties
  end subroutine user_setup

  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    case_ready = .true.
  end subroutine initialize

  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type (field_t), pointer :: rho, u, v, w, p
    integer :: i
    real(kind=rp) :: x, mux

    rho => fields%get_by_name("fluid_rho")
    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")
    p => fields%get_by_name("p")


    ! rho_L = 1, u_L = (0, 0), p_L = 1
    ! rho_R = 0.125, u_R = (0, 0), p_R = 0.1
    ! end time = 0.2
    mux = 0.5
    do i = 1, rho%dof%size()
       x = rho%dof%x(i,1,1,1)

       u%x(i,1,1,1) = 0.0
       v%x(i,1,1,1) = 0.0
       w%x(i,1,1,1) = 0.0

       if (x < mux) then
          rho%x(i,1,1,1) = 1.0
          p%x(i,1,1,1) = 1.0
       else
          rho%x(i,1,1,1) = 0.125
          p%x(i,1,1,1) = 0.1
       end if
    end do
  end subroutine initial_conditions

  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time

    integer :: i, n, nx
    real(kind=rp) :: mu

    if (scheme_name .ne. "fluid") return

    if (.not. case_ready) return

    select type (fluid => neko_user_access%case%fluid)
    class is (fluid_scheme_compressible_t)
       if (.not. associated(fluid%artificial_visc)) return

       n = fluid%artificial_visc%dof%size()
       nx = fluid%msh%glb_nelv
       mu = 1.0_rp / (nx * 30)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(fluid%artificial_visc%x, &
               fluid%artificial_visc%x_d, n, DEVICE_TO_HOST, sync = .true.)
       end if

       do i = 1, n
          fluid%artificial_visc%x(i,1,1,1) = &
               max(fluid%artificial_visc%x(i,1,1,1), mu)
       end do

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(fluid%artificial_visc%x, &
               fluid%artificial_visc%x_d, n, HOST_TO_DEVICE, sync = .false.)
       end if
    end select
  end subroutine material_properties

end module user
