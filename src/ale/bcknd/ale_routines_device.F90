module ale_routines_device
  use num_types, only : rp
  use field, only : field_t
  use coefs, only : coef_t
  use field_series, only : field_series_t
  use time_state, only : time_state_t
  use ab_time_scheme, only : ab_time_scheme_t
  use mesh, only : mesh_t
  use utils, only : neko_error
  use device_math, only : device_add2s2
  use math, only : rzero
  use ale_rigid_kinematics
  implicit none
  private

  public :: compute_stiffness_ale_device
  public :: compute_cheap_dist_device
  public :: add_kinematics_to_mesh_velocity_device
  public :: update_ale_mesh_device

contains

  subroutine compute_stiffness_ale_device(coef, params)
    type(coef_t), intent(inout) :: coef
    type(ale_config_t), intent(in) :: params
    call neko_error("ALE: compute_stiffness_ale_device not implemented yet")
  end subroutine compute_stiffness_ale_device

  subroutine compute_cheap_dist_device(d, coef, msh, zone_indices)
    real(kind=rp), intent(inout), target :: d(:)
    type(coef_t), intent(in) :: coef
    type(mesh_t), intent(in) :: msh
    integer, intent(in) :: zone_indices(:)
    call neko_error("ALE: compute_cheap_dist_device not implemented yet")
  end subroutine compute_cheap_dist_device

  subroutine add_kinematics_to_mesh_velocity_device(wx, wy, wz, &
       x_ref, y_ref, z_ref, phi, coef, kinematics, rot_mat, inital_pivot_loc)
    type(field_t), intent(inout) :: wx, wy, wz
    type(field_t), intent(in) :: x_ref, y_ref, z_ref
    type(field_t), intent(in) :: phi
    type(coef_t), intent(in) :: coef
    type(body_kinematics_t), intent(in) :: kinematics
    real(kind=rp), intent(in) :: inital_pivot_loc(3)
    real(kind=rp), intent(in) :: rot_mat(3,3)
    call neko_error("ALE: add_kinematics_to_mesh_velocity_device not implemented yet")
  end subroutine add_kinematics_to_mesh_velocity_device


  subroutine update_ale_mesh_device(c_Xh, wm_x, wm_y, wm_z, wm_x_lag, &
       wm_y_lag, wm_z_lag, time, nadv, scheme_type)

    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: wm_x, wm_y, wm_z
    type(field_series_t), intent(in) :: wm_x_lag, wm_y_lag, wm_z_lag
    type(time_state_t), intent(in) :: time
    type(ab_time_scheme_t) :: ab_scheme_obj
    integer, intent(in) :: nadv
    integer :: j, n
    character(len=*), intent(in) :: scheme_type
    real(kind=rp) :: ab_coeffs(4), dt_history(10), factor

    call rzero(ab_coeffs, 4)
    if (trim(scheme_type) .eq. 'ab') then
       dt_history(1) = time%dt
       dt_history(2) = time%dtlag(1)
       dt_history(3) = time%dtlag(2)
       call ab_scheme_obj%compute_coeffs(ab_coeffs, dt_history, nadv)
    else
       call neko_error("ALE: Unknown mesh time-integration scheme")
    end if

    n = c_Xh%dof%size()

    factor = time%dt * ab_coeffs(1)
    ! for now we leave it like this. Can be replaced by a single routine
    ! that does all three components at once
    call device_add2s2(c_Xh%dof%x_d, wm_x%x_d, factor, n)
    call device_add2s2(c_Xh%dof%y_d, wm_y%x_d, factor, n)
    call device_add2s2(c_Xh%dof%z_d, wm_z%x_d, factor, n)

    ! History Terms
    do j = 2, nadv
       factor = time%dt * ab_coeffs(j)
       call device_add2s2(c_Xh%dof%x_d, wm_x_lag%lf(j - 1)%x_d, factor, n)
       call device_add2s2(c_Xh%dof%y_d, wm_y_lag%lf(j - 1)%x_d, factor, n)
       call device_add2s2(c_Xh%dof%z_d, wm_z_lag%lf(j - 1)%x_d, factor, n)
    end do

  end subroutine update_ale_mesh_device

end module ale_routines_device
