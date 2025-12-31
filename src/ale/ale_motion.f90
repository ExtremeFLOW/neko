!> ALE Motion(s) Logic
module ale_motion
  use num_types, only : rp
  use field, only : field_t
  use coefs, only : coef_t
  use math, only : pi, rzero
  use field_series, only : field_series_t
  use time_state, only : time_state_t
  use ab_time_scheme, only : ab_time_scheme_t
  use utils, only : neko_error
  use neko_config, only : NEKO_BCKND_DEVICE
  use comm, only : pe_rank

  implicit none
  private

  public :: update_ale_mesh_velocity
  public :: init_ale_pivot
  public :: update_ale_mesh
  public :: update_ale_mass_history
  public :: ale_config_t
  public :: pivot_state_t
  
  !> Container for ALE parameters (read from JSON)
  type, public :: ale_config_t
     !> Stiffness Control (for ale_solver)
     character(len=32) :: stiffness_type = 'built-in'
     real(kind=rp)     :: stiffness_gain  = 9.0_rp
     real(kind=rp)     :: stiffness_decay = 2.0_rp
     logical           :: if_output_phi  = .true.
     !> Oscillation (x, y, z)
     real(kind=rp) :: osc_amp(3)  = 0.0_rp
     real(kind=rp) :: osc_freq(3) = 0.0_rp

     !> Rotation Control
     !> "fixed" = Pivot is fixed in space over time (e.g. hinge)
     !> "relativ_sin" = Pivot moves with the oscillation (e.g. center of mass)
     !> In this mode, pivot position follows a sinisudal, based on translational velocity.
     !> "relative = Pivot moves with the oscillation of the object,
     !> but the exact motion of the object is knot know a priori, 
     !> so we use time-integration to update the pivot location.
     !> Used when the motion is not known a priori, like FSI problems.
     character(len=32) :: rotation_center_type = 'fixed'       
     real(kind=rp) :: rot_amp_degree(3) = 0.0_rp
     real(kind=rp) :: rot_freq(3)       = 0.0_rp
     real(kind=rp) :: rot_center(3)     = 0.0_rp
  end type ale_config_t

  type :: ale_kinematics_t
     real(kind=rp) :: vel_trans(3) = 0.0_rp !< Linear Velocity 
     real(kind=rp) :: vel_ang(3)   = 0.0_rp !< Angular Velocity 
     real(kind=rp) :: center(3)    = 0.0_rp !< Current Center of Rotation
  end type ale_kinematics_t

  type, public :: pivot_state_t
      real(kind=rp) :: pos(3)          = 0.0_rp
      real(kind=rp) :: vel_lag(3, 3)   = 0.0_rp
  end type pivot_state_t

contains

  subroutine init_ale_pivot(pivot, config)
    type(pivot_state_t), intent(out) :: pivot
    type(ale_config_t), intent(in) :: config
    pivot%pos = config%rot_center
    pivot%vel_lag = 0.0_rp
  end subroutine init_ale_pivot

  !> Todo: GPU compatibility 
  !> Todo: User motion types (compute_user_kinematics or sth like this)
  !> Todo: compute_fsi_kinematics interface
  subroutine update_ale_mesh_velocity (wx, wy, wz, phi, coef, time, config, pivot, nadv)
    type(field_t), intent(inout) :: wx, wy, wz !< mesh velocities
    type(field_t), intent(in)    :: phi
    type(coef_t),  intent(in)    :: coef
    type(time_state_t), intent(in)    :: time
    type(ale_config_t), intent(in) :: config 
    type(ale_kinematics_t) :: solid_state
    type(pivot_state_t), intent(inout) :: pivot 
    integer, intent(in)  :: nadv

    ! Compute the solid object kinematics based on prescribed motion
    ! No update on the mesh velocities yet! Just compute the kinematics.
    call compute_prescribed_kinematics(solid_state, config, time)

    call update_pivot_location(pivot, solid_state%center, &
                               solid_state%vel_trans, time, nadv, config)

    ! Calculate the mesh velocities based on the solid_state
    call apply_solid_motion(wx, wy, wz, phi, coef, solid_state)

  end subroutine update_ale_mesh_velocity


  !> Calculates kinematics based on Prescribed Sinusoidal Motion
  subroutine compute_prescribed_kinematics(state, config, time)
    type(ale_kinematics_t), intent(out) :: state
    type(ale_config_t), intent(in)      :: config
    type(time_state_t), intent(in)      :: time
    
    real(kind=rp) :: omega_osc(3), omega_rot(3)
    real(kind=rp) :: rad_amp(3), trans_disp(3)
    
    ! --- Oscillation ---
    ! X = A * sin(wt)
    ! V = A * w * cos(wt)
    omega_osc = config%osc_freq * 2.0_rp * pi
    
    state%vel_trans(1) = config%osc_amp(1) * omega_osc(1) * cos(omega_osc(1)*time%t)
    state%vel_trans(2) = config%osc_amp(2) * omega_osc(2) * cos(omega_osc(2)*time%t)
    state%vel_trans(3) = config%osc_amp(3) * omega_osc(3) * cos(omega_osc(3)*time%t)

    ! --- Rotation ---
    ! X_rot = A_rad * sin(wt)
    ! w_rot = A_rad * w * cos(wt)
    ! input is in Degrees, convert to Radians
    rad_amp   = config%rot_amp_degree * (pi / 180.0_rp)
    omega_rot = config%rot_freq * 2.0_rp * pi
    
    state%vel_ang(1) = omega_rot(1) * rad_amp(1) * cos(omega_rot(1)*time%t)
    state%vel_ang(2) = omega_rot(2) * rad_amp(2) * cos(omega_rot(2)*time%t)
    state%vel_ang(3) = omega_rot(3) * rad_amp(3) * cos(omega_rot(3)*time%t)

  end subroutine compute_prescribed_kinematics

  subroutine update_pivot_location(pivot, pivot_loc, vel_trans, time, nadv, config)  
    type(pivot_state_t), intent(inout) :: pivot 
    real(kind=rp), intent(out)     :: pivot_loc(3)
    real(kind=rp), intent(in)      :: vel_trans(3)
    type(time_state_t), intent(in) :: time
    integer, intent(in)            :: nadv
    type(ale_config_t), intent(in) :: config
    real(kind=rp) :: omega_osc(3), trans_disp(3)
    type(ab_time_scheme_t) :: ab_scheme_obj
    real(kind=rp) :: beta(4), dt_history(10)
    integer :: j

   ! --- Center of Rotation (Pivot) ---
   if (trim(config%rotation_center_type) == 'relative') then
      ! Here, we update the pivot point coord based on translational
      ! velocity of the object. For more general use (like FSI in future)
      ! motion of the object may not be explicitly known. So we need this also.
      
      if (time%tstep .gt. 0) then
          call rzero(beta, 4)
          dt_history(1) = time%dt
          dt_history(2) = time%dtlag(1)
          dt_history(3) = time%dtlag(2)
          call ab_scheme_obj%compute_coeffs(beta, dt_history, nadv)
          ! Integrate
          do j = 1, nadv
             pivot%pos = pivot%pos + time%dt * beta(j) * pivot%vel_lag(:, j)
          end do
      end if

      ! Shift History 
      do j = 3, 2, -1
         pivot%vel_lag(:, j) = pivot%vel_lag(:, j-1)
      end do
      pivot%vel_lag(:, 1) = vel_trans
      
      pivot_loc = pivot%pos

   else if (trim(config%rotation_center_type) == 'relative_sin') then
       ! Calculate Displacement of the pivot point: d = A * sin(wt)
       ! Here, we calculate to position of the pivot point (rotation center)
       ! analytically based on the oscillation of the object.
       ! Can be extended later for other motion types, if needed.
       omega_osc = config%osc_freq * 2.0_rp * pi

       trans_disp(1) = config%osc_amp(1) * sin(omega_osc(1)*time%t)
       trans_disp(2) = config%osc_amp(2) * sin(omega_osc(2)*time%t)
       trans_disp(3) = config%osc_amp(3) * sin(omega_osc(3)*time%t)
       
       ! Pivot moves WITH the object translation
       pivot_loc = config%rot_center + trans_disp

   else if (trim(config%rotation_center_type) == 'fixed') then
       ! Pivot is fixed in space
       pivot_loc = config%rot_center

   end if
   
  end subroutine update_pivot_location

  !> Applies Rigid Body Motion formula to the mesh fields
  !> V_mesh = ( V_trans + Omega x (x - center) ) * Phi
  !> Todo: extend to other motion types
  !> Todo: GPU compatibility
  subroutine apply_solid_motion(wx, wy, wz, phi, coef, state)
    type(field_t), intent(inout) :: wx, wy, wz
    type(field_t), intent(in)    :: phi
    type(coef_t),  intent(in)    :: coef
    type(ale_kinematics_t), intent(in) :: state
    
    integer :: i, n
    real(kind=rp) :: dx, dy, dz
    real(kind=rp) :: v_rot_x, v_rot_y, v_rot_z

    n = phi%dof%size()
    
    wx%x = 0.0_rp
    wy%x = 0.0_rp
    wz%x = 0.0_rp

    associate ( x => coef%dof%x, &
                y => coef%dof%y, &
                z => coef%dof%z )

    do concurrent (i = 1:n)
       ! Calculate position vector r relative to the CURRENT center
       dx = x(i,1,1,1) - state%center(1)
       dy = y(i,1,1,1) - state%center(2)
       dz = z(i,1,1,1) - state%center(3)
       
       ! Calculate Cross Product: Omega \cross r
       v_rot_x = state%vel_ang(2)*dz - state%vel_ang(3)*dy
       v_rot_y = state%vel_ang(3)*dx - state%vel_ang(1)*dz
       v_rot_z = state%vel_ang(1)*dy - state%vel_ang(2)*dx

       ! Total Velocity = (Translation + Rotation) * Phi
       wx%x(i,1,1,1) = (state%vel_trans(1) + v_rot_x) * phi%x(i,1,1,1)
       wy%x(i,1,1,1) = (state%vel_trans(2) + v_rot_y) * phi%x(i,1,1,1)
       wz%x(i,1,1,1) = (state%vel_trans(3) + v_rot_z) * phi%x(i,1,1,1)
    end do
    end associate

  end subroutine apply_solid_motion


  !> Save the previous mass matrix history for ALE computations
  subroutine update_ale_mass_history(B, Blag, Blaglag, n)
    integer, intent(in)           :: n       
    real(kind=rp), intent(in)     :: B(n)    
    real(kind=rp), intent(inout)  :: Blag(n) 
    real(kind=rp), intent(inout)  :: Blaglag(n)

    Blaglag = Blag
    Blag = B    
    
  end subroutine update_ale_mass_history

  !> Integrates the mesh position in time (AB Scheme)
  !> Todo: GPU compatibility
  subroutine update_ale_mesh(c_Xh, wm_x, wm_y, wm_z, wm_lag_x, wm_lag_y, wm_lag_z, &
                             time, nadv, scheme_type)
    
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: wm_x, wm_y, wm_z
    type(field_series_t), intent(in) :: wm_lag_x, wm_lag_y, wm_lag_z
    type(time_state_t), intent(in) :: time
    type(ab_time_scheme_t) :: ab_scheme_obj

    integer, intent(in) :: nadv
    character(len=*), intent(in) :: scheme_type
    
    real(kind=rp) :: beta(4)
    real(kind=rp) :: dt_history(10)
    integer :: j, n
    
    if (NEKO_BCKND_DEVICE .eq. 1) call neko_error("ALE Mesh Update not on GPU yet")

    n = c_Xh%dof%size()
    call rzero(beta, 4)

    if (trim(scheme_type) .eq. 'ab') then
       dt_history(1) = time%dt
       dt_history(2) = time%dtlag(1)
       dt_history(3) = time%dtlag(2)
       call ab_scheme_obj%compute_coeffs(beta, dt_history, nadv)       
    else
       call neko_error("ALE: Unknown mesh time-integration scheme '"//trim(scheme_type)//"'")
    end if

    ! Update Mesh Position
    ! Current Time Step
    c_Xh%dof%x = c_Xh%dof%x + time%dt * beta(1) * wm_x%x
    c_Xh%dof%y = c_Xh%dof%y + time%dt * beta(1) * wm_y%x
    c_Xh%dof%z = c_Xh%dof%z + time%dt * beta(1) * wm_z%x

    ! Lagged Time Steps
    do j = 2, nadv       
       c_Xh%dof%x = c_Xh%dof%x + time%dt * beta(j) * wm_lag_x%lf(j-1)%x
       c_Xh%dof%y = c_Xh%dof%y + time%dt * beta(j) * wm_lag_y%lf(j-1)%x
       c_Xh%dof%z = c_Xh%dof%z + time%dt * beta(j) * wm_lag_z%lf(j-1)%x
    end do

  end subroutine update_ale_mesh


end module ale_motion