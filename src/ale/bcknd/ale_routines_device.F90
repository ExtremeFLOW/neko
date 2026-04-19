! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!

module ale_routines_device
  use num_types, only : rp, c_rp
  use field, only : field_t
  use coefs, only : coef_t
  use field_series, only : field_series_t
  use time_state, only : time_state_t
  use ab_time_scheme, only : ab_time_scheme_t
  use mesh, only : mesh_t
  use utils, only : neko_error
  use device_math, only : device_add2s2
  use zero_dirichlet, only : zero_dirichlet_t
  use math, only : rzero, glimax, cfill
  use gather_scatter, only : GS_OP_MIN
  use ale_rigid_kinematics, only : ale_config_t, body_kinematics_t
  use device, only : device_map, device_memcpy, device_unmap, &
       HOST_TO_DEVICE, DEVICE_TO_HOST
  use comm, only : NEKO_COMM
  use mpi_f08, only : MPI_WTIME, MPI_Barrier
  use logger, only : neko_log
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int, C_NULL_PTR

  implicit none
  private

  public :: compute_stiffness_ale_device
  public :: add_kinematics_to_mesh_velocity_device
  public :: update_ale_mesh_device

  type, bind(c) :: kinematics_params_t
    real(c_rp) :: cx, cy, cz
    real(c_rp) :: vtx, vty, vtz
    real(c_rp) :: vax, vay, vaz
    real(c_rp) :: px, py, pz
    real(c_rp) :: r11, r12, r13
    real(c_rp) :: r21, r22, r23
    real(c_rp) :: r31, r32, r33
  end type kinematics_params_t

#ifdef HAVE_HIP
  interface
    subroutine add_kinematics_to_mesh_velocity_hip(wx, wy, wz, &
         x_ref, y_ref, z_ref, phi, x, y, z, &
         kin_params, n) bind(c, name="add_kinematics_to_mesh_velocity_hip")
      use, intrinsic :: iso_c_binding
      import :: kinematics_params_t
      type(c_ptr), value :: wx, wy, wz, x_ref, y_ref, z_ref, phi, x, y, z
      type(kinematics_params_t), value :: kin_params
      integer(c_int), value :: n
    end subroutine add_kinematics_to_mesh_velocity_hip

    subroutine compute_cheap_dist_hip(d_d, x_d, y_d, z_d, lx, ly, lz, nel, &
         local_iters, nchange_d) bind(c, name="compute_cheap_dist_hip")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: d_d, x_d, y_d, z_d, nchange_d
      integer(c_int), value :: lx, ly, lz, nel, local_iters
    end subroutine compute_cheap_dist_hip
  end interface
#elif HAVE_CUDA
  interface
    subroutine add_kinematics_to_mesh_velocity_cuda(wx, wy, wz, &
         x_ref, y_ref, z_ref, phi, x, y, z, &
         kin_params, n) bind(c, name="add_kinematics_to_mesh_velocity_cuda")
      use, intrinsic :: iso_c_binding
      import :: kinematics_params_t
      type(c_ptr), value :: wx, wy, wz, x_ref, y_ref, z_ref, phi, x, y, z
      type(kinematics_params_t), value :: kin_params
      integer(c_int), value :: n
    end subroutine add_kinematics_to_mesh_velocity_cuda

    subroutine compute_cheap_dist_cuda(d_d, x_d, y_d, z_d, lx, ly, lz, nel, &
         local_iters, nchange_d) bind(c, name="compute_cheap_dist_cuda")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: d_d, x_d, y_d, z_d, nchange_d
      integer(c_int), value :: lx, ly, lz, nel, local_iters
    end subroutine compute_cheap_dist_cuda
  end interface

#endif

contains

  !> Compute mesh stiffness with per-body gain/decay from stiff_geom
  subroutine compute_stiffness_ale_device(coef, params)
    type(coef_t), intent(inout) :: coef
    type(ale_config_t), intent(in) :: params
    integer :: i, n, b, ierr
    integer, allocatable :: cheap_map(:)
    integer :: n_cheap, map_idx
    real(kind=rp) :: x, y, z
    real(kind=rp) :: raw_dist, body_stiff_val, max_added_stiff
    real(kind=rp) :: cx, cy, cz
    real(kind=rp) :: arg, decay, gain, norm_dist
    real(kind=rp) :: sample_start_time, sample_end_time, sample_time
    real(kind=rp), allocatable :: dist_fields(:,:)
    character(len=128) :: log_buf

    n = coef%dof%size()

    ! Check how many bodies need cheap_dist and create Map
    allocate(cheap_map(params%nbodies))
    cheap_map = 0
    n_cheap = 0

    do b = 1, params%nbodies
       if (trim(params%bodies(b)%stiff_geom%type) == 'cheap_dist') then
          n_cheap = n_cheap + 1
          cheap_map(b) = n_cheap
       end if
    end do

    ! Allocate and Compute cheap_dist only for required bodies
    if (n_cheap > 0) then
       allocate(dist_fields(n, n_cheap))
       dist_fields = huge(0.0_rp)

       do b = 1, params%nbodies
          map_idx = cheap_map(b)
          if (map_idx > 0) then
             call neko_log%message(' ')
             call neko_log%message(" Start: cheap dist calculation " // &
                  "for body '" // trim(params%bodies(b)%name) // "'")
             call MPI_Barrier(NEKO_COMM, ierr)
             sample_start_time = MPI_WTIME()

             ! Compute into the specific slot for this body

             call compute_cheap_dist_device(dist_fields(:, map_idx), coef, &
                  coef%msh, params%bodies(b)%zone_indices)

             call MPI_Barrier(NEKO_COMM, ierr)
             sample_end_time = MPI_WTIME()
             sample_time = sample_end_time - sample_start_time

             write(log_buf, '(A, A, A, ES11.4, A)') "   cheap dist for '", &
                  trim(params%bodies(b)%name), "' took ", sample_time, " (s)"
             call neko_log%message(log_buf)
          end if
       end do
    end if
    call neko_log%message(' ')

    ! Build stiffness field
    select case (trim(params%stiffness_type))
    case ('built-in')

       do concurrent (i = 1:n)
          x = coef%dof%x(i, 1, 1, 1)
          y = coef%dof%y(i, 1, 1, 1)
          z = coef%dof%z(i, 1, 1, 1)

          max_added_stiff = 0.0_rp

          ! Loop over bodies, calculate local contribution
          do b = 1, params%nbodies
             gain = params%bodies(b)%stiff_geom%gain
             if (trim(params%bodies(b)%stiff_geom%type) == 'cheap_dist') then
                decay = params%bodies(b)%stiff_geom%stiff_dist
             else
                decay = params%bodies(b)%stiff_geom%radius
             end if

             ! Geometry Center
             cx = params%bodies(b)%stiff_geom%center(1)
             cy = params%bodies(b)%stiff_geom%center(2)
             cz = params%bodies(b)%stiff_geom%center(3)

             raw_dist = huge(0.0_rp)

             ! Calculate Distance
             select case (trim(params%bodies(b)%stiff_geom%type))
             case ('sphere')

                raw_dist = sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)

             case ('cylinder')

                ! Distance to Z-axis centered at (cx, cy)
                raw_dist = sqrt((x - cx)**2 + (y - cy)**2)

             case ('box')

                ! ToDO

             case ('cheap_dist')

                map_idx = cheap_map(b)
                if (map_idx > 0) then
                   raw_dist = dist_fields(i, map_idx)
                end if

             end select

             ! Apply Profile
             body_stiff_val = 0.0_rp
             select case (trim(params%bodies(b)%stiff_geom%decay_profile))
             case ('gaussian')

                ! exp( -(r/decay)^2 )
                arg = -(raw_dist**2) / (decay**2)
                arg = arg * params%bodies(b)%stiff_geom%cutoff_coef
                body_stiff_val = gain * exp(arg)

             case ('tanh')

                ! Tanh profile
                norm_dist = (raw_dist / decay)
                norm_dist = norm_dist * &
                     params%bodies(b)%stiff_geom%cutoff_coef
                body_stiff_val = gain * (1.0_rp - tanh(norm_dist))

             end select

             if (body_stiff_val > max_added_stiff) then
                max_added_stiff = body_stiff_val
             end if

          end do

          coef%h1(i, 1, 1, 1) = 1.0_rp + max_added_stiff
          coef%h2(i, 1, 1, 1) = 0.0_rp
       end do

    case default
       call neko_error("ALE Manager: Unknown stiffness type")
    end select

    coef%ifh2 = .false.
    if (allocated(dist_fields)) deallocate(dist_fields)
    if (allocated(cheap_map)) deallocate(cheap_map)
  end subroutine compute_stiffness_ale_device

!> Cheap dist device implementation
  subroutine compute_cheap_dist_device(d, coef, msh, zone_indices)
    real(kind=rp), intent(inout), target :: d(:)
    type(coef_t), intent(in) :: coef
    type(mesh_t), intent(in) :: msh
    type(zero_dirichlet_t) :: bc_wall
    integer, intent(in) :: zone_indices(:)
    integer :: i, k, n, m, idx
    integer :: ipass, max_pass, local_iters
    integer :: lx, ly, lz, nel, z_idx
    integer, target :: change_vec(1) 
    logical :: done
    character(len=128) :: log_buf
    type(c_ptr) :: d_d, nchange_d

    d_d = C_NULL_PTR
    nchange_d = C_NULL_PTR

    lx = coef%dof%Xh%lx
    ly = coef%dof%Xh%ly
    lz = coef%dof%Xh%lz
    nel = msh%nelv
    n = size(d)
    max_pass = 10000

    ! Limit for worst case scenario such that all nodes can propagate
    ! their values across the element before triggering an MPI call.
    local_iters = lx + ly + lz

    ! Initialize array on host
    call cfill(d, huge(0.0_rp), n)

    if (size(zone_indices) > 0) then
       call bc_wall%init_from_components(coef)
       do k = 1, size(zone_indices)
          z_idx = zone_indices(k)
          call bc_wall%mark_zone(msh%labeled_zones(z_idx))
       end do
       call bc_wall%finalize()
       m = bc_wall%msk(0)
       do i = 1, m
          idx = bc_wall%msk(i)
          d(idx) = 0.0_rp
       end do
       call bc_wall%free()
    end if

    call device_map(d, d_d, n)
    call device_map(change_vec, nchange_d, 1) 
    call device_memcpy(d, d_d, n, HOST_TO_DEVICE, .true.)
    ipass = 1
    done = .false.
    do while (ipass <= max_pass .and. .not. done)
       
       change_vec(1) = 0
       call device_memcpy(change_vec, nchange_d, 1, HOST_TO_DEVICE, .true.)

#ifdef HAVE_HIP
       call compute_cheap_dist_hip(d_d, coef%dof%x_d, coef%dof%y_d, coef%dof%z_d, &
            lx, ly, lz, nel, local_iters, nchange_d)
#elif HAVE_CUDA
       call compute_cheap_dist_cuda(d_d, coef%dof%x_d, coef%dof%y_d, coef%dof%z_d, &
            lx, ly, lz, nel, local_iters, nchange_d)
#else
      call neko_error("ALE: compute_cheap_dist_device supports only HIP or CUDA backends currently")
#endif

       call device_memcpy(change_vec, nchange_d, 1, DEVICE_TO_HOST, .true.)
       call coef%gs_h%gs_op_vector(d, n, GS_OP_MIN)

       if (glimax(change_vec, 1) == 0) done = .true.
       ipass = ipass + 1
    end do

    call device_memcpy(d, d_d, n, DEVICE_TO_HOST, .true.)

    call device_unmap(change_vec, nchange_d)
    call device_unmap(d, d_d)
         
    write(log_buf, '(A, I0, A)') "   converged in: ", ipass, " passes"
    call neko_log%message(log_buf)
  end subroutine compute_cheap_dist_device


  !> Add Kinematics to Mesh Velocity
  subroutine add_kinematics_to_mesh_velocity_device(wx, wy, wz, &
        x_ref, y_ref, z_ref, phi, coef, kinematics, rot_mat, inital_pivot_loc)
    type(field_t), intent(inout) :: wx, wy, wz
    type(field_t), intent(in) :: x_ref, y_ref, z_ref
    type(field_t), intent(in) :: phi
    type(coef_t), intent(in) :: coef
    type(body_kinematics_t), intent(in), target :: kinematics
    real(kind=rp), intent(in), target :: inital_pivot_loc(3)
    real(kind=rp), intent(in), target :: rot_mat(3,3)
    
    integer(c_int), target :: n
    type(kinematics_params_t) :: kin_params

    n = phi%dof%size()

    kin_params%cx = kinematics%center(1)
    kin_params%cy = kinematics%center(2)
    kin_params%cz = kinematics%center(3)
    kin_params%vtx = kinematics%vel_trans(1)
    kin_params%vty = kinematics%vel_trans(2)
    kin_params%vtz = kinematics%vel_trans(3)
    kin_params%vax = kinematics%vel_ang(1)
    kin_params%vay = kinematics%vel_ang(2)
    kin_params%vaz = kinematics%vel_ang(3)
    kin_params%px = inital_pivot_loc(1)
    kin_params%py = inital_pivot_loc(2)
    kin_params%pz = inital_pivot_loc(3)
    kin_params%r11 = rot_mat(1,1)
    kin_params%r12 = rot_mat(1,2)
    kin_params%r13 = rot_mat(1,3)
    kin_params%r21 = rot_mat(2,1)
    kin_params%r22 = rot_mat(2,2)
    kin_params%r23 = rot_mat(2,3)
    kin_params%r31 = rot_mat(3,1)
    kin_params%r32 = rot_mat(3,2)
    kin_params%r33 = rot_mat(3,3)

#ifdef HAVE_HIP
    call add_kinematics_to_mesh_velocity_hip( &
         wx%x_d, wy%x_d, wz%x_d, &
         x_ref%x_d, y_ref%x_d, z_ref%x_d, &
         phi%x_d, coef%dof%x_d, coef%dof%y_d, coef%dof%z_d, & 
         kin_params, n)
#elif HAVE_CUDA
    call add_kinematics_to_mesh_velocity_cuda( &
         wx%x_d, wy%x_d, wz%x_d, &
         x_ref%x_d, y_ref%x_d, z_ref%x_d, &
         phi%x_d, coef%dof%x_d, coef%dof%y_d, coef%dof%z_d, & 
         kin_params, n)
#else
    call neko_error("ALE: add_kinematics_to_mesh_velocity_device " // &
         "supports only HIP or CUDA backends")
#endif

  end subroutine add_kinematics_to_mesh_velocity_device


  !> Update ALE Mesh
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

    ! Current timestep update
    factor = time%dt * ab_coeffs(1)
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
