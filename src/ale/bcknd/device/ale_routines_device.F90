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
  use logger, only : neko_log, LOG_SIZE
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int, C_NULL_PTR

  implicit none
  private

  public :: add_kinematics_to_mesh_velocity_device
  public :: update_ale_mesh_device
  public :: compute_cheap_dist_device

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

!> Cheap dist device implementation
  subroutine compute_cheap_dist_device(dist_field, coef, msh, zone_indices, copy_to_host)
    type(field_t), intent(inout) :: dist_field
    type(coef_t), intent(in) :: coef
    type(mesh_t), intent(in) :: msh
    type(zero_dirichlet_t) :: bc_wall
    integer, intent(in) :: zone_indices(:)
    logical, intent(in) :: copy_to_host
    integer :: i, k, n, m, idx
    integer :: ipass, max_pass, local_iters
    integer :: lx, ly, lz, nel, z_idx
    integer, target :: change_vec(1)
    logical :: done
    character(len=LOG_SIZE) :: log_buf
    type(c_ptr) :: nchange_d

    nchange_d = C_NULL_PTR
    lx = coef%dof%Xh%lx
    ly = coef%dof%Xh%ly
    lz = coef%dof%Xh%lz
    nel = msh%nelv
    n = coef%dof%size()
    max_pass = 10000

    ! Limit for worst case scenario such that all nodes can propagate
    ! their values across the element before triggering an MPI call.
    local_iters = lx + ly + lz

    call cfill(dist_field%x, huge(0.0_rp), n)

    if (size(zone_indices) .gt. 0) then
       call bc_wall%init_from_components(coef)
       do k = 1, size(zone_indices)
          z_idx = zone_indices(k)
          call bc_wall%mark_zone(msh%labeled_zones(z_idx))
       end do
       call bc_wall%finalize()
       m = bc_wall%msk(0)
       do i = 1, m
          idx = bc_wall%msk(i)
          dist_field%x(idx, 1, 1, 1) = 0.0_rp
       end do
       call bc_wall%free()
    end if

    call device_map(change_vec, nchange_d, 1)
    call dist_field%copy_from(HOST_TO_DEVICE, sync = .true.)

    ipass = 1
    done = .false.

    do while ((ipass .le. max_pass) .and. .not. done)

       change_vec(1) = 0
       call device_memcpy(change_vec, nchange_d, 1, HOST_TO_DEVICE, .true.)

#ifdef HAVE_HIP
       call compute_cheap_dist_hip(dist_field%x_d, &
            coef%dof%x_d, coef%dof%y_d, coef%dof%z_d, &
            lx, ly, lz, nel, local_iters, nchange_d)
#elif HAVE_CUDA
       call compute_cheap_dist_cuda(dist_field%x_d, &
            coef%dof%x_d, coef%dof%y_d, coef%dof%z_d, &
            lx, ly, lz, nel, local_iters, nchange_d)
#endif

       call device_memcpy(change_vec, nchange_d, 1, DEVICE_TO_HOST, .true.)

       call coef%gs_h%gs_op_vector(dist_field%x, n, GS_OP_MIN)

       if (glimax(change_vec, 1) .eq. 0) done = .true.
       ipass = ipass + 1
    end do

    call device_unmap(change_vec, nchange_d)

    if (copy_to_host) then
       call dist_field%copy_from(DEVICE_TO_HOST, sync = .true.)
    end if

    if (done) then
       write(log_buf, '(A, I0, A)') "    converged in: ", ipass, " passes"
       call neko_log%message(log_buf)
    else
       write(log_buf, '(A, I0, A)') "    reached max passes: ", ipass, &
            " without convergence"
       call neko_log%message(log_buf)
    end if
  end subroutine compute_cheap_dist_device


  !> Add Kinematics to Mesh Velocity
  subroutine add_kinematics_to_mesh_velocity_device(wx, wy, wz, &
        x_ref, y_ref, z_ref, phi, coef, kinematics, rot_mat, inital_pivot_loc)
    type(field_t), intent(inout) :: wx, wy, wz
    type(field_t), intent(in) :: x_ref, y_ref, z_ref
    type(field_t), intent(in) :: phi
    type(coef_t), intent(in) :: coef
    type(body_kinematics_t), intent(in) :: kinematics
    real(kind=rp), intent(in) :: inital_pivot_loc(3)
    real(kind=rp), intent(in) :: rot_mat(3,3)
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
