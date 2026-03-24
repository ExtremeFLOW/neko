! Copyright (c) 2025-2026 The Neko Authors
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
! Routines for expensive ALE calculations on CPU
module ale_routines_cpu
  use num_types, only : rp
  use field, only : field_t
  use coefs, only : coef_t
  use math, only : add2s2, cfill, glimax, rzero
  use field_series, only : field_series_t
  use time_state, only : time_state_t
  use ab_time_scheme, only : ab_time_scheme_t
  use utils, only : neko_error
  use comm, only : NEKO_COMM
  use zero_dirichlet, only : zero_dirichlet_t
  use mesh, only : mesh_t
  use gather_scatter, only : GS_OP_MIN
  use mpi_f08, only : MPI_WTIME, MPI_Barrier
  use logger, only : neko_log
  use ale_rigid_kinematics, only : ale_config_t, body_kinematics_t
  implicit none
  private

  public :: compute_stiffness_ale_cpu
  public :: add_kinematics_to_mesh_velocity_cpu
  public :: update_ale_mesh_cpu

contains

  !> Compute mesh stiffness with per-body gain/decay from stiff_geom
  subroutine compute_stiffness_ale_cpu(coef, params)
    type(coef_t), intent(inout) :: coef
    type(ale_config_t), intent(in) :: params
    integer :: i, n, b, ierr
    integer, allocatable :: cheap_map(:)
    integer :: n_cheap, map_idx
    real(kind=rp) :: x, y, z
    real(kind=rp) :: raw_dist, body_stiff_val, max_added_stiff
    real(kind=rp) :: cx, cy, cz, dx, dy, dz
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

             ! Nek5000 algorithm
             ! call compute_cheap_dist_cpu(dist_fields(:, map_idx), coef, &
             !     coef%msh, params%bodies(b)%zone_indices)

             call compute_cheap_dist_v2_cpu(dist_fields(:, map_idx), coef, &
                  coef%msh, params%bodies(b)%zone_indices)

             call MPI_Barrier(NEKO_COMM, ierr)
             sample_end_time = MPI_WTIME()
             sample_time = sample_end_time - sample_start_time

             write(log_buf, '(A, A, A, F10.4, A)') "   cheap dist for '", &
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
  end subroutine compute_stiffness_ale_cpu

  !> Implementation of cheap_dist in Nek5000 (CPU)
  subroutine compute_cheap_dist_cpu(d, coef, msh, zone_indices)
    real(kind=rp), intent(inout), target :: d(:)
    type(coef_t), intent(in) :: coef
    type(mesh_t), intent(in) :: msh
    type(zero_dirichlet_t) :: bc_wall
    integer, intent(in) :: zone_indices(:)
    real(kind=rp), pointer :: d4(:, :, :, :)
    integer :: i, j, k, e, n
    integer :: ipass, nchange, max_pass
    integer :: ii, jj, kk, i0, i1, j0, j1, k0, k1
    integer :: lx, ly, lz, nel, z_idx
    integer :: m, idx
    real(kind=rp) :: dtmp, x1, y1, z1, x2, y2, z2
    integer :: change_vec(1)
    logical :: done

    lx = coef%dof%Xh%lx
    ly = coef%dof%Xh%ly
    lz = coef%dof%Xh%lz
    nel = msh%nelv
    n = size(d)
    d4(1:lx, 1:ly, 1:lz, 1:nel) => d
    max_pass = 10000

    !d = huge(0.0_rp)
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

    ipass = 1
    done = .false.
    do while (ipass <= max_pass .and. .not. done)
       nchange = 0
       do e = 1, nel
          do k = 1, lz
             do j = 1, ly
                do i = 1, lx
                   x1 = coef%dof%x(i, j, k, e)
                   y1 = coef%dof%y(i, j, k, e)
                   z1 = coef%dof%z(i, j, k, e)
                   i0 = max(1, i - 1)
                   i1 = min(lx, i + 1)
                   j0 = max(1, j - 1)
                   j1 = min(ly, j + 1)
                   k0 = max(1, k - 1)
                   k1 = min(lz, k + 1)
                   do kk = k0, k1
                      do jj = j0, j1
                         do ii = i0, i1
                            if (ii == i .and. jj == j .and. kk == k) cycle
                            x2 = coef%dof%x(ii, jj, kk, e)
                            y2 = coef%dof%y(ii, jj, kk, e)
                            z2 = coef%dof%z(ii, jj, kk, e)
                            dtmp = d4(ii, jj, kk, e) + &
                                 sqrt((x1 - x2)**2 + &
                                 (y1 - y2)**2 + (z1 - z2)**2)
                            if (dtmp < d4(i, j, k, e)) then
                               d4(i, j, k, e) = dtmp
                               nchange = nchange + 1
                            end if
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
       call coef%gs_h%gs_op_vector(d, n, GS_OP_MIN)
       change_vec(1) = nchange
       if (glimax(change_vec, 1) == 0) done = .true.
       ipass = ipass + 1
    end do
  end subroutine compute_cheap_dist_cpu

  !> Compute cheap_dist field by passing distance information
  !> throughout an entire local element
  !> before doing a global MPI synchronization.
  subroutine compute_cheap_dist_v2_cpu(d, coef, msh, zone_indices)
    real(kind=rp), intent(inout), target :: d(:)
    type(coef_t), intent(in) :: coef
    type(mesh_t), intent(in) :: msh
    type(zero_dirichlet_t) :: bc_wall
    integer, intent(in) :: zone_indices(:)
    real(kind=rp), pointer :: d4(:, :, :, :)
    integer :: i, j, k, e, n
    integer :: ipass, nchange, max_pass
    integer :: ii, jj, kk, i0, i1, j0, j1, k0, k1
    integer :: lx, ly, lz, nel, z_idx
    integer :: m, idx, iter, local_iters
    real(kind=rp) :: dtmp, x1, y1, z1, x2, y2, z2
    integer :: change_vec(1)
    logical :: done, changed_local, element_changed_ever
    character(len=128) :: log_buf

    lx = coef%dof%Xh%lx
    ly = coef%dof%Xh%ly
    lz = coef%dof%Xh%lz
    nel = msh%nelv
    n = size(d)
    d4(1:lx, 1:ly, 1:lz, 1:nel) => d
    max_pass = 10000

!   Limit for worst case scenario such that all nodes can propagate
!   their values across the element before triggering an MPI call.
    local_iters = lx + ly + lz

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

    ipass = 1
    done = .false.
    do while (ipass <= max_pass .and. .not. done)
       nchange = 0

       do e = 1, nel
          iter = 1
          element_changed_ever = .false.
          changed_local = .true.
          do while (changed_local .and. iter <= local_iters)

             changed_local = .false.
             do k = 1, lz
                do j = 1, ly
                   do i = 1, lx
                      x1 = coef%dof%x(i, j, k, e)
                      y1 = coef%dof%y(i, j, k, e)
                      z1 = coef%dof%z(i, j, k, e)
                      i0 = max(1, i - 1)
                      i1 = min(lx, i + 1)
                      j0 = max(1, j - 1)
                      j1 = min(ly, j + 1)
                      k0 = max(1, k - 1)
                      k1 = min(lz, k + 1)
                      do kk = k0, k1
                         do jj = j0, j1
                            do ii = i0, i1
                               if (ii == i .and. jj == j .and. kk == k) cycle

                               x2 = coef%dof%x(ii, jj, kk, e)
                               y2 = coef%dof%y(ii, jj, kk, e)
                               z2 = coef%dof%z(ii, jj, kk, e)

                               dtmp = d4(ii, jj, kk, e) + &
                               sqrt((x1 - x2)**2 + &
                               (y1 - y2)**2 + (z1 - z2)**2)

                               if (dtmp < d4(i, j, k, e)) then
                                  d4(i, j, k, e) = dtmp
                                  changed_local = .true.
                               end if

                            end do
                         end do
                      end do
                   end do
                end do
             end do
             if (changed_local) element_changed_ever = .true.
             iter = iter + 1
          end do

          if (element_changed_ever) nchange = nchange + 1
       end do

       call coef%gs_h%gs_op_vector(d, n, GS_OP_MIN)
       change_vec(1) = nchange

       if (glimax(change_vec, 1) == 0) done = .true.
       ipass = ipass + 1
    end do

    write(log_buf, '(A, I0, A)') "   converged in: ", ipass, " passes"
    call neko_log%message(log_buf)
  end subroutine compute_cheap_dist_v2_cpu


  !> Adds kinematics to mesh velocity (CPU)
  subroutine add_kinematics_to_mesh_velocity_cpu(wx, wy, wz, &
       x_ref, y_ref, z_ref, phi, coef, kinematics, rot_mat, inital_pivot_loc)
    type(field_t), intent(inout) :: wx, wy, wz
    type(field_t), intent(in) :: x_ref, y_ref, z_ref
    type(field_t), intent(in) :: phi
    type(coef_t), intent(in) :: coef
    type(body_kinematics_t), intent(in) :: kinematics
    real(kind=rp), intent(in) :: inital_pivot_loc(3)
    real(kind=rp), intent(in) :: rot_mat(3,3)
    integer :: i, n
    real(kind=rp) :: rx, ry, rz
    real(kind=rp) :: v_tan_x, v_tan_y, v_tan_z
    real(kind=rp) :: dx_ref, dy_ref, dz_ref
    real(kind=rp) :: rx_target, ry_target, rz_target
    real(kind=rp) :: p_val
    n = phi%dof%size()
    associate (x => coef%dof%x, y => coef%dof%y, z => coef%dof%z)
      do concurrent (i = 1:n)
         ! for points on the wall (phi=1) we do this to avoid any
         ! numeric error due to computation of rotation matrix.
         ! It ensures, the walls are always where they need to be!
         if ( abs(phi%x(i, 1, 1, 1) - 1.0_rp) < 1e-6_rp ) then

            rx = x(i, 1, 1, 1) - kinematics%center(1)
            ry = y(i, 1, 1, 1) - kinematics%center(2)
            rz = z(i, 1, 1, 1) - kinematics%center(3)
            v_tan_x = kinematics%vel_ang(2) * rz - kinematics%vel_ang(3) * ry
            v_tan_y = kinematics%vel_ang(3) * rx - kinematics%vel_ang(1) * rz
            v_tan_z = kinematics%vel_ang(1) * ry - kinematics%vel_ang(2) * rx
            wx%x(i, 1, 1, 1) = wx%x(i, 1, 1, 1) + &
                  (kinematics%vel_trans(1) + v_tan_x) * phi%x(i, 1, 1, 1)
            wy%x(i, 1, 1, 1) = wy%x(i, 1, 1, 1) + &
                  (kinematics%vel_trans(2) + v_tan_y) * phi%x(i, 1, 1, 1)
            wz%x(i, 1, 1, 1) = wz%x(i, 1, 1, 1) + &
                  (kinematics%vel_trans(3) + v_tan_z) * phi%x(i, 1, 1, 1)
         else
            ! For other points we do this to avoid the time-dependnent
            ! drift in some special cases, which happens due to the nature of
            ! the ODE that we integrate above!
            dx_ref = x_ref%x(i, 1, 1, 1) - inital_pivot_loc(1)
            dy_ref = y_ref%x(i, 1, 1, 1) - inital_pivot_loc(2)
            dz_ref = z_ref%x(i, 1, 1, 1) - inital_pivot_loc(3)

            ! Rotate to find the "ghost" target vector
            ! Apply the rotation matrix R(t) to the reference vector
            rx_target = rot_mat(1,1)*dx_ref + rot_mat(1,2)*dy_ref + &
                 rot_mat(1,3)*dz_ref
            ry_target = rot_mat(2,1)*dx_ref + rot_mat(2,2)*dy_ref + &
                 rot_mat(2,3)*dz_ref
            rz_target = rot_mat(3,1)*dx_ref + rot_mat(3,2)*dy_ref + &
                 rot_mat(3,3)*dz_ref

            ! Calculate tangential velocity at the ghost target
            ! v_tan = Omega \corss R_target
            v_tan_x = kinematics%vel_ang(2) * rz_target - &
                 kinematics%vel_ang(3) * ry_target
            v_tan_y = kinematics%vel_ang(3) * rx_target - &
                 kinematics%vel_ang(1) * rz_target
            v_tan_z = kinematics%vel_ang(1) * ry_target - &
                 kinematics%vel_ang(2) * rx_target

            p_val = phi%x(i, 1, 1, 1)

            ! Total Mesh Velocity
            wx%x(i, 1, 1, 1) = wx%x(i, 1, 1, 1) + &
                 (kinematics%vel_trans(1) + v_tan_x) * p_val
            wy%x(i, 1, 1, 1) = wy%x(i, 1, 1, 1) + &
                 (kinematics%vel_trans(2) + v_tan_y) * p_val
            wz%x(i, 1, 1, 1) = wz%x(i, 1, 1, 1) + &
                 (kinematics%vel_trans(3) + v_tan_z) * p_val

         end if
      end do
    end associate
  end subroutine add_kinematics_to_mesh_velocity_cpu

  !> Updates mesh position by integrating mesh velocity in time using AB (CPU)
  subroutine update_ale_mesh_cpu(c_Xh, wm_x, wm_y, wm_z, wm_x_lag, wm_y_lag, &
       wm_z_lag, time, nadv, scheme_type)
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: wm_x, wm_y, wm_z
    type(field_series_t), intent(in) :: wm_x_lag, wm_y_lag, wm_z_lag
    type(time_state_t), intent(in) :: time
    type(ab_time_scheme_t) :: ab_scheme_obj
    integer, intent(in) :: nadv
    integer :: i, j, n
    character(len=*), intent(in) :: scheme_type
    real(kind=rp) :: ab_coeffs(4), dt_history(10)

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

    do concurrent (i = 1:n)
       c_Xh%dof%x(i, 1, 1, 1) = c_Xh%dof%x(i, 1, 1, 1) + &
            time%dt * ab_coeffs(1) * wm_x%x(i, 1, 1, 1)
       c_Xh%dof%y(i, 1, 1, 1) = c_Xh%dof%y(i, 1, 1, 1) + &
            time%dt * ab_coeffs(1) * wm_y%x(i, 1, 1, 1)
       c_Xh%dof%z(i, 1, 1, 1) = c_Xh%dof%z(i, 1, 1, 1) + &
            time%dt * ab_coeffs(1) * wm_z%x(i, 1, 1, 1)

       do j = 2, nadv
          c_Xh%dof%x(i, 1, 1, 1) = c_Xh%dof%x(i, 1, 1, 1) + &
               time%dt * ab_coeffs(j) * wm_x_lag%lf(j - 1)%x(i, 1, 1, 1)
          c_Xh%dof%y(i, 1, 1, 1) = c_Xh%dof%y(i, 1, 1, 1) + &
               time%dt * ab_coeffs(j) * wm_y_lag%lf(j - 1)%x(i, 1, 1, 1)
          c_Xh%dof%z(i, 1, 1, 1) = c_Xh%dof%z(i, 1, 1, 1) + &
               time%dt * ab_coeffs(j) * wm_z_lag%lf(j - 1)%x(i, 1, 1, 1)
       end do

    end do
  end subroutine update_ale_mesh_cpu

end module ale_routines_cpu
