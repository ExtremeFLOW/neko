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
!> Moving-mesh (ALE) support for the probes simulation component.
!
module probes_moving_mesh
  use num_types, only : rp, dp
  use vector, only : vector_t
  use global_interpolation, only : global_interpolation_t
  use dofmap, only : dofmap_t
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use logger, only : neko_log
  use utils, only : neko_error
  use time_state, only : time_state_t
  use ale_manager, only : neko_ale
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_MAX, MPI_SUM, &
       MPI_INTEGER, MPI_Wtime
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_memcpy, DEVICE_TO_HOST
  use math, only : copy
  implicit none
  private

  character(len=*), parameter :: PROBES_MM_FIXED = 'fixed'
  character(len=*), parameter :: PROBES_MM_BODY = 'body_attached'

  !> Keeps the probe interpolation mapping valid on a moving mesh.
  type, public :: probes_moving_mesh_t
     !> Whether moving-mesh handling is active.
     logical :: enabled = .false.
     !> Whether the mesh moves at all.
     logical :: mesh_has_moved = .false.
     !> Probes move rigidly with a body (otherwise lab-frame fixed).
     logical :: body_attached = .false.
     !> Current probe coordinates are written alongside the data.
     logical :: output_coords = .false.
     !> Boundary zone id linking the ALE body.
     integer :: zone_id = -1
     !> Index of the linked body in the ALE manager arrays.
     integer :: body_idx = -1
     !> Name of the linked body.
     character(len=32) :: body_name = 'NOT_LINKED'
     !> Initial rotation center of the linked body.
     real(kind=rp) :: pivot0(3) = 0.0_rp
     !> Number of local probes (fixed after probes init).
     integer :: n_points = 0
     !> Whether set_points has been called.
     logical :: points_set = .false.
     !> Reference coordinates (3, n_points), in the initial frame.
     real(kind=rp), allocatable :: xyz0(:,:)
     !> Current target coordinates (3, n_points).
     real(kind=rp), allocatable :: xyz_target(:,:)
     !> Work vectors for the validity check.
     type(vector_t) :: chk_x, chk_y, chk_z
     !> Newton tolerance and padding to re-initialize the
     !! interpolation with.
     real(kind=dp) :: tolerance
     real(kind=dp) :: padding
     !> Max accepted distance between mapped and target position
     !! before a re-find is triggered.
     real(kind=rp) :: check_tolerance = 1e-8_rp
     !> Global number of points currently not found in the mesh.
     integer :: n_lost = 0
   contains
     procedure, pass(this) :: init => probes_moving_mesh_init
     procedure, pass(this) :: set_points => probes_moving_mesh_set_points
     procedure, pass(this) :: free => probes_moving_mesh_free
     procedure, pass(this) :: update => probes_moving_mesh_update
     procedure, private, pass(this) :: compute_targets => &
          probes_moving_mesh_compute_targets
     procedure, private, pass(this) :: mapping_error => &
          probes_moving_mesh_mapping_error
     procedure, private, pass(this) :: refind => probes_moving_mesh_refind
  end type probes_moving_mesh_t

contains

  !> Parse the configuration and link the ALE body.
  !! @param json The `moving_mesh` JSON subdict.
  subroutine probes_moving_mesh_init(this, json)
    class(probes_moving_mesh_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    character(len=:), allocatable :: mode
    character(len=256) :: log_buf
    integer :: i, j
    logical :: body_found

    call this%free()

    if (associated(neko_ale)) then
       if (neko_ale%active) this%mesh_has_moved = .true.
    end if

    ! Return if the mesh is fixed.
    if (.not. this%mesh_has_moved) return

    call json_get(json, 'mode', mode)

    ! The default scales with the working precision.
    call json_get_or_default(json, 'check_tolerance', &
         this%check_tolerance, max(1e-8_rp, 100.0_rp * epsilon(1.0_rp)))
    if (this%check_tolerance .le. 0.0_rp) then
       call neko_error('probes moving_mesh: check_tolerance must be > 0.')
    end if

    select case (trim(mode))
    case (PROBES_MM_FIXED)
       this%body_attached = .false.
    case (PROBES_MM_BODY)
       this%body_attached = .true.
    case default
       call neko_error("probes moving_mesh: mode must be '" // &
            PROBES_MM_FIXED // "' or '" // PROBES_MM_BODY // "'.")
    end select

    if (this%body_attached) then
       if (.not. neko_ale%has_moving_boundary) then
          call neko_error('probes moving_mesh: body_attached requires ' // &
               'a moving boundary')
       end if

       ! Link to the body through one of its zone ids.
       call json_get(json, 'zone_id', this%zone_id)
       body_found = .false.
       do i = 1, neko_ale%config%nbodies
          if (allocated(neko_ale%config%bodies(i)%zone_indices)) then
             do j = 1, size(neko_ale%config%bodies(i)%zone_indices)
                if (neko_ale%config%bodies(i)%zone_indices(j) .eq. &
                     this%zone_id) then
                   this%body_idx = i
                   this%body_name = neko_ale%config%bodies(i)%name
                   this%pivot0 = neko_ale%config%bodies(i)%rot_center
                   body_found = .true.
                end if
             end do
          end if
       end do
       if (.not. body_found) then
          write(log_buf, '(A,I0,A)') 'probes moving_mesh: zone ', &
               this%zone_id, ' is not registered as an ALE body.'
          call neko_error(trim(log_buf))
       end if
    end if

    this%output_coords = this%body_attached
    this%enabled = .true.

    call neko_log%section('Probes moving mesh')
    if (this%body_attached) then
       call neko_log%message('Mode            : ' // PROBES_MM_BODY)
       call neko_log%message('Linked ALE body : ' // trim(this%body_name))
       write(log_buf, '(A,I0)') 'Zone id         : ', this%zone_id
       call neko_log%message(log_buf)
    else
       call neko_log%message('Mode            : ' // PROBES_MM_FIXED)
    end if
    write(log_buf, '(A,E15.7)') 'Check tolerance : ', this%check_tolerance
    call neko_log%message(log_buf)
    call neko_log%end_section()

  end subroutine probes_moving_mesh_init

  !> Capture the reference coordinates.
  !! @param interp The global interpolation object of the probes.
  !! @param xyz Local probe coordinates (3, n_points).
  !! @param n_points Number of local probes.
  subroutine probes_moving_mesh_set_points(this, interp, xyz, n_points)
    class(probes_moving_mesh_t), intent(inout) :: this
    type(global_interpolation_t), intent(in) :: interp
    real(kind=rp), intent(in) :: xyz(:,:)
    integer, intent(in) :: n_points

    if (.not. this%enabled) return

    ! Re-initializing the interpolation resets these to defaults, so
    ! remember the values the probes were set up with.
    this%tolerance = interp%tolerance
    this%padding = interp%padding

    this%n_points = n_points
    allocate(this%xyz0(3, n_points))
    allocate(this%xyz_target(3, n_points))
    call copy(this%xyz0, xyz, 3 * n_points)
    call copy(this%xyz_target, xyz, 3 * n_points)
    call this%chk_x%init(n_points)
    call this%chk_y%init(n_points)
    call this%chk_z%init(n_points)

    this%n_lost = 0
    this%points_set = .true.

  end subroutine probes_moving_mesh_set_points

  !> Destructor.
  subroutine probes_moving_mesh_free(this)
    class(probes_moving_mesh_t), intent(inout) :: this

    if (allocated(this%xyz0)) deallocate(this%xyz0)
    if (allocated(this%xyz_target)) deallocate(this%xyz_target)
    call this%chk_x%free()
    call this%chk_y%free()
    call this%chk_z%free()

    this%enabled = .false.
    this%mesh_has_moved = .false.
    this%body_attached = .false.
    this%output_coords = .false.
    this%zone_id = -1
    this%body_idx = -1
    this%body_name = 'NOT_LINKED'
    this%pivot0 = 0.0_rp
    this%n_points = 0
    this%points_set = .false.
    this%n_lost = 0

  end subroutine probes_moving_mesh_free

  !> Update the probe coordinates, validate the mapping and re-find the
  !! points if needed.
  !! @param interp The global interpolation object of the probes.
  !! @param dof The dofmap holding the current mesh coordinates.
  !! @param xyz The probe coordinate array of the probes component,
  !! kept in sync with the current target coordinates.
  !! @param n_points Number of local probes.
  !! @param time Current time state (for logging).
  subroutine probes_moving_mesh_update(this, interp, dof, xyz, n_points, time)
    class(probes_moving_mesh_t), intent(inout) :: this
    type(global_interpolation_t), intent(inout) :: interp
    type(dofmap_t), intent(inout) :: dof
    real(kind=rp), intent(inout) :: xyz(:,:)
    integer, intent(in) :: n_points
    type(time_state_t), intent(in) :: time

    character(len=256) :: log_buf
    real(kind=rp) :: err
    real(kind=rp) :: t_start

    if (.not. this%enabled) return

    if (.not. this%points_set) then
       call neko_error('probes moving_mesh: set_points was not called.')
    end if
    if (n_points .ne. this%n_points) then
       call neko_error('probes moving_mesh: number of local probes ' // &
            'changed after initialization.')
    end if

    call this%compute_targets()
    call copy(xyz, this%xyz_target, 3 * this%n_points)

    t_start = MPI_Wtime()
    err = this%mapping_error(interp, dof)
    write(log_buf, '(A,E12.5)') &
         'Probes moving mesh timing (s): check    ', MPI_Wtime() - t_start
    call neko_log%message(log_buf)

    ! Re-find when the mapping no longer samples within tolerance of
    ! the targets, or while any point is lost.
    if (err .gt. this%check_tolerance .or. this%n_lost .gt. 0) then
       write(log_buf, '(A,I0,A,E15.7)') &
            'Probes moving mesh: re-finding points at step ', time%tstep, &
            ', mapping error: ', err
       call neko_log%message(log_buf)

       t_start = MPI_Wtime()
       call this%refind(interp, dof)
       write(log_buf, '(A,E12.5)') &
            'Probes moving mesh timing (s): re-find  ', MPI_Wtime() - t_start
       call neko_log%message(log_buf)

       if (this%n_lost .gt. 0) then
          write(log_buf, '(A,I0,A)') 'Probes moving mesh: ', this%n_lost, &
               ' point(s) not found in the mesh, evaluating to 0.'
          call neko_log%message(log_buf)
       end if
    end if

  end subroutine probes_moving_mesh_update

  !> Re-find all points in the current mesh.
  subroutine probes_moving_mesh_refind(this, interp, dof)
    class(probes_moving_mesh_t), intent(inout) :: this
    type(global_interpolation_t), intent(inout) :: interp
    type(dofmap_t), intent(inout) :: dof

    character(len=256) :: log_buf
    real(kind=rp) :: t_start
    integer :: i, n, ierr

    if (NEKO_BCKND_DEVICE .eq. 1) then
       t_start = MPI_Wtime()
       n = dof%size()
       call device_memcpy(dof%x, dof%x_d, n, DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(dof%y, dof%y_d, n, DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(dof%z, dof%z_d, n, DEVICE_TO_HOST, sync = .true.)
       write(log_buf, '(A,E12.5)') &
            'Probes moving mesh timing (s): dof sync D2H  ', &
            MPI_Wtime() - t_start
       call neko_log%message(log_buf)
    end if

    call interp%init(dof, tol = this%tolerance, pad = this%padding)

    ! Only find points, no redistribution.
    call interp%find_points(this%xyz_target, this%n_points)

    this%n_lost = 0
    ! A point is lost if no rank finds a containing element.
    do i = 1, this%n_points
       if ( (interp%pe_owner(i) .eq. -1) .or. &
            (interp%el_owner0(i) .eq. -1) .or. &
            (maxval(abs(interp%rst(:, i))) .gt. 1.1_rp) ) then
         this%n_lost = this%n_lost + 1
       end if
    end do
    call MPI_Allreduce(MPI_IN_PLACE, this%n_lost, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)

  end subroutine probes_moving_mesh_refind

  !> Compute the current target coordinates of all local probes.
  subroutine probes_moving_mesh_compute_targets(this)
    class(probes_moving_mesh_t), intent(inout) :: this
    real(kind=rp) :: R(3,3), P(3)
    integer :: i

    if (this%body_attached) then
       R = neko_ale%body_rot_matrices(:, :, this%body_idx)
       P = neko_ale%ale_pivot(this%body_idx)%pos
       do i = 1, this%n_points
          this%xyz_target(:, i) = P + &
               matmul(R, this%xyz0(:, i) - this%pivot0)
       end do
    end if

  end subroutine probes_moving_mesh_compute_targets

  !> Largest distance, over all ranks, between the position the stored
  !! mapping samples in the current mesh and the target position. Lost
  !! points are excluded.
  function probes_moving_mesh_mapping_error(this, interp, dof) result(err)
    class(probes_moving_mesh_t), intent(inout) :: this
    type(global_interpolation_t), intent(inout) :: interp
    type(dofmap_t), intent(inout) :: dof
    real(kind=rp) :: err

    integer :: i, ierr

    ! Interpolating the current coordinate fields through the stored
    ! mapping gives the point the probe would sample now.
    call interp%evaluate(this%chk_x%x, dof%x, .false.)
    call interp%evaluate(this%chk_y%x, dof%y, .false.)
    call interp%evaluate(this%chk_z%x, dof%z, .false.)
    call this%chk_x%copy_from(DEVICE_TO_HOST, .false.)
    call this%chk_y%copy_from(DEVICE_TO_HOST, .false.)
    call this%chk_z%copy_from(DEVICE_TO_HOST, .true.)

    err = 0.0_rp
    do i = 1, this%n_points
       ! Excluding lost points, until they're found again.
       if ( (interp%pe_owner(i) .eq. -1) .or. &
            (interp%el_owner0(i) .eq. -1) .or. &
            (maxval(abs(interp%rst(:, i))) .gt. 1.1_rp) ) cycle

       err = max(err, norm2([this%chk_x%x(i) - this%xyz_target(1, i), &
            this%chk_y%x(i) - this%xyz_target(2, i), &
            this%chk_z%x(i) - this%xyz_target(3, i)]))
    end do

    call MPI_Allreduce(MPI_IN_PLACE, err, 1, MPI_REAL_PRECISION, &
         MPI_MAX, NEKO_COMM, ierr)

  end function probes_moving_mesh_mapping_error

end module probes_moving_mesh
