! Copyright (c) 2020-2025, The Neko Authors
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
!> Dirichlet condition applied in the facet normal direction
module facet_normal
  use num_types, only : rp
  use neko_config, only : NEKO_BCKND_DEVICE
  use math, only : cfill_mask
  use device_math, only : device_col2, device_masked_gather_copy_0, &
       device_masked_scatter_copy_0
  use vector, only : vector_t
  use coefs, only : coef_t
  use bc, only : bc_t
  use utils, only : neko_error, nonlinear_index
  use json_module, only : json_file
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_associated
  use htable, only : htable_i4_t
  use device, only : device_map, device_memcpy, device_unmap, &
       HOST_TO_DEVICE, DEVICE_TO_HOST, glb_cmd_queue
  use time_state, only : time_state_t
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use amr_reconstruct, only : amr_reconstruct_t
  implicit none
  private

  !> Dirichlet condition in facet normal direction
  type, public, extends(bc_t) :: facet_normal_t
     integer, allocatable :: unique_mask(:)
     type(c_ptr) :: unique_mask_d = c_null_ptr
     type(vector_t) :: nx, ny, nz, work
   contains
     procedure, pass(this) :: apply_scalar => facet_normal_apply_scalar
     procedure, pass(this) :: apply_scalar_dev => facet_normal_apply_scalar_dev
     procedure, pass(this) :: apply_vector => facet_normal_apply_vector
     procedure, pass(this) :: apply_vector_dev => facet_normal_apply_vector_dev
     procedure, pass(this) :: apply_surfvec => facet_normal_apply_surfvec
     procedure, pass(this) :: apply_surfvec_dev => &
          facet_normal_apply_surfvec_dev
     !> Constructor.
     procedure, pass(this) :: init => facet_normal_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          facet_normal_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => facet_normal_free
     !> Finalize.
     procedure, pass(this) :: finalize => facet_normal_finalize
     !> AMR restart
     procedure, pass(this) :: amr_restart => facet_normal_amr_restart
  end type facet_normal_t

contains

  !> Constructor.
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine facet_normal_init(this, coef, json)
    class(facet_normal_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json

    call this%init_from_components(coef)
  end subroutine facet_normal_init

  !> Constructor from components.
  !! @param[in] coef The SEM coefficients.
  subroutine facet_normal_init_from_components(this, coef)
    class(facet_normal_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef

    call this%init_base(coef)
  end subroutine facet_normal_init_from_components

  !> No-op scalar apply
  subroutine facet_normal_apply_scalar(this, x, n, time, strong)
    class(facet_normal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
  end subroutine facet_normal_apply_scalar

  !> No-op scalar apply on device
  subroutine facet_normal_apply_scalar_dev(this, x_d, time, strong, strm)
    class(facet_normal_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

  end subroutine facet_normal_apply_scalar_dev

  !> No-op vector apply on device
  subroutine facet_normal_apply_vector_dev(this, x_d, y_d, z_d, time, &
       strong, strm)
    class(facet_normal_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

  end subroutine facet_normal_apply_vector_dev

  !> No-op vector apply
  subroutine facet_normal_apply_vector(this, x, y, z, n, time, strong)
    class(facet_normal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
  end subroutine facet_normal_apply_vector

  !> Apply in facet normal direction (vector valued)
  subroutine facet_normal_apply_surfvec(this, x, y, z, u, v, w, n, time)
    class(facet_normal_t), intent(in) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    real(kind=rp), intent(inout), dimension(n) :: u
    real(kind=rp), intent(inout), dimension(n) :: v
    real(kind=rp), intent(inout), dimension(n) :: w
    type(time_state_t), intent(in), optional :: time
    integer :: i, m, k, idx(4), facet
    real(kind=rp) :: normal(3), area

    m = this%unique_mask(0)

    do i = 1, m
       k = this%unique_mask(i)
       x(k) = u(k) * this%nx%x(i)
       y(k) = v(k) * this%ny%x(i)
       z(k) = w(k) * this%nz%x(i)
    end do

  end subroutine facet_normal_apply_surfvec

  !> Apply in facet normal direction (vector valued, device version)
  subroutine facet_normal_apply_surfvec_dev(this, x_d, y_d, z_d, &
       u_d, v_d, w_d, time, strm)
    class(facet_normal_t), intent(in), target :: this
    type(c_ptr) :: x_d, y_d, z_d, u_d, v_d, w_d
    type(time_state_t), intent(in), optional :: time
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_
    integer :: n, m

    n = this%coef%dof%size()
    m = this%unique_mask(0)

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

    if (m .gt. 0) then
       call device_masked_gather_copy_0(this%work%x_d, u_d, &
            this%unique_mask_d, n, m, strm_)
       call device_col2(this%work%x_d, this%nx%x_d, m, strm_)
       call device_masked_scatter_copy_0(x_d, this%work%x_d, &
            this%unique_mask_d, n, m, strm_)
       call device_masked_gather_copy_0(this%work%x_d, v_d, &
            this%unique_mask_d, n, m, strm_)
       call device_col2(this%work%x_d, this%ny%x_d, m, strm_)
       call device_masked_scatter_copy_0(y_d, this%work%x_d, &
            this%unique_mask_d, n, m, strm_)
       call device_masked_gather_copy_0(this%work%x_d, w_d, &
            this%unique_mask_d, n, m, strm_)
       call device_col2(this%work%x_d, this%nz%x_d, m, strm)
       call device_masked_scatter_copy_0(z_d, this%work%x_d, &
            this%unique_mask_d, n, m, strm_)
    end if

  end subroutine facet_normal_apply_surfvec_dev

  !> Destructor
  subroutine facet_normal_free(this)
    class(facet_normal_t), target, intent(inout) :: this

    call this%free_base()
    if (allocated(this%unique_mask)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_unmap(this%unique_mask, this%unique_mask_d)
       end if
       deallocate(this%unique_mask)
    end if

    call this%nx%free()
    call this%ny%free()
    call this%nz%free()
    call this%work%free()

    call this%free_amr_base()

  end subroutine facet_normal_free

  !> Finalize
  subroutine facet_normal_finalize(this, only_facets)
    class(facet_normal_t), target, intent(inout) :: this
    logical, optional, intent(in) :: only_facets
    logical :: only_facets_
    type(htable_i4_t) :: unique_point_idx
    integer :: htable_data, rcode, i, j, idx(4), facet
    real(kind=rp) :: area, normal(3)

    if (present(only_facets)) then
       if (.not. only_facets) then
          call neko_error("For facet_normal_t, only_facets has to be true.")
       end if
    end if

    call this%finalize_base(.true.)
    ! This part is purely needed to ensure that contributions
    ! for all faces a point is on is properly summed up.
    ! If one simply uses the original mask, if a point is on a corner
    ! where both faces are on the boundary
    ! one will only get the contribution from one face, not both
    ! We solve this by adding up the normals of both faces for these points
    ! and storing this sum in this%nx, this%ny, this%nz.
    ! As both contrbutions are added already,
    ! we also ensure that we only visit each point once
    ! and create a new mask with only unique points (this%unique_mask).
    if (allocated(this%unique_mask)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_unmap(this%unique_mask, this%unique_mask_d)
       end if
       deallocate(this%unique_mask)
    end if

    call unique_point_idx%init(this%msk(0), htable_data)
    j = 0
    do i = 1, this%msk(0)
       if (unique_point_idx%get(this%msk(i), htable_data) .ne. 0) then
          j = j + 1
          htable_data = j
          call unique_point_idx%set(this%msk(i), j)
       end if
    end do

    ! Only allocate work vectors if size is non-zero
    if (unique_point_idx%num_entries() .gt. 0 ) then
       call this%nx%init(unique_point_idx%num_entries())
       call this%ny%init(unique_point_idx%num_entries())
       call this%nz%init(unique_point_idx%num_entries())
       call this%work%init(unique_point_idx%num_entries())
    end if
    allocate(this%unique_mask(0:unique_point_idx%num_entries()))

    this%unique_mask(0) = unique_point_idx%num_entries()
    do i = 1, this%unique_mask(0)
       this%unique_mask(i) = 0
    end do


    do i = 1, this%msk(0)
       rcode = unique_point_idx%get(this%msk(i), htable_data)
       if (rcode .ne. 0) call neko_error("Facet normal: htable get failed.")
       this%unique_mask(htable_data) = this%msk(i)
       facet = this%facet(i)

       idx = nonlinear_index(this%msk(i), this%Xh%lx, this%Xh%lx, this%Xh%lx)
       normal = this%coef%get_normal(idx(1), idx(2), idx(3), idx(4), facet)
       area = this%coef%get_area(idx(1), idx(2), idx(3), idx(4), facet)
       normal = normal * area !Scale normal by area
       this%nx%x(htable_data) = this%nx%x(htable_data) + normal(1)
       this%ny%x(htable_data) = this%ny%x(htable_data) + normal(2)
       this%nz%x(htable_data) = this%nz%x(htable_data) + normal(3)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1 .and. &
         (unique_point_idx%num_entries() .gt. 0 )) then
       call device_map(this%unique_mask, this%unique_mask_d, &
            size(this%unique_mask))
       call device_memcpy(this%unique_mask, this%unique_mask_d, &
            size(this%unique_mask), HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%nx%x, this%nx%x_d, &
            this%nx%size(), HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%ny%x, this%ny%x_d, &
            this%ny%size(), HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%nz%x, this%nz%x_d, &
            this%nz%size(), HOST_TO_DEVICE, sync = .true.)
    end if

    call unique_point_idx%free()

  end subroutine facet_normal_finalize

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     tstep         time step
  subroutine facet_normal_amr_restart(this, reconstruct, counter, tstep)
    class(facet_normal_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter, tstep
    character(len=LOG_SIZE) :: log_buf
    integer :: il

    ! Was this component already restarted?
    if (this%counter .eq. counter) return

    this%counter = counter

    ! For defined zone indices perform full reconstruction including
    ! finalisation. If zones are missing just prepare for collecting
    ! facets. Do not forget to finalise those bc later.
    if (allocated(this%zone_indices)) then
       if (allocated(this%type)) then
          log_buf = 'Reconstruct Facet normal: '//trim(this%type)
       else
          log_buf = 'Reconstruct Facet normal'
       end if
       call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

       ! reconstruct dofmap; No problem, as AMR restart prevents recursive
       ! reconstructions
       if (associated(this%dof)) call this%dof%amr_restart(reconstruct, &
            counter, tstep)
       ! reconstruct coef; No problem, as AMR restart prevents recursive
       ! reconstructions
       if (associated(this%coef)) call this%coef%amr_restart(reconstruct, &
            counter, tstep)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          ! added utils module; could be removed
          call neko_error('Facet normal:: Nothing done for device.')
       end if

       ! free space
       if (allocated(this%msk)) deallocate(this%msk)
       if (allocated(this%facet)) deallocate(this%facet)
!       call this%marked_facet%free()
!       call this%marked_facet%init()
       call this%marked_facet%clear()

       this%iffinalised = .false.

       ! get zones
       do il = 1, size(this%zone_indices)
          call this%mark_zone(this%coef%msh%labeled_zones(&
               this%zone_indices(il)))
       end do
       call this%finalize()

    else
       if (allocated(this%type)) then
          log_buf = 'Clean Facet normal: '//trim(this%type)
       else
          log_buf = 'Clean Facet normal'
       end if
       call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

       ! reconstruct dofmap; No problem, as AMR restart prevents recursive
       ! reconstructions
       if (associated(this%dof)) call this%dof%amr_restart(reconstruct, &
            counter, tstep)
       ! reconstruct coef; No problem, as AMR restart prevents recursive
       ! reconstructions
       if (associated(this%coef)) call this%coef%amr_restart(reconstruct, &
            counter, tstep)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          ! added utils module; could be removed
          call neko_error('Facet normal:: Nothing done for device.')
       end if

       ! free space
       if (allocated(this%msk)) deallocate(this%msk)
       if (allocated(this%facet)) deallocate(this%facet)
!       call this%marked_facet%free()
!       call this%marked_facet%init()
       call this%marked_facet%clear()

       if (allocated(this%unique_mask)) then
          deallocate(this%unique_mask)
       end if

       call this%nx%free()
       call this%ny%free()
       call this%nz%free()
       call this%work%free()

       this%iffinalised = .false.

    end if

  end subroutine facet_normal_amr_restart

end module facet_normal
