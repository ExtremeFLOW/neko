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
!> Defines overset interface scalar boundary conditions
module overset_interface
  use comm, only : NEKO_GLOBAL_COMM
  use neko_config, only : NEKO_BCKND_DEVICE
  use registry, only : neko_registry
  use num_types, only : rp
  use coefs, only : coef_t
  use global_interpolation, only : global_interpolation_t, &
       global_interpolation_settings_t
  use mask, only : mask_t
  use bc, only : bc_t, BC_DIRICHLET
  use field_list, only : field_list_t
  use math, only : masked_copy_0, copy
  use device_math, only : device_masked_copy_0, device_copy
  use vector, only : vector_t
  use vector_series, only : vector_series_t
  use vector_list, only : vector_list_t
  use vector_math, only : vector_masked_gather_copy, &
       vector_masked_scatter_copy, vector_add2s2, vector_cmult, vector_cmult2, &
       vector_glsc2
  use device, only : DEVICE_TO_HOST
  use field_dirichlet, only : field_dirichlet_t
  use iextm_time_scheme, only : iextm_time_scheme_t
  use utils, only : neko_error, nonlinear_index, linear_index
  use stack, only : stack_i4_t
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use field, only : field_t
  use field_series, only : field_series_t
  use logger, only : neko_log, LOG_SIZE
  use scratch_registry, only : neko_scratch_registry
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_SUM
  use, intrinsic :: iso_c_binding, only : c_ptr
  use time_state, only : time_state_t
  implicit none
  private

  !> Overset interface BC for a scalar field.
  type, public, extends(bc_t) :: overset_interface_t
     !> Underlying scalar field Dirichlet bc.
     type(field_dirichlet_t) :: bc_s
     !> Single-field list for compatibility with field-based update patterns.
     type(field_list_t) :: field_list
     !> Name of the scalar field to interpolate from the registry.
     character(len=:), allocatable :: field_name
     !> Interpolator
     type(global_interpolation_t) :: interface_interpolator
     !> Mask for overset interface points.
     type(mask_t) :: interface_dof_mask
     type(mask_t) :: domain_element_mask
     !> Vectors holding dof coordinates for cases where mesh moves.
     type(vector_t) :: x_dof, y_dof, z_dof
     type(vector_t) :: x_interface_dof, y_interface_dof, z_interface_dof
     !> Interpolated scalar values on the interface.
     type(vector_t) :: s_interface
     type(vector_series_t) :: s_interface_lag
     integer :: iextm_order = 1
     !> Under-relaxation factor for updated interface values.
     real(kind=rp) :: relaxation = 1.0_rp
     integer :: last_tstep = -1
     type(vector_list_t) :: interface_dof, interface_field
     !> Interpolation settings.
     type(global_interpolation_settings_t) :: interpolation_settings
     integer :: n_int_tot = 0
     logical :: find_interface = .false.
     logical :: setup = .false.
     logical :: log = .false.
     !> Skip donor interpolation on the first update after restoring history.
     logical :: restart_pending = .false.

     !> Function pointer to the user routine performing the update of the values
     !! of the boundary fields.
     procedure(morph_overset_interface), nopass, pointer :: &
          morph_interface => null()

   contains
     !> Constructor.
     procedure, pass(this) :: init => overset_interface_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          overset_interface_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => overset_interface_free
     !> Finalize.
     procedure, pass(this) :: finalize => overset_interface_finalize
     !> Apply scalar by performing a masked copy.
     procedure, pass(this) :: apply_scalar => overset_interface_apply_scalar
     !> (No-op) Apply vector.
     procedure, pass(this) :: apply_vector => overset_interface_apply_vector
     !> (No-op) Apply vector (device).
     procedure, pass(this) :: apply_vector_dev => &
          overset_interface_apply_vector_dev
     !> Apply scalar (device).
     procedure, pass(this) :: apply_scalar_dev => &
          overset_interface_apply_scalar_dev
     procedure, pass(this) :: update => overset_interface_update
     !> Restore scalar interface history from accepted solution fields.
     procedure, pass(this) :: restart_scalar => overset_interface_restart

     !> Build domain masks for the overset interface.
     procedure, pass(this), private :: build_masks_ => build_masks_
     !> Gather the dofs at the interface.
     procedure, pass(this), private :: gather_interface_dofs_ => &
          gather_interface_dofs_
     !> Set up the interpolator.
     procedure, pass(this), private :: setup_interpolator_ => &
          setup_interpolator_
     !> Log interface interpolation error diagnostics.
     procedure, pass(this), private :: log_interface_error_ => &
          log_interface_error_
     !> Under-relax the new interface value with the previously applied one.
     procedure, pass(this), private :: relax_interface_value_ => &
          relax_interface_value_
  end type overset_interface_t

  abstract interface

     !> User callback for overset-interface morphing and boundary-value updates.
     !!
     !! Implementations may update interface coordinates and/or prescribed
     !! interface field values before interpolation is evaluated.
     !!
     !! @param[inout] interface_dof Interface coordinates as a vector list
     !!               (x, y, z). To be updated in the routine
     !! @param[inout] interface_field Interface boundary field values
     !! @param[in] interface_mask Mask describing active interface degrees of
     !!                freedom.
     !! @param[in] time Current simulation time state.
     !! @param[in] bc_name Name of the boundary condition invoking the callback.
     !! @param[inout] find_interface Set to .true. when interpolation
     !! points must be rediscovered after coordinate changes.
     subroutine morph_overset_interface(interface_dof, interface_field, &
          interface_mask, time, bc_name, &
          find_interface)
       import vector_list_t, mask_t, time_state_t
       type(vector_list_t), intent(inout) :: interface_dof
       type(vector_list_t), intent(inout) :: interface_field
       type(mask_t), intent(in) :: interface_mask
       type(time_state_t), intent(in) :: time
       character(len=*), intent(in) :: bc_name
       logical, intent(inout) :: find_interface
     end subroutine morph_overset_interface
  end interface

  public :: morph_overset_interface

contains

  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine overset_interface_init(this, coef, json)
    class(overset_interface_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: field_name
    real(kind=rp) :: tol, pad, relaxation
    logical :: log

    call json_get(json, "field_name", field_name)
    call json_get_or_default(json, "interpolation.tolerance", tol, -1.0_rp)
    call json_get_or_default(json, "interpolation.padding", pad, -1.0_rp)
    call json_get_or_default(json, "order", this%iextm_order, 1)
    if (this%iextm_order .lt. 1 .or. this%iextm_order .gt. 3) then
       call neko_error("The order of the IEXTm time scheme must be 1 to 3.")
    end if
    call json_get_or_default(json, "relaxation", relaxation, 1.0_rp)
    if (relaxation .le. 0.0_rp .or. relaxation .gt. 1.0_rp) then
       call neko_error("The overset relaxation factor must be in (0, 1].")
    end if
    call json_get_or_default(json, "log", log, .false.)

    call this%init_from_components(coef, field_name, tol, pad, log, relaxation)
    if (allocated(field_name)) deallocate(field_name)

  end subroutine overset_interface_init

  !> Constructor from components
  !! @param[in] coef The SEM coefficients.
  !! @param[in] relaxation Under-relaxation factor for interface updates.
  subroutine overset_interface_init_from_components(this, coef, field_name, &
       tol, pad, log, relaxation)
    class(overset_interface_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    character(len=*), intent(in) :: field_name
    real(kind=rp), intent(in), optional :: tol, pad, relaxation
    logical, intent(in), optional :: log
    character(len=256) :: log_buf

    call this%init_base(coef)
    this%bc_type = BC_DIRICHLET
    this%relaxation = 1.0_rp
    this%last_tstep = -1
    this%restart_pending = .false.

    if (present(tol)) then
       if (tol .gt. 0.0_rp) then
          this%interpolation_settings%tolerance = tol
       end if
    end if

    if (present(pad)) then
       if (pad .gt. 0.0_rp) then
          this%interpolation_settings%padding = pad
       end if
    end if

    if (present(log)) then
       this%log = log
    end if

    if (present(relaxation)) then
       if (relaxation .le. 0.0_rp .or. relaxation .gt. 1.0_rp) then
          call neko_error("The overset relaxation factor must be in (0, 1].")
       end if
       this%relaxation = relaxation
    end if

    this%field_name = field_name
    write (log_buf, '(A,A)') "Coupling overset interface for: ", &
         trim(this%field_name)
    call neko_log%message(log_buf)

    call this%bc_s%init_from_components(coef, this%field_name)
    call this%field_list%init(1)
    call this%field_list%assign_to_field(1, this%bc_s%field_bc)

    call this%x_dof%init(this%dof%size(), 'x')
    call this%y_dof%init(this%dof%size(), 'y')
    call this%z_dof%init(this%dof%size(), 'z')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%x_dof%x_d, this%dof%x_d, this%dof%size())
       call device_copy(this%y_dof%x_d, this%dof%y_d, this%dof%size())
       call device_copy(this%z_dof%x_d, this%dof%z_d, this%dof%size())

       call this%x_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
       call this%y_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
       call this%z_dof%copy_from(DEVICE_TO_HOST, sync = .true.)
    else
       call copy(this%x_dof%x, this%dof%x, this%dof%size())
       call copy(this%y_dof%x, this%dof%y, this%dof%size())
       call copy(this%z_dof%x, this%dof%z, this%dof%size())
    end if

  end subroutine overset_interface_init_from_components

  !> Destructor.
  subroutine overset_interface_free(this)
    class(overset_interface_t), target, intent(inout) :: this

    call this%bc_s%free()
    call this%field_list%free()
    call this%interface_dof%free()
    call this%interface_field%free()

    call this%x_dof%free()
    call this%y_dof%free()
    call this%z_dof%free()

    call this%x_interface_dof%free()
    call this%y_interface_dof%free()
    call this%z_interface_dof%free()
    call this%s_interface%free()
    call this%s_interface_lag%free()

    if (allocated(this%field_name)) then
       deallocate(this%field_name)
    end if

    call this%interface_interpolator%free()

    call this%interface_dof_mask%free()
    call this%domain_element_mask%free()

    call this%free_base()
    this%restart_pending = .false.
  end subroutine overset_interface_free

  !> Restore scalar interface history from an accepted solution field.
  !! The lag series is present for transported scalars and absent for
  !! pressure. Historical values are inserted oldest-to-newest, while the
  !! current accepted value remains in the base interface vector for the first
  !! regular update to shift into lag position one.
  subroutine overset_interface_restart(this, s, slag)
    class(overset_interface_t), intent(inout) :: this
    type(field_t), intent(in) :: s
    type(field_series_t), intent(in), optional :: slag
    integer :: i, n_previous

    call this%s_interface_lag%reset()

    n_previous = 0
    if (present(slag)) then
       n_previous = min(this%iextm_order - 1, slag%size())
       do i = n_previous, 1, -1
          call vector_masked_gather_copy(this%s_interface, &
               slag%lf(i)%x(:,1,1,1), this%interface_dof_mask, &
               slag%lf(i)%dof%size())
          call this%s_interface_lag%update()
       end do
    end if

    call vector_masked_gather_copy(this%s_interface, &
         s%x(:,1,1,1), this%interface_dof_mask, s%dof%size())

    this%restart_pending = .true.
  end subroutine overset_interface_restart

  !> Apply scalar.
  !! @param x Field onto which to copy the values.
  !! @param n Size of the array `x`.
  !! @param time The current time state.
  subroutine overset_interface_apply_scalar(this, x, n, time, strong)
    class(overset_interface_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then
       if (.not. this%updated) then
          call this%update(time)
          this%updated = .true.
       end if

       call masked_copy_0(x, this%bc_s%field_bc%x, this%msk, n, this%msk(0))
    end if

  end subroutine overset_interface_apply_scalar

  !> Apply scalar (device).
  !! @param x_d Device pointer to the field onto which to copy the values.
  !! @param time The current time state.
  !! @param strm Device stream
  subroutine overset_interface_apply_scalar_dev(this, x_d, time, strong, strm)
    class(overset_interface_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then
       if (.not. this%updated) then
          call this%update(time)
          this%updated = .true.
       end if

       if (this%msk(0) .gt. 0) then
          call device_masked_copy_0(x_d, this%bc_s%field_bc%x_d, &
               this%bc_s%msk_d, this%bc_s%dof%size(), this%msk(0), strm)
       end if
    end if

  end subroutine overset_interface_apply_scalar_dev

  !> (No-op) Apply vector.
  subroutine overset_interface_apply_vector(this, x, y, z, n, time, strong)
    class(overset_interface_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong

    call neko_error("overset_interface cannot apply vector BCs.&
    & Use overset_interface_vector instead!")

  end subroutine overset_interface_apply_vector

  !> (No-op) Apply vector (device).
  subroutine overset_interface_apply_vector_dev(this, x_d, y_d, z_d, time, &
       strong, strm)
    class(overset_interface_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

    call neko_error("overset_interface cannot apply vector BCs.&
    & Use overset_interface_vector instead!")

  end subroutine overset_interface_apply_vector_dev

  !> Finalize by building the mask arrays and preparing interpolation data.
  subroutine overset_interface_finalize(this)
    class(overset_interface_t), target, intent(inout) :: this

    call this%finalize_base()

    call this%bc_s%mark_facets(this%marked_facet)
    call this%bc_s%finalize()

    call this%build_masks_()

    call this%x_interface_dof%init(this%interface_dof_mask%size(), &
         'x_interface')
    call this%y_interface_dof%init(this%interface_dof_mask%size(), &
         'y_interface')
    call this%z_interface_dof%init(this%interface_dof_mask%size(), &
         'z_interface')
    call this%gather_interface_dofs_()

    call this%setup_interpolator_()

    call this%s_interface%init(this%interface_dof_mask%size(), 's_interface')

    call this%interface_dof%init(3)
    call this%interface_dof%assign_to_vector(1, this%x_interface_dof)
    call this%interface_dof%assign_to_vector(2, this%y_interface_dof)
    call this%interface_dof%assign_to_vector(3, this%z_interface_dof)

    call this%interface_field%init(1)
    call this%interface_field%assign_to_vector(1, this%s_interface)

    call this%s_interface_lag%init(this%s_interface, this%iextm_order)

    call MPI_Allreduce(this%s_interface%size(), this%n_int_tot, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_GLOBAL_COMM)

  end subroutine overset_interface_finalize

  !> Update values at the overset interface.
  subroutine overset_interface_update(this, time)
    class(overset_interface_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: s
    type(iextm_time_scheme_t) :: time_scheme
    integer :: nhist, ihist
    real(kind=rp) :: iextm_coeffs(4)
    logical :: new_tstep

    !> Change the coordinates of the interface if set up by the user
    call this%morph_interface(this%interface_dof, this%interface_field, &
         this%interface_dof_mask, time, this%name, &
         this%find_interface)

    !> Find points if needed - later make sure only in first substep
    if (this%find_interface) then

       ! sync
       call this%x_interface_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
       call this%y_interface_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
       call this%z_interface_dof%copy_from(DEVICE_TO_HOST, sync = .true.)

       call this%interface_interpolator%find_points(this%x_interface_dof%x, &
            this%y_interface_dof%x, this%z_interface_dof%x, &
            this%x_interface_dof%size())
       this%find_interface = .false.

    end if

    s => neko_registry%get_field(trim(this%field_name))

    ! The current accepted interface value is restored locally from the
    ! solution field, so cross-domain interpolation is unnecessary once.
    if (.not. this%restart_pending) then
       call this%interface_interpolator%evaluate_masked(this%s_interface%x, &
            s%x, this%domain_element_mask, .false.)

       if (this%log) then
          call this%log_interface_error_(s)
       end if
    end if

    new_tstep = time%tstep .ne. this%last_tstep

    if (new_tstep) then
       this%last_tstep = time%tstep

       call this%s_interface_lag%update()

       nhist = min(this%s_interface_lag%filled_size(), this%iextm_order)
       call time_scheme%compute_coeffs(iextm_coeffs, &
            real(time%dtlag, kind=rp), nhist)

       call vector_cmult2(this%s_interface, this%s_interface_lag%lv(1), &
            iextm_coeffs(1))
       do ihist = 2, nhist
          call vector_add2s2(this%s_interface, this%s_interface_lag%lv(ihist), &
               iextm_coeffs(ihist))
       end do

       this%restart_pending = .false.
    end if

    ! Preserve the IEXT prediction on the first pass of every physical
    ! timestep. Relax only subsequent Schwarz corrections at the same tstep.
    if (.not. new_tstep) call this%relax_interface_value_()

    call vector_masked_scatter_copy(this%bc_s%field_bc%x(:,1,1,1), &
         this%s_interface, this%interface_dof_mask, this%bc_s%dof%size())

    nullify(s)

  end subroutine overset_interface_update

  !> Under-relax a Schwarz correction using the previously applied interface.
  !! Blend the new donor value `g_new` with the preceding Schwarz iterate
  !! `g_old` as
  !! `g = relaxation * g_new + (1 - relaxation) * g_old`.
  !! The caller skips this routine on the first pass of every physical
  !! timestep, preserving the temporal accuracy of the IEXT prediction.
  subroutine relax_interface_value_(this)
    class(overset_interface_t), intent(inout) :: this
    type(vector_t), pointer :: previous
    integer :: ind(1)
    logical :: clear_scratch = .false.

    ! A factor of one recovers the original, unrelaxed Schwarz iteration.
    if (this%relaxation .ge. 1.0_rp) return

    ! Gather the preceding scalar iterate into reusable temporary storage.
    call neko_scratch_registry%request_vector(previous, ind(1), &
         this%s_interface%size(), clear_scratch)
    call vector_masked_gather_copy(previous, this%bc_s%field_bc%x(:,1,1,1), &
         this%interface_dof_mask, this%bc_s%dof%size())

    ! Relax the new scalar donor data against the preceding iterate.
    call vector_cmult(this%s_interface, this%relaxation)
    call vector_add2s2(this%s_interface, previous, &
         1.0_rp - this%relaxation)

    ! Return the temporary storage to the scratch registry.
    call neko_scratch_registry%relinquish(ind)

  end subroutine relax_interface_value_

  !> Log interface RMSE for the scalar field.
  subroutine log_interface_error_(this, s)
    class(overset_interface_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: s
    real(kind=rp) :: s_int_norm
    type(vector_t), pointer :: error
    integer :: ind(1)
    logical :: clear_scratch = .false.
    character(len=256) :: log_buf

    call neko_scratch_registry%request_vector(error, ind(1), this%s_interface%size(), &
         clear_scratch)
    call vector_masked_gather_copy(error, s%x(:,1,1,1), this%interface_dof_mask, &
         this%dof%size())
    call vector_add2s2(error, this%s_interface, -1.0_rp)
    s_int_norm = sqrt(vector_glsc2(error, error)) / sqrt(real(this%n_int_tot, kind=rp))
    call neko_scratch_registry%relinquish(ind)

    write(log_buf, '(A12,A3,A10,1x,E15.7)') 'Interface BC', ' | ', &
         'L2 Error: ', s_int_norm
    call neko_log%message(log_buf)

  end subroutine log_interface_error_

  !===================
  ! Helper subroutines
  !===================

  !> Build masks.
  subroutine build_masks_(this)
    class(overset_interface_t), intent(inout) :: this
    type(mask_t) :: temp_mask
    logical, allocatable :: found(:)
    integer :: i, j, k, e, nelems
    integer :: lx, ly, lz
    integer :: nonlinear_idx(4), linear_idx
    type(stack_i4_t) :: idx_stack

    call this%interface_dof_mask%init(this%msk(1:this%msk(0)), this%msk(0))

    lx = this%Xh%lx
    ly = this%Xh%ly
    lz = this%Xh%lz

    allocate(found(this%msh%nelv))
    found = .false.

    do i = 1, this%msk(0)
       linear_idx = this%msk(i)
       nonlinear_idx = nonlinear_index(linear_idx, lx, ly, lz)
       found(nonlinear_idx(4)) = .true.
    end do

    nelems = 0
    call idx_stack%init()
    do e = 1, this%msh%nelv
       if (found(e)) then
          nelems = nelems + 1
          do k = 1, this%Xh%lz
             do j = 1, this%Xh%ly
                do i = 1, this%Xh%lx
                   linear_idx = linear_index(i, j, k, e, lx, ly, lz)
                   call idx_stack%push(linear_idx)
                end do
             end do
          end do
       end if
    end do

    deallocate(found)

    call temp_mask%init(idx_stack%array(), idx_stack%size())
    call idx_stack%free()

    call this%domain_element_mask%invert_mask(temp_mask, this%dof%size())
    call temp_mask%free()

  end subroutine build_masks_

  !> Gather interface dofs.
  subroutine gather_interface_dofs_(this)
    class(overset_interface_t), intent(inout) :: this

    call vector_masked_gather_copy(this%x_interface_dof, this%dof%x(:,1,1,1), &
         this%interface_dof_mask, this%dof%size())
    call vector_masked_gather_copy(this%y_interface_dof, this%dof%y(:,1,1,1), &
         this%interface_dof_mask, this%dof%size())
    call vector_masked_gather_copy(this%z_interface_dof, this%dof%z(:,1,1,1), &
         this%interface_dof_mask, this%dof%size())

    call this%x_interface_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
    call this%y_interface_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
    call this%z_interface_dof%copy_from(DEVICE_TO_HOST, sync = .true.)

  end subroutine gather_interface_dofs_

  !> Set up the global interpolator.
  subroutine setup_interpolator_(this)
    class(overset_interface_t), intent(inout) :: this

    call this%interface_interpolator%init(this%dof, &
         comm = NEKO_GLOBAL_COMM, &
         tol = this%interpolation_settings%tolerance, &
         pad = this%interpolation_settings%padding, &
         mask = this%domain_element_mask)

    call this%interface_interpolator%find_points(this%x_interface_dof%x, &
         this%y_interface_dof%x, this%z_interface_dof%x, &
         this%x_interface_dof%size())

  end subroutine setup_interpolator_

end module overset_interface
