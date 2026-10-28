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
!
!> Implements `duprat_t`.
module duprat
  use field, only : field_t
  use num_types, only : rp
  use json_module, only : json_file
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use wall_model, only : wall_model_t
  use wall_sampler, only : wall_sampler_t
  use wall_sampler_fctry, only : wall_sampler_factory
  use user_intf, only : user_t
  use registry, only : neko_registry
  use json_utils, only : json_get_or_lookup, json_get_or_lookup_or_default, &
       json_get_or_default, json_get_subdict_or_empty
  use duprat_cpu, only : duprat_compute_cpu, duprat_update_dpds_cpu
  use duprat_device, only : duprat_compute_device, duprat_update_dpds_device
  use operators, only : grad
  use field_math, only : field_invcol3
  use vector, only : vector_t
  use math, only : masked_gather_copy_0, cfill
  use device_math, only : device_masked_gather_copy_0, &
       device_masked_scatter_copy_0, device_cfill
  use scratch_registry, only : neko_scratch_registry
  use logger, only : LOG_SIZE, neko_log
  use utils, only : neko_error, neko_warning
  use comm, only : pe_rank

  implicit none
  private

  !> Wall model of Duprat et al. (2011), which accounts for the streamwise
  !! pressure gradient through the extended velocity scale
  !! \f$ u_{\tau p} = \sqrt{u_\tau^2 + u_p^2} \f$, with the pressure velocity
  !! \f$ u_p = |(\nu / \rho) \partial p / \partial s|^{1/3} \f$, and a
  !! pressure-gradient-dependent mixing-length eddy viscosity. The velocity
  !! profile is obtained by integrating the simplified thin-boundary-layer
  !! equation across the wall layer.
  !! Reference: https://doi.org/10.1063/1.3529358
  !! The pressure gradient is either prescribed as a constant, or evaluated
  !! locally from the pressure field and filtered in time.
  !! @note The device kernels have not yet been validated on GPU hardware.
  type, public, extends(wall_model_t) :: duprat_t
     !> The von Karman coefficient.
     real(kind=rp) :: kappa = 0.41_rp
     !> The exponent of the pressure-gradient term in the eddy viscosity.
     real(kind=rp) :: beta = 0.78_rp
     !> The damping constant of the eddy viscosity.
     real(kind=rp) :: A = 17.0_rp
     !> Whether the pressure gradient is evaluated locally from the pressure
     !! field (true) or prescribed as a constant (false).
     logical :: local_dpds = .false.
     !> The prescribed wall-tangential pressure gradient, positive when
     !! adverse to the local flow direction.
     real(kind=rp) :: dpds_constant = 0.0_rp
     !> Time constant of the low-pass filter of the local pressure gradient.
     real(kind=rp) :: filter_time = 0.0_rp
     !> Time from which the local pressure gradient is taken into account.
     real(kind=rp) :: start_time = 0.0_rp
     !> Time of the previous filter update.
     real(kind=rp) :: t_prev = 0.0_rp
     !> Time step of the previous filter update, -1 before the first one.
     integer :: tstep_prev = -1
     !> The kinematic viscosity.
     type(vector_t) :: nu
     !> The fluid density at the boundary.
     type(vector_t) :: rho_w
     !> The wall-tangential pressure gradient at the wall nodes, filtered in
     !! time in the local mode.
     type(vector_t) :: dpds
     !> The ratio alpha = u_tau^2 / u_tau_p^2 at the wall nodes.
     type(vector_t) :: alpha
     !> Velocity sampled away from the wall.
     type(vector_t) :: u_s, v_s, w_s
     !> Pressure gradient sampled away from the wall (local mode).
     type(vector_t) :: dpx_s, dpy_s, dpz_s
     !> Registry field with dp/ds at the wall nodes, for output.
     type(field_t), pointer :: dpds_field => null()
     !> Registry field with alpha at the wall nodes, for output.
     type(field_t), pointer :: alpha_field => null()
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => duprat_init
     !> Partial constructor from JSON, meant to work as the first stage of
     !! initialization before the `finalize` call.
     procedure, pass(this) :: partial_init => duprat_partial_init
     !> Finalize the construction using the mask and facet arrays of the bc.
     procedure, pass(this) :: finalize => duprat_finalize
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          duprat_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => duprat_free
     !> Compute the kinematic viscosity at the wall.
     procedure, pass(this) :: compute_nu => duprat_compute_nu
     !> Update the filtered local pressure gradient.
     procedure, pass(this) :: update_dpds => duprat_update_dpds
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => duprat_compute
  end type duprat_t

contains
  !> Constructor from JSON.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param json A dictionary with parameters.
  subroutine duprat_init(this, scheme_name, coef, msk, facet, json)
    class(duprat_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, beta, A, dpds_constant, filter_time, start_time
    logical :: local_dpds
    class(wall_sampler_t), allocatable :: sampler

    call duprat_read_json(json, kappa, beta, A, local_dpds, dpds_constant, &
         filter_time, start_time)

    call wall_sampler_factory(sampler, json)
    call this%init_from_components(scheme_name, coef, msk, facet, sampler, &
         kappa, beta, A, local_dpds, dpds_constant, filter_time, start_time)
  end subroutine duprat_init

  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param json A dictionary with parameters.
  subroutine duprat_partial_init(this, coef, scheme_name, json)
    class(duprat_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    character(len=*), intent(in) :: scheme_name
    type(json_file), intent(inout) :: json
    character(len=LOG_SIZE) :: log_buf

    call duprat_device_warning()

    call this%partial_init_base(coef, scheme_name, json)
    call duprat_read_json(json, this%kappa, this%beta, this%A, &
         this%local_dpds, this%dpds_constant, this%filter_time, &
         this%start_time)

    call neko_log%section('Wall model')
    write(log_buf, '(A)') 'Model : Duprat'
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'kappa : ', this%kappa
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'beta : ', this%beta
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'A : ', this%A
    call neko_log%message(log_buf)
    if (this%local_dpds) then
       write(log_buf, '(A)') 'dp/ds : local'
       call neko_log%message(log_buf)
       write(log_buf, '(A, E15.7)') 'filter time : ', this%filter_time
       call neko_log%message(log_buf)
       write(log_buf, '(A, E15.7)') 'start time : ', this%start_time
       call neko_log%message(log_buf)
    else
       write(log_buf, '(A, E15.7)') 'dp/ds (constant) : ', this%dpds_constant
       call neko_log%message(log_buf)
    end if
    call neko_log%end_section()

  end subroutine duprat_partial_init

  !> Read the model parameters from JSON.
  !! @param json A dictionary with parameters.
  !! @param kappa The von Karman coefficient.
  !! @param beta The exponent of the pressure-gradient term.
  !! @param A The damping constant.
  !! @param local_dpds Whether the pressure gradient is evaluated locally.
  !! @param dpds_constant The prescribed wall-tangential pressure gradient.
  !! @param filter_time The time constant of the pressure-gradient filter.
  !! @param start_time The time from which the local gradient is used.
  subroutine duprat_read_json(json, kappa, beta, A, local_dpds, &
       dpds_constant, filter_time, start_time)
    type(json_file), intent(inout) :: json
    real(kind=rp), intent(out) :: kappa, beta, A, dpds_constant
    real(kind=rp), intent(out) :: filter_time, start_time
    logical, intent(out) :: local_dpds
    type(json_file) :: pg_json
    character(len=:), allocatable :: pg_type

    call json_get_or_lookup(json, "kappa", kappa)
    call json_get_or_lookup_or_default(json, "beta", beta, 0.78_rp)
    call json_get_or_lookup_or_default(json, "A", A, 17.0_rp)

    ! The pressure gradient settings. Without them, the model runs with a
    ! zero pressure gradient, i.e. in its equilibrium limit.
    call json_get_subdict_or_empty(json, "pressure_gradient", pg_json)
    call json_get_or_default(pg_json, "type", pg_type, "constant")
    dpds_constant = 0.0_rp
    filter_time = 0.0_rp
    start_time = 0.0_rp
    select case (trim(pg_type))
    case ("constant")
       local_dpds = .false.
       call json_get_or_lookup_or_default(pg_json, "value", dpds_constant, &
            0.0_rp)
    case ("local")
       local_dpds = .true.
       call json_get_or_lookup(pg_json, "filter_time", filter_time)
       call json_get_or_lookup_or_default(pg_json, "start_time", start_time, &
            0.0_rp)
       if (filter_time .lt. 0.0_rp) then
          call neko_error("The pressure_gradient filter_time of the " // &
               "duprat wall model must not be negative")
       end if
    case default
       call neko_error("Unknown pressure_gradient type '" // trim(pg_type) &
            // "' for the duprat wall model. Supported: constant, local")
    end select
    call pg_json%destroy()

  end subroutine duprat_read_json

  !> Finalize the construction using the mask and facet arrays of the bc.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param bc_name The name of the boundary condition.
  !! @param user The user interface.
  subroutine duprat_finalize(this, msk, facet, bc_name, user)
    class(duprat_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    character(len=*), optional, intent(in) :: bc_name
    type(user_t), target, optional, intent(in) :: user

    call this%finalize_base(msk, facet, bc_name, user)
    call duprat_init_vectors(this)
  end subroutine duprat_finalize

  !> Constructor from components.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param sampler The sampling strategy. Ownership is transferred.
  !! @param kappa The von Karman coefficient.
  !! @param beta The exponent of the pressure-gradient term.
  !! @param A The damping constant.
  !! @param local_dpds Whether the pressure gradient is evaluated locally
  !! from the pressure field instead of being prescribed.
  !! @param dpds_constant The prescribed wall-tangential pressure gradient,
  !! positive when adverse to the local flow direction.
  !! @param filter_time The time constant of the low-pass filter of the local
  !! pressure gradient, 0 for no filtering.
  !! @param start_time The time from which the local pressure gradient is
  !! used. Before, the model runs with a zero pressure gradient.
  subroutine duprat_init_from_components(this, scheme_name, coef, msk, &
       facet, sampler, kappa, beta, A, local_dpds, dpds_constant, &
       filter_time, start_time)
    class(duprat_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    class(wall_sampler_t), allocatable, intent(inout) :: sampler
    real(kind=rp), intent(in) :: kappa, beta, A, dpds_constant
    real(kind=rp), intent(in) :: filter_time, start_time
    logical, intent(in) :: local_dpds

    call duprat_device_warning()

    call this%free()
    call this%init_base(scheme_name, coef, msk, facet, sampler)

    this%kappa = kappa
    this%beta = beta
    this%A = A
    this%local_dpds = local_dpds
    this%dpds_constant = dpds_constant
    this%filter_time = filter_time
    this%start_time = start_time

    call duprat_init_vectors(this)
  end subroutine duprat_init_from_components

  !> Allocate the work vectors on the wall nodes.
  subroutine duprat_init_vectors(this)
    class(duprat_t), intent(inout) :: this

    call this%nu%init(this%n_nodes)
    call this%rho_w%init(this%n_nodes)
    call this%dpds%init(this%n_nodes)
    call this%alpha%init(this%n_nodes)
    call this%validate_single_sample()
    call this%u_s%init(this%n_nodes)
    call this%v_s%init(this%n_nodes)
    call this%w_s%init(this%n_nodes)

    ! The filtered local gradient starts from zero, so that the pressure
    ! gradient enters smoothly over about one filter time constant.
    if (this%local_dpds) then
       call duprat_fill(this%dpds, 0.0_rp)
       call this%dpx_s%init(this%n_nodes)
       call this%dpy_s%init(this%n_nodes)
       call this%dpz_s%init(this%n_nodes)
    else
       call duprat_fill(this%dpds, this%dpds_constant)
    end if
    call duprat_fill(this%alpha, 1.0_rp)
    this%tstep_prev = -1

    ! Fields for output of dp/ds and alpha at the wall, shared by all
    ! boundaries using this model.
    call neko_registry%add_field(this%dof, "duprat_dpds", &
         ignore_existing = .true.)
    this%dpds_field => neko_registry%get_field("duprat_dpds")
    call neko_registry%add_field(this%dof, "duprat_alpha", &
         ignore_existing = .true.)
    this%alpha_field => neko_registry%get_field("duprat_alpha")
  end subroutine duprat_init_vectors

  !> Fill a vector on the wall nodes with a constant, on host or device.
  !! @param vec The vector.
  !! @param val The value.
  subroutine duprat_fill(vec, val)
    type(vector_t), intent(inout) :: vec
    real(kind=rp), intent(in) :: val

    if (vec%size() .lt. 1) return
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(vec%x_d, val, vec%size())
    else
       call cfill(vec%x, val, vec%size())
    end if
  end subroutine duprat_fill

  !> Warn that the device kernels are not yet validated on GPU hardware.
  subroutine duprat_device_warning()
    if (NEKO_BCKND_DEVICE .eq. 1 .and. pe_rank .eq. 0) then
       call neko_warning("The duprat wall model GPU implementation " // &
            "is partially untested so far")
    end if
  end subroutine duprat_device_warning

  !> Compute the kinematic viscosity vector.
  subroutine duprat_compute_nu(this)
    class(duprat_t), intent(inout) :: this
    type(field_t), pointer :: temp
    integer :: idx

    call neko_scratch_registry%request_field(temp, idx, .false.)
    call field_invcol3(temp, this%mu, this%rho)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_masked_gather_copy_0(this%nu%x_d, temp%x_d, this%msk_d, &
            temp%size(), this%nu%size())
       call device_masked_gather_copy_0(this%rho_w%x_d, this%rho%x_d, &
            this%msk_d, this%rho%size(), this%rho_w%size())
    else
       call masked_gather_copy_0(this%nu%x, temp%x, this%msk, temp%size(), &
            this%nu%size())
       call masked_gather_copy_0(this%rho_w%x, this%rho%x, this%msk, &
            this%rho%size(), this%rho_w%size())
    end if

    call neko_scratch_registry%relinquish_field(idx)
  end subroutine duprat_compute_nu

  !> Destructor for the duprat_t class.
  subroutine duprat_free(this)
    class(duprat_t), intent(inout) :: this

    call this%nu%free()
    call this%rho_w%free()
    call this%dpds%free()
    call this%alpha%free()
    call this%dpx_s%free()
    call this%dpy_s%free()
    call this%dpz_s%free()
    nullify(this%dpds_field)
    nullify(this%alpha_field)
    call this%u_s%free()
    call this%v_s%free()
    call this%w_s%free()
    call this%free_base()

  end subroutine duprat_free

  !> Compute the wall shear stress.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine duprat_compute(this, t, tstep)
    class(duprat_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w
    integer :: i

    call this%compute_nu()

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    call this%sampler%sample(u, this%u_s)
    call this%sampler%sample(v, this%v_s)
    call this%sampler%sample(w, this%w_s)

    if (this%local_dpds) call this%update_dpds(t, tstep)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call duprat_compute_device(this%u_s%x_d, this%v_s%x_d, this%w_s%x_d, &
            this%n_x%x_d, this%n_y%x_d, this%n_z%x_d, &
            this%nu%x_d, this%rho_w%x_d, this%sampler%h%x_d, this%dpds%x_d, &
            this%tau_x%x_d, this%tau_y%x_d, this%tau_z%x_d, this%alpha%x_d, &
            this%n_nodes, this%kappa, this%beta, this%A, tstep)

       ! Copy dp/ds and alpha to the output fields
       call device_masked_scatter_copy_0(this%dpds_field%x_d, &
            this%dpds%x_d, this%msk_d, this%dpds_field%size(), this%n_nodes)
       call device_masked_scatter_copy_0(this%alpha_field%x_d, &
            this%alpha%x_d, this%msk_d, this%alpha_field%size(), &
            this%n_nodes)
    else
       call duprat_compute_cpu(this%u_s%x, this%v_s%x, this%w_s%x, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%rho_w%x, this%sampler%h%x, this%dpds%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, this%alpha%x, &
            this%n_nodes, this%kappa, this%beta, this%A, tstep)

       ! Copy dp/ds and alpha to the output fields
       do i = 1, this%n_nodes
          this%dpds_field%x(this%msk(i), 1, 1, 1) = this%dpds%x(i)
          this%alpha_field%x(this%msk(i), 1, 1, 1) = this%alpha%x(i)
       end do
    end if

    nullify(u, v, w)

  end subroutine duprat_compute

  !> Update the filtered wall-tangential pressure gradient from the current
  !! pressure field.
  !!
  !! The wall model is evaluated at the beginning of a time step, before the
  !! pressure and velocity solves. The pressure field is therefore the one of
  !! the latest completed step, i.e. the same time level as the sampled
  !! velocity. Its gradient is sampled at the same points as the velocity,
  !! projected on the local flow direction and filtered in time with an
  !! exponential moving average,
  !! (dp/ds)_filt^n = (1 - eps) (dp/ds)_filt^(n-1) + eps (dp/ds)^n,
  !! with eps = 1 - exp(-dt / filter_time); see `duprat_update_dpds_cpu`.
  !! The update is done once per time step. Before `start_time`, the
  !! filtered gradient stays zero.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine duprat_update_dpds(this, t, tstep)
    class(duprat_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: p, dpdx, dpdy, dpdz
    real(kind=rp) :: dt, eps
    integer :: idx(3)

    ! Update once per time step only
    if (tstep .eq. this%tstep_prev) return

    ! No time increment is known at the first call, so only record the time
    if (this%tstep_prev .lt. 0) then
       dt = 0.0_rp
    else
       dt = t - this%t_prev
    end if
    this%tstep_prev = tstep
    this%t_prev = t

    if (t .lt. this%start_time .or. dt .le. 0.0_rp) return

    ! Weight of the current gradient in the exponential moving average,
    ! eps = 1 - exp(-dt / T). This is the exact discrete form of a
    ! first-order low-pass filter with time constant T, also for variable
    ! time steps. T = 0 gives eps = 1, i.e. no filtering.
    if (this%filter_time .gt. 0.0_rp) then
       eps = 1.0_rp - exp(-dt / this%filter_time)
    else
       eps = 1.0_rp
    end if

    ! Physical pressure gradient, sampled at the velocity sampling points
    p => neko_registry%get_field("p")
    call neko_scratch_registry%request_field(dpdx, idx(1), .false.)
    call neko_scratch_registry%request_field(dpdy, idx(2), .false.)
    call neko_scratch_registry%request_field(dpdz, idx(3), .false.)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call grad(dpdx%x_d, dpdy%x_d, dpdz%x_d, p%x_d, this%coef)
    else
       call grad(dpdx%x, dpdy%x, dpdz%x, p%x, this%coef)
    end if
    call this%sampler%sample(dpdx, this%dpx_s)
    call this%sampler%sample(dpdy, this%dpy_s)
    call this%sampler%sample(dpdz, this%dpz_s)

    call neko_scratch_registry%relinquish_field(idx)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call duprat_update_dpds_device(this%dpx_s%x_d, this%dpy_s%x_d, &
            this%dpz_s%x_d, this%u_s%x_d, this%v_s%x_d, this%w_s%x_d, &
            this%n_x%x_d, this%n_y%x_d, this%n_z%x_d, &
            this%nu%x_d, this%rho_w%x_d, this%sampler%h%x_d, eps, &
            this%dpds%x_d, this%n_nodes)
    else
       call duprat_update_dpds_cpu(this%dpx_s%x, this%dpy_s%x, &
            this%dpz_s%x, this%u_s%x, this%v_s%x, this%w_s%x, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%rho_w%x, this%sampler%h%x, eps, this%dpds%x, &
            this%n_nodes)
    end if

    nullify(p, dpdx, dpdy, dpdz)

  end subroutine duprat_update_dpds
end module duprat
