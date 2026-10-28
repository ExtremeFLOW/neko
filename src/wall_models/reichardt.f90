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
!> Implements `reichardt_t`.
module reichardt
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
  use json_utils, only : json_get_or_lookup, json_get_or_lookup_or_default
  use reichardt_cpu, only : reichardt_compute_cpu
  use reichardt_device, only : reichardt_compute_device
  use field_math, only : field_invcol3
  use vector, only : vector_t
  use math, only : masked_gather_copy_0
  use device_math, only : device_masked_gather_copy_0
  use scratch_registry, only : neko_scratch_registry
  use logger, only : LOG_SIZE, neko_log
  use utils, only : neko_warning
  use comm, only : pe_rank

  implicit none
  private

  !> Wall model based on Reichardt's law of the wall,
  !! \f$ u^+ = \frac{1}{\kappa} \ln(1 + \kappa y^+) + C \left[1 -
  !! e^{-y^+/B_1} - \frac{y^+}{B_1} e^{-y^+/B_2} \right] \f$,
  !! with the original constants \f$ C = 7.8 \f$, \f$ B_1 = 11 \f$ and
  !! \f$ B_2 = 3 \f$ as defaults.
  !! Reference: https://doi.org/10.1002/zamm.19510310704
  !! @note The device kernels have not yet been validated on GPU hardware.
  type, public, extends(wall_model_t) :: reichardt_t
     !> The von Karman coefficient.
     real(kind=rp) :: kappa = 0.41_rp
     !> The amplitude of the exponential correction.
     real(kind=rp) :: C = 7.8_rp
     !> The damping length scale, in wall units.
     real(kind=rp) :: B1 = 11.0_rp
     !> The decay length scale of the second exponential term, in wall units.
     real(kind=rp) :: B2 = 3.0_rp
     !> The kinematic viscosity.
     type(vector_t) :: nu
     !> The fluid density at the boundary.
     type(vector_t) :: rho_w
     !> Velocity sampled away from the wall.
     type(vector_t) :: u_s, v_s, w_s
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => reichardt_init
     !> Partial constructor from JSON, meant to work as the first stage of
     !! initialization before the `finalize` call.
     procedure, pass(this) :: partial_init => reichardt_partial_init
     !> Finalize the construction using the mask and facet arrays of the bc.
     procedure, pass(this) :: finalize => reichardt_finalize
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          reichardt_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => reichardt_free
     !> Compute the kinematic viscosity at the wall.
     procedure, pass(this) :: compute_nu => reichardt_compute_nu
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => reichardt_compute
  end type reichardt_t

contains
  !> Constructor from JSON.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param json A dictionary with parameters.
  subroutine reichardt_init(this, scheme_name, coef, msk, facet, json)
    class(reichardt_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, C, B1, B2
    class(wall_sampler_t), allocatable :: sampler

    call json_get_or_lookup(json, "kappa", kappa)
    call json_get_or_lookup_or_default(json, "C", C, 7.8_rp)
    call json_get_or_lookup_or_default(json, "B1", B1, 11.0_rp)
    call json_get_or_lookup_or_default(json, "B2", B2, 3.0_rp)

    call wall_sampler_factory(sampler, json)
    call this%init_from_components(scheme_name, coef, msk, facet, sampler, &
         kappa, C, B1, B2)
  end subroutine reichardt_init

  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param json A dictionary with parameters.
  subroutine reichardt_partial_init(this, coef, scheme_name, json)
    class(reichardt_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    character(len=*), intent(in) :: scheme_name
    type(json_file), intent(inout) :: json
    character(len=LOG_SIZE) :: log_buf

    call reichardt_device_warning()

    call this%partial_init_base(coef, scheme_name, json)
    call json_get_or_lookup(json, "kappa", this%kappa)
    call json_get_or_lookup_or_default(json, "C", this%C, 7.8_rp)
    call json_get_or_lookup_or_default(json, "B1", this%B1, 11.0_rp)
    call json_get_or_lookup_or_default(json, "B2", this%B2, 3.0_rp)

    call neko_log%section('Wall model')
    write(log_buf, '(A)') 'Model : Reichardt'
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'kappa : ', this%kappa
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'C : ', this%C
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'B1 : ', this%B1
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'B2 : ', this%B2
    call neko_log%message(log_buf)
    call neko_log%end_section()

  end subroutine reichardt_partial_init

  !> Finalize the construction using the mask and facet arrays of the bc.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param bc_name The name of the boundary condition.
  !! @param user The user interface.
  subroutine reichardt_finalize(this, msk, facet, bc_name, user)
    class(reichardt_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    character(len=*), optional, intent(in) :: bc_name
    type(user_t), target, optional, intent(in) :: user

    call this%finalize_base(msk, facet, bc_name, user)
    call this%nu%init(this%n_nodes)
    call this%rho_w%init(this%n_nodes)
    call this%validate_single_sample()
    call this%u_s%init(this%n_nodes)
    call this%v_s%init(this%n_nodes)
    call this%w_s%init(this%n_nodes)
  end subroutine reichardt_finalize

  !> Constructor from components.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param sampler The sampling strategy. Ownership is transferred.
  !! @param kappa The von Karman coefficient.
  !! @param C The amplitude of the exponential correction.
  !! @param B1 The damping length scale, in wall units.
  !! @param B2 The decay length scale of the second exponential term, in
  !! wall units.
  subroutine reichardt_init_from_components(this, scheme_name, coef, msk, &
       facet, sampler, kappa, C, B1, B2)
    class(reichardt_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    class(wall_sampler_t), allocatable, intent(inout) :: sampler
    real(kind=rp), intent(in) :: kappa, C, B1, B2

    call reichardt_device_warning()

    call this%free()
    call this%init_base(scheme_name, coef, msk, facet, sampler)

    this%kappa = kappa
    this%C = C
    this%B1 = B1
    this%B2 = B2

    call this%nu%init(this%n_nodes)
    call this%rho_w%init(this%n_nodes)
    call this%validate_single_sample()
    call this%u_s%init(this%n_nodes)
    call this%v_s%init(this%n_nodes)
    call this%w_s%init(this%n_nodes)
  end subroutine reichardt_init_from_components

  !> Warn that the device kernels are not yet validated on GPU hardware.
  subroutine reichardt_device_warning()
    if (NEKO_BCKND_DEVICE .eq. 1 .and. pe_rank .eq. 0) then
       call neko_warning("The reichardt wall model GPU implementation " // &
            "is partially untested so far")
    end if
  end subroutine reichardt_device_warning

  !> Compute the kinematic viscosity vector.
  subroutine reichardt_compute_nu(this)
    class(reichardt_t), intent(inout) :: this
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
  end subroutine reichardt_compute_nu

  !> Destructor for the reichardt_t class.
  subroutine reichardt_free(this)
    class(reichardt_t), intent(inout) :: this

    call this%nu%free()
    call this%rho_w%free()
    call this%u_s%free()
    call this%v_s%free()
    call this%w_s%free()
    call this%free_base()

  end subroutine reichardt_free

  !> Compute the wall shear stress.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine reichardt_compute(this, t, tstep)
    class(reichardt_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w

    call this%compute_nu()

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    call this%sampler%sample(u, this%u_s)
    call this%sampler%sample(v, this%v_s)
    call this%sampler%sample(w, this%w_s)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call reichardt_compute_device(this%u_s%x_d, this%v_s%x_d, &
            this%w_s%x_d, &
            this%n_x%x_d, this%n_y%x_d, this%n_z%x_d, &
            this%nu%x_d, this%rho_w%x_d, this%sampler%h%x_d, &
            this%tau_x%x_d, this%tau_y%x_d, this%tau_z%x_d, &
            this%n_nodes, this%kappa, this%C, this%B1, this%B2, tstep)
    else
       call reichardt_compute_cpu(this%u_s%x, this%v_s%x, this%w_s%x, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%rho_w%x, this%sampler%h%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, this%kappa, this%C, this%B1, this%B2, tstep)
    end if

    nullify(u, v, w)

  end subroutine reichardt_compute
end module reichardt
