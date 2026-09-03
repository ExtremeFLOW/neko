! Copyright (c) 2024, The Neko Authors
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
!> Implements `rough_log_law_t`.
module rough_log_law
  use field, only: field_t
  use num_types, only : rp
  use json_module, only : json_file
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use vector, only : vector_t
  use wall_model, only : wall_model_t
  use wall_sampler, only : wall_sampler_t
  use wall_sampler_fctry, only : wall_sampler_factory
  use user_intf, only : user_t
  use utils, only : neko_error
  use registry, only : neko_registry
  use json_utils, only : json_get_or_lookup, &
       json_get_or_lookup_or_default
  use rough_log_law_device, only : rough_log_law_compute_device
  use rough_log_law_cpu, only : rough_log_law_compute_cpu
  use scratch_registry, only : neko_scratch_registry
  use math, only: masked_gather_copy_0
  use device_math, only: device_masked_gather_copy_0
  use logger, only : LOG_SIZE, neko_log
  implicit none
  private

  !> Wall model based on the log-law for a rough wall.
  !! The formula defining the law is \f$ u^+ = log(z/z_0)/\kappa + B \f$.
  !! Here, \f$ z \f$ is the wall-normal distance, as per tradition in
  !! atmospheric sciences, where this law is often used.
  type, public, extends(wall_model_t) :: rough_log_law_t

     !> The von Karman coefficient.
     real(kind=rp) :: kappa
     !> The log-law intercept
     real(kind=rp) :: B
     !> The roughness height
     real(kind=rp) :: z0
     ! The fluid density at the boundary
     type(vector_t) :: rho_w
     !> The sampled velocity.
     type(vector_t) :: u_s, v_s, w_s
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => rough_log_law_init
     !> Partial constructor from JSON, meant to work as the first stage of
     !! initialization before the `finalize` call.
     procedure, pass(this) :: partial_init => rough_log_law_partial_init
     !> Finalize the construction using the mask and facet arrays of the bc.
     procedure, pass(this) :: finalize => rough_log_law_finalize
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          rough_log_law_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => rough_log_law_free
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => rough_log_law_compute
     ! Extract fluid properties at the wall (rho)
     procedure, pass(this) :: extract_properties => &
          rough_log_law_extract_properties
  end type rough_log_law_t

contains
  !> Constructor from JSON.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param json A dictionary with parameters.
  subroutine rough_log_law_init(this, scheme_name, coef, msk, facet, &
       json)
    class(rough_log_law_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, B, z0
    class(wall_sampler_t), allocatable :: sampler

    call json_get_or_lookup_or_default(json, "kappa", kappa, 0.4_rp)
    call json_get_or_lookup_or_default(json, "B", B, 0.0_rp)
    call json_get_or_lookup(json, "z0", z0)

    call wall_sampler_factory(sampler, json)
    call this%init_from_components(scheme_name, coef, msk, facet, sampler, &
         kappa, B, z0)
  end subroutine rough_log_law_init

  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  subroutine rough_log_law_partial_init(this, coef, scheme_name, json)
    class(rough_log_law_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    character(len=*), intent(in) :: scheme_name
    type(json_file), intent(inout) :: json
    character(len=LOG_SIZE) :: log_buf

    call this%partial_init_base(coef, scheme_name, json)
    call json_get_or_lookup_or_default(json, "kappa", this%kappa, 0.4_rp)
    call json_get_or_lookup_or_default(json, "B", this%B, 0.0_rp)
    call json_get_or_lookup(json, "z0", this%z0)

    call neko_log%section('Wall model')
    write(log_buf, '(A)') 'Model : Rough log law'
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'kappa : ', this%kappa
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'B : ', this%B
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'z0 : ', this%z0
    call neko_log%message(log_buf)
    call neko_log%end_section()

  end subroutine rough_log_law_partial_init

  !> Finalize the construction using the mask and facet arrays of the bc.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  subroutine rough_log_law_finalize(this, msk, facet, bc_name, user)
    class(rough_log_law_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    character(len=*), optional, intent(in) :: bc_name
    type(user_t), target, optional, intent(in) :: user

    call this%finalize_base(msk, facet, bc_name, user)

    call this%rho_w%init(this%n_nodes)
    call this%validate_single_sample()
    call this%u_s%init(this%n_nodes)
    call this%v_s%init(this%n_nodes)
    call this%w_s%init(this%n_nodes)
    call rough_log_law_validate_sampling_height(this)

  end subroutine rough_log_law_finalize

  !> Extract the values of rho at the boundary.
  subroutine rough_log_law_extract_properties(this)
    class(rough_log_law_t), intent(inout) :: this

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_masked_gather_copy_0(this%rho_w%x_d, this%rho%x_d, &
            this%msk_d, &
            this%rho%size(), this%rho_w%size())
    else
       call masked_gather_copy_0(this%rho_w%x, this%rho%x, this%msk, &
            this%rho%size(), this%rho_w%size())
    end if
  end subroutine rough_log_law_extract_properties

  !> Constructor from components.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param sampler The sampling strategy. Ownership is transferred.
  !! @param kappa The von Karman coefficient.
  !! @param B The log-law intercept.
  !! @param z0 The roughness height.
  subroutine rough_log_law_init_from_components(this, scheme_name, coef, msk, &
       facet, sampler, kappa, B, z0)
    class(rough_log_law_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    class(wall_sampler_t), allocatable, intent(inout) :: sampler
    real(kind=rp), intent(in) :: kappa, B, z0

    call this%free()
    call this%init_base(scheme_name, coef, msk, facet, sampler)

    this%kappa = kappa
    this%B = B
    this%z0 = z0

    call this%rho_w%init(this%n_nodes)
    call this%validate_single_sample()
    call this%u_s%init(this%n_nodes)
    call this%v_s%init(this%n_nodes)
    call this%w_s%init(this%n_nodes)

    call rough_log_law_validate_sampling_height(this)

  end subroutine rough_log_law_init_from_components

  subroutine rough_log_law_validate_sampling_height(this)
    class(rough_log_law_t), intent(in) :: this

    if (any(this%sampler%h%x(1:this%n_nodes) .le. this%z0)) then
       call neko_error("Roughlog WM: Sampling height h must be greater " // &
            "than roughness z0. Increase the sampling height or decrease z0.")
    else if (this%z0 .eq. 0.0_rp) then
       call neko_error("Roughlog WM: Roughness z0 must be greater than 0.")
    end if
  end subroutine rough_log_law_validate_sampling_height

  !> Destructor for the rough_log_law_t (base) class.
  subroutine rough_log_law_free(this)
    class(rough_log_law_t), intent(inout) :: this

    call this%rho_w%free()
    call this%u_s%free()
    call this%v_s%free()
    call this%w_s%free()
    call this%free_base()

  end subroutine rough_log_law_free

  !> Compute the wall shear stress.
  !> @param t The time value.
  !> @param tstep The time iteration.
  subroutine rough_log_law_compute(this, t, tstep)
    class(rough_log_law_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w

    ! Extract boundary values for rho
    call this%extract_properties()

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    call this%sampler%sample(u, this%u_s)
    call this%sampler%sample(v, this%v_s)
    call this%sampler%sample(w, this%w_s)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call rough_log_law_compute_device(this%u_s%x_d, this%v_s%x_d, &
            this%w_s%x_d, &
            this%n_x%x_d, this%n_y%x_d, this%n_z%x_d, &
            this%sampler%h%x_d, this%tau_x%x_d, this%tau_y%x_d, &
            this%tau_z%x_d, this%n_nodes, this%kappa, &
            this%rho_w%x_d, this%B, this%z0, tstep)
    else
       call rough_log_law_compute_cpu(this%u_s%x, this%v_s%x, this%w_s%x, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%sampler%h%x, this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, this%kappa, &
            this%rho_w%x, this%B, this%z0, tstep)
    end if

    nullify(u, v, w)

  end subroutine rough_log_law_compute


end module rough_log_law
