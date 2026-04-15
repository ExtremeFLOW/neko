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
  use wall_model, only : wall_model_t
  use utils, only : neko_error
  use registry, only : neko_registry
  use json_utils, only : json_get_or_lookup, json_get_or_lookup_or_default
  use rough_log_law_device, only : rough_log_law_compute_device
  use rough_log_law_cpu, only : rough_log_law_compute_cpu
  use scratch_registry, only : neko_scratch_registry
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
     !> The fluid density
     real(kind=rp) :: rho_val
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
  end type rough_log_law_t

contains
  !> Constructor from JSON.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param json A dictionary with parameters.
  subroutine rough_log_law_init(this, scheme_name, coef, msk, facet, &
       h_index, json)
    class(rough_log_law_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, B, z0, rho_val

    call json_get_or_lookup(json, "kappa", kappa)
    call json_get_or_lookup(json, "B", B)
    call json_get_or_lookup(json, "z0", z0)
    call json_get_or_lookup(json, "rho", rho_val)

    call this%init_from_components(scheme_name, coef, msk, facet, h_index, &
         kappa, rho_val, B, z0)
  end subroutine rough_log_law_init

  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  subroutine rough_log_law_partial_init(this, coef, json)
    class(rough_log_law_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json
    character(len=LOG_SIZE) :: log_buf

    call this%partial_init_base(coef, json)
    call json_get_or_lookup(json, "kappa", this%kappa)
    call json_get_or_lookup(json, "B", this%B)
    call json_get_or_lookup(json, "z0", this%z0)
    call json_get_or_lookup_or_default(json, "rho", this%rho_val, 1.0_rp)

    call neko_log%section('Wall model')
    write(log_buf, '(A)') 'Model : Rough log law'
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'kappa : ', this%kappa
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'B : ', this%B
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'z0 : ', this%z0
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'rho : ', this%rho_val
    call neko_log%message(log_buf)
    call neko_log%end_section()

  end subroutine rough_log_law_partial_init

  !> Finalize the construction using the mask and facet arrays of the bc.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  subroutine rough_log_law_finalize(this, msk, facet)
    class(rough_log_law_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)

    call this%finalize_base(msk, facet)

  end subroutine rough_log_law_finalize

  !> Constructor from components.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param kappa The von Karman coefficient.
  !! @param rho_val fluid density
  !! @param B The log-law intercept.
  !! @param z0 The roughness height.
  subroutine rough_log_law_init_from_components(this, scheme_name, coef, msk, &
       facet, h_index, kappa, rho_val, B, z0)
    class(rough_log_law_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    real(kind=rp), intent(in) :: kappa, B, z0, rho_val

    call this%init_base(scheme_name, coef, msk, facet, h_index)

    this%kappa = kappa
    this%rho_val = rho_val
    this%B = B
    this%z0 = z0

    ! Check that the sampling height is above the roughness length
    if (any(this%h%x(1:this%n_nodes) .le. this%z0)) then
       call neko_error("Roughlog WM: Sampling height h must be greater than roughness z0. " // &
            "Increase h_index or decrease z0.")
    else if (this%z0 .eq. 0.0_rp) then
       call neko_error("Roughlog WM: Roughness z0 must be greater than 0.")
    end if

  end subroutine rough_log_law_init_from_components

  !> Destructor for the rough_log_law_t (base) class.
  subroutine rough_log_law_free(this)
    class(rough_log_law_t), intent(inout) :: this

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

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call rough_log_law_compute_device(u%x_d, v%x_d, w%x_d, this%ind_r_d, &
            this%ind_s_d, this%ind_t_d, this%ind_e_d, &
            this%n_x%x_d, this%n_y%x_d, this%n_z%x_d, &
            this%h%x_d, this%tau_x%x_d, this%tau_y%x_d, &
            this%tau_z%x_d, this%n_nodes, u%Xh%lx, this%kappa, &
            this%rho_val, this%B, this%z0, tstep)
    else
       call rough_log_law_compute_cpu(u%x, v%x, w%x, this%ind_r, this%ind_s, &
            this%ind_t, this%ind_e, this%n_x%x, this%n_y%x, this%n_z%x, &
            this%h%x, this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, this%kappa, &
            this%rho_val, this%B, this%z0, tstep)
    end if

  end subroutine rough_log_law_compute


end module rough_log_law
