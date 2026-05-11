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
!> Implements [cai_sagaut_model_ii_t](#cai_sagaut_model_ii::cai_sagaut_model_ii_t).
module cai_sagaut_model_ii
  use field, only : field_t
  use num_types, only : rp
  use json_module, only : json_file
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use wall_model, only : wall_model_t
  use registry, only : neko_registry
  use json_utils, only : json_get_or_lookup, json_get_or_lookup_or_default
  use cai_sagaut_model_ii_cpu, only : cai_sagaut_model_ii_compute_cpu
  use cai_sagaut_model_ii_device, only : cai_sagaut_model_ii_compute_device
  use field_math, only : field_invcol3
  use vector, only : vector_t
  use math, only : masked_gather_copy_0
  use device_math, only : device_masked_gather_copy_0
  use scratch_registry, only : neko_scratch_registry

  implicit none
  private

  !> Explicit wall model based on Model-II from Cai and Sagaut (2021).
  !! @details Reference: https://doi.org/10.1063/5.0048563
  type, public, extends(wall_model_t) :: cai_sagaut_model_ii_t
     !> The von Karman coefficient.
     real(kind=rp) :: kappa = 0.41_rp
     !> The log-law intercept.
     real(kind=rp) :: B = 5.2_rp
     !> Blending exponent.
     real(kind=rp) :: p = 1.138_rp
     !> Blending scale.
     real(kind=rp) :: s = 217.8_rp
     !> The kinematic viscosity.
     type(vector_t) :: nu
     !> The fluid density at the boundary.
     type(vector_t) :: rho_w
   contains
     !> Initialise the wall model from the case file.
     procedure, pass(this) :: init => cai_sagaut_model_ii_init
     !> Partially initialise the wall model from the case file.
     procedure, pass(this) :: partial_init => cai_sagaut_model_ii_partial_init
     !> Finalise allocation of derived data structures.
     procedure, pass(this) :: finalize => cai_sagaut_model_ii_finalize
     !> Initialise the wall model from explicit components.
     procedure, pass(this) :: init_from_components => &
          cai_sagaut_model_ii_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => cai_sagaut_model_ii_free
     !> Compute wall-node viscosity and density samples.
     procedure, pass(this) :: compute_nu_and_rho => &
          cai_sagaut_model_ii_compute_nu_and_rho
     !> Evaluate wall shear stresses.
     procedure, pass(this) :: compute => cai_sagaut_model_ii_compute
  end type cai_sagaut_model_ii_t

contains
  !> Initialise the wall model from the case file.
  !! @param scheme_name The solver scheme name.
  !! @param coef The SEM coefficients.
  !! @param msk The wall-point mask.
  !! @param facet The wall-point facet indices.
  !! @param h_index The GLL index used for wall-model sampling.
  !! @param json The case-file parameters for this wall model.
  subroutine cai_sagaut_model_ii_init(this, scheme_name, coef, msk, facet, &
       h_index, json)
    class(cai_sagaut_model_ii_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, B, p, s

    call json_get_or_lookup(json, "kappa", kappa)
    call json_get_or_lookup(json, "B", B)
    call json_get_or_lookup_or_default(json, "p", p, 1.138_rp)
    call json_get_or_lookup_or_default(json, "s", s, 217.8_rp)

    call this%init_from_components(scheme_name, coef, msk, facet, h_index, &
         kappa, B, p, s)
  end subroutine cai_sagaut_model_ii_init

  !> Partially initialise the wall model from the case file.
  !! @param coef The SEM coefficients.
  !! @param json The case-file parameters for this wall model.
  subroutine cai_sagaut_model_ii_partial_init(this, coef, json)
    class(cai_sagaut_model_ii_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json
    call this%partial_init_base(coef, json)
    call json_get_or_lookup(json, "kappa", this%kappa)
    call json_get_or_lookup(json, "B", this%B)
    call json_get_or_lookup_or_default(json, "p", this%p, 1.138_rp)
    call json_get_or_lookup_or_default(json, "s", this%s, 217.8_rp)
  end subroutine cai_sagaut_model_ii_partial_init

  !> Finalise allocation of derived data structures.
  !! @param msk The wall-point mask.
  !! @param facet The wall-point facet indices.
  subroutine cai_sagaut_model_ii_finalize(this, msk, facet)
    class(cai_sagaut_model_ii_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)

    call this%finalize_base(msk, facet)
    call this%nu%init(this%n_nodes)
    call this%rho_w%init(this%n_nodes)
  end subroutine cai_sagaut_model_ii_finalize

  !> Initialise the wall model from explicit components.
  !! @param scheme_name The solver scheme name.
  !! @param coef The SEM coefficients.
  !! @param msk The wall-point mask.
  !! @param facet The wall-point facet indices.
  !! @param h_index The GLL index used for wall-model sampling.
  !! @param kappa The von Karman coefficient.
  !! @param B The log-law intercept.
  !! @param p The blending exponent.
  !! @param s The blending scale.
  subroutine cai_sagaut_model_ii_init_from_components(this, scheme_name, coef, &
       msk, facet, h_index, kappa, B, p, s)
    class(cai_sagaut_model_ii_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    real(kind=rp), intent(in) :: kappa, B, p, s

    call this%free()
    call this%init_base(scheme_name, coef, msk, facet, h_index)

    this%kappa = kappa
    this%B = B
    this%p = p
    this%s = s

    call this%nu%init(this%n_nodes)
    call this%rho_w%init(this%n_nodes)
  end subroutine cai_sagaut_model_ii_init_from_components

  !> Gather viscosity and density values at the wall-model points.
  subroutine cai_sagaut_model_ii_compute_nu_and_rho(this)
    class(cai_sagaut_model_ii_t), intent(inout) :: this
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
  end subroutine cai_sagaut_model_ii_compute_nu_and_rho

  !> Destructor.
  subroutine cai_sagaut_model_ii_free(this)
    class(cai_sagaut_model_ii_t), intent(inout) :: this

    call this%rho_w%free()
    call this%nu%free()
    call this%free_base()
  end subroutine cai_sagaut_model_ii_free

  !> Evaluate wall shear stresses with the Cai-Sagaut Model-II closure.
  !! @param t The current physical time.
  !! @param tstep The current time-step number.
  subroutine cai_sagaut_model_ii_compute(this, t, tstep)
    class(cai_sagaut_model_ii_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w

    call this%compute_nu_and_rho()

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call cai_sagaut_model_ii_compute_device(u%x_d, v%x_d, w%x_d, &
            this%ind_r_d, this%ind_s_d, this%ind_t_d, this%ind_e_d, &
            this%n_x%x_d, this%n_y%x_d, this%n_z%x_d, this%nu%x_d, &
            this%rho_w%x_d, this%h%x_d, this%tau_x%x_d, this%tau_y%x_d, &
            this%tau_z%x_d, this%n_nodes, u%Xh%lx, this%kappa, this%B, &
            this%p, this%s)
    else
       call cai_sagaut_model_ii_compute_cpu(u%x, v%x, w%x, this%ind_r, &
            this%ind_s, this%ind_t, this%ind_e, this%n_x%x, this%n_y%x, &
            this%n_z%x, this%nu%x, this%rho_w%x, this%h%x, this%tau_x%x, &
            this%tau_y%x, this%tau_z%x, this%n_nodes, u%Xh%lx, u%msh%nelv, &
            this%kappa, this%B, this%p, this%s)
    end if
  end subroutine cai_sagaut_model_ii_compute
end module cai_sagaut_model_ii
