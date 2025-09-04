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
!> Implements `spalding_t`.
module spalding
  use field, only: field_t
  use num_types, only : rp
  use json_module, only : json_file
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use wall_model, only : wall_model_t
  use field_registry, only : neko_field_registry
  use json_utils, only : json_get_or_default
  use spalding_cpu, only : spalding_compute_cpu
  use spalding_device, only : spalding_compute_device
  use field_math, only: field_invcol3
  use vector, only : vector_t
  use math, only: masked_gather_copy_0
  use device_math, only: device_masked_gather_copy_0
  use scratch_registry, only : neko_scratch_registry

  implicit none
  private

  !> Wall model based on Spalding's law of the wall.
  !! Reference: http://dx.doi.org/10.1115/1.3641728
  type, public, extends(wall_model_t) :: spalding_t
     !> The von Karman coefficient.
     real(kind=rp) :: kappa = 0.41_rp
     !> The log-law intercept.
     real(kind=rp) :: B = 5.2_rp
     !> The kinematic viscosity.
     type(vector_t) :: nu
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => spalding_init
     !> Partial constructor from JSON, meant to work as the first stage of
     !! initialization before the `finalize` call.
     procedure, pass(this) :: partial_init => spalding_partial_init
     !> Finalize the construction using the mask and facet arrays of the bc.
     procedure, pass(this) :: finalize => spalding_finalize
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          spalding_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => spalding_free
     !> Compute the kinematic viscosity at the wall.
     procedure, pass(this) :: compute_nu => spalding_compute_nu
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => spalding_compute
  end type spalding_t

contains
  !> Constructor from JSON.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param json A dictionary with parameters.
  subroutine spalding_init(this, scheme_name, coef, msk, facet, h_index, json)
    class(spalding_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, B

    call json_get_or_default(json, "kappa", kappa, 0.41_rp)
    call json_get_or_default(json, "B", B, 5.2_rp)

    call this%init_from_components(scheme_name, coef, msk, facet, h_index, &
         kappa, B)
  end subroutine spalding_init

  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  subroutine spalding_partial_init(this, coef, json)
    class(spalding_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json

    call this%partial_init_base(coef, json)
    call json_get_or_default(json, "kappa", this%kappa, 0.41_rp)
    call json_get_or_default(json, "B", this%B, 5.2_rp)

  end subroutine spalding_partial_init

  !> Finalize the construction using the mask and facet arrays of the bc.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  subroutine spalding_finalize(this, msk, facet)
    class(spalding_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)

    call this%finalize_base(msk, facet)
    call this%nu%init(this%n_nodes)
  end subroutine spalding_finalize

  !> Constructor from components.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param kappa The von Karman coefficient.
  !! @param B The log-law intercept.
  subroutine spalding_init_from_components(this, scheme_name, coef, msk, &
       facet, h_index, kappa, B)
    class(spalding_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    real(kind=rp), intent(in) :: kappa
    real(kind=rp), intent(in) :: B

    call this%free()
    call this%init_base(scheme_name, coef, msk, facet, h_index)

    this%kappa = kappa
    this%B = B

    call this%nu%init(this%n_nodes)
  end subroutine spalding_init_from_components

  !> Compute the kinematic viscosity vector.
  subroutine spalding_compute_nu(this)
    class(spalding_t), intent(inout) :: this
    type(field_t), pointer :: temp
    integer :: idx

    call neko_scratch_registry%request_field(temp, idx)
    call field_invcol3(temp, this%mu, this%rho)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_masked_gather_copy_0(this%nu%x_d, temp%x_d, this%msk_d, &
            temp%size(), this%nu%size())
    else
       call masked_gather_copy_0(this%nu%x, temp%x, this%msk, temp%size(), &
            this%nu%size())
    end if

    call neko_scratch_registry%relinquish_field(idx)
  end subroutine spalding_compute_nu

  !> Destructor for the spalding_t (base) class.
  subroutine spalding_free(this)
    class(spalding_t), intent(inout) :: this

    call this%free_base()

  end subroutine spalding_free

  !> Compute the wall shear stress.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine spalding_compute(this, t, tstep)
    class(spalding_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu, guess

    call this%compute_nu()

    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call spalding_compute_device(u%x_d, v%x_d, w%x_d, this%ind_r_d, &
            this%ind_s_d, this%ind_t_d, this%ind_e_d, &
            this%n_x%x_d, this%n_y%x_d, this%n_z%x_d, &
            this%nu%x_d, this%h%x_d, &
            this%tau_x%x_d, this%tau_y%x_d, this%tau_z%x_d, &
            this%n_nodes, u%Xh%lx, &
            this%kappa, this%B, tstep)
    else
       call spalding_compute_cpu(u%x, v%x, w%x, &
            this%ind_r, this%ind_s, this%ind_t, this%ind_e, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%h%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, &
            this%kappa, this%B, tstep)
    end if

  end subroutine spalding_compute
end module spalding
