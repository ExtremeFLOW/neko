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
!> Implements `most_t`.
module most
  use field, only: field_t
  use num_types, only : rp
  use json_module, only : json_file
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use wall_model, only : wall_model_t
  use registry, only : neko_registry, neko_const_registry
  use json_utils, only : json_get_or_default, json_get
  use most_device, only : most_compute_device
  use most_cpu, only : most_compute_cpu
  use scratch_registry, only : neko_scratch_registry
  use utils, only : neko_error
  implicit none
  private

  !> Wall model based on the Monin-Obukhov Similarity Theory for atmospheric
  !! boundary layer flows. Automatically switches between stable, unstable and
  !! neutral layer formulations based on the Richardson number.
  !!
  type, public, extends(wall_model_t) :: most_t
     !> The von Karman coefficient.
     real(kind=rp) :: kappa = 0.41_rp
     !> The roughness height
     real(kind=rp) :: z0 = 0.1_rp
     !> The thermal roughness height
     real(kind=rp) :: z0h_in = 0.1_rp
     !> The type of temperature boundary condition set in the case file
     character(len=:), allocatable :: bc_type
     !> the heat flux or temperature value set in the case file
     real(kind=rp) :: bc_value
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => most_init
     !> Partial constructor from JSON, meant to work as the first stage of
     !! initialization before the `finalize` call.
     procedure, pass(this) :: partial_init => most_partial_init
     !> Finalize the construction using the mask and facet arrays of the bc.
     procedure, pass(this) :: finalize => most_finalize
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          most_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => most_free
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => most_compute
  end type most_t

contains
  !> Constructor from JSON.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param json A dictionary with parameters.
  subroutine most_init(this, scheme_name, coef, msk, facet, &
       h_index, json)
    class(most_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, z0, z0h_in
    character(len=:), allocatable :: bc_type
    real(kind=rp) :: bc_value
    real(kind=rp), pointer :: bc_value

    call json_get_or_default(json, "kappa", kappa, 0.4_rp)
    call json_get_or_default(json, "z0", z0, 0.1_rp)
    call json_get_or_default(json, "z0h", z0h_in, -10.0_rp)
    ! if z0h not specified, assign negative
    ! and compute automatically with Zilitinkevich
    call json_get(json, "type_of_temp_bc", bc_type)

    call neko_const_registry%add_real_scalar(this%bc_value, "bc_value")

    bc_value => neko_const_registry%get_real_scalar("bc_value")

    this%bc_value = bc_type

    call this%init_from_components(scheme_name, coef, msk, facet, h_index, &
         kappa, z0, z0h_in, bc_type, bc_value)
  end subroutine most_init

  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  subroutine most_partial_init(this, coef, json)
    class(most_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json
    real(kind=rp), pointer :: bc_value

    call this%partial_init_base(coef, json)
    call json_get_or_default(json, "kappa", this%kappa, 0.4_rp)
    call json_get_or_default(json, "z0", this%z0, 0.1_rp)
    call json_get_or_default(json, "z0h", this%z0h_in, -10.0_rp)
    call json_get(json, "type_of_temp_bc", this%bc_type)
    ! call json_get(json, "bottom_bc_flux_or_temp", this%bc_value)

    call neko_const_registry%add_real_scalar(this%bc_value, "bc_value")

    bc_value => neko_const_registry%get_real_scalar("bc_value")

    this%bc_value = bc_value

  end subroutine most_partial_init

  !> Finalize the construction using the mask and facet arrays of the bc.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  subroutine most_finalize(this, msk, facet)
    class(most_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)

    call this%finalize_base(msk, facet)

  end subroutine most_finalize

  !> Constructor from components.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param kappa The von Karman coefficient.
  !! @param z0 The roughness height.
  !! @param z0h_in The thermal roughness height. If negative, set automatically from Zilitinkevich, 1995.
  !! @param bc_type The type of bc set for temperature in the case file.
  !! @param bc_value The heat flux at the surface boundary condition.
  subroutine most_init_from_components(this, scheme_name, coef, msk, &
       facet, h_index, kappa, z0, z0h_in, bc_type, bc_value)
    class(most_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    character(len=*), intent(in) :: bc_type
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    real(kind=rp), intent(in) :: kappa
    real(kind=rp), intent(in) :: z0, z0h_in, bc_value

    call this%init_base(scheme_name, coef, msk, facet, h_index)

    this%kappa = kappa
    this%z0 = z0
    this%z0h_in = z0h_in
    this%bc_type = bc_type
    this%bc_value = bc_value

    ! Check that the sampling height is above the roughness length
    if (any(this%h%x(1:this%n_nodes) .le. this%z0)) then
       call neko_error("MOST WM: Sampling height h must be greater than roughness z0. " // &
                       "Increase h_index or decrease z0.")
    end if

  end subroutine most_init_from_components

  !> Destructor for the most_t (base) class.
  subroutine most_free(this)
    class(most_t), intent(inout) :: thi

    if (allocated(this%bc_type)) then
      deallocate(this%bc_type)
    end if

    call this%free_base()

  end subroutine most_free

  !> Compute the wall shear stress.
  !> @param t The time value.
  !> @param tstep The time iteration.
  subroutine most_compute(this, t, tstep)
    class(most_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w
    type(field_t), pointer :: temp
    real(kind=rp), pointer :: bc_value

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")
    temp => neko_registry%get_field("temperature")
    bc_value => neko_const_registry%get_real_scalar("bc_value")
    this%bc_value = bc_value

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call most_compute_device(u%x_d, v%x_d, w%x_d, temp%x_d, this%ind_r_d, &
            this%ind_s_d, this%ind_t_d, this%ind_e_d, &
            this%n_x%x_d, this%n_y%x_d, this%n_z%x_d, &
            this%h%x_d, this%tau_x%x_d, this%tau_y%x_d, &
            this%tau_z%x_d, this%n_nodes, u%Xh%lx, this%kappa, &
            this%z0, this%z0h_in, this%bc_type, &
            this%bc_value, tstep)
    else
       call most_compute_cpu(u%x, v%x, w%x, temp%x, this%ind_r, this%ind_s, &
            this%ind_t, this%ind_e, this%n_x%x, this%n_y%x, this%n_z%x, &
            this%h%x, this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, this%kappa, &
            this%z0, this%z0h_in, this%bc_type, &
            this%bc_value, tstep)
    end if

  end subroutine most_compute

end module most
