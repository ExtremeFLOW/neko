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
!> Common configuration and diagnostics for the Euler IDP solver.
module euler_idp
  use json_module, only : json_file
  use json_utils, only : json_get_or_default
  use num_types, only : rp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, public, parameter :: EULER_IDP_NCOMP = 5

  !> Configuration of the opt-in Euler IDP path.
  type, public :: euler_idp_config_t
     logical :: enabled = .false.
     logical :: low_order_only = .false.
     logical :: relax_density_bounds = .false.
     logical :: limit_internal_energy = .true.
     logical :: limit_entropy = .true.
     real(kind=rp) :: density_bound_relaxation_factor = 1.0_rp
     real(kind=rp) :: internal_energy_floor = 1.0e-12_rp
   contains
     procedure, pass(this) :: init => euler_idp_config_init
     procedure, pass(this) :: valid => euler_idp_config_valid
  end type euler_idp_config_t

  !> Diagnostics shared by the Euler IDP implementation modules.
  type, public :: euler_idp_diagnostics_t
     integer :: stage = 0
     real(kind=rp) :: min_density = huge(1.0_rp)
     real(kind=rp) :: min_internal_energy = huge(1.0_rp)
     real(kind=rp) :: min_pressure = huge(1.0_rp)
     real(kind=rp) :: min_specific_entropy = huge(1.0_rp)
     real(kind=rp) :: max_graph_cfl = 0.0_rp
     real(kind=rp) :: min_convex_weight = 1.0_rp
     real(kind=rp) :: maximum_floor_timestep = huge(1.0_rp)
     logical :: entropy_viscosity_enabled = .false.
     real(kind=rp) :: max_entropy_viscosity = 0.0_rp
     real(kind=rp) :: mean_entropy_viscosity = 0.0_rp
     real(kind=rp) :: entropy_viscosity_conservation(EULER_IDP_NCOMP) = &
          0.0_rp
     real(kind=rp) :: entropy_viscosity_element_compatibility( &
          EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: high_order_conservation(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: low_order_conservation(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: limited_conservation(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: correction_assembly_error(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: correction_global_compatibility(EULER_IDP_NCOMP) = &
          0.0_rp
     real(kind=rp) :: correction_element_compatibility(EULER_IDP_NCOMP) = &
          0.0_rp
     real(kind=rp) :: correction_face_mismatch(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: element_compatibility(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: graph_compatibility(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: reconstruction_residual(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: forward_euler_time = 0.0_rp
     real(kind=rp) :: reconstruction_time = 0.0_rp
     real(kind=rp) :: graph_gs_time = 0.0_rp
     real(kind=rp) :: graph_reduction_time = 0.0_rp
     real(kind=rp) :: max_correction_flux = 0.0_rp
     real(kind=rp) :: rms_correction_flux = 0.0_rp
     real(kind=rp) :: min_limiter = 1.0_rp
     real(kind=rp) :: mean_limiter = 1.0_rp
     real(kind=rp) :: max_limiter = 1.0_rp
     real(kind=rp) :: limited_edge_fraction = 0.0_rp
     real(kind=rp) :: limiter_weight_error = 0.0_rp
     real(kind=rp) :: min_density_lower_bound = huge(1.0_rp)
     real(kind=rp) :: max_density_upper_bound = -huge(1.0_rp)
     real(kind=rp) :: max_density_bound_relaxation = 0.0_rp
     real(kind=rp) :: max_density_lower_violation = 0.0_rp
     real(kind=rp) :: max_density_upper_violation = 0.0_rp
     real(kind=rp) :: max_entropy_lower_violation = 0.0_rp
     integer :: density_limited_edges = 0
     integer :: internal_energy_limited_edges = 0
     integer :: entropy_limited_edges = 0
   contains
     procedure, pass(this) :: reset => euler_idp_diagnostics_reset
  end type euler_idp_diagnostics_t

contains

  !> Read the Euler IDP configuration from a case file.
  subroutine euler_idp_config_init(this, params)
    class(euler_idp_config_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    character(len=*), parameter :: root = "case.numerics.euler_idp."

    this%enabled = .false.
    this%low_order_only = .false.
    this%relax_density_bounds = .false.
    this%limit_internal_energy = .true.
    this%limit_entropy = .true.
    this%density_bound_relaxation_factor = 1.0_rp
    this%internal_energy_floor = 1.0e-12_rp

    if (.not. params%valid_path("case.numerics.euler_idp")) return

    call json_get_or_default(params, root // "enabled", this%enabled, .false.)
    call json_get_or_default(params, root // "low_order_only", &
         this%low_order_only, .false.)
    call json_get_or_default(params, root // "relax_density_bounds", &
         this%relax_density_bounds, .false.)
    call json_get_or_default(params, root // "limit_internal_energy", &
         this%limit_internal_energy, .true.)
    call json_get_or_default(params, root // "limit_entropy", &
         this%limit_entropy, .true.)
    call json_get_or_default(params, &
         root // "density_bound_relaxation_factor", &
         this%density_bound_relaxation_factor, 1.0_rp)
    call json_get_or_default(params, root // "internal_energy_floor", &
         this%internal_energy_floor, 1.0e-12_rp)
  end subroutine euler_idp_config_init

  !> Validate scalar Euler IDP configuration values.
  logical function euler_idp_config_valid(this, message) result(valid)
    class(euler_idp_config_t), intent(in) :: this
    character(len=*), intent(out) :: message

    valid = .false.
    message = ""

    if (.not. ieee_is_finite(this%internal_energy_floor) .or. &
         this%internal_energy_floor .le. 0.0_rp) then
       message = "internal_energy_floor must be finite and positive"
       return
    end if
    if (.not. ieee_is_finite(this%density_bound_relaxation_factor) .or. &
         this%density_bound_relaxation_factor .le. 0.0_rp) then
       message = "density_bound_relaxation_factor must be finite and positive"
       return
    end if
    valid = .true.
  end function euler_idp_config_valid

  !> Reset diagnostics before a Forward Euler map or SSPRK stage.
  subroutine euler_idp_diagnostics_reset(this)
    class(euler_idp_diagnostics_t), intent(inout) :: this

    this%stage = 0
    this%min_density = huge(1.0_rp)
    this%min_internal_energy = huge(1.0_rp)
    this%min_pressure = huge(1.0_rp)
    this%min_specific_entropy = huge(1.0_rp)
    this%max_graph_cfl = 0.0_rp
    this%min_convex_weight = 1.0_rp
    this%maximum_floor_timestep = huge(1.0_rp)
    this%entropy_viscosity_enabled = .false.
    this%max_entropy_viscosity = 0.0_rp
    this%mean_entropy_viscosity = 0.0_rp
    this%entropy_viscosity_conservation = 0.0_rp
    this%entropy_viscosity_element_compatibility = 0.0_rp
    this%high_order_conservation = 0.0_rp
    this%low_order_conservation = 0.0_rp
    this%limited_conservation = 0.0_rp
    this%correction_assembly_error = 0.0_rp
    this%correction_global_compatibility = 0.0_rp
    this%correction_element_compatibility = 0.0_rp
    this%correction_face_mismatch = 0.0_rp
    this%element_compatibility = 0.0_rp
    this%graph_compatibility = 0.0_rp
    this%reconstruction_residual = 0.0_rp
    this%forward_euler_time = 0.0_rp
    this%reconstruction_time = 0.0_rp
    this%graph_gs_time = 0.0_rp
    this%graph_reduction_time = 0.0_rp
    this%max_correction_flux = 0.0_rp
    this%rms_correction_flux = 0.0_rp
    this%min_limiter = 1.0_rp
    this%mean_limiter = 1.0_rp
    this%max_limiter = 1.0_rp
    this%limited_edge_fraction = 0.0_rp
    this%limiter_weight_error = 0.0_rp
    this%min_density_lower_bound = huge(1.0_rp)
    this%max_density_upper_bound = -huge(1.0_rp)
    this%max_density_bound_relaxation = 0.0_rp
    this%max_density_lower_violation = 0.0_rp
    this%max_density_upper_violation = 0.0_rp
    this%max_entropy_lower_violation = 0.0_rp
    this%density_limited_edges = 0
    this%internal_energy_limited_edges = 0
    this%entropy_limited_edges = 0
  end subroutine euler_idp_diagnostics_reset

end module euler_idp
