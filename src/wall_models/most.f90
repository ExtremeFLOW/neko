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
  use registry, only : neko_registry
  use json_utils, only : json_get, json_get_or_default, &
       json_get_or_lookup, json_get_or_lookup_or_default
  use most_device, only : most_compute_device
  use most_cpu, only : most_compute_cpu
  use scratch_registry, only : neko_scratch_registry
  use utils, only : neko_error, neko_warning
  use logger, only : LOG_SIZE, neko_log
  use vector, only : vector_t
  use vector_math, only : vector_glsum, vector_glmin, vector_glmax
  implicit none
  private

  !> Wall model based on the Monin-Obukhov Similarity Theory for atmospheric
  !! boundary layer flows. Automatically switches between stable, unstable and
  !! neutral layer formulations based on the Richardson number.
  !!
  type, public, extends(wall_model_t) :: most_t
     !> The von Karman coefficient.
     real(kind=rp) :: kappa
     !> The roughness height
     real(kind=rp) :: z0
     !> The thermal roughness height
     real(kind=rp) :: z0h_in
     !> The gravity vector
     real(kind=rp) :: g(3)
     !> The fluid density
     real(kind=rp) :: rho_val
     !> The fluid dynamic viscosity
     real(kind=rp) :: mu_val
     !> The type of temperature boundary condition set in the case file
     character(len=:), allocatable :: bc_type
     !> The heat flux or temperature value set in the case file
     real(kind=rp) :: bc_value
     !> The name of the temperature variable
     character(len=:), allocatable :: scalar_name
     !> Diagnostics
     type(vector_t) :: Ri_b
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
    real(kind=rp) :: mu_val, rho_val
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, z0, z0h_in
    character(len=:), allocatable :: bc_type
    character(len=:), allocatable :: scalar_name
    real(kind=rp) :: bc_value
    real(kind=rp), allocatable :: g_tmp(:)
    real(kind=rp) :: g(3)

    call json_get_or_lookup_or_default(json, "kappa", kappa, 0.4_rp)
    call json_get_or_lookup_or_default(json, "z0", z0, 0.1_rp)
    ! If z0h is specified and positive, z0h will be constant and equal to
    ! what's specified in the case file.
    ! If z0h is specified and negative, the Zilitinkevich 1995 formulation
    ! is used, with the specified value acting as -C_Zil.
    ! If z0h is to specified, assign default value of -0.8,
    ! corresponding to the Zilitinkevich constant value used in Zilitinkevich 1995.
    call json_get_or_lookup_or_default(json, "z0h", z0h_in, -0.8_rp)
    call json_get_or_lookup_or_default(json, "mu", mu_val, 1e-10_rp)
    call json_get_or_lookup_or_default(json, "rho", rho_val, 1.0_rp)
    call json_get(json, "type_of_temp_bc", bc_type)
    call json_get(json, "scalar_field", scalar_name)
    call json_get_or_lookup(json, "bottom_bc_flux_or_temp", bc_value)

    call json_get_or_lookup(json, "g", g_tmp)
    if (size(g_tmp) == 3) then
       g = g_tmp
    else
       call neko_error("MOST WM: The gravity vector should have exactly 3 components")
    end if
    deallocate(g_tmp)

    call this%init_from_components(scheme_name, scalar_name, coef, msk, facet, h_index, &
         kappa, mu_val, rho_val, g, z0, z0h_in, bc_type, bc_value)
    deallocate(bc_type)
    deallocate(scalar_name)
  end subroutine most_init

  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  subroutine most_partial_init(this, coef, json)
    class(most_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json
    real(kind=rp), allocatable :: g_tmp(:)
    character(len=LOG_SIZE) :: log_buf

    call this%partial_init_base(coef, json)
    call json_get_or_lookup_or_default(json, "kappa", this%kappa, 0.4_rp)
    call json_get_or_lookup_or_default(json, "z0", this%z0, 0.1_rp)
    call json_get_or_lookup_or_default(json, "z0h", this%z0h_in, -0.8_rp)
    call json_get_or_lookup_or_default(json, "mu", this%mu_val, 1e-10_rp)
    call json_get_or_lookup_or_default(json, "rho", this%rho_val, 1.0_rp)
    call json_get(json, "type_of_temp_bc", this%bc_type)
    call json_get(json, "scalar_field", this%scalar_name)
    call json_get_or_lookup(json, "bottom_bc_flux_or_temp", this%bc_value)

    call json_get_or_lookup(json, "g", g_tmp)
    if (size(g_tmp) == 3) then
       this%g = g_tmp
    else
       call neko_error("MOST WM: The gravity vector should have exactly 3 components")
    end if
    deallocate(g_tmp)

    call neko_log%section('Wall model')
    write(log_buf, '(A, A)') 'Model : MOST'
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'scalar_name : ', trim(this%scalar_name)
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'bc_type : ', trim(this%bc_type)
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'bc_value : ', this%bc_value
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'kappa : ', this%kappa
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'z0 : ', this%z0
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'z0h : ', this%z0h_in
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'rho : ', this%rho_val
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'mu : ', this%mu_val
    call neko_log%message(log_buf)
    write(log_buf, '(A, 3(E15.7,1X))') 'g : ', this%g
    call neko_log%message(log_buf)
    call neko_log%end_section()

  end subroutine most_partial_init

  !> Finalize the construction using the mask and facet arrays of the bc.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  subroutine most_finalize(this, msk, facet)
    class(most_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)

    call this%finalize_base(msk, facet)

    call this%Ri_b%init(this%n_nodes)

  end subroutine most_finalize

  !> Constructor from components.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param kappa The von Karman coefficient.
  !! @param rho_val fluid density
  !! @param mu_val fluid dynamic viscosity
  !! @param g The gravity vector.
  !! @param z0 The roughness height.
  !! @param z0h_in The thermal roughness height. If negative, set automatically from Zilitinkevich, 1995.
  !! @param bc_type The type of bc set for temperature in the case file.
  !! @param scalar_name The name of the scalar field (temperature) for MOST.
  !! @param bc_value The heat flux at the surface boundary condition.
  subroutine most_init_from_components(this, scheme_name, scalar_name, coef, msk, &
       facet, h_index, kappa, mu_val, rho_val, g, z0, z0h_in, bc_type, bc_value)
    class(most_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    character(len=*), intent(in) :: bc_type
    character(len=*), intent(in) :: scalar_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    real(kind=rp), intent(in) :: g(3)
    real(kind=rp) :: g_mag, g_dot_n, cos_alpha, max_ang
    integer :: i
    real(kind=rp), intent(in) :: kappa, mu_val, rho_val
    real(kind=rp), intent(in) :: z0, z0h_in, bc_value
    character(len=LOG_SIZE) :: log_buf

    call this%init_base(scheme_name, coef, msk, facet, h_index)

    this%kappa = kappa
    this%g = g
    this%mu_val = mu_val
    this%rho_val = rho_val
    this%z0 = z0
    this%z0h_in = z0h_in
    this%bc_type = bc_type
    this%bc_value = bc_value
    this%scalar_name = scalar_name

    !> Check magnitude of g
    g_mag = sqrt(sum(g**2))
    if (g_mag < 1.0e-6_rp) then
       call neko_error("MOST WM: Gravity magnitude is zero. Check your input configuration.")
    end if

    !> Check alignment across all nodes (handling hills/slopes)
    max_ang = 0.0_rp
    do i = 1, this%n_nodes
       g_dot_n = abs(g(1)*this%n_x%x(i) + g(2)*this%n_y%x(i) + g(3)*this%n_z%x(i))
       cos_alpha = g_dot_n / g_mag
       max_ang = max(max_ang, acos(min(1.0_rp, cos_alpha)))
    end do
    max_ang = max_ang * 180.0_rp / (4.0_rp * atan(1.0_rp))
    if (max_ang > 8.0_rp) then
       write(log_buf, '(A, F6.2, A)') "MOST WM: Significant gravity-normal misalignment (max ", &
            max_ang, " deg). Stability corrections will use projected gravity."
       call neko_warning(trim(log_buf))
    end if

    !> Check sampling height
    if (any(this%h%x(1:this%n_nodes) .le. this%z0)) then
       call neko_error("MOST WM: Sampling height h must be greater than roughness z0.")
    end if

  end subroutine most_init_from_components

  !> Destructor for the most_t (base) class.
  subroutine most_free(this)
    class(most_t), intent(inout) :: this

    if (allocated(this%bc_type)) then
       deallocate(this%bc_type)
    end if

    if (allocated(this%scalar_name)) then
       deallocate(this%scalar_name)
    end if

    call this%free_base()

    call this%Ri_b%free()

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

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")
    temp => neko_registry%get_field(this%scalar_name)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call most_compute_device(u%x_d, v%x_d, w%x_d, temp%x_d, this%ind_r_d, &
            this%ind_s_d, this%ind_t_d, this%ind_e_d, &
            this%n_x%x_d, this%n_y%x_d, this%n_z%x_d, &
            this%h%x_d, this%tau_x%x_d, this%tau_y%x_d, &
            this%tau_z%x_d, this%n_nodes, u%Xh%lx, this%kappa, &
            this%mu_val, this%rho_val, this%g, this%z0, this%z0h_in, this%bc_type, &
            this%bc_value, tstep, this%Ri_b%x_d)
    else
       call most_compute_cpu(u%x, v%x, w%x, temp%x, this%ind_r, this%ind_s, &
            this%ind_t, this%ind_e, this%n_x%x, this%n_y%x, this%n_z%x, &
            this%h%x, this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, this%kappa, &
            this%mu_val, this%rho_val, this%g, this%z0, this%z0h_in, this%bc_type, &
            this%bc_value, tstep, this%Ri_b%x)
    end if

    call most_log_diagnostics(this%Ri_b, this%n_nodes)

  end subroutine most_compute

  subroutine most_log_diagnostics(Ri_b, n_nodes)
    use logger, only : neko_log, LOG_SIZE
    use math, only : glsum, glmin, glmax
    character(len=LOG_SIZE) :: log_buf
    integer, intent(in) :: n_nodes
    type(vector_t), intent(in) :: Ri_b

    call neko_log%section("Wall model diagnostics (MOST)")
    write(log_buf, '(A)') 'mean min max'
    write(log_buf,'(A,3E15.7)') "Ri_b: ",&
    vector_glsum(Ri_b, n_nodes) / n_nodes, &
    vector_glmin(Ri_b, n_nodes), vector_glmax(Ri_b, n_nodes)
    call neko_log%message(trim(log_buf))

    call neko_log%end_section()

end subroutine most_log_diagnostics


end module most
