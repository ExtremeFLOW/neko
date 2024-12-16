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
!> Defines a shear stress boundary condition for a vector field.
!! Maintainer: Timofey Mukha.
module shear_stress
  use num_types, only : rp
  use bc, only : bc_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  use utils, only : neko_error
  use coefs, only : coef_t
  use symmetry, only : symmetry_t
  use neumann, only : neumann_t
  use json_module, only : json_file
  use json_utils, only : json_get
  implicit none
  private

  !> A shear stress boundary condition.
  !! @warning Currently strictly for axis-aligned boundaries.
  type, public, extends(bc_t) :: shear_stress_t
     ! This bc takes care of setting the wall-normal component to zero.
     ! It can be passed to associated bc lists, which take care of masking
     ! changes in residuals and solution increments.
     type(symmetry_t) :: symmetry

     !> Neumann condition for the x direction.
     type(neumann_t) :: neumann_x
     !> Neumann condition for the y direction.
     type(neumann_t) :: neumann_y
     !> Neumann condition for the z direction.
     type(neumann_t) :: neumann_z
   contains
     procedure, pass(this) :: apply_scalar => shear_stress_apply_scalar
     procedure, pass(this) :: apply_vector => shear_stress_apply_vector
     procedure, pass(this) :: apply_scalar_dev => shear_stress_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => shear_stress_apply_vector_dev
     procedure, pass(this) :: init => shear_stress_init
     procedure, pass(this) :: set_stress_scalar => &
          shear_stress_set_stress_scalar
     procedure, pass(this) :: set_stress_array => &
          shear_stress_set_stress_array
     generic :: set_stress => set_stress_scalar, set_stress_array
     procedure, pass(this) :: free => shear_stress_free
     procedure, pass(this) :: finalize => shear_stress_finalize
  end type shear_stress_t

contains

  !> Apply shear stress for a scalar field @a x.
  subroutine shear_stress_apply_scalar(this, x, n, t, tstep)
    class(shear_stress_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k, facet
    ! Store non-linear index
    integer :: idx(4)

    call neko_error("The shear stress bc is not applicable to scalar fields.")

  end subroutine shear_stress_apply_scalar

  !> Boundary condition apply for a generic shear_stress condition
  !! to vectors @a x, @a y and @a z
  subroutine shear_stress_apply_vector(this, x, y, z, n, t, tstep)
    class(shear_stress_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep


    call this%neumann_x%apply_scalar(x, n, t, tstep)
    call this%neumann_y%apply_scalar(y, n, t, tstep)
    call this%neumann_z%apply_scalar(z, n, t, tstep)

  end subroutine shear_stress_apply_vector

  !> Boundary condition apply for a generic shear_stress condition
  !! to a vector @a x (device version)
  subroutine shear_stress_apply_scalar_dev(this, x_d, t, tstep)
    class(shear_stress_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call neko_error("shear_stress bc not implemented on the device")

  end subroutine shear_stress_apply_scalar_dev

  !> Boundary condition apply for a generic shear_stress condition
  !! to vectors @a x, @a y and @a z (device version)
  subroutine shear_stress_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
    class(shear_stress_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call neko_error("shear_stress bc not implemented on the device")

  end subroutine shear_stress_apply_vector_dev

  !> Constructor.
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine shear_stress_init(this, coef, json)
    class(shear_stress_t), target, intent(inout) :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) ::json
    real(kind=rp), allocatable :: x(:)

    call this%init_base(coef)

    ! This bc is set to weak, because the vector_apply pertains to the Neumann
    ! part of the bc. The strong part is currently done by the nested
    ! symmetry bc, which is added separately  to the list of velocity bcs.
    this%strong = .false.

    call json_get(json, 'value', x)

    if (size(x) .ne. 3) then
       call neko_error ("The shear stress vector provided for the shear stress &
            & boundary condition should have 3 components.")
    end if

    call this%symmetry%free()
    call this%symmetry%init_from_components(this%coef)

    call this%neumann_x%free()
    call this%neumann_y%free()
    call this%neumann_z%free()

    call this%neumann_x%init_from_components(this%coef, x(1))
    call this%neumann_y%init_from_components(this%coef, x(2))
    call this%neumann_z%init_from_components(this%coef, x(3))

  end subroutine shear_stress_init

  subroutine shear_stress_finalize(this)
    class(shear_stress_t), target, intent(inout) :: this

    call this%finalize_base()

    call this%symmetry%mark_facets(this%marked_facet)
    call this%symmetry%finalize()


    call this%neumann_x%mark_facets(this%marked_facet)
    call this%neumann_y%mark_facets(this%marked_facet)
    call this%neumann_z%mark_facets(this%marked_facet)

    call this%neumann_x%finalize()
    call this%neumann_y%finalize()
    call this%neumann_z%finalize()

  end subroutine shear_stress_finalize

  !> Set the value of the shear stress vector using 3 scalars.
  subroutine shear_stress_set_stress_scalar(this, tau_x, tau_y, tau_z)
    class(shear_stress_t), intent(inout) :: this
    real(kind=rp), intent(in) :: tau_x
    real(kind=rp), intent(in) :: tau_y
    real(kind=rp), intent(in) :: tau_z

    ! Calls finalize and allocates the flux arrays
    call this%neumann_x%set_flux(tau_x)
    call this%neumann_y%set_flux(tau_y)
    call this%neumann_z%set_flux(tau_z)


  end subroutine shear_stress_set_stress_scalar

  !> Set the shear stress components.
  subroutine shear_stress_set_stress_array(this, tau_x, tau_y, tau_z)
    class(shear_stress_t), intent(inout) :: this
    real(kind=rp), intent(in) :: tau_x(this%msk(0))
    real(kind=rp), intent(in) :: tau_y(this%msk(0))
    real(kind=rp), intent(in) :: tau_z(this%msk(0))

    call this%neumann_x%set_flux(tau_x)
    call this%neumann_y%set_flux(tau_y)
    call this%neumann_z%set_flux(tau_z)

  end subroutine shear_stress_set_stress_array

  !> Destructor.
  subroutine shear_stress_free(this)
    class(shear_stress_t), target, intent(inout) :: this
    call this%free_base
    call this%symmetry%free

    call this%neumann_x%free
    call this%neumann_y%free
    call this%neumann_z%free

  end subroutine shear_stress_free
end module shear_stress
