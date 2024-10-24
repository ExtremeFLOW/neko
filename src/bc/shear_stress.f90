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
module shear_stress
  use num_types
  use bc, only : bc_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  use utils, only : neko_error, nonlinear_index
  use coefs, only : coef_t
  use dirichlet, only : dirichlet_t
  use symmetry, only : symmetry_t
  use neumann, only : neumann_t
  use math, only : copy
  use device_math, only : device_copy
  use tuple, only : tuple_i4_t
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
     procedure, pass(this) :: init_shear_stress => &
       shear_stress_init_shear_stress
     procedure, pass(this) :: set_stress => shear_stress_set_stress
     procedure, pass(this) :: shear_stress_finalize
     procedure, pass(this) :: free => shear_stress_free
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
  !> @param coef The SEM coefficients.
  subroutine shear_stress_init_shear_stress(this, coef)
    class(shear_stress_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    integer :: i, j, l
    type(tuple_i4_t), pointer :: bfp(:)
    real(kind=rp) :: sx, sy, sz
    real(kind=rp), parameter :: TOL = 1d-3
    type(tuple_i4_t) :: bc_facet
    integer :: facet, el

    this%coef => coef

    call this%symmetry%free()
    call this%symmetry%init_base(coef)
    call this%symmetry%mark_facets(this%marked_facet)
    call this%symmetry%finalize()
    call this%symmetry%init(coef)

    call this%neumann_x%init_base(coef)
    call this%neumann_y%init_base(coef)
    call this%neumann_z%init_base(coef)

    call this%neumann_x%mark_facets(this%marked_facet)
    call this%neumann_y%mark_facets(this%marked_facet)
    call this%neumann_z%mark_facets(this%marked_facet)


  end subroutine shear_stress_init_shear_stress

  !> Finalize by allocating the stress arrays and marking the facets for
  !! the bc components.
  subroutine shear_stress_finalize(this, tau_x, tau_y, tau_z)
    class(shear_stress_t), intent(inout) :: this
    real(kind=rp), intent(in) :: tau_x
    real(kind=rp), intent(in) :: tau_y
    real(kind=rp), intent(in) :: tau_z

    ! Calls finalize and allocates the flux arrays
    call this%neumann_x%finalize_neumann(tau_x)
    call this%neumann_y%finalize_neumann(tau_y)
    call this%neumann_z%finalize_neumann(tau_z)


  end subroutine shear_stress_finalize

  !> Set the shear stress components.
  subroutine shear_stress_set_stress(this, tau_x, tau_y, tau_z)
    class(shear_stress_t), intent(inout) :: this
    real(kind=rp), intent(in) :: tau_x(this%msk(0))
    real(kind=rp), intent(in) :: tau_y(this%msk(0))
    real(kind=rp), intent(in) :: tau_z(this%msk(0))

    call this%neumann_x%set_flux(tau_x)
    call this%neumann_y%set_flux(tau_x)
    call this%neumann_z%set_flux(tau_x)

  end subroutine shear_stress_set_stress

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
