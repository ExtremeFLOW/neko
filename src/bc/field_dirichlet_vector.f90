! Copyright (c) 2020-2024, The Neko Authors
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
!> Defines inflow dirichlet conditions
module field_dirichlet_vector
  use num_types, only: rp
  use coefs, only: coef_t
  use dirichlet, only: dirichlet_t
  use bc, only: bc_list_t, bc_t, bc_list_free
  use device, only: c_ptr, c_size_t
  use utils, only: split_string
  use field, only : field_t
  use field_list, only : field_list_t
  use math, only: masked_copy
  use device_math, only: device_masked_copy
  use dofmap, only : dofmap_t
  use field_dirichlet, only: field_dirichlet_t, field_dirichlet_update
  use utils, only: neko_error
  use field_list, only : field_list_t
  implicit none
  private

  !> Extension of the user defined dirichlet condition `field_dirichlet`
  ! for the application on a vector field.
  type, public, extends(bc_t) :: field_dirichlet_vector_t
     ! The bc for the first compoent.
     type(field_dirichlet_t) :: bc_u
     ! The bc for the second compoent.
     type(field_dirichlet_t) :: bc_v
     ! The bc for the third compoent.
     type(field_dirichlet_t) :: bc_w
     !> A field list to store the bcs for passing to various subroutines.
     type(field_list_t) :: field_list
     !> A bc list to store the bcs for passing to various subroutines.
     type(bc_list_t) :: bc_list
     !> Function pointer to the user routine performing the update of the values
     !! of the boundary fields.
     procedure(field_dirichlet_update), nopass, pointer :: update => null()
   contains
     !> Initializes this%field_bc.
     procedure, pass(this) :: init_field => field_dirichlet_vector_init
     !> Destructor
     procedure, pass(this) :: free => field_dirichlet_vector_free
     !> Apply scalar by performing a masked copy.
     procedure, pass(this) :: apply_scalar => field_dirichlet_vector_apply_scalar
     !> (No-op) Apply vector.
     procedure, pass(this) :: apply_vector => field_dirichlet_vector_apply_vector
     !> (No-op) Apply vector (device).
     procedure, pass(this) :: apply_vector_dev => &
          field_dirichlet_vector_apply_vector_dev
     !> Apply scalar (device).
     procedure, pass(this) :: apply_scalar_dev => &
          field_dirichlet_vector_apply_scalar_dev
  end type field_dirichlet_vector_t

contains

  !> Initializes this%field_bc.
  subroutine field_dirichlet_vector_init(this, bc_name)
    class(field_dirichlet_vector_t), intent(inout) :: this
    character(len=*), intent(in) :: bc_name

    call neko_error("Fields must be initialized individually!")

  end subroutine field_dirichlet_vector_init

  !> Destructor. Currently unused as is, all field_dirichlet attributes
  !! are freed in `fluid_scheme::free`.
  subroutine field_dirichlet_vector_free(this)
    class(field_dirichlet_vector_t), target, intent(inout) :: this

    call this%bc_u%free()
    call this%bc_v%free()
    call this%bc_w%free()

    call this%field_list%free()
    call bc_list_free(this%bc_list)

    if (associated(this%update)) then
       nullify(this%update)
    end if


  end subroutine field_dirichlet_vector_free

  !> Apply scalar by performing a masked copy.
  !! @param x Field onto which to copy the values (e.g. u,v,w,p or s).
  !! @param n Size of the array `x`.
  !! @param t Time.
  !! @param tstep Time step.
  subroutine field_dirichlet_vector_apply_scalar(this, x, n, t, tstep)
    class(field_dirichlet_vector_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call neko_error("field_dirichlet_vector cannot apply scalar BCs.&
&Use field_dirichlet instead!")

  end subroutine field_dirichlet_vector_apply_scalar

  !> Apply scalar (device).
  !! @param x_d Device pointer to the field onto which to copy the values.
  !! @param t Time.
  !! @param tstep Time step.
  subroutine field_dirichlet_vector_apply_scalar_dev(this, x_d, t, tstep)
    class(field_dirichlet_vector_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call neko_error("field_dirichlet_vector cannot apply scalar BCs.&
&Use field_dirichlet instead!")

  end subroutine field_dirichlet_vector_apply_scalar_dev

  !> (No-op) Apply vector.
  !! @param x x-component of the field onto which to apply the values.
  !! @param y y-component of the field onto which to apply the values.
  !! @param z z-component of the field onto which to apply the values.
  !! @param n Size of the `x`, `y` and `z` arrays.
  !! @param t Time.
  !! @param tstep Time step.
  subroutine field_dirichlet_vector_apply_vector(this, x, y, z, n, t, tstep)
    class(field_dirichlet_vector_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    if (present(t) .and. present(tstep)) then
       call this%bc_u%apply_scalar(x, n, t, tstep)
       call this%bc_v%apply_scalar(y, n, t, tstep)
       call this%bc_w%apply_scalar(z, n, t, tstep)
    else if (present(t)) then
       call this%bc_u%apply_scalar(x, n, t=t)
       call this%bc_v%apply_scalar(y, n, t=t)
       call this%bc_w%apply_scalar(z, n, t=t)
    else if (present(tstep)) then
       call this%bc_u%apply_scalar(x, n, tstep=tstep)
       call this%bc_v%apply_scalar(y, n, tstep=tstep)
       call this%bc_w%apply_scalar(z, n, tstep=tstep)
    end if

  end subroutine field_dirichlet_vector_apply_vector

  !> (No-op) Apply vector (device).
  !! @param x x-component of the field onto which to apply the values.
  !! @param y y-component of the field onto which to apply the values.
  !! @param z z-component of the field onto which to apply the values.
  !! @param t Time.
  !! @param tstep Time step.
  subroutine field_dirichlet_vector_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
    class(field_dirichlet_vector_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    if (present(t) .and. present(tstep)) then
       call this%bc_u%apply_scalar_dev(x_d, t, tstep)
       call this%bc_v%apply_scalar_dev(y_d, t, tstep)
       call this%bc_w%apply_scalar_dev(z_d, t, tstep)
    else if (present(t)) then
       call this%bc_u%apply_scalar_dev(x_d, t=t)
       call this%bc_v%apply_scalar_dev(y_d, t=t)
       call this%bc_w%apply_scalar_dev(z_d, t=t)
    else if (present(tstep)) then
       call this%bc_u%apply_scalar_dev(x_d, tstep=tstep)
       call this%bc_v%apply_scalar_dev(y_d, tstep=tstep)
       call this%bc_w%apply_scalar_dev(z_d, tstep=tstep)
    end if

   end subroutine field_dirichlet_vector_apply_vector_dev

end module field_dirichlet_vector
