! Copyright (c) 2020-2025, The Neko Authors
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
!> Defines user dirichlet condition for a scalar field.
module field_dirichlet
  use num_types, only: rp
  use coefs, only: coef_t
  use dirichlet, only: dirichlet_t
  use bc, only: bc_t
  use bc_list, only : bc_list_t
  use utils, only: split_string
  use field, only : field_t
  use field_list, only : field_list_t
  use math, only: masked_copy_0
  use device_math, only: device_masked_copy_0
  use dofmap, only : dofmap_t
  use utils, only: neko_error
  use json_module, only : json_file
  use field_list, only : field_list_t
  use json_utils, only : json_get
  use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
  use time_state, only : time_state_t
  implicit none
  private

  !> User defined dirichlet condition, for which the user can work
  !! with an entire field.
  !! The type stores a separate dummy field `field_bc`, which is passed
  !! to the user routine and can be populated with arbitrary values. The
  !! boundary condition then copy-pastes these values to the actual solution
  !! field using the mask of the boundary condition. So, in the end, only the
  !! relevant boundary values are updated.
  !! @note Would be neat to add another class that contains all three
  !! dirichlet bcs for the velocity, this bc would then implement
  !! apply_vector.
  type, public, extends(bc_t) :: field_dirichlet_t
     !> A dummy field which can be manipulated by the user to set the boundary
     !! values.
     type(field_t) :: field_bc
     !> A field list, which just stores `field_bc`, for convenience.
     type(field_list_t) :: field_list
     !> Function pointer to the user routine performing the update of the values
     !! of the boundary fields.
     procedure(field_dirichlet_update), nopass, pointer :: update => null()
   contains
     !> Constructor.
     procedure, pass(this) :: init => field_dirichlet_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          field_dirichlet_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => field_dirichlet_free
     !> Finalize.
     procedure, pass(this) :: finalize => field_dirichlet_finalize
     !> Apply scalar by performing a masked copy.
     procedure, pass(this) :: apply_scalar => field_dirichlet_apply_scalar
     !> (No-op) Apply vector.
     procedure, pass(this) :: apply_vector => field_dirichlet_apply_vector
     !> (No-op) Apply vector (device).
     procedure, pass(this) :: apply_vector_dev => &
          field_dirichlet_apply_vector_dev
     !> Apply scalar (device).
     procedure, pass(this) :: apply_scalar_dev => &
          field_dirichlet_apply_scalar_dev

  end type field_dirichlet_t

  !> Abstract interface defining a dirichlet condition on a list of fields.
  !! @param fields List of fields that are used to extract values for
  !! field_dirichlet.
  !! @param bc  The bc for which the condition is applied.
  !! @param coef Coef object.
  !! @param time Current time state.
  !! @param which_solver Indicates wether the fields provided come from "fluid"
  !! or "scalar".
  abstract interface
     subroutine field_dirichlet_update(fields, bc, time)
       import rp, field_list_t, bc_t, field_dirichlet_t, time_state_t
       type(field_list_t), intent(inout) :: fields
       type(field_dirichlet_t), intent(in) :: bc
       type(time_state_t), intent(in) :: time
     end subroutine field_dirichlet_update
  end interface

  public :: field_dirichlet_update

contains
  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine field_dirichlet_init(this, coef, json)
    class(field_dirichlet_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) ::json
    character(len=:), allocatable :: field_name

    call json_get(json, "field_name", field_name)
    call this%init_from_components(coef, field_name)

  end subroutine field_dirichlet_init

  !> Constructor from components.
  !! @param[in] coef The SEM coefficients.
  subroutine field_dirichlet_init_from_components(this, coef, field_name)
    class(field_dirichlet_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    character(len=*), intent(in) :: field_name

    call this%init_base(coef)

    call this%field_bc%init(this%dof, field_name)
    call this%field_list%init(1)
    call this%field_list%assign_to_field(1, this%field_bc)
  end subroutine field_dirichlet_init_from_components

  !> Destructor. Currently this%field_bc is being freed in
  !! `fluid_scheme_incompressible::free`
  subroutine field_dirichlet_free(this)
    class(field_dirichlet_t), target, intent(inout) :: this

    call this%free_base
    call this%field_bc%free()
    call this%field_list%free()

    if (associated(this%update)) then
       this%update => null()
    end if

  end subroutine field_dirichlet_free

  !> Apply scalar by performing a masked copy.
  !! @param x Field onto which to copy the values (e.g. u,v,w,p or s).
  !! @param n Size of the array `x`.
  !! @param time The current time state.
  subroutine field_dirichlet_apply_scalar(this, x, n, time, strong)
    class(field_dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then

       if (.not. this%updated) then
          call this%update(this%field_list, this, time)
          this%updated = .true.
       end if

       call masked_copy_0(x, this%field_bc%x, this%msk, n, this%msk(0))
    end if

  end subroutine field_dirichlet_apply_scalar

  !> Apply scalar (device).
  !! @param x_d Device pointer to the field onto which to copy the values.
  !! @param time The current time state.
  !! @param strm Device stream
  subroutine field_dirichlet_apply_scalar_dev(this, x_d, time, strong, strm)
    class(field_dirichlet_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then
       if (.not. this%updated) then
          call this%update(this%field_list, this, time)
          this%updated = .true.
       end if

       if (this%msk(0) .gt. 0) then
          call device_masked_copy_0(x_d, this%field_bc%x_d, this%msk_d, &
               this%field_bc%dof%size(), this%msk(0), strm)
       end if
    end if

  end subroutine field_dirichlet_apply_scalar_dev

  !> (No-op) Apply vector.
  !! @param x x-component of the field onto which to apply the values.
  !! @param y y-component of the field onto which to apply the values.
  !! @param z z-component of the field onto which to apply the values.
  !! @param n Size of the `x`, `y` and `z` arrays.
  !! @param time The current time state.
  subroutine field_dirichlet_apply_vector(this, x, y, z, n, time, strong)
    class(field_dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong

    call neko_error("field_dirichlet cannot apply vector BCs.&
    & Use field_dirichlet_vector instead!")

  end subroutine field_dirichlet_apply_vector

  !> (No-op) Apply vector (device).
  !! @param x x-component of the field onto which to apply the values.
  !! @param y y-component of the field onto which to apply the values.
  !! @param z z-component of the field onto which to apply the values.
  !! @param time The current time state.
  !! @param strm Device stream
  subroutine field_dirichlet_apply_vector_dev(this, x_d, y_d, z_d, time, &
       strong, strm)
    class(field_dirichlet_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    call neko_error("field_dirichlet cannot apply vector BCs.&
    & Use field_dirichlet_vector instead!")

  end subroutine field_dirichlet_apply_vector_dev

  !> Finalize
  subroutine field_dirichlet_finalize(this, only_facets)
    class(field_dirichlet_t), target, intent(inout) :: this
    logical, optional, intent(in) :: only_facets
    logical :: only_facets_

    if (present(only_facets)) then
       only_facets_ = only_facets
    else
       only_facets_ = .false.
    end if

    call this%finalize_base(only_facets_)
  end subroutine field_dirichlet_finalize
end module field_dirichlet
