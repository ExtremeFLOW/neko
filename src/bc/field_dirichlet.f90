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
module field_dirichlet
  use num_types, only: rp
  use coefs, only: coef_t
  use dirichlet, only: dirichlet_t
  use bc, only: bc_list_t, bc_t
  use device, only: c_ptr, c_size_t
  use utils, only: split_string
  use field, only : field_t
  use field_list, only : field_list_t
  use math, only: masked_copy
  use device_math, only: device_masked_copy
  use dofmap, only : dofmap_t
  use utils, only: neko_error
  use field_list, only : field_list_t
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
     procedure, pass(this) :: init_field => field_dirichlet_init
     !> Apply scalar by performing a masked copy.
     procedure, pass(this) :: apply_scalar => field_dirichlet_apply_scalar
     !> (No-op) Apply vector.
     procedure, pass(this) :: apply_vector => field_dirichlet_apply_vector
     !> (No-op) Apply vector (device).
     procedure, pass(this) :: apply_vector_dev => field_dirichlet_apply_vector_dev
     !> Apply scalar (device).
     procedure, pass(this) :: apply_scalar_dev => field_dirichlet_apply_scalar_dev
     !> Destructor
     procedure, pass(this) :: free => field_dirichlet_free

  end type field_dirichlet_t

  !> Abstract interface defining a dirichlet condition on a list of fields.
  !! @param field_bc_list List of fields that are used to extract values for
  !! field_dirichlet.
  !! @param dirichlet_bc_list List of BCs containing field_dirichlet_t BCs only.
  !! @param coef Coef object.
  !! @param t Current time.
  !! @param tstep Current time step.
  !! @param which_solver Indicates wether the fields provided come from "fluid"
  !! or "scalar".
  abstract interface
     subroutine field_dirichlet_update(dirichlet_field_list, dirichlet_bc_list, coef, t, tstep, which_solver)
       import rp
       import field_list_t
       import bc_list_t
       import coef_t
       type(field_list_t), intent(inout) :: dirichlet_field_list
       type(bc_list_t), intent(inout) :: dirichlet_bc_list
       type(coef_t), intent(inout) :: coef
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
       character(len=*), intent(in) :: which_solver
     end subroutine field_dirichlet_update
  end interface

  public :: field_dirichlet_update

contains

  !> Initializes this%field_bc.
  !! @param bc_name Name of this%field_bc
  subroutine field_dirichlet_init(this, bc_name)
    class(field_dirichlet_t), target, intent(inout) :: this
    character(len=*), intent(in) :: bc_name

    call this%field_bc%init(this%dof, bc_name)
    call this%field_list%init(1)
    call this%field_list%assign_to_field(1, this%field_bc)

  end subroutine field_dirichlet_init

  !> Destructor. Currently this%field_bc is being freed in `fluid_scheme::free`
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
  !! @param t Time.
  !! @param tstep Time step.
  subroutine field_dirichlet_apply_scalar(this, x, n, t, tstep)
    class(field_dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    if (this%msk(0) .gt. 0) then
       call masked_copy(x, this%field_bc%x, this%msk, n, this%msk(0))
    end if

  end subroutine field_dirichlet_apply_scalar

  !> Apply scalar (device).
  !! @param x_d Device pointer to the field onto which to copy the values.
  !! @param t Time.
  !! @param tstep Time step.
  subroutine field_dirichlet_apply_scalar_dev(this, x_d, t, tstep)
    class(field_dirichlet_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    if (this%msk(0) .gt. 0) then
       call device_masked_copy(x_d, this%field_bc%x_d, this%msk_d, &
            this%field_bc%dof%size(), this%msk(0))
    end if

  end subroutine field_dirichlet_apply_scalar_dev

  !> (No-op) Apply vector.
  !! @param x x-component of the field onto which to apply the values.
  !! @param y y-component of the field onto which to apply the values.
  !! @param z z-component of the field onto which to apply the values.
  !! @param n Size of the `x`, `y` and `z` arrays.
  !! @param t Time.
  !! @param tstep Time step.
  subroutine field_dirichlet_apply_vector(this, x, y, z, n, t, tstep)
    class(field_dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call neko_error("field_dirichlet cannot apply vector BCs.&
&Use field_dirichlet_vector instead!")

  end subroutine field_dirichlet_apply_vector

  !> (No-op) Apply vector (device).
  !! @param x x-component of the field onto which to apply the values.
  !! @param y y-component of the field onto which to apply the values.
  !! @param z z-component of the field onto which to apply the values.
  !! @param t Time.
  !! @param tstep Time step.
  subroutine field_dirichlet_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
    class(field_dirichlet_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call neko_error("field_dirichlet cannot apply vector BCs.&
&Use field_dirichlet_vector instead!")

  end subroutine field_dirichlet_apply_vector_dev

end module field_dirichlet
