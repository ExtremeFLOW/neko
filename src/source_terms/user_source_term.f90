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
!> Implements the `user_source_term_t` type.
module user_source_term
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use utils, only : neko_error
  use source_term, only : source_term_t
  use json_module, only : json_file
  use field_list, only : field_list_t
  use coefs, only : coef_t
  use device, only : device_map, device_free
  use device_math, only : device_add2
  use field_math, only : field_add2
  use dofmap, only : dofmap_t
  use user_intf, only : user_source_term_intf
  use time_state, only : time_state_t
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> A source term wrapping the user source term routine.
  !! Stores fields that are passed to the user routine so tha the user routine
  !! never touches the actual RHS fields directly.
  !! @warning
  !! The user source term does not support init from JSON and should instead be
  !! directly initialized from components.
  type, public, extends(source_term_t) :: user_source_term_t
     !> The name of the scheme that owns this source term.
     character(len=:), allocatable :: scheme_name
     !> Pointer to the dofmap of the right-hand-side fields.
     type(dofmap_t), pointer :: dof
     !> Field list passed to the user source term routine. The values are then
     !! added to this%fields, i.e. the actual RHS.
     type(field_list_t) :: user_fields
     !> Compute the source term for the entire boundary
     procedure(user_source_term_intf), nopass, pointer :: compute_user_ &
          => null()
   contains
     !> Constructor from JSON (will throw!).
     procedure, pass(this) :: init => user_source_term_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          user_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => user_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => user_source_term_compute
  end type user_source_term_t

contains

  !> Costructor from JSON.
  !! @details
  !! This will throw, as the user source term should be initialized directly
  !! from components.
  subroutine user_source_term_init(this, json, fields, coef, &
       variable_name)
    class(user_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name

    call neko_error("The user fluid source term should be init from components")

  end subroutine user_source_term_init

  !> Costructor from components.
  !! @param fields A list of 3 fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param user_proc The procedure user procedure to compute the source term.
  !! @param scheme_name The name of the scheme that owns this source term.
  subroutine user_source_term_init_from_components(this, fields, coef, &
       user_proc, scheme_name)
    class(user_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    procedure(user_source_term_intf) :: user_proc
    character(len=*), intent(in) :: scheme_name
    integer :: i

    call this%free()
    call this%init_base(fields, coef, 0.0_rp, huge(0.0_rp))

    this%scheme_name = scheme_name
    this%dof => fields%dof(1)

    call this%user_fields%init(3)

    do i = 1, this%fields%size()
       allocate(this%user_fields%items(i)%ptr)
       call this%user_fields%items(i)%ptr%init(this%dof)
    end do

    this%compute_user_ => user_proc
  end subroutine user_source_term_init_from_components

  !> Destructor.
  subroutine user_source_term_free(this)
    class(user_source_term_t), intent(inout) :: this

    call this%user_fields%free()

    if (allocated(this%scheme_name)) deallocate(this%scheme_name)

    nullify(this%compute_user_)
    nullify(this%dof)

    call this%free_base()
  end subroutine user_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time state.
  subroutine user_source_term_compute(this, time)
    class(user_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i

    if (time%t .ge. this%start_time .and. time%t .le. this%end_time) then
       call this%compute_user_(this%scheme_name, this%fields, time)

       do i = 1, this%fields%size()
          call this%user_fields%items(i)%ptr%init(this%dof)
          call field_add2(this%fields%items(i)%ptr, &
               this%user_fields%items(i)%ptr)
       end do
    end if
  end subroutine user_source_term_compute

end module user_source_term
