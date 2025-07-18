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
!> Implements the `fluid_user_source_term_t` type.
module fluid_user_source_term
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use utils, only : neko_error
  use source_term, only : source_term_t
  use json_module, only : json_file
  use field_list, only : field_list_t
  use coefs, only : coef_t
  use device, only : device_map, device_free
  use device_math, only : device_add2
  use math, only : add2
  use dofmap, only : dofmap_t
  use user_intf, only : user_source_term
  use time_state, only : time_state_t
  use, intrinsic :: iso_c_binding
  implicit none
  private

  public :: fluid_source_compute_vector

  !> A source-term for the fluid, with procedure pointers pointing to the
  !! actual implementation in the user file.
  !! @details The user source term can be applied either pointiwse or acting
  !! on the whole array in a single call, which is referred to as "vector"
  !! application.
  !! @warning
  !! The user source term does not support init from JSON and should instead be
  !! directly initialized from components.
  type, public, extends(source_term_t) :: fluid_user_source_term_t
     !> The name of the scheme that owns this source term.
     character(len=:), allocatable :: scheme_name
     !> Pointer to the dofmap of the right-hand-side fields.
     type(dofmap_t), pointer :: dm
     !> x-component of source term.
     real(kind=rp), allocatable :: u(:, :, :, :)
     !> y-component of source term.
     real(kind=rp), allocatable :: v(:, :, :, :)
     !> z-component of source term.
     real(kind=rp), allocatable :: w(:, :, :, :)

     !> Device pointer for `u`.
     type(c_ptr) :: u_d = C_NULL_PTR
     !> Device pointer for `v`.
     type(c_ptr) :: v_d = C_NULL_PTR
     !> Device pointer for `w`.
     type(c_ptr) :: w_d = C_NULL_PTR
     !> Compute the source term for the entire boundary
     procedure(user_source_term), nopass, pointer :: compute_vector_&
          => null()
   contains
     !> Constructor from JSON (will throw!).
     procedure, pass(this) :: init => fluid_user_source_term_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          fluid_user_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => fluid_user_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => fluid_user_source_term_compute
  end type fluid_user_source_term_t

contains

  !> Costructor from JSON.
  !! @details
  !! This will throw, as the user source term should be initialized directly
  !! from components.
  subroutine fluid_user_source_term_init(this, json, fields, coef, &
       variable_name)
    class(fluid_user_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name

    call neko_error("The user fluid source term should be init from components")

  end subroutine fluid_user_source_term_init

  !> Costructor from components.
  !! @param fields A list of 3 fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param sourc_termtype The type of the user source term, "user_vector" or
  !! "user_pointwise".
  !! @param eval_vector The procedure to vector-compute the source term.
  subroutine fluid_user_source_term_init_from_components(this, fields, coef, &
       source_term_type, eval_vector, scheme_name)
    class(fluid_user_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*) :: source_term_type
    procedure(user_source_term), optional :: eval_vector
    character(len=*), intent(in) :: scheme_name

    call this%free()
    call this%init_base(fields, coef, 0.0_rp, huge(0.0_rp))

    this%scheme_name = scheme_name

    this%dm => fields%dof(1)

    allocate(this%u(this%dm%Xh%lx, this%dm%Xh%ly, this%dm%Xh%lz, &
         this%dm%msh%nelv))
    allocate(this%v(this%dm%Xh%lx, this%dm%Xh%ly, this%dm%Xh%lz, &
         this%dm%msh%nelv))
    allocate(this%w(this%dm%Xh%lx, this%dm%Xh%ly, this%dm%Xh%lz, &
         this%dm%msh%nelv))

    this%u = 0d0
    this%v = 0d0
    this%w = 0d0

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%u, this%u_d, this%dm%size())
       call device_map(this%v, this%v_d, this%dm%size())
       call device_map(this%w, this%w_d, this%dm%size())
    end if

    if (trim(source_term_type) .eq. 'user_vector' .and. &
         present(eval_vector)) then
       this%compute_vector_ => eval_vector
    else
       call neko_error('Invalid fluid source term '//source_term_type)
    end if
  end subroutine fluid_user_source_term_init_from_components

  !> Destructor.
  subroutine fluid_user_source_term_free(this)
    class(fluid_user_source_term_t), intent(inout) :: this

    if (allocated(this%u)) deallocate(this%u)
    if (allocated(this%v)) deallocate(this%v)
    if (allocated(this%w)) deallocate(this%w)

    if (c_associated(this%u_d)) call device_free(this%u_d)
    if (c_associated(this%v_d)) call device_free(this%v_d)
    if (c_associated(this%w_d)) call device_free(this%w_d)

    nullify(this%compute_vector_)
    nullify(this%compute_pw_)
    nullify(this%dm)

    call this%free_base()
  end subroutine fluid_user_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time state.
  subroutine fluid_user_source_term_compute(this, time)
    class(fluid_user_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: n

    if (time.t .ge. this%start_time .and. time.t .le. this%end_time) then
       call this%compute_vector_(this%scheme_name, this%fields, time)
       n = this%fields%item_size(1)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_add2(this%fields%x_d(1), this%u_d, n)
          call device_add2(this%fields%x_d(2), this%v_d, n)
          call device_add2(this%fields%x_d(3), this%w_d, n)
       else
          call add2(this%fields%items(1)%ptr%x, this%u, n)
          call add2(this%fields%items(2)%ptr%x, this%v, n)
          call add2(this%fields%items(3)%ptr%x, this%w, n)
       end if
    end if
  end subroutine fluid_user_source_term_compute

end module fluid_user_source_term
