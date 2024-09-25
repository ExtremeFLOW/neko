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
!> Implements the `scalar_user_source_term_t` type.
module scalar_user_source_term
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
  use, intrinsic :: iso_c_binding
  implicit none
  private

  public :: scalar_source_compute_pointwise, scalar_source_compute_vector

  !> A source-term for the scalar, with procedure pointers pointing to the
  !! actual implementation in the user file.
  !! @details The user source term can be applied either pointiwse or acting
  !! on the whole array in a single call, which is referred to as "vector"
  !! application.
  !! @warning
  !! The user source term does not support init from JSON and should instead be
  !! directly initialized from components.
  type, public, extends(source_term_t) :: scalar_user_source_term_t
     !> Pointer to the dofmap of the right-hand-side fields.
     type(dofmap_t), pointer :: dm
     !> The source term.
     real(kind=rp), allocatable :: s(:, :, :, :)

     !> Device pointer for `s`.
     type(c_ptr) :: s_d = C_NULL_PTR
     !> Compute the source term for a single point
     procedure(scalar_source_compute_pointwise), nopass, pointer :: &
          compute_pw_ => null()
     !> Compute the source term for the entire boundary
     procedure(scalar_source_compute_vector), nopass, pointer :: &
          compute_vector_ => null()
   contains
     !> Constructor from JSON (will throw!).
     procedure, pass(this) :: init => scalar_user_source_term_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          scalar_user_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => scalar_user_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => scalar_user_source_term_compute
  end type scalar_user_source_term_t

  abstract interface
     !> Computes the source term and adds the result to `fields`.
     !! @param t The time value.
     !! @param tstep The current time-step.
     subroutine scalar_source_compute_vector(this, t)
       import scalar_user_source_term_t, rp
       class(scalar_user_source_term_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t
     end subroutine scalar_source_compute_vector
  end interface

  abstract interface
     !> Computes the source term at a single point.
     !! @param s The source value.
     !! @param j The x-index of GLL point.
     !! @param k The y-index of GLL point.
     !! @param l The z-index of GLL point.
     !! @param e The index of element.
     !! @param t The time value.
     subroutine scalar_source_compute_pointwise(s, j, k, l, e, t)
       import rp
       real(kind=rp), intent(inout) :: s
       integer, intent(in) :: j
       integer, intent(in) :: k
       integer, intent(in) :: l
       integer, intent(in) :: e
       real(kind=rp), intent(in) :: t
     end subroutine scalar_source_compute_pointwise
  end interface

contains

  !> Constructor from JSON.
  !! @details
  !! This will throw, as the user source term should be initialized directly
  !! from components.
  subroutine scalar_user_source_term_init(this, json, fields, coef)
    class(scalar_user_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout), target :: coef

    call neko_error("The user scalar source term &
         &should be init from components")

  end subroutine scalar_user_source_term_init

  !> Constructor from components.
  !! @param fields A list with 1 field for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param sourc_term_type The type of the user source term, "user_vector" or
  !! "user_pointwise".
  !! @param eval_vector The procedure to vector-compute the source term.
  !! @param eval_pointwise The procedure to pointwise-compute the source term.
  subroutine scalar_user_source_term_init_from_components(this, fields, coef, &
       source_term_type, eval_vector, eval_pointwise)
    class(scalar_user_source_term_t), intent(inout) :: this
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout) :: coef
    character(len=*) :: source_term_type
    procedure(scalar_source_compute_vector), optional :: eval_vector
    procedure(scalar_source_compute_pointwise), optional :: eval_pointwise

    call this%free()
    call this%init_base(fields, coef, 0.0_rp, huge(0.0_rp))

    this%dm => fields%dof(1)

    allocate(this%s(this%dm%Xh%lx, this%dm%Xh%ly, this%dm%Xh%lz, &
         this%dm%msh%nelv))

    this%s = 0d0

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%s, this%s_d, this%dm%size())
    end if


    if (trim(source_term_type) .eq. 'user_pointwise' .and. &
         present(eval_pointwise)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error('Pointwise source terms not &
               &supported on accelerators')
       end if
       this%compute_vector_ => pointwise_eval_driver
       this%compute_pw_ => eval_pointwise
    else if (trim(source_term_type) .eq. 'user_vector' .and. &
         present(eval_vector)) then
       this%compute_vector_ => eval_vector
    else
       call neko_error('Invalid fluid source term '//source_term_type)
    end if
  end subroutine scalar_user_source_term_init_from_components

  !> Destructor.
  subroutine scalar_user_source_term_free(this)
    class(scalar_user_source_term_t), intent(inout) :: this

    if (allocated(this%s)) deallocate(this%s)

    if (c_associated(this%s_d)) call device_free(this%s_d)

    nullify(this%compute_vector_)
    nullify(this%compute_pw_)
    nullify(this%dm)

    call this%free_base()
  end subroutine scalar_user_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine scalar_user_source_term_compute(this, t, tstep)
    class(scalar_user_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: n

    if (t .ge. this%start_time .and. t .le. this%end_time) then
       call this%compute_vector_(this, t)
       n = this%fields%item_size(1)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_add2(this%fields%x_d(1), this%s_d, n)
       else
          call add2(this%fields%items(1)%ptr%x, this%s, n)
       end if
    end if
  end subroutine scalar_user_source_term_compute

  !> Driver for all pointwise source term evaluatons.
  !! @param t The time value.
  subroutine pointwise_eval_driver(this, t)
    class(scalar_user_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: j, k, l, e
    integer :: jj, kk, ll, ee

    select type (this)
      type is (scalar_user_source_term_t)
       do e = 1, size(this%s, 4)
          ee = e
          do l = 1, size(this%s, 3)
             ll = l
             do k = 1, size(this%s, 2)
                kk = k
                do j = 1, size(this%s, 1)
                   jj = j
                   call this%compute_pw_(this%s(j,k,l,e), jj, kk, ll, ee, t)
                end do
             end do
          end do
       end do
      class default
       call neko_error('Incorrect source type in pointwise eval driver!')
    end select

  end subroutine pointwise_eval_driver

end module scalar_user_source_term
