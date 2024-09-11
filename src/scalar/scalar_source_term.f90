
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
!> Implements the `scalar_source_term_t` type.
module scalar_source_term
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use scalar_user_source_term, only: scalar_user_source_term_t
  use source_term, only : source_term_wrapper_t, source_term_t, &
       source_term_factory
  use field, only : field_t
  use field_list, only : field_list_t
  use json_utils, only : json_get
  use json_module, only : json_file, json_core, json_value
  use coefs, only : coef_t
  use user_intf, only : user_t
  use utils, only : neko_warning
  implicit none
  private

  !> Wrapper contaning and executing the scalar source terms.
  !! @details
  !! Exists mainly to keep the `scalar_scheme_t` type smaller and also as
  !! placeholder for future optimizations.
  type, public :: scalar_source_term_t
     !> Array of ordinary source terms.
     class(source_term_wrapper_t), allocatable :: source_terms(:)
     !> The right-hand side.
     type(field_t), pointer :: f => null()
   contains
     !> Constructor.
     procedure, pass(this) :: init => scalar_source_term_init
     !> Destructor.
     procedure, pass(this) :: free => scalar_source_term_free
     !> Add all the source terms to the passed right-hand side fields.
     procedure, pass(this) :: compute => scalar_source_term_compute
     !> Append a new source term to the source_terms array.
     procedure, pass(this) :: add_source_term => &
          scalar_source_term_add_source_term
     !> Initialize the user source term.
     procedure, nopass, private :: init_user_source

  end type scalar_source_term_t

contains

  !> Constructor.
  subroutine scalar_source_term_init(this, json, f, coef, user)
    class(scalar_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_t), pointer, intent(in) :: f
    type(coef_t), intent(inout), target :: coef
    type(user_t), intent(in) :: user

    type(field_list_t) :: rhs_fields
    ! Json low-level manipulator.
    type(json_core) :: core
    ! Pointer to the source_terms JSON object and the individual sources.
    type(json_value), pointer :: source_object, source_pointer
    ! Buffer for serializing the json.
    character(len=:), allocatable :: buffer
    ! A single source term as its own json_file.
    type(json_file) :: source_subdict
    ! Source type
    character(len=:), allocatable :: type
    logical :: found
    integer :: n_sources, i

    call this%free()

    this%f => f

    if (json%valid_path('case.scalar.source_terms')) then
       ! We package the fields for the source term to operate on in a field list
       call rhs_fields%init(1)
       call rhs_fields%assign(1, f)

       call json%get_core(core)
       call json%get('case.scalar.source_terms', source_object, found)

       n_sources = core%count(source_object)
       allocate(this%source_terms(n_sources))

       do i = 1, n_sources
          ! Create a new json containing just the subdict for this source.
          call core%get_child(source_object, i, source_pointer, found)
          call core%print_to_string(source_pointer, buffer)
          call source_subdict%load_from_string(buffer)
          call json_get(source_subdict, "type", type)

          ! The user source is treated separately
          if ((trim(type) .eq. "user_vector") .or. &
               (trim(type) .eq. "user_pointwise")) then
             if (source_subdict%valid_path("start_time") .or. &
                  source_subdict%valid_path("end_time")) then
                call neko_warning("The start_time and end_time parameters have&
                     & no effect on the scalar user source term")
             end if

             call init_user_source(this%source_terms(i)%source_term, &
                  rhs_fields, coef, type, user)
          else

             call source_term_factory(this%source_terms(i)%source_term, &
                  source_subdict, rhs_fields, coef)
          end if
       end do
    end if

  end subroutine scalar_source_term_init

  !> Initialize the user source term.
  !! @param source_term The allocatable source term to be initialized to a user.
  !! @param rhs_fields The field list with the right-hand-side.
  !! @param coef The SEM coefs.
  !! @param type The type of the user source term, "user_vector" or
  !! "user_poinwise".
  !! @param user The user type containing the user source term routines.
  subroutine init_user_source(source_term, rhs_fields, coef, type, user)
    class(source_term_t), allocatable, intent(inout) :: source_term
    type(field_list_t) :: rhs_fields
    type(coef_t), intent(inout) :: coef
    character(len=*) :: type
    type(user_t), intent(in) :: user

    allocate(scalar_user_source_term_t::source_term)

    select type (source_term)
      type is (scalar_user_source_term_t)
       call source_term%init_from_components(rhs_fields, coef, type, &
            user%scalar_user_f_vector, &
            user%scalar_user_f)
    end select
  end subroutine init_user_source

  !> Destructor.
  subroutine scalar_source_term_free(this)
    class(scalar_source_term_t), intent(inout) :: this
    integer :: i

    nullify(this%f)

    if (allocated(this%source_terms)) then
       do i = 1, size(this%source_terms)
          call this%source_terms(i)%free()
       end do
       deallocate(this%source_terms)
    end if

  end subroutine scalar_source_term_free

  !> Add new sourceterm to the list.
  !! @param source_term The source term to be added.
  subroutine scalar_source_term_add_source_term(this, source_term)
    class(scalar_source_term_t), intent(inout) :: this
    class(source_term_t), intent(in) :: source_term

    integer :: n_sources

    n_sources = size(this%source_terms)
    allocate(this%source_terms(n_sources + 1))
    this%source_terms(n_sources + 1)%source_term = source_term

  end subroutine scalar_source_term_add_source_term

  !> Add all the source term to the passed right-hand side fields.
  !! @param t The time value.
  !! @param tstep The current time step.
  subroutine scalar_source_term_compute(this, t, tstep)
    class(scalar_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: i

    this%f = 0.0_rp

    ! Add contribution from all source terms.
    if (allocated(this%source_terms)) then
       do i = 1, size(this%source_terms)
          call this%source_terms(i)%source_term%compute(t, tstep)
       end do
    end if

  end subroutine scalar_source_term_compute
end module scalar_source_term
