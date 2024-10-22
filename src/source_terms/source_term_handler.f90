
! Copyright (c) 2023, The Neko Authors
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
!> Implements the `source_term_handler_t` type.
module source_term_handler
  use neko_config, only: NEKO_BCKND_DEVICE
  use num_types, only: rp
  use source_term, only: source_term_wrapper_t, source_term_t, &
       source_term_factory
  use field, only: field_t
  use field_list, only: field_list_t
  use json_utils, only: json_get, json_extract_item, json_get_or_default
  use json_module, only: json_file
  use coefs, only: coef_t
  use user_intf, only: user_t
  use utils, only: neko_warning
  use field_math, only: field_rzero
  use math, only : col2
  use device_math, only : device_col2
  implicit none
  private

  !> Abstract class for handling source terms.
  !!
  !! @details
  !! This class is responsible for managing the source terms and adding their
  !! contributions to the right-hand side fields. This serve as a common
  !! interface for fluid, scalar and other source terms.
  !!
  !! In general, a derived class should implement the `init_user_source` method
  !! to handle user-defined source terms and should manually implement an
  !! initializer which calls `this%init_base`.
  type, abstract, public :: source_term_handler_t
     !> Array of ordinary source terms.
     class(source_term_wrapper_t), allocatable :: source_terms(:)
     !> The right-hand side.
     type(field_list_t) :: rhs_fields
     !> The coefficients of the (space, mesh) pair.
     type(coef_t), pointer :: coef
     !> The user object.
     type(user_t), pointer :: user

   contains
     !> Constructor.
     procedure, pass(this) :: init_base => source_term_handler_init_base
     !> Destructor.
     procedure, pass(this) :: free => source_term_handler_free
     !> Add all the source terms to the passed right-hand side fields.
     procedure, pass(this) :: compute => source_term_handler_compute
     !> Generic interface to add a source term to the list.
     generic :: add => add_source_term, add_json_source_terms
     !> Append a new source term to the source_terms array.
     procedure, pass(this) :: add_source_term => &
          source_term_handler_add_source_term
     !> Read from the json file and initialize the source terms.
     procedure, pass(this) :: add_json_source_terms => &
          source_term_handler_add_json_source_terms
     !> Initialize the user source term.
     procedure(source_term_handler_init_user_source), &
          nopass, deferred :: init_user_source
  end type source_term_handler_t

  abstract interface
     subroutine source_term_handler_init_user_source(source_term, rhs_fields, &
          coef, type, user)
       import :: source_term_t, field_list_t, coef_t, user_t
       class(source_term_t), allocatable, intent(inout) :: source_term
       type(field_list_t) :: rhs_fields
       type(coef_t), intent(inout) :: coef
       character(len=*) :: type
       type(user_t), intent(in) :: user
     end subroutine source_term_handler_init_user_source
  end interface

contains

  !> Constructor.
  subroutine source_term_handler_init_base(this, rhs_fields, coef, user)
    class(source_term_handler_t), intent(inout) :: this
    type(field_list_t), intent(in) :: rhs_fields
    type(coef_t), target, intent(inout) :: coef
    type(user_t), target, intent(in) :: user

    call this%free()

    ! We package the fields for the source term to operate on in a field list.
    this%rhs_fields = rhs_fields
    this%coef => coef
    this%user => user

  end subroutine source_term_handler_init_base


  !> Destructor.
  subroutine source_term_handler_free(this)
    class(source_term_handler_t), intent(inout) :: this
    integer :: i

    call this%rhs_fields%free()

    if (allocated(this%source_terms)) then
       do i = 1, size(this%source_terms)
          call this%source_terms(i)%free()
       end do
       deallocate(this%source_terms)
    end if

  end subroutine source_term_handler_free

  !> Add all the source term to the passed right-hand side fields.
  !! @param t The time value.
  !! @param tstep The current time step.
  subroutine source_term_handler_compute(this, t, tstep)
    class(source_term_handler_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: i
    type(field_t), pointer :: f

    do i = 1, this%rhs_fields%size()
       f => this%rhs_fields%get(i)
       call field_rzero(f)
    end do

    ! Add contribution from all source terms. If time permits.
    if (allocated(this%source_terms)) then

       do i = 1, size(this%source_terms)
          call this%source_terms(i)%source_term%compute(t, tstep)
       end do

       ! Multiply by mass matrix
       do i = 1, this%rhs_fields%size()
          f => this%rhs_fields%get(i)
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call device_col2(f%x_d, this%coef%B_d, f%size())
          else
             call col2(f%x, this%coef%B, f%size())
          end if
       end do

    end if

  end subroutine source_term_handler_compute

  !> Read from the json file and initialize the source terms.
  subroutine source_term_handler_add_json_source_terms(this, json, name)
    class(source_term_handler_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name

    class(source_term_wrapper_t), dimension(:), allocatable :: temp

    ! A single source term as its own json_file.
    type(json_file) :: source_subdict
    character(len=:), allocatable :: type
    integer :: n_sources, i, i0

    if (json%valid_path(name)) then
       ! Get the number of source terms.
       call json%info(name, n_children = n_sources)

       if (allocated(this%source_terms)) then
          i0 = size(this%source_terms)
          call move_alloc(this%source_terms, temp)
          allocate(this%source_terms(i0 + n_sources))
          if (allocated(temp)) then
             do i = 1, i0
                call move_alloc(temp(i)%source_term, this%source_terms(i)%source_term)
             end do
          end if
       else
          i0 = 0
          allocate(this%source_terms(n_sources))
       end if

       do i = 1, n_sources
          ! Create a new json containing just the subdict for this source.
          call json_extract_item(json, name, i, source_subdict)
          call json_get(source_subdict, "type", type)

          ! The user source is treated separately
          if ((trim(type) .eq. "user_vector") .or. &
               (trim(type) .eq. "user_pointwise")) then

             call this%init_user_source(this%source_terms(i+ i0)%source_term, &
                  this%rhs_fields, this%coef, type, this%user)

             call json_get_or_default(source_subdict, "start_time", &
                  this%source_terms(i + i0)%source_term%start_time, 0.0_rp)
             call json_get_or_default(source_subdict, "end_time", &
                  this%source_terms(i + i0)%source_term%end_time, huge(0.0_rp))
          else

             call source_term_factory(this%source_terms(i + i0)%source_term, &
                  source_subdict, this%rhs_fields, this%coef)
          end if
       end do
    end if

  end subroutine source_term_handler_add_json_source_terms

  !> Add new source term to the list.
  !! @param source_term The source term to be added.
  subroutine source_term_handler_add_source_term(this, source_term)
    class(source_term_handler_t), intent(inout) :: this
    class(source_term_t), intent(in) :: source_term
    class(source_term_wrapper_t), dimension(:), allocatable :: temp

    integer :: n_sources, i

    if (allocated(this%source_terms)) then
       n_sources = size(this%source_terms)
    else
       n_sources = 0
    end if

    call move_alloc(this%source_terms, temp)
    allocate(this%source_terms(n_sources + 1))

    if (allocated(temp)) then
       do i = 1, n_sources
          call move_alloc(temp(i)%source_term, &
               this%source_terms(i)%source_term)
       end do
    end if

    this%source_terms(n_sources + 1)%source_term = source_term

  end subroutine source_term_handler_add_source_term
end module source_term_handler
