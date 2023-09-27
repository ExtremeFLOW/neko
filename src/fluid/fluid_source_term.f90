
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
!> Implements the `fluid_source_term_t` type.
module fluid_source_term
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use fluid_user_source_term, only: fluid_user_source_term_t
  use source_term, only : source_term_wrapper_t, source_term_t
  use source_term_fctry, only : source_term_factory
  use field, only : field_t
  use field_list, only : field_list_t
  use json_utils, only : json_get
  use json_module, only : json_file, json_core, json_value
  use coefs, only : coef_t
  use user_intf, only : user_t 
  implicit none
  private

  !> Wrapper contaning and executing the fluid source terms.
  !! @details
  !! Exists mainly to keep the `fluid_scheme_t` type smaller and also as
  !! placeholder for future optimizations.
  type, public :: fluid_source_term_t
     !> Array of ordinary source terms.
     class(source_term_wrapper_t), allocatable :: source_terms(:)
     !> X-component of the right-hand side.
     type(field_t), pointer :: f_x => null()
     !> Y-component of the right-hand side.
     type(field_t), pointer :: f_y => null()
     !> Z-component of the right-hand side.
     type(field_t), pointer :: f_z => null()
   contains
     !> Constructor.
     procedure, pass(this) :: init => fluid_source_term_init
     !> Destructor.
     procedure, pass(this) :: free => fluid_source_term_free
     !> Add all the source term to the passed right-hand side fields.
     procedure, pass(this) :: compute => fluid_source_term_compute
     !> Initialize the user source term.
     procedure, nopass, private :: init_user_source

  end type fluid_source_term_t

contains

  !> Costructor.
  subroutine fluid_source_term_init(this, json, f_x, f_y, f_z, coef, user)
    class(fluid_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_t), pointer, intent(in) :: f_x, f_y, f_z
    type(coef_t), intent(inout) :: coef
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
    ! Dummy source strenth values
    real(kind=rp) :: values(3)
    logical :: found
    integer :: n_sources, i

    call this%free()

    this%f_x => f_x
    this%f_y => f_y
    this%f_z => f_z


    if (json%valid_path('case.fluid.source_terms')) then
      ! We package the fields for the source term to operate on in a field list.
       allocate(rhs_fields%fields(3))
       rhs_fields%fields(1)%f => f_x
       rhs_fields%fields(2)%f => f_y
       rhs_fields%fields(3)%f => f_z

       call json%get_core(core)
       call json%get('case.fluid.source_terms', source_object, found)

       n_sources = core%count(source_object)
       allocate(this%source_terms(n_sources))
       write(*,*) "HI"


       do i=1, n_sources
         ! Create a new json containing just the subdict for this source.
          call core%get_child(source_object, i, source_pointer, found)
          call core%print_to_string(source_pointer, buffer)
          call source_subdict%load_from_string(buffer)
          call json_get(source_subdict, "type", type)

          ! The user source is treated separately
          if ((trim(type) .eq. "user_vector") .or. &
              (trim(type) .eq. "user_pointwise")) then
              
              call init_user_source(this%source_terms(i)%source_term, &
                                    rhs_fields, coef, type, user)
          else 
              
              call source_term_factory(this%source_terms(i)%source_term, &
                                       source_subdict, rhs_fields, coef)
          end if
      end do 
    end if
    
  end subroutine fluid_source_term_init

  !> Initialize the user source term.
  !! @param source_term The allocatable source term to be initialized to a user.
  !! @param rhs_fields The field list with the 3 right-hand-side components.
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

    allocate(fluid_user_source_term_t::source_term)

    select type (source_term)
    type is (fluid_user_source_term_t)
      call source_term%init_from_components(rhs_fields, coef, type, &
                                            user%fluid_user_f_vector, &
                                            user%fluid_user_f)
    end select
  end subroutine

  !> Destructor.
  subroutine fluid_source_term_free(this)
    class(fluid_source_term_t), intent(inout) :: this
    integer :: i

    nullify(this%f_x)
    nullify(this%f_y)
    nullify(this%f_z)

    if (allocated(this%source_terms)) then
       do i=1, size(this%source_terms)
          call this%source_terms(i)%free()
       end do
       deallocate(this%source_terms)
    end if

  end subroutine fluid_source_term_free

  !> Add all the source term to the passed right-hand side fields.
  !! @param t The time value.
  !! @param tstep The current time step.
  subroutine fluid_source_term_compute(this, t, tstep)
    class(fluid_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: i, n

    this%f_x = 0.0_rp
    this%f_y = 0.0_rp
    this%f_z = 0.0_rp

    ! Add contribution from all source terms.
    if (allocated(this%source_terms)) then
       do i=1, size(this%source_terms)
          call this%source_terms(i)%source_term%compute(t, tstep)
       end do
    end if

  end subroutine fluid_source_term_compute
end module fluid_source_term
