! Copyright (c) 2020-2023, The Neko Authors
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
  use dofmap, only : dofmap_t
  use utils, only : neko_error
  use device_math, only : device_rzero
  use device, only : device_free, device_map
  use source_term
  use field, only : field_t
  use json_module, only : json_file
  use field_list, only : field_list_t
  use coefs, only : coef_t
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> A source-term for the fluid, with procedure pointers pointing to the
  !! actual implementation in the user file.
  !! @details The user source term can be applied either pointiwse or acting
  !! on the whole array in a single call, which is referred to as "vector"
  !! application. 
  type, public, extends(source_term_t) :: fluid_user_source_term_t
     !> x-component of source term.
     type(field_t), pointer :: u
     !> y-component of source term.
     type(field_t), pointer :: v
     !> z-component of source term.
     type(field_t), pointer :: w
     !>
     procedure(source_term_compute_pointwise), nopass, pointer :: eval_pw => null()
     procedure(source_term_compute), nopass, pointer :: eval => null()
   contains
     !> Constructor from JSON (will throw!).
     procedure, pass(this) :: init => fluid_user_source_term_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_compenents => &
       fluid_user_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => fluid_user_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute => fluid_user_source_term_compute
  end type fluid_user_source_term_t

contains

  !> Costructor from JSON.
  !! @details
  !! This will throw, as the user source term should be initialized directly
  !! from components.
  subroutine fluid_user_source_term_init(this, json, fields, coef)
    class(fluid_user_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout) :: coef

    call neko_error("The user fluid source term should be init from components")

  end subroutine fluid_user_source_term_init

  !> Costructor from components.
  subroutine fluid_user_source_term_init_from_components(this, fields, coef, &
    source_term_type, eval, eval_pointwise)
    class(fluid_user_source_term_t), intent(inout) :: this
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout) :: coef
    character(len=*) :: source_term_type
    procedure(source_term_compute), optional :: eval
    procedure(source_term_compute_pointwise), optional :: eval_pointwise

    call this%free()
    call this%init_base(fields, coef)

    this%u => fields%fields(1)%f
    this%v => fields%fields(2)%f
    this%w => fields%fields(3)%f


    if (trim(source_term_type) .eq. 'user_poinwise' .and. &
              present(eval_pointwise)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error('Pointwise source terms not supported on accelerators')
       end if
       this%eval => pointwise_eval_driver
       this%eval_pw => eval_pointwise
    else if (trim(source_term_type) .eq. 'user_vector' .and. &
             present(eval)) then
       this%eval => eval
    else
       call neko_error('Invalid fluid source term '//source_term_type)
    end if
  end subroutine fluid_user_source_term_init_from_components

  !> Destructor.
  subroutine fluid_user_source_term_free(this)
    class(fluid_user_source_term_t), intent(inout) :: this

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%eval)
    nullify(this%eval_pw)

    call this%free_base()
  end subroutine fluid_user_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  !> Computes the source term and adds the result to `fields`.
  subroutine fluid_user_source_term_compute(this, t, tstep)
    class(fluid_user_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call this%eval(this, t, tstep)

  end subroutine fluid_user_source_term_compute

  !> Driver for all pointwise source term evaluatons
  !! @param f The source term.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine pointwise_eval_driver(this, t, tstep)
    class(source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: j, k, l, e
    integer :: jj, kk, ll, ee

    select type (this)
    type is (fluid_user_source_term_t)
      do e = 1, this%u%dof%msh%nelv
         ee = e
         do l = 1, this%u%dof%Xh%lz
            ll = l
            do k = 1, this%u%dof%Xh%ly
               kk = k
               do j = 1, this%u%dof%Xh%lx
                  jj =j
                  call this%eval_pw(this%u%x(j,k,l,e), &
                                    this%v%x(j,k,l,e), &
                                    this%w%x(j,k,l,e), &
                                    jj, kk, ll, ee, t)
               end do
            end do
         end do
      end do
    class default
       call neko_error('Incorrect source type in pointwise eval driver!')
    end select
    
  end subroutine pointwise_eval_driver
  
end module fluid_user_source_term
