! Copyright (c) 2025, The Neko Authors
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
module regularization
  use num_types, only : rp
  use json_module, only : json_file
  use field, only: field_t
  use coefs, only: coef_t
  use dofmap, only: dofmap_t
  use time_state, only: time_state_t
  implicit none
  private

  type, abstract, public :: regularization_t
     type(field_t), pointer :: reg_coeff => null()
     type(coef_t), pointer :: coef => null()
     type(dofmap_t), pointer :: dof => null()
   contains
     procedure, pass(this) :: init_base => regularization_init_base
     procedure, pass(this) :: free_base => regularization_free_base
     procedure(reg_init), pass(this), deferred :: init
     procedure(reg_free), pass(this), deferred :: free
     procedure(reg_compute), pass(this), deferred :: compute
  end type regularization_t

  abstract interface
     subroutine reg_init(this, json, coef, dof, reg_coeff)
       import regularization_t, json_file, coef_t, dofmap_t, field_t
       class(regularization_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       type(coef_t), intent(in), target :: coef
       type(dofmap_t), intent(in), target :: dof
       type(field_t), intent(in), target :: reg_coeff
     end subroutine reg_init
  end interface

  abstract interface
     subroutine reg_free(this)
       import regularization_t
       class(regularization_t), intent(inout) :: this
     end subroutine reg_free
  end interface

  abstract interface
     subroutine reg_compute(this, time, tstep, dt)
       import regularization_t, rp, time_state_t
       class(regularization_t), intent(inout) :: this
       type(time_state_t), intent(in) :: time
       integer, intent(in) :: tstep
       real(kind=rp), intent(in) :: dt
     end subroutine reg_compute
  end interface

  interface
     module subroutine regularization_factory(object, type_name, json, &
                                              coef, dof, reg_coeff)
       class(regularization_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: type_name
       type(json_file), intent(inout) :: json
       type(coef_t), intent(in), target :: coef
       type(dofmap_t), intent(in), target :: dof
       type(field_t), intent(in), target :: reg_coeff
     end subroutine regularization_factory
  end interface

  public :: regularization_factory

contains

  subroutine regularization_init_base(this, json, coef, dof, reg_coeff)
    class(regularization_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in), target :: coef
    type(dofmap_t), intent(in), target :: dof
    type(field_t), intent(in), target :: reg_coeff

    this%coef => coef
    this%dof => dof
    this%reg_coeff => reg_coeff

  end subroutine regularization_init_base

  subroutine regularization_free_base(this)
    class(regularization_t), intent(inout) :: this

    nullify(this%coef)
    nullify(this%dof)
    nullify(this%reg_coeff)

  end subroutine regularization_free_base

end module regularization

