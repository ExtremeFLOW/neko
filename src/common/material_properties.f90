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
!> Implements material_properties_t type.
module material_properties
  use num_types, only: rp
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file, json_core, json_value
  use logger, only : neko_log, LOG_SIZE
  use user_intf, only : user_t, dummy_user_material_properties, &
                        user_material_properties
  use utils, only : neko_warning, neko_error
  use comm, only : pe_rank
  implicit none
  private
  
  !> Contains all the material properties necessary in the simulation.
  type, public :: material_properties_t
     !> Density, f$\rho \f$.
     real(kind=rp) :: rho
     !> Dynamic viscosity, f$\mu \f$.
     real(kind=rp) :: mu
     !> Scalar conductivity, f$\lambda \f$.
     real(kind=rp) :: lambda
     !> Scalar specific heat capacity, f$c_p \f$.
     real(kind=rp) :: cp
   contains
     !> Constructor.
     procedure, pass(this) :: init => material_properties_init
     !> Write final dimensional values to the log.
     procedure, private, pass(this) :: write_to_log
   end type material_properties_t

contains

  !> Constructor.
  !! @param params Case file.
  subroutine material_properties_init(this, params, user)
    class(material_properties_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(user_t), target, intent(in) :: user
    character(len=LOG_SIZE) :: log_buf
    ! A local pointer that is needed to make Intel happy
    procedure(user_material_properties),  pointer :: dummy_mp_ptr
    logical :: nondimensional

    call neko_log%section('Material properties')
    dummy_mp_ptr => dummy_user_material_properties
    nondimensional = .false.

    if (.not. associated(user%material_properties, dummy_mp_ptr)) then

          write(log_buf, '(A)') "Material properties must be set in the user&
                              & file!"
          call neko_log%message(log_buf)
          call user%material_properties(0.0_rp, 0, this%rho, this%mu, &
                                        this%cp, this%lambda, params)
    else

       !
       ! Fluid
       !

       ! Incorrect user input
       if (params%valid_path('case.fluid.Re') .and. &
           (params%valid_path('case.fluid.mu') .or. &
            params%valid_path('case.fluid.rho'))) then
           call neko_error("To set the material properties for the fluid,&
              & either provide Re OR mu and rho in the case file.")

       ! Non-dimensional case
       else if (params%valid_path('case.fluid.Re')) then
          nondimensional = .true.

          write(log_buf, '(A)') 'Non-dimensional fluid material properties &
                              & input.'
          call neko_log%message(log_buf, lvl=2)
          write(log_buf, '(A)') 'Density will be set to 1, dynamic viscosity to&
                              & 1/Re.'
          call neko_log%message(log_buf, lvl=2)

          ! Read Re into mu for further manipulation.
          call json_get(params, 'case.fluid.Re', this%mu)
          write(log_buf, '(A)') 'Read non-dimensional values:'
          call neko_log%message(log_buf)
          write(log_buf, '(A,ES13.6)') 'Re         :',  this%mu
          call neko_log%message(log_buf)

          ! Set rho to 1 since the setup is non-dimensional.
          this%rho = 1.0_rp
          ! Invert the Re to get viscosity.
          this%mu = 1.0_rp/this%mu
       ! Dimensional case
       else 
          call json_get(params, 'case.fluid.mu', this%mu)
          call json_get(params, 'case.fluid.rho', this%rho)
       end if

       !
       ! Scalar
       !
       if (.not. params%valid_path('case.scalar')) then
         ! Set dummy values
         this%cp = 1.0_rp
         this%lambda = 1.0_rp
         call this%write_to_log(.false.)
         return
       end if

       ! Incorrect user input
       if (nondimensional .and. &
           (params%valid_path('case.scalar.lambda') .or. &
            params%valid_path('case.scalar.cp'))) then
           call neko_error("For non-dimensional setup set the Pe number for&
                         & the scalar")
       else if (.not. nondimensional .and. &
                params%valid_path('case.scalar.Pe')) then
           call neko_error("Dimensional material properties input detected,&
                         & because you set rho and mu for the fluid. &
                         & Please set cp and lambda for the scalar.")

       ! Non-dimensional case
       else if (nondimensional) then
          write(log_buf, '(A)') 'Non-dimensional scalar material properties &
                              & input.'
          call neko_log%message(log_buf, lvl=2)
          write(log_buf, '(A)') 'Specific heat capacity will be set to 1, &
                              & conductivity to 1/Pe.'
          call neko_log%message(log_buf, lvl=2)

          ! Read Pe into lambda for further manipulation.
          call json_get(params, 'case.scalar.Pe', this%lambda)
          write(log_buf, '(A,ES13.6)') 'Pe         :',  this%lambda
          call neko_log%message(log_buf)

          ! Set cp and rho to 1 since the setup is non-dimensional.
          this%cp = 1.0_rp
          this%rho = 1.0_rp
          ! Invert the Pe to get conductivity
          this%lambda = 1.0_rp/this%lambda
       ! Dimensional case
       else 
          call json_get(params, 'case.scalar.lambda', this%lambda)
          call json_get(params, 'case.scalar.cp', this%cp)
       end if
    end if

    call this%write_to_log(.true.)

  end subroutine material_properties_init

  !> Write final dimensional values to the log.
  subroutine write_to_log(this, scalar)
    class(material_properties_t), intent(inout) :: this
    logical, intent(in) :: scalar
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, '(A)') 'Set dimensional values:'
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'rho        :',  this%rho
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'mu         :',  this%mu
    call neko_log%message(log_buf)
    if (scalar) then
       write(log_buf, '(A,ES13.6)') 'cp         :',  this%cp
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES13.6)') 'lambda     :',  this%lambda
       call neko_log%message(log_buf)
    end if
  end subroutine write_to_log

end module material_properties
