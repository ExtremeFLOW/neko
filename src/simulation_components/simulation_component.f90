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
!
!> Simulation components are objects that encapsulate functionality that can be
!! fit to a particular execution pattern.
!! @note
!! The canonical way to abbreviate simulation_component is simcomp.
module simulation_component
  use num_types, only : rp
  use json_module, only : json_file
  use case, only : case_t
  implicit none
  private
  
  !> Base abstract class for simulation components.
  type, abstract, public :: simulation_component_t
     !> Pointer to the simulation case.
     class(case_t), pointer :: case
   contains
     !> Constructor for the simulation_component_t (base) class.
     procedure, pass(this) :: init_base => simulation_component_init_base
     !> The common constructor using a JSON dictionary.
     procedure(simulation_component_init), pass(this), deferred :: init
     !> Destructor.
     procedure(simulation_component_free), pass(this), deferred :: free
     !> The main function to be executed during the run.
     procedure(simulation_component_compute), pass(this), deferred :: compute
  end type simulation_component_t
  
  !> A helper type that is needed to have an array of polymorphic objects
  type, public :: simulation_component_wrapper_t
    class(simulation_component_t), allocatable :: simcomp
  end type simulation_component_wrapper_t

  
  abstract interface
     !> The common constructor using a JSON dictionary.
     !! @param json The JSON with properties.
     !! @param case The case_t object. 
     subroutine simulation_component_init(this, json, case)  
       import simulation_component_t, json_file, case_t
       class(simulation_component_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       class(case_t), intent(inout), target :: case
     end subroutine
  end interface

  abstract interface
     !> Destructor.
     subroutine simulation_component_free(this)  
       import simulation_component_t
       class(simulation_component_t), intent(inout) :: this
     end subroutine
  end interface

  abstract interface
     !> The main function to be executed during the run.
     !! @param t The time value.
     !! @param tstep The current time-step
     subroutine simulation_component_compute(this, t, tstep)  
       import simulation_component_t, rp
       class(simulation_component_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
     end subroutine
  end interface

contains
  !> Constructor for the `simulation_component_t` (base) class.
  subroutine simulation_component_init_base(this, json, case)  
    class(simulation_component_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case

    this%case => case

  end subroutine simulation_component_init_base

  
end module simulation_component
