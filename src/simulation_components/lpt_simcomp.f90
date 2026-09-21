! Copyright (c) 2026, The Neko Authors
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
!> Simulation-component handler for Lagrangian particle tracking.
module lpt_simcomp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use case, only : case_t
  use time_state, only : time_state_t
  use num_types, only : dp
  use lpt, only : lpt_t
  implicit none
  private

  !> A simulation component that drives passive Lagrangian particle tracking.
  type, public, extends(simulation_component_t) :: lpt_simcomp_t
     !> Lagrangian particle tracking solver.
     type(lpt_t) :: lpt
   contains
     !> Constructor from json.
     procedure, pass(this) :: init => lpt_simcomp_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => lpt_simcomp_free
     !> Preprocess hook.
     procedure, pass(this) :: preprocess_ => lpt_simcomp_preprocess
     !> Main compute hook.
     procedure, pass(this) :: compute_ => lpt_simcomp_compute
     !> The time until the next scheduled execution, output included.
     procedure, pass(this) :: time_to_next => lpt_simcomp_time_to_next
  end type lpt_simcomp_t

contains

  !> Construct from JSON.
  subroutine lpt_simcomp_init_from_json(this, json, case)
    class(lpt_simcomp_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case

    call this%free()
    call this%init_base(json, case)
    this%preprocess_controller = this%compute_controller
    call this%lpt%init(json, case)
    this%lpt%output_controller = this%output_controller
    this%name = this%lpt%name
  end subroutine lpt_simcomp_init_from_json

  !> The time until the next scheduled execution, including the output of
  !! the particle tracker, which is driven by a controller of its own.
  !! @param time The current time.
  !! @param dt The time step about to be taken, before it is shortened to
  !! land on a scheduled time.
  function lpt_simcomp_time_to_next(this, time, dt) result(t)
    class(lpt_simcomp_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    real(kind=dp), intent(in) :: dt
    real(kind=dp) :: t

    t = min(this%preprocess_controller%time_to_next(time, dt), &
         this%compute_controller%time_to_next(time, dt), &
         this%output_controller%time_to_next(time, dt), &
         this%lpt%output_controller%time_to_next(time, dt))

  end function lpt_simcomp_time_to_next

  !> Free the component.
  subroutine lpt_simcomp_free(this)
    class(lpt_simcomp_t), intent(inout) :: this

    call this%free_base()
    call this%lpt%free()
  end subroutine lpt_simcomp_free

  !> Update LPT before the fluid solve.
  subroutine lpt_simcomp_preprocess(this, time)
    class(lpt_simcomp_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call this%lpt%preprocess(time)
  end subroutine lpt_simcomp_preprocess

  !> Update LPT after the fluid solve.
  subroutine lpt_simcomp_compute(this, time)
    class(lpt_simcomp_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call this%lpt%compute(time)
  end subroutine lpt_simcomp_compute

end module lpt_simcomp
