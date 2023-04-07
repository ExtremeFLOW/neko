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
!! fit to a particular software pattern.
!! The canonical way to abbreviate simulation_component is simcomp
module simulation_component
  use num_types
  use math
  implicit none
  private
  
  !> Base abstract class for simulation components
  type, abstract, public :: simulation_component_t
   contains
     !> Constructor
     procedure(simulation_component_init), pass(this), deferred :: init
     !> Destructor
     procedure(simulation_component_free), pass(this), deferred :: free
     !> The main function to be called at each timestep
     procedure(simulation_component_compute), pass(this), deferred :: compute
  end type simulation_component_t
  
  !> A helper type that is needed to have an array of polymorphic objects
  type simulation_component_wrapper_t
    class(simulation_component_t), allocatable :: simcomp
  end type simulation_component_wrapper_t

  class(simulation_component_wrapper_t), public, allocatable :: neko_simcomps(:)
  
  abstract interface
     subroutine simulation_component_init(this, json_dict)  
       import simulation_component_t
       class(simulation_component_t), intent(inout) :: this
       !! stub for a json dict. The idea is that the construct will plock the
       !! necessary stuff for each concrete component from the dict, which is
       !! extracted from the case file.
       integer, intent(in) :: json_dict
     end subroutine
  end interface

  abstract interface
     subroutine simulation_component_free(this)  
       import simulation_component_t
       class(simulation_component_t), intent(inout) :: this
     end subroutine
  end interface

  abstract interface
     subroutine simulation_component_compute(this)  
       import simulation_component_t
       class(simulation_component_t), intent(inout) :: this
     end subroutine
  end interface
  
end module simulation_component
