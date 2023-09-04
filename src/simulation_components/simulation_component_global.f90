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
!> Contains the global simulation component array and its init for the given 
!! case.
module simulation_component_global
  use simulation_component, only : simulation_component_wrapper_t
  use simulation_component_fctry, only : simulation_component_factory
  use json_module, only : json_file, json_core, json_value
  use case, only : case_t
  implicit none
  private
  
  class(simulation_component_wrapper_t), public, allocatable :: neko_simcomps(:)
  public :: simcomps_global_init

contains

  !> Populates the global `neko_simcomps` array with objects based on the json
  !! configuration and the simulation case.
  !! @param case Simulation case.
  subroutine simcomps_global_init(case)
    type(case_t), intent(inout) :: case
    integer :: n_simcomps, i
    type(json_core) :: core
    type(json_value), pointer :: simcomp_object, comp_pointer 
    type(json_file) :: comp_subdict
    character(len=:), allocatable :: buffer
    logical :: found

    if (case%params%valid_path('case.simulation_components')) then

      call case%params%info('case.simulation_components', n_children=n_simcomps)
      allocate(neko_simcomps(n_simcomps))

      call case%params%get_core(core)
      call case%params%get('case.simulation_components', simcomp_object, found)
      print *, n_simcomps
      do i=1, n_simcomps
         ! Create a new json containing just the subdict for this simcomp
         call core%get_child(simcomp_object, i, comp_pointer, found)
         call core%print_to_string(comp_pointer, buffer)
         call comp_subdict%load_from_string(buffer)
         call simulation_component_factory(neko_simcomps(i)%simcomp, &
                                           comp_subdict, case) 
      end do
    end if
  end subroutine

end module simulation_component_global
