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
!> Defines a factory subroutine for simulation components.
module simulation_component_fctry
  use simulation_component, only : simulation_component_t
  use vorticity, only : vorticity_t
  use lambda2, only : lambda2_t
  use probes, only : probes_t
  use json_module, only : json_file
  use case, only : case_t
  use json_utils, only : json_get
  implicit none
  private
  
  public :: simulation_component_factory
  
  contains

  !> Simulation component factory. Both constructs and initializes the object.
  !! @param json JSON object initializing the simulation component.
  subroutine simulation_component_factory(simcomp, json, case)
       class(simulation_component_t), allocatable, intent(inout) :: simcomp
       type(json_file), intent(inout) :: json
       class(case_t), intent(inout), target :: case
       character(len=:), allocatable :: simcomp_type
              
       call json_get(json, "type", simcomp_type)

       if (trim(simcomp_type) .eq. "vorticity") then 
          allocate(vorticity_t::simcomp)
       end if
       if (trim(simcomp_type) .eq. "lambda2") then 
          allocate(lambda2_t::simcomp)
       end if
       if (trim(simcomp_type) .eq. "probes") then 
          allocate(probes_t::simcomp)
       end if
       
       ! Initialize
       call simcomp%init(json, case)

  end subroutine simulation_component_factory

end module simulation_component_fctry
