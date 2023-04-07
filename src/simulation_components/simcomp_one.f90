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
module simcomp_one
  use num_types
  use math
  use simulation_component
  implicit none
  private
  
  type, public, extends(simulation_component_t) :: simcomp_one_t
     integer :: foo
     contains
      !> Constructor from json
      procedure, pass(this) :: init => simcomp_one_init_from_json
      !> Actual constructor
      procedure, pass(this) :: init_from_int => simcomp_one_init_from_int
      !> Destructor
      procedure, pass(this) :: free => simcomp_one_free
      !> The main function to be called at each timestep
     procedure, pass(this) :: compute => simcomp_one_compute
  end type simcomp_one_t
  
  contains
  
  subroutine simcomp_one_init_from_json(this, json_dict)
       class(simcomp_one_t), intent(inout) :: this
       !! stub for a json dict. The idea is that the construct will plock the
       !! necessary stuff for each concrete component from the dict, which is
       !! extracted from the case file.
       integer, intent(in) :: json_dict
       integer :: foo
       
       ! foo = json_dict.simulation_components.simcomp_one.foo
       foo = 2
       call simcomp_one_init_from_int(this, foo)
  end subroutine simcomp_one_init_from_json

  subroutine simcomp_one_init_from_int(this, foo)
       class(simcomp_one_t), intent(inout) :: this
       integer, intent(in) :: foo

       this%foo = foo
  end subroutine simcomp_one_init_from_int

  subroutine simcomp_one_free(this)
       class(simcomp_one_t), intent(inout) :: this
  end subroutine simcomp_one_free

  subroutine simcomp_one_compute(this)
       class(simcomp_one_t), intent(inout) :: this
    write(*,*) "Computing simcomp_one"
  end subroutine simcomp_one_compute
  

  
end module simcomp_one
