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
!> Contains the `simcomp_executor_t` type.
module simcomp_executor
  use num_types, only : rp
  use simulation_component, only : simulation_component_wrapper_t
  use simulation_component_fctry, only : simulation_component_factory
  use json_module, only : json_file, json_core, json_value
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use case, only : case_t
  implicit none
  private

  !> Singleton type that serves as a driver for the simulation components.
  !! Stores all the components in the case and provides an interface matching
  !! that of a single simcomp, which executes the corresponding routines for
  !! each stored simcomp.
  !! The execution order is based on the order property of each simcomp.
  !! By default, the order is by the order of apparence in the case file.
  type, public :: simcomp_executor_t
     !> The simcomps.
     class(simulation_component_wrapper_t), allocatable :: simcomps(:)
     !> Index array defining the order of execution, i.e. simcomps(order(1)) is
     !! first to execute, and so on.
     integer, allocatable ::  order(:)
   contains
     !> Constructor.
     procedure, pass(this) :: init => simcomp_executor_init
     !> Destructor.
     procedure, pass(this) :: free => simcomp_executor_free
     !> Execute compute_ for all simcomps.
     procedure, pass(this) :: compute => simcomp_executor_compute
     !> Execute restart for all simcomps.
     procedure, pass(this) :: restart=> simcomp_executor_restart

  end type simcomp_executor_t

  !> Global variable for the simulation component driver.
  type(simcomp_executor_t), target, public :: neko_simcomps

contains

  !> Constructor.
  subroutine simcomp_executor_init(this, case)
    class(simcomp_executor_t), intent(inout) :: this
    type(case_t), target, intent(inout) :: case
    integer :: n_simcomps, i
    type(json_core) :: core
    type(json_value), pointer :: simcomp_object
    type(json_file) :: comp_subdict
    logical :: found
    ! Help array for finding minimal values
    logical, allocatable :: mask(:)
    ! The order value for each simcomp in order of appearance in the case file.
    integer, allocatable :: read_order(:)
    ! Location of the min value
    integer :: loc(1)

    call this%free()

    if (case%params%valid_path('case.simulation_components')) then

       call case%params%info('case.simulation_components', n_children=n_simcomps)
       allocate(this%simcomps(n_simcomps))
       allocate(this%order(n_simcomps))
       allocate(read_order(n_simcomps))
       allocate(mask(n_simcomps))
       mask = .true.

       call case%params%get_core(core)
       call case%params%get('case.simulation_components', simcomp_object, found)

       ! We need a separate loop to figure out the order, so that we can
       ! apply the order to the initialization as well.
       do i=1, n_simcomps
          ! Create a new json containing just the subdict for this simcomp
          call json_extract_item(core, simcomp_object, i, comp_subdict)
          call json_get_or_default(comp_subdict, "order", read_order(i), i)
       end do

       ! Figure out the execution order using a poor man's argsort.
       ! Searches for the location of the min value, each time masking out the
       ! found location prior to the next search.
       do i= 1, n_simcomps
         loc = minloc(read_order, mask=mask)
         this%order(i) = loc(1)
         mask(loc) = .false.
       end do

       ! Init in the determined order.
       do i=1, n_simcomps
          call json_extract_item(core, simcomp_object, this%order(i),&
                                    comp_subdict)
          ! Have to add, the simcomp constructor expects it.
          if (.not. comp_subdict%valid_path("order")) then
             call comp_subdict%add("order", this%order(i))
          end if
          call simulation_component_factory(this%simcomps(i)%simcomp, &
                                            comp_subdict, case)
       end do
    end if
  end subroutine simcomp_executor_init

  !> Destructor.
  subroutine simcomp_executor_free(this)
    class(simcomp_executor_t), intent(inout) :: this
    integer :: i

    if (allocated(this%order)) deallocate(this%order)

    if (allocated(this%simcomps)) then
       do i=1, size(this%simcomps)
          call this%simcomps(i)%simcomp%free
       end do
       deallocate(this%simcomps)
    end if
  end subroutine simcomp_executor_free

  !> Execute compute_ for all simcomps.
  !! @param t The time value.
  !! @param tstep The timestep number.
  subroutine simcomp_executor_compute(this, t, tstep)
    class(simcomp_executor_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: i

    if (allocated(this%simcomps)) then
       do i=1, size(this%simcomps)
          call this%simcomps(this%order(i))%simcomp%compute(t, tstep)
       end do
    end if

  end subroutine simcomp_executor_compute

  !> Execute restart for all simcomps.
  !! @param t The time value.
  subroutine simcomp_executor_restart(this, t)
    class(simcomp_executor_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i

    if (allocated(this%simcomps)) then
       do i=1, size(this%simcomps)
          call this%simcomps(this%order(i))%simcomp%restart(t)
       end do
    end if

  end subroutine simcomp_executor_restart

end module simcomp_executor
