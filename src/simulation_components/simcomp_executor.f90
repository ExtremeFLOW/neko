! Copyright (c) 2024, The Neko Authors
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
  use simulation_component, only : simulation_component_t, simulation_component_wrapper_t
  use simulation_component_fctry, only : simulation_component_factory
  use json_module, only : json_file, json_core, json_value
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use case, only : case_t
  use utils, only : neko_error
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
     !> Number of simcomps
     integer :: n_simcomps
   contains
     !> Constructor.
     procedure, pass(this) :: init => simcomp_executor_init
     !> Destructor.
     procedure, pass(this) :: free => simcomp_executor_free
     !> Appending a new simcomp to the executor.
     procedure, pass(this) :: add => simcomp_executor_add
     !> Execute compute_ for all simcomps.
     procedure, pass(this) :: compute => simcomp_executor_compute
     !> Execute restart for all simcomps.
     procedure, pass(this) :: restart=> simcomp_executor_restart
     !> Finalize the initialization.
     procedure, pass(this) :: finalize => simcomp_executor_finalize
  end type simcomp_executor_t

  !> Global variable for the simulation component driver.
  type(simcomp_executor_t), public :: neko_simcomps

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

    call this%free()

    if (case%params%valid_path('case.simulation_components')) then

       call case%params%info('case.simulation_components', n_children=n_simcomps)

       this%n_simcomps = n_simcomps
       allocate(this%simcomps(n_simcomps))

       call case%params%get_core(core)
       call case%params%get('case.simulation_components', simcomp_object, found)

       ! Init in the determined order.
       do i = 1, n_simcomps
          call json_extract_item(core, simcomp_object, i, comp_subdict)
          call simulation_component_factory(this%simcomps(i)%simcomp, &
                                            comp_subdict, case)
       end do
    end if
  end subroutine simcomp_executor_init

  !> Destructor.
  subroutine simcomp_executor_free(this)
    class(simcomp_executor_t), intent(inout) :: this
    integer :: i

    if (allocated(this%simcomps)) then
       do i = 1, this%n_simcomps
          call this%simcomps(i)%simcomp%free
       end do
       deallocate(this%simcomps)
    end if
  end subroutine simcomp_executor_free

  !> Appending a new simcomp to the executor.
  !! @param simcomp The simcomp to append.
  subroutine simcomp_executor_add(this, simcomp)
    class(simcomp_executor_t), intent(inout) :: this
    class(simulation_component_t), intent(in) :: simcomp

    class(simulation_component_wrapper_t), allocatable :: tmp_simcomps(:)
    integer :: i

    if (allocated(this%simcomps)) then
       call move_alloc(this%simcomps, tmp_simcomps)
    end if

    ! Insert the simulation component into the list
    allocate(this%simcomps(this%n_simcomps+1))
    do i = 1, this%n_simcomps
       call move_alloc(tmp_simcomps(i)%simcomp, this%simcomps(i)%simcomp)
    end do

    this%n_simcomps = this%n_simcomps + 1
    this%simcomps(this%n_simcomps)%simcomp = simcomp

    if (allocated(tmp_simcomps)) then
       deallocate(tmp_simcomps)
    end if

  end subroutine simcomp_executor_add

  !> Finalize the initialization.
  !! Sorts the simcomps based on the order property.
  !! Additionally we check that the order is unique, contiguous, starts at 1 and
  !! within bounds.
  subroutine simcomp_executor_finalize(this)
    class(simcomp_executor_t), intent(inout) :: this
    integer :: i, order, max_order
    logical :: order_found, previous_found

    class(simulation_component_wrapper_t), allocatable :: tmp_simcomps(:)
    integer, allocatable :: order_list(:)

    ! Check that the order is unique and contiguous
    previous_found = .true.
    do order = 1, this%n_simcomps
       order_found = .false.
       do i = 1, this%n_simcomps
          if (this%simcomps(i)%simcomp%order == order .and. order_found) then
             call neko_error("Simulation component order must be unique.")
          else if (this%simcomps(i)%simcomp%order == order) then
             order_found = .true.
          end if
       end do
       if (order_found .and. .not. previous_found) then
          call neko_error("Simulation component order must be contiguous &
            &starting at 1.")
       end if
       previous_found = order_found
    end do

    allocate(order_list(this%n_simcomps))
    max_order = 0
    do i = 1, this%n_simcomps
       order_list(i) = this%simcomps(i)%simcomp%order
       if (order_list(i) .gt. max_order) then
          max_order = order_list(i)
       end if
    end do

    do i = 1, this%n_simcomps
       if (order_list(i) .eq. -1) then
          order_list(i) = max_order + 1
          max_order = max_order + 1
       end if
    end do

    ! Check that the order is within bounds
    do i = 1, this%n_simcomps
       if (order_list(i) .gt. this%n_simcomps) then
          call neko_error("Simulation component order is out of bounds.")
       end if
    end do

    ! Reorder the simcomps based on the order specified
    call move_alloc(this%simcomps, tmp_simcomps)
    allocate(this%simcomps(this%n_simcomps))
    do i = 1, this%n_simcomps
       order = order_list(i)
       call move_alloc(tmp_simcomps(i)%simcomp, this%simcomps(order)%simcomp)
    end do

    deallocate(tmp_simcomps)
    deallocate(order_list)

  end subroutine simcomp_executor_finalize

  !> Execute compute_ for all simcomps.
  !! @param t The time value.
  !! @param tstep The timestep number.
  subroutine simcomp_executor_compute(this, t, tstep)
    class(simcomp_executor_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: i

    if (allocated(this%simcomps)) then
       do i = 1, this%n_simcomps
          call this%simcomps(i)%simcomp%compute(t, tstep)
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
       do i = 1, this%n_simcomps
          call this%simcomps(i)%simcomp%restart(t)
       end do
    end if

  end subroutine simcomp_executor_restart

end module simcomp_executor
