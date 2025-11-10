! Copyright (c) 2024-2025, The Neko Authors
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
!> Contains the `simcomp_executor_t` type.
module simcomp_executor
  use num_types, only : rp
  use simulation_component, only : simulation_component_t, &
       simulation_component_wrapper_t, simulation_component_factory
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use case, only : case_t
  use time_state, only : time_state_t
  use utils, only : neko_error
  use logger, only : neko_log
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
     integer, private :: n_simcomps
     !> The case
     type(case_t), pointer :: case
     !> Flag to indicate if the simcomp executor has been finalized.
     logical, private :: finalized = .false.
   contains
     !> Constructor.
     procedure, pass(this) :: init => simcomp_executor_init
     !> Destructor.
     procedure, pass(this) :: free => simcomp_executor_free
     !> Appending a new simcomp to the executor.
     procedure, pass(this) :: add => simcomp_executor_add
     !> Execute preprocess_ for all simcomps.
     procedure, pass(this) :: preprocess => simcomp_executor_preprocess
     !> Execute compute_ for all simcomps.
     procedure, pass(this) :: compute => simcomp_executor_compute
     !> Execute restart for all simcomps.
     procedure, pass(this) :: restart=> simcomp_executor_restart
     !> Finalize the initialization.
     procedure, private, pass(this) :: finalize => simcomp_executor_finalize
     !> Get the number of simcomps.
     procedure, pass(this) :: get_n => simcomp_executor_get_n
  end type simcomp_executor_t

  !> Global variable for the simulation component driver.
  type(simcomp_executor_t), target, public :: neko_simcomps

contains

  !> Constructor.
  !! @param case The case.
  !! @param simcomp_root The root name of the simulation components in the case.
  !! If not provided, the default is 'case.simulation_components'.
  subroutine simcomp_executor_init(this, case, simcomp_root)
    class(simcomp_executor_t), intent(inout) :: this
    type(case_t), target, intent(inout) :: case
    character(len=*), optional, intent(in) :: simcomp_root
    integer :: n_simcomps, i
    type(json_file) :: comp_subdict
    logical :: found
    ! Help array for finding minimal values
    logical, allocatable :: mask(:)
    ! The order value for each simcomp in order of appearance in the case file.
    integer, allocatable :: read_order(:), order(:)
    ! Location of the min value
    integer :: loc(1)
    integer :: max_order
    character(len=:), allocatable :: root_name, comp_type

    call this%free()
    this%case => case

    ! Get the root name of the simulation components if specified
    if (present(simcomp_root)) then
       root_name = simcomp_root
    else
       root_name = 'case.simulation_components'
    end if

    ! Get the core json object and the simulation components object
    if (.not. (root_name .in. case%params)) return
    call neko_log%section('Initialize simcomp')

    ! Set the number of simcomps and allocate the arrays
    call case%params%info(root_name, n_children = n_simcomps)
    this%n_simcomps = n_simcomps
    allocate(this%simcomps(n_simcomps))
    allocate(order(n_simcomps))
    allocate(read_order(n_simcomps))
    allocate(mask(n_simcomps), source = .true.)

    ! We need a separate loop to figure out the order, so that we can
    ! apply the order to the initialization as well.
    max_order = 0
    do i = 1, n_simcomps
       ! Create a new json containing just the subdict for this simcomp
       call json_extract_item(case%params, root_name, i, comp_subdict)

       call json_get_or_default(comp_subdict, "order", read_order(i), -1)
       if (read_order(i) .gt. max_order) then
          max_order = read_order(i)
       end if
    end do

    ! If the order was not specified, we use the order of appearance in the
    ! case file.
    do i = 1, n_simcomps
       if (read_order(i) == -1) then
          max_order = max_order + 1
          read_order(i) = max_order
       end if
    end do

    ! Figure out the execution order using a poor man's argsort.
    ! Searches for the location of the min value, each time masking out the
    ! found location prior to the next search.
    do i = 1, n_simcomps
       loc = minloc(read_order, mask = mask)
       order(i) = loc(1)
       mask(loc) = .false.
    end do

    ! Init in the determined order.
    do i = 1, n_simcomps
       call json_extract_item(case%params, root_name, order(i), comp_subdict)

       ! Log the component type
       call json_get(comp_subdict, "type", comp_type)
       call neko_log%message('- ' // trim(comp_type))

       call simulation_component_factory(this%simcomps(i)%simcomp, &
            comp_subdict, case)
    end do

    ! Cleanup
    deallocate(order)
    deallocate(read_order)
    deallocate(mask)

    call neko_log%end_section()
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

    nullify(this%case)
    this%finalized = .false.
  end subroutine simcomp_executor_free

  !> Appending a new simcomp to the executor.
  !! @param new_object The simcomp to append.
  !! @param settings The settings for the simcomp.
  subroutine simcomp_executor_add(this, object, settings)
    class(simcomp_executor_t), intent(inout) :: this
    class(simulation_component_t), intent(in) :: object
    type(json_file), intent(inout), optional :: settings

    class(simulation_component_wrapper_t), allocatable :: tmp_simcomps(:)
    integer :: i, position

    ! Find the first empty position
    position = 0
    do i = 1, this%n_simcomps
       if (.not. allocated(this%simcomps(i)%simcomp)) then
          position = i
          exit
       end if
    end do

    ! If no empty position was found, append to the end
    if (position .eq. 0) then
       if (this%n_simcomps .gt. 0) call move_alloc(this%simcomps, tmp_simcomps)
       allocate(this%simcomps(this%n_simcomps + 1))

       if (allocated(tmp_simcomps)) then
          do i = 1, this%n_simcomps
             call move_alloc(tmp_simcomps(i)%simcomp, this%simcomps(i)%simcomp)
          end do
       end if

       this%n_simcomps = this%n_simcomps + 1
       position = this%n_simcomps

       if (allocated(tmp_simcomps)) deallocate(tmp_simcomps)
    end if

    this%simcomps(position)%simcomp = object
    if (present(settings)) then
       call this%simcomps(position)%simcomp%init(settings, this%case)
    end if

    this%finalized = .false.

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

    ! Check that all components are initialized
    do i = 1, this%n_simcomps
       if (.not. allocated(this%simcomps(i)%simcomp)) then
          call neko_error("Simulation component not initialized.")
       end if
    end do

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
          call neko_error("Simulation component order must be contiguous " // &
               "starting at 1.")
       end if
       previous_found = order_found
    end do

    allocate(order_list(this%n_simcomps))
    order_list = 0
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
          deallocate(order_list)
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

    if (allocated(tmp_simcomps)) then
       deallocate(tmp_simcomps)
    end if
    if (allocated(order_list)) then
       deallocate(order_list)
    end if

    this%finalized = .true.
  end subroutine simcomp_executor_finalize

  !> Execute preprocess_ for all simcomps.
  !! @param time The current time
  subroutine simcomp_executor_preprocess(this, time)
    class(simcomp_executor_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i

    if (.not. this%finalized) call this%finalize()

    if (allocated(this%simcomps)) then
       do i = 1, size(this%simcomps)
          call this%simcomps(i)%simcomp%preprocess(time)
       end do
    end if

  end subroutine simcomp_executor_preprocess

  !> Execute compute_ for all simcomps.
  !! @param time The current time
  subroutine simcomp_executor_compute(this, time)
    class(simcomp_executor_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i

    if (.not. this%finalized) call this%finalize()

    if (allocated(this%simcomps)) then
       do i = 1, this%n_simcomps
          call this%simcomps(i)%simcomp%compute(time)
       end do
    end if

  end subroutine simcomp_executor_compute

  !> Execute restart for all simcomps.
  !! @param time The current time
  subroutine simcomp_executor_restart(this, time)
    class(simcomp_executor_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i

    if (allocated(this%simcomps)) then
       do i = 1, this%n_simcomps
          call this%simcomps(i)%simcomp%restart(time)
       end do
    end if

  end subroutine simcomp_executor_restart

  !> Get the number of simcomps.
  pure function simcomp_executor_get_n(this) result(n)
    class(simcomp_executor_t), intent(in) :: this
    integer :: n

    n = this%n_simcomps
  end function simcomp_executor_get_n

end module simcomp_executor
