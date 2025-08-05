! Copyright (c) 2025, The Neko Authors
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
!> Implements the `mass_flux_t` type.

module mass_flux
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t, simulation_component_allocate, &
      register_simulation_component
  use time_state, only : time_state_t
  use field_registry, only : neko_field_registry
  use case, only : case_t
  use json_utils, only : json_get, json_get_or_default
  use time_based_controller, only : time_based_controller_t
  use field, only : field_t
  use utils, only : neko_error
  use mpi_f08, only: MPI_SUM, MPI_IN_PLACE, MPI_Allreduce
  use comm, only: NEKO_COMM, MPI_REAL_PRECISION, MPI_EXTRA_PRECISION, pe_rank
  use math, only : glsum
  implicit none
  private

  !> A simulation component that computes the mass flux through a zone.
  type, public, extends(simulation_component_t) :: mass_flux_t
     !> X velocity component.
     type(field_t), pointer :: u => null()
     !> X velocity component.
     type(field_t), pointer :: v => null()
     !> X velocity component.
     type(field_t), pointer :: w => null()
     !> Zone indices over which the flux is computed
     integer, allocatable :: zone_indices(:)

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => mass_flux_init_from_json
     !> Common part of both constructors.
     procedure, private, pass(this) :: init_common => mass_flux_init_common
     !> Destructor.
     procedure, pass(this) :: free => mass_flux_free
     !> Here to compy with the interface, does nothing.
     procedure, pass(this) :: compute_ => mass_flux_compute
  end type mass_flux_t

  public :: mass_flux_register_types

contains

  !> Constructor from json.
  !> @param json JSON object with the parameters.
  !! @param case The case object.
  subroutine mass_flux_init_from_json(this, json, case)
    class(mass_flux_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    integer, allocatable :: zone_indices(:)
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: precision
    character(len=20), allocatable :: velocity(:)

    call this%init_base(json, case)

    if (json%valid_path("velocity_components")) then
       call json_get(json, "velocity_components", velocity)
    else
       velocity = ["u", "v", "w"]
    end if

    call json_get(json, "zone_indices", zone_indices)

    this%zone_indices = zone_indices

    if (size(velocity) .ne. 3) then
       call neko_error("The velocity_keyword should contain 3 items.")
    end if

    this%u => neko_field_registry%get_field_by_name(velocity(1))
    this%v => neko_field_registry%get_field_by_name(velocity(2))
    this%w => neko_field_registry%get_field_by_name(velocity(3))

    !call this%init_common(velocity)
  end subroutine mass_flux_init_from_json

  !> Common part of both constructors.
  !! @param fields Array of field names to be sampled.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! provided, fields are added to the main output file.
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine mass_flux_init_common(this, velocity, filename, precision)
    class(mass_flux_t), intent(inout) :: this
    character(len=20), intent(in) :: velocity(:)
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision


  end subroutine mass_flux_init_common

  !> Destructor.
  subroutine mass_flux_free(this)
    class(mass_flux_t), intent(inout) :: this
    call this%free_base()

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    deallocate(this%zone_indices)
  end subroutine mass_flux_free

  !> Here to comply with the interface, does nothing.
  subroutine mass_flux_compute(this, time)
    class(mass_flux_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp) :: flux_i, flux_sum
    real(kind=rp), allocatable :: flux(:)
    real(kind=rp) :: normal(3)
    integer :: zone_i, zi, i, j, k, l, face_i, elem_i, ierr
    integer :: lx, ly, lz


    allocate(flux(size(this%zone_indices)))
    flux = 0.0_rp

    lx = this%case%fluid%Xh%lx
    ly = this%case%fluid%Xh%ly
    lz = this%case%fluid%Xh%lz


    associate(coef => this%case%fluid%c_Xh)
      do zi = 1, size(this%zone_indices)
         zone_i = this%zone_indices(zi)
         do i = 1, this%case%msh%labeled_zones(zone_i)%size
            face_i = this%case%msh%labeled_zones(zone_i)%facet_el(i)%x(1)
            elem_i = this%case%msh%labeled_zones(zone_i)%facet_el(i)%x(2)

            select case (face_i)
            case (1)
               do l = 1, lz
                  do k = 1, ly
                     normal = -coef%get_normal(1, k, l, elem_i, face_i)
                     flux_i = &
                          this%u%x(1, k, l, elem_i) * normal(1) + &
                          this%v%x(1, k, l, elem_i) * normal(2) + &
                          this%w%x(1, k, l, elem_i) * normal(3)
                     flux(zi) = flux(zi) + &
                          flux_i * coef%area(k, l, face_i, elem_i)
                  end do
               end do
            case (2)
               do l = 1, lz
                  do k = 1, ly
                     normal = -coef%get_normal(lx, k, l, elem_i, face_i)
!                     write(*, *), "NORMAL 2", normal
                     flux_i = &
                          this%u%x(lx, k, l, elem_i) * normal(1) + &
                          this%v%x(lx, k, l, elem_i) * normal(2) + &
                          this%w%x(lx, k, l, elem_i) * normal(3)
                     flux(zi) = flux(zi) + &
                          flux_i * coef%area(k, l, face_i, elem_i)
                  end do
               end do
            case (3)
               do l = 1, lz
                  do j = 1, lx
                     normal = -coef%get_normal(j, 1, l, elem_i, face_i)
                     flux_i = &
                          this%u%x(j, 1, l, elem_i) * normal(1) + &
                          this%v%x(j, 1, l, elem_i) * normal(2) + &
                          this%w%x(j, 1, l, elem_i) * normal(3)
                     flux(zi) = flux(zi) + &
                          flux_i * coef%area(j, l, face_i, elem_i)
                  end do
               end do
            case (4)
               do l = 1, lz
                  do j = 1, lx
                     normal = -coef%get_normal(j, ly, l, elem_i, face_i)
                     flux_i = &
                          this%u%x(j, ly, l, elem_i) * normal(1) + &
                          this%v%x(j, ly, l, elem_i) * normal(2) + &
                          this%w%x(j, ly, l, elem_i) * normal(3)
                     flux(zi) = flux(zi) + &
                          flux_i * coef%area(j, l, face_i, elem_i)
                  end do
               end do
            case (5)
               do k = 1, ly
                  do j = 1, lx
                     normal = -coef%get_normal(j, k, 1, elem_i, face_i)
                     flux_i = &
                          this%u%x(j, k, 1, elem_i) * normal(1) + &
                          this%v%x(j, k, 1, elem_i) * normal(2) + &
                          this%w%x(j, k, 1, elem_i) * normal(3)
                     flux(zi) = flux(zi) + &
                          flux_i * coef%area(j, k, face_i, elem_i)
                  end do
               end do
            case (6)
               do k = 1, ly
                  do j = 1, lx
                     normal = -coef%get_normal(j, k, lz, elem_i, face_i)
                     flux_i = &
                          this%u%x(j, k, lz, elem_i) * normal(1) + &
                          this%v%x(j, k, lz, elem_i) * normal(2) + &
                          this%w%x(j, k, lz, elem_i) * normal(3)
                     flux(zi) = flux(zi) + &
                          flux_i * coef%area(j, k, face_i, elem_i)
                  end do
               end do
            end select
         end do
      end do
      call MPI_Allreduce(MPI_IN_PLACE, flux, size(flux), &
           MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)

      flux_sum = glsum(flux, size(flux))
      if (pe_rank .eq. 0) then
          write(*,*) "Flux across zone", this%zone_indices,  flux_sum
      end if
    end associate


  end subroutine mass_flux_compute

  subroutine my_alloctor(obj)
    class(simulation_component_t), allocatable, intent(inout) :: obj
    allocate(mass_flux_t::obj)
  end subroutine


  subroutine mass_flux_register_types()
      procedure(simulation_component_allocate), pointer :: allocator

      allocator => my_alloctor

      call register_simulation_component("check_flux", allocator)
  end subroutine

end module mass_flux
