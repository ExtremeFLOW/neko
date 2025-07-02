! Copyright (c) 2020-2023, The Neko Authors
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
!> Defines an output for a fluid
module fluid_output
  use num_types, only : rp
  use fluid_scheme_incompressible, only : fluid_scheme_incompressible_t
  use fluid_scheme_base, only : fluid_scheme_base_t
  use scalar_scheme, only : scalar_scheme_t
  use field_list, only : field_list_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use device
  use output, only : output_t
  use scalars, only : scalars_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  implicit none
  private

  !> Fluid output
  type, public, extends(output_t) :: fluid_output_t
     type(field_list_t) :: fluid
   contains
     procedure, pass(this) :: init => fluid_output_init
     procedure, pass(this) :: sample => fluid_output_sample
     procedure, pass(this) :: free => fluid_output_free
  end type fluid_output_t

contains

  subroutine fluid_output_init(this, precision, fluid, scalar_fields, name, path, &
       fmt, layout)
    class(fluid_output_t), intent(inout) :: this
    integer, intent(inout) :: precision
    class(fluid_scheme_base_t), intent(in), target :: fluid
    class(scalars_t), intent(in), optional, target :: scalar_fields
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    character(len=*), intent(in), optional :: fmt
    integer, intent(in), optional :: layout
    character(len=1024) :: fname
    integer :: i, j, n_scalars
    character(len=10) :: suffix
    logical :: has_max_wave_speed, has_density
    type(field_t), pointer :: max_wave_speed_field

    suffix = '.fld'
    if (present(fmt)) then
       if (fmt .eq. 'adios2') then
          suffix = '.bp'
       end if
    end if

    if (present(name) .and. present(path)) then
       fname = trim(path) // trim(name) // trim(suffix)
    else if (present(name)) then
       fname = trim(name) // trim(suffix)
    else if (present(path)) then
       fname = trim(path) // 'field' // trim(suffix)
    else
       fname = 'field' // trim(suffix)
    end if

    if (present(layout)) then
       call this%init_base(fname, precision, layout)
    else
       call this%init_base(fname, precision)
    end if

    ! Calculate total number of fields
    n_scalars = 0
    if (present(scalar_fields)) then
       n_scalars = size(scalar_fields%scalar_fields)
    end if

    ! Check if max_wave_speed field exists (for compressible flows)
    has_max_wave_speed = neko_field_registry%field_exists("max_wave_speed")

    ! Check if density field exists (for compressible flows)
    has_density = associated(fluid%rho)

    ! Initialize field list with appropriate size
    ! Standard fields: p, u, v, w (4)
    ! Scalar fields: n_scalars
    ! Compressible fields: density + max_wave_speed (2 additional)
    i = 4

    if (has_density) then
       i = i + 1
    end if

    if (has_max_wave_speed) then
       i = i + 1
    end if

    call this%fluid%init(i + n_scalars)

    call this%fluid%assign(1, fluid%p)
    call this%fluid%assign(2, fluid%u)
    call this%fluid%assign(3, fluid%v)
    call this%fluid%assign(4, fluid%w)

    ! Assign all scalar fields first
    i = 4
    if (present(scalar_fields)) then
       do j = 1, n_scalars
          i = i + 1
          call this%fluid%assign(i, scalar_fields%scalar_fields(j)%s)
       end do
    end if

    ! Add density field if it exists (for compressible flows)
    if (has_density) then
       i = i + 1
       call this%fluid%assign(i, fluid%rho)
    end if

    ! Add max_wave_speed field if it exists (for compressible flows)
    if (has_max_wave_speed) then
       i = i + 1
       max_wave_speed_field => neko_field_registry%get_field("max_wave_speed")
       call this%fluid%assign(i, max_wave_speed_field)
    end if

  end subroutine fluid_output_init

  !> Destroy a fluid output list
  subroutine fluid_output_free(this)
    class(fluid_output_t), intent(inout) :: this

    call this%fluid%free()

  end subroutine fluid_output_free

  !> Sample a fluid solution at time @a t
  subroutine fluid_output_sample(this, t)
    class(fluid_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i

    if (NEKO_BCKND_DEVICE .eq. 1) then

       associate(fields => this%fluid%items)
         do i = 1, size(fields)
            call device_memcpy(fields(i)%ptr%x, fields(i)%ptr%x_d, &
                 fields(i)%ptr%dof%size(), DEVICE_TO_HOST, &
                 sync = (i .eq. size(fields))) ! Sync on the last field
         end do
       end associate

    end if

    call this%file_%write(this%fluid, t)

  end subroutine fluid_output_sample

end module fluid_output
