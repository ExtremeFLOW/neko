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
!     copyright notice, this list of conditions and the following disclaimer
!     in the documentation and/or other materials provided with the
!     distribution.
!
!   * Neither the name of the authors nor the names of its contributors may
!     be used to endorse or promote products derived from this software
!     without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Data and filter construction for spectral vanishing viscosity.
module spectral_vanishing_viscosity
  use num_types, only : rp
  use elementwise_filter, only : elementwise_filter_t
  use registry, only : neko_registry
  use field, only : field_t
  use utils, only : neko_error, neko_type_error
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default, json_get_or_lookup
  use coefs, only : coef_t
  use math, only : cfill, copy, rzero, col2
  use device_math, only : device_rzero, device_cfill, device_copy, device_col2
  use device, only : device_map, device_unmap, device_memcpy, HOST_TO_DEVICE
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  use neko_config, only : NEKO_BCKND_DEVICE
  implicit none
  private

  character(len=3), parameter :: KNOWN_DIRECTIONS(7) = [character(len=3) :: &
       "rst", "rs", "rt", "st", "r", "s", "t"]
  character(len=13), parameter :: KNOWN_FORMULATIONS(1) = [character(len=13) :: &
       "Kirby-Sherwin"]
  character(len=5), parameter :: KNOWN_NU_TYPES(2) = [character(len=5) :: &
       "value", "field"]
  character(len=5), parameter :: KNOWN_KERNEL_TYPES(1) = [character(len=5) :: &
       "power"]


  !> Spectral vanishing viscosity configuration and coefficients.
  type, public :: svv_t
     !> Element-local modal low-pass filter.
     type(elementwise_filter_t) :: filter
     !> Reference directions in which the filter is applied.
     character(len=:), allocatable :: direction
     !> Kernel type defining the modal transfer function.
     character(len=:), allocatable :: kernel_type
     !> SEM coefficients associated with this SVV object.
     type(coef_t), pointer :: coef => null()
     !> Density-weighted SVV viscosity.
     real(kind=rp), allocatable :: h1(:,:,:,:)
     type(c_ptr) :: h1_d = C_NULL_PTR
     !> Identity matrix used to disable filtering in selected directions.
     real(kind=rp), allocatable :: ident(:,:)
     type(c_ptr) :: ident_d = C_NULL_PTR
     !> Name and pointer for a field-valued SVV viscosity.
     character(len=:), allocatable :: nue_field_name
     type(field_t), pointer :: nue => null()
     !> Whether the viscosity field is updated every time step.
     logical :: tvar_h1 = .false.
   contains
     procedure, pass(this) :: init => svv_init_from_json
     procedure, pass(this) :: free => svv_free
     procedure, pass(this) :: update => svv_update_h1
  end type svv_t

contains

  !> Construct an SVV object from case parameters.
  !! @param this SVV object.
  !! @param json JSON object containing an `svv` dictionary.
  !! @param coef SEM coefficients.
  !! @param rho Density field.
  subroutine svv_init_from_json(this, json, coef, rho)
    class(svv_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in), target :: coef
    type(field_t), intent(in) :: rho
    real(kind=rp), allocatable :: transfer(:)
    real(kind=rp) :: nu_val, power_coef
    character(len=:), allocatable :: nu_type, direction, formulation
    integer :: i, lx


    call this%free()
    this%coef => coef
    lx = coef%Xh%lx

    call json_get_or_default(json, "svv.formulation", formulation, &
         "Kirby-Sherwin")
    if (trim(formulation) .ne. "Kirby-Sherwin") then
       call neko_error("This SVV operator only supports the Kirby-Sherwin " // &
            "formulation")
    end if

    call json_get_or_default(json, "svv.direction", direction, "rst")
    this%direction = trim(direction)
    if (all(trim(this%direction) .ne. KNOWN_DIRECTIONS)) then
       call neko_type_error("The SVV direction", this%direction, &
            KNOWN_DIRECTIONS)
    end if

    allocate(transfer(lx))

    call json_get(json, "svv.kernel_type", this%kernel_type)
    select case (trim(this%kernel_type))
    case ("power")
       call json_get_or_lookup(json, "svv.power_coefficient", power_coef)
       if (power_coef .eq. 0.0_rp) then
          transfer = 0.0_rp
       else
          do i = 1, lx
             transfer(i) = 1.0_rp - ((i - 1.0_rp) / (lx - 1.0_rp)) ** &
                  ((lx - 1.0_rp) * power_coef)
          end do
       end if
    case default
       call neko_type_error("The SVV kernel type `", &
            trim(this%kernel_type), KNOWN_KERNEL_TYPES)
    end select

    call this%filter%init_from_components(coef, "nonBoyd", transfer)
    deallocate(transfer)

    allocate(this%ident(lx, lx))
    this%ident = 0.0_rp
    do i = 1, lx
       this%ident(i, i) = 1.0_rp
    end do

    allocate(this%h1(lx, lx, lx, coef%msh%nelv))
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%ident, this%ident_d, lx * lx)
       call device_memcpy(this%ident, this%ident_d, lx * lx, &
            HOST_TO_DEVICE, sync = .false.)
       call device_map(this%h1, this%h1_d, coef%dof%size())
       call device_rzero(this%h1_d, coef%dof%size())
    else
       call rzero(this%h1, coef%dof%size())
    end if

    call json_get(json, "svv.nu.type", nu_type)
    select case (trim(nu_type))
    case ("value")
       call json_get_or_lookup(json, "svv.nu.value", nu_val)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_cfill(this%h1_d, nu_val, coef%dof%size())
          call device_col2(this%h1_d, rho%x_d, coef%dof%size())
       else
          call cfill(this%h1, nu_val, coef%dof%size())
          call col2(this%h1, rho%x, coef%dof%size())
       end if
    case ("field")
       call json_get_or_default(json, "svv.nu.time_variable", &
            this%tvar_h1, .true.)
       call json_get(json, "svv.nu.field_name", this%nue_field_name)
       if (neko_registry%field_exists(this%nue_field_name)) then
          this%nue => neko_registry%get_field(this%nue_field_name)
          call svv_copy_viscosity(this, rho)
       else if (.not. this%tvar_h1) then
          call neko_error("The static SVV viscosity field `" // &
               this%nue_field_name // "` does not exist")
       end if
    case default
       call neko_error("Invalid svv.nu.type: " // trim(nu_type))
    end select

    if (allocated(nu_type)) deallocate(nu_type)
    if (allocated(direction)) deallocate(direction)
    if (allocated(formulation)) deallocate(formulation)
  end subroutine svv_init_from_json

  !> Update a time-varying, field-valued SVV viscosity.
  !! @param this SVV object.
  !! @param rho Density field.
  !! @param tstep Current time-step index.
  subroutine svv_update_h1(this, rho, tstep)
    class(svv_t), intent(inout) :: this
    type(field_t), intent(in) :: rho
    integer, intent(in) :: tstep

    if (.not. this%tvar_h1) return
    if (.not. associated(this%nue)) then
       if (.not. neko_registry%field_exists(this%nue_field_name)) then
          call neko_error("The SVV viscosity field `" // &
               this%nue_field_name // "` does not exist at time step " // &
               trim(adjustl(to_string(tstep))))
       end if
       this%nue => neko_registry%get_field(this%nue_field_name)
    end if
    call svv_copy_viscosity(this, rho)
  end subroutine svv_update_h1

  !> Copy the configured viscosity field and multiply it by density.
  !! @param this SVV object.
  !! @param rho Density field.
  subroutine svv_copy_viscosity(this, rho)
    class(svv_t), intent(inout) :: this
    type(field_t), intent(in) :: rho

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%h1_d, this%nue%x_d, this%coef%dof%size())
       call device_col2(this%h1_d, rho%x_d, this%coef%dof%size())
    else
       call copy(this%h1, this%nue%x, this%coef%dof%size())
       call col2(this%h1, rho%x, this%coef%dof%size())
    end if
  end subroutine svv_copy_viscosity

  !> Release all resources owned by an SVV object.
  !! @param this SVV object.
  subroutine svv_free(this)
    class(svv_t), intent(inout) :: this

    call this%filter%free()
    if (allocated(this%h1)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_unmap(this%h1, this%h1_d)
       end if
       deallocate(this%h1)
    end if
    if (allocated(this%ident)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_unmap(this%ident, this%ident_d)
       end if
       deallocate(this%ident)
    end if
    if (allocated(this%direction)) deallocate(this%direction)
    if (allocated(this%kernel_type)) deallocate(this%kernel_type)
    if (allocated(this%nue_field_name)) deallocate(this%nue_field_name)
    nullify(this%coef)
    nullify(this%nue)
    this%h1_d = C_NULL_PTR
    this%ident_d = C_NULL_PTR
    this%tvar_h1 = .false.
  end subroutine svv_free

  !> Convert an integer to a compact string.
  function to_string(value) result(string)
    integer, intent(in) :: value
    character(len=32) :: string
    write(string, '(I0)') value
  end function to_string

end module spectral_vanishing_viscosity
