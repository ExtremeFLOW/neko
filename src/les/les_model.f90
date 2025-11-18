! Copyright (c) 2023-2024, The Neko Authors
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
!> Implements `les_model_t`.
module les_model
  use num_types, only : rp
  use fluid_scheme_base, only : fluid_scheme_base_t
  use fluid_pnpn, only : fluid_pnpn_t
  use time_scheme_controller, only : time_scheme_controller_t
  use rhs_maker, only : rhs_maker_sumab_t, rhs_maker_sumab_fctry
  use field, only : field_t
  use field_series, only : field_series_t
  use json_module, only : json_file
  use field_registry, only : neko_field_registry
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use gs_ops, only : GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_memcpy, HOST_TO_DEVICE
  use math, only : col2
  use device_math, only : device_col2
  use utils, only : neko_type_error, neko_error, neko_warning
  use comm, only : pe_rank
  implicit none
  private

  ! List of all possible types created by the factory routine
  character(len=20) :: DELTA_KNOWN_TYPES(3) = [character(len=20) :: &
       "pointwise", &
       "elementwise_average", &
       "elementwise_max"]

  !> Base abstract type for LES models based on the Boussinesq approximation.
  type, abstract, public :: les_model_t
     !> Pointer to the extrapolation scheme.
     type(time_scheme_controller_t), pointer :: ext_bdf => null()
     !> Pointer to the lag list of the velocities
     type(field_series_t), pointer :: ulag => null()
     type(field_series_t), pointer :: vlag => null()
     type(field_series_t), pointer :: wlag => null()
     !> Summation of AB/BDF contributions to extrapolate the field
     class(rhs_maker_sumab_t), allocatable :: sumab
     !> Logical variable for extrapolation
     logical :: if_ext = .false.
     !> Subgrid kinematic viscosity.
     type(field_t), pointer :: nut => null()
     !> LES lengthscale type
     character(len=:), allocatable :: delta_type
     !> LES lengthscale.
     type(field_t), pointer :: delta => null()
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
   contains
     !> Constructor for the les_model_t (base) class.
     procedure, pass(this) :: init_base => les_model_init_base
     !> Destructor for the les_model_t (base) class.
     procedure, pass(this) :: free_base => les_model_free_base
     !> Compute the LES length-scale
     procedure, pass(this) :: compute_delta => les_model_compute_delta
     !> The common constructor.
     procedure(les_model_init), pass(this), deferred :: init
     !> Destructor.
     procedure(les_model_free), pass(this), deferred :: free
     !> Compute eddy viscosity.
     procedure(les_model_compute), pass(this), deferred :: compute
  end type les_model_t

  abstract interface
     !> Compute eddy viscosity.
     !! @param t The time value.
     !! @param tstep The current time-step.
     subroutine les_model_compute(this, t, tstep)
       import les_model_t, rp
       class(les_model_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
     end subroutine les_model_compute
  end interface

  abstract interface
     !> Common constructor.
     !! @param fluid The fluid_scheme_t object.
     !! @param json A dictionary with parameters.
     subroutine les_model_init(this, fluid, json)
       import les_model_t, json_file, fluid_scheme_base_t
       class(les_model_t), intent(inout) :: this
       class(fluid_scheme_base_t), intent(inout), target :: fluid
       type(json_file), intent(inout) :: json
     end subroutine les_model_init
  end interface

  abstract interface
     !> Destructor.
     subroutine les_model_free(this)
       import les_model_t
       class(les_model_t), intent(inout) :: this
     end subroutine les_model_free
  end interface

  interface
     !> LES model allocator.
     !! @param object The object to be allocated.
     !! @param type_name The name of the LES model.
     module subroutine les_model_allocator(object, type_name)
       class(les_model_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: type_name
     end subroutine les_model_allocator
  end interface

  interface
     !> LES model factory. Both allocates and initializes the object.
     !! @param object The object to be allocated.
     !! @param type_name The name of the LES model.
     !! @param fluid The fluid scheme base type pointer.
     !! @param dofmap SEM map of degrees of freedom.
     !! @param coef SEM coefficients.
     !! @param json A dictionary with parameters.
     module subroutine les_model_factory(object, type_name, fluid, json)
       class(les_model_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: type_name
       class(fluid_scheme_base_t), intent(inout) :: fluid
       type(json_file), intent(inout) :: json
     end subroutine les_model_factory
  end interface

  !
  ! Machinery for injecting user-defined types
  !

  !> Interface for an object allocator.
  !! Implemented in the user modules, should allocate the `obj` to the custom
  !! user type.
  abstract interface
     subroutine les_model_allocate(obj)
       import les_model_t
       class(les_model_t), allocatable, intent(inout) :: obj
     end subroutine les_model_allocate
  end interface

  interface
     !> Called in user modules to add an allocator for custom types.
     module subroutine register_les_model(type_name, allocator)
       character(len=*), intent(in) :: type_name
       procedure(les_model_allocate), pointer, intent(in) :: allocator
     end subroutine register_les_model
  end interface

  ! A name-allocator pair for user-defined types. A helper type to define a
  ! registry of custom allocators.
  type allocator_entry
     character(len=20) :: type_name
     procedure(les_model_allocate), pointer, nopass :: allocator
  end type allocator_entry

  !> Registry of LES model allocators for user-defined types
  type(allocator_entry), allocatable :: les_model_registry(:)

  !> The size of the `les_model_registry`
  integer :: les_model_registry_size = 0

  public :: les_model_factory, les_model_allocator, register_les_model, &
       les_model_allocate


contains
  !> Constructor for the les_model_t (base) class.
  !! @param fluid The fluid_scheme_t object.
  !! @param nu_name The name of the turbulent viscosity field.
  !! @param delta_type The type of filter size.
  !! @param if_ext Whether trapolate the velocity.
  subroutine les_model_init_base(this, fluid, nut_name, delta_type, if_ext)
    class(les_model_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    character(len=*), intent(in) :: nut_name
    character(len=*), intent(in) :: delta_type
    logical, intent(in) :: if_ext

    associate(dofmap => fluid%dm_Xh, &
         coef => fluid%c_Xh)

      call neko_field_registry%add_field(dofmap, trim(nut_name), .true.)
      call neko_field_registry%add_field(dofmap, "les_delta", .true.)
      this%nut => neko_field_registry%get_field(trim(nut_name))
      this%delta => neko_field_registry%get_field("les_delta")
      this%coef => fluid%c_Xh
      this%delta_type = delta_type
      this%if_ext = if_ext

      if (pe_rank .eq. 0) then
         if (if_ext .eqv. .true.) then
            call neko_warning("Extrapolation of the velocity in eddy &
            &viscosity estimation might be unstable.")
         else
            call neko_warning("The time integration for eddy viscosity &
            &estimation is only first-order accurate")
         end if
      end if

      call this%compute_delta()

      select type (fluid)
      type is (fluid_pnpn_t)
         this%ulag => fluid%ulag
         this%vlag => fluid%vlag
         this%wlag => fluid%wlag
         ! Setup backend dependent summation of AB/BDF
         this%ext_bdf => fluid%ext_bdf
         call rhs_maker_sumab_fctry(this%sumab)
      class default
         if (this%if_ext .eqv. .true.) then
            call neko_error("Fluid scheme does not support &
            &velocity extrapolation")
         end if
      end select

    end associate
  end subroutine les_model_init_base

  !> Destructor for the les_model_t (base) class.
  subroutine les_model_free_base(this)
    class(les_model_t), intent(inout) :: this

    nullify(this%nut)
    nullify(this%delta)
    nullify(this%coef)
    nullify(this%ulag)
    nullify(this%vlag)
    nullify(this%wlag)
    nullify(this%ext_bdf)

    if (allocated(this%delta_type)) then
       deallocate(this%delta_type)
    end if

    if (allocated(this%sumab)) then
       deallocate(this%sumab)
    end if
  end subroutine les_model_free_base

  !> Compute the LES lengthscale.
  !! For each GLL point, we take the distance between its neighbours in all 3
  !! directions divided by 2 with the exception of face nodes, where only one
  !! neighbour exists. To form the lengthscale the distances along 3 directions
  !! are multiplied, and a cubic root is extracted from the result. This
  !! roughly corresponds to a cubic root of the cell volume in FVM computations.
  subroutine les_model_compute_delta(this)
    class(les_model_t), intent(inout) :: this
    integer :: e, i, j, k
    integer :: im, ip, jm, jp, km, kp
    real(kind=rp) :: di, dj, dk, ndim_inv, volume_element
    integer :: lx_half, ly_half, lz_half
    character(len=:), allocatable :: type_string

    lx_half = this%coef%Xh%lx / 2
    ly_half = this%coef%Xh%ly / 2
    lz_half = this%coef%Xh%lz / 2

    if (this%delta_type .eq. "elementwise_max") then
       ! use a same length scale throughout an entire element
       ! the length scale is based on maximum GLL spacing
       do e = 1, this%coef%msh%nelv
          di = (this%coef%dof%x(lx_half, 1, 1, e) &
               - this%coef%dof%x(lx_half + 1, 1, 1, e))**2 &
               + (this%coef%dof%y(lx_half, 1, 1, e) &
               - this%coef%dof%y(lx_half + 1, 1, 1, e))**2 &
               + (this%coef%dof%z(lx_half, 1, 1, e) &
               - this%coef%dof%z(lx_half + 1, 1, 1, e))**2

          dj = (this%coef%dof%x(1, ly_half, 1, e) &
               - this%coef%dof%x(1, ly_half + 1, 1, e))**2 &
               + (this%coef%dof%y(1, ly_half, 1, e) &
               - this%coef%dof%y(1, ly_half + 1, 1, e))**2 &
               + (this%coef%dof%z(1, ly_half, 1, e) &
               - this%coef%dof%z(1, ly_half + 1, 1, e))**2

          dk = (this%coef%dof%x(1, 1, lz_half, e) &
               - this%coef%dof%x(1, 1, lz_half + 1, e))**2 &
               + (this%coef%dof%y(1, 1, lz_half, e) &
               - this%coef%dof%y(1, 1, lz_half + 1, e))**2 &
               + (this%coef%dof%z(1, 1, lz_half, e) &
               - this%coef%dof%z(1, 1, lz_half + 1, e))**2
          di = sqrt(di)
          dj = sqrt(dj)
          dk = sqrt(dk)
          this%delta%x(:,:,:,e) = (di * dj * dk)**(1.0_rp / 3.0_rp)
       end do
    else if (this%delta_type .eq. "elementwise_average") then
       ! use a same length scale throughout an entire element
       ! the length scale is based on (volume)^(1/3)/(N+1)
       do e = 1, this%coef%msh%nelv
          volume_element = 0.0_rp
          do k = 1, this%coef%Xh%lx * this%coef%Xh%ly * this%coef%Xh%lz
             volume_element = volume_element + this%coef%B(k, 1, 1, e)
          end do
          this%delta%x(:,:,:,e) = (volume_element / &
               (this%coef%Xh%lx - 1.0_rp) / &
               (this%coef%Xh%ly - 1.0_rp) / &
               (this%coef%Xh%lz - 1.0_rp) ) ** (1.0_rp / 3.0_rp)
       end do
    else if (this%delta_type .eq. "pointwise") then
       do e = 1, this%coef%msh%nelv
          do k = 1, this%coef%Xh%lz
             km = max(1, k-1)
             kp = min(this%coef%Xh%lz, k+1)

             do j = 1, this%coef%Xh%ly
                jm = max(1, j-1)
                jp = min(this%coef%Xh%ly, j+1)

                do i = 1, this%coef%Xh%lx
                   im = max(1, i-1)
                   ip = min(this%coef%Xh%lx, i+1)

                   di = (this%coef%dof%x(ip, j, k, e) - &
                        this%coef%dof%x(im, j, k, e))**2 &
                        + (this%coef%dof%y(ip, j, k, e) - &
                        this%coef%dof%y(im, j, k, e))**2 &
                        + (this%coef%dof%z(ip, j, k, e) - &
                        this%coef%dof%z(im, j, k, e))**2

                   dj = (this%coef%dof%x(i, jp, k, e) - &
                        this%coef%dof%x(i, jm, k, e))**2 &
                        + (this%coef%dof%y(i, jp, k, e) - &
                        this%coef%dof%y(i, jm, k, e))**2 &
                        + (this%coef%dof%z(i, jp, k, e) - &
                        this%coef%dof%z(i, jm, k, e))**2

                   dk = (this%coef%dof%x(i, j, kp, e) - &
                        this%coef%dof%x(i, j, km, e))**2 &
                        + (this%coef%dof%y(i, j, kp, e) - &
                        this%coef%dof%y(i, j, km, e))**2 &
                        + (this%coef%dof%z(i, j, kp, e) - &
                        this%coef%dof%z(i, j, km, e))**2

                   di = sqrt(di) / (ip - im)
                   dj = sqrt(dj) / (jp - jm)
                   dk = sqrt(dk) / (kp - km)
                   this%delta%x(i,j,k,e) = (di * dj * dk)**(1.0_rp / 3.0_rp)

                end do
             end do
          end do
       end do
    else
       call neko_type_error("delta_type for LES model", &
            this%delta_type, DELTA_KNOWN_TYPES)
       stop
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%delta%x, this%delta%x_d, this%delta%dof%size(),&
            HOST_TO_DEVICE, sync = .false.)
       call this%coef%gs_h%op(this%delta%x, this%delta%dof%size(), GS_OP_ADD)
       call device_col2(this%delta%x_d, this%coef%mult_d, this%delta%dof%size())
    else
       call this%coef%gs_h%op(this%delta%x, this%delta%dof%size(), GS_OP_ADD)
       call col2(this%delta%x, this%coef%mult, this%delta%dof%size())
    end if

  end subroutine les_model_compute_delta

end module les_model
