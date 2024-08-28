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
  use field, only : field_t, field_ptr_t
  use json_module, only : json_file
  use field_registry, only : neko_field_registry
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use gs_ops, only : GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_memcpy, HOST_TO_DEVICE
  use math, only : col2
  use device_math, only : device_col2
  implicit none
  private

  !> Base abstract type for LES models based on the Boussinesq approximation.
  type, abstract, public :: les_model_t
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
     !! @param dofmap SEM map of degrees of freedom.
     !! @param coef SEM coefficients.
     !! @param json A dictionary with parameters.
     subroutine les_model_init(this, dofmap, coef, json)
       import les_model_t, json_file, dofmap_t, coef_t
       class(les_model_t), intent(inout) :: this
       type(coef_t), intent(in) :: coef
       type(dofmap_t), intent(in) :: dofmap
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
     !> LES model factory. Both constructs and initializes the object.
     !! @param object The object to be allocated.
     !! @param type_name The name of the LES model.
     !! @param dofmap SEM map of degrees of freedom.
     !! @param coef SEM coefficients.
     !! @param json A dictionary with parameters.
     module subroutine les_model_factory(object, type_name, dofmap, coef, json)
       class(les_model_t), allocatable, target, intent(inout) :: object
       character(len=*), intent(in) :: type_name
       type(dofmap_t), intent(in) :: dofmap
       type(coef_t), intent(in) :: coef
       type(json_file), intent(inout) :: json
     end subroutine les_model_factory
  end interface

  public :: les_model_factory
  
contains
  !> Constructor for the les_model_t (base) class.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  !! @param nu_name The name of the turbulent viscosity field.
  subroutine les_model_init_base(this, dofmap, coef, nut_name, delta_type)
    class(les_model_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), target, intent(in) :: coef
    character(len=*), intent(in) :: nut_name
    character(len=*), intent(in) :: delta_type

    if (.not. neko_field_registry%field_exists(trim(nut_name))) then
       call neko_field_registry%add_field(dofmap, trim(nut_name))
    end if
    if (.not. neko_field_registry%field_exists("les_delta")) then
       call neko_field_registry%add_field(dofmap, "les_delta")
    end if
    this%nut => neko_field_registry%get_field(trim(nut_name))
    this%delta => neko_field_registry%get_field("les_delta")
    this%coef => coef
    this%delta_type = delta_type

    call this%compute_delta()
  end subroutine les_model_init_base

  !> Destructor for the les_model_t (base) class.
  subroutine les_model_free_base(this)
    class(les_model_t), intent(inout) :: this

    nullify(this%nut)
    nullify(this%delta)
    nullify(this%coef)
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
    real(kind=rp) :: di, dj, dk, ndim_inv
    integer :: lx_half, ly_half, lz_half

    lx_half = this%coef%Xh%lx / 2
    ly_half = this%coef%Xh%ly / 2
    lz_half = this%coef%Xh%lz / 2

    if (this%delta_type .eq. "elementwise") then
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
