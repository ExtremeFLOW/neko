! Copyright (c) 2021-2024, The Neko Authors
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
!> Subroutines to add advection terms to the RHS of a transport equation.
module adv_no_dealias
  use advection, only : advection_t
  use num_types, only : rp
  use math, only : subcol3, rzero, addcol3
  use field_math, only : field_col3
  use space, only : space_t
  use field, only : field_t
  use coefs, only : coef_t
  use device_math, only : device_subcol3, device_rzero, device_addcol3
  use neko_config, only : NEKO_BCKND_DEVICE
  use operators, only : conv1, div
  use device, only : device_free, device_map, device_get_ptr
  use scratch_registry, only : neko_scratch_registry
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated
  implicit none
  private

  !> Type encapsulating advection routines with no dealiasing applied
  type, public, extends(advection_t) :: adv_no_dealias_t
   contains
     !> Constructor
     procedure, pass(this) :: init => init_no_dealias
     !> Destructor
     procedure, pass(this) :: free => free_no_dealias
     !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
     !! the RHS
     procedure, pass(this) :: compute => compute_advection_no_dealias
     !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to
     !! the RHS
     procedure, pass(this) :: compute_scalar => &
          compute_scalar_advection_no_dealias
     !> add the advection term for ALE, i.e. \f$(u - w_m) \cdot \nabla s \f$, to
     !! the RHS
     procedure, pass(this) :: compute_ale => compute_ale_advection_no_dealias
     !> Update any metrics needed for the advection computation in ALE.
     procedure, pass(this) :: recompute_metrics => recompute_metrics_no_dealias
  end type adv_no_dealias_t

contains

  !> Constructor
  !! @param coef The coefficients of the (space, mesh) pair.
  subroutine init_no_dealias(this, coef)
    class(adv_no_dealias_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

  end subroutine init_no_dealias

  !> Destructor
  subroutine free_no_dealias(this)
    class(adv_no_dealias_t), intent(inout) :: this

  end subroutine free_no_dealias

  !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$ to the
  !! RHS.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param fx The x component of source term.
  !! @param fy The y component of source term.
  !! @param fz The z component of source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step, not required for this method.
  subroutine compute_advection_no_dealias(this, vx, vy, vz, fx, fy, fz, Xh, &
       coef, n, dt)
    class(adv_no_dealias_t), intent(inout) :: this
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: fx, fy, fz
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt

    type(field_t), pointer :: temp
    integer :: id_temp

    call neko_scratch_registry%request_field(temp, id_temp, .false.)

    if (NEKO_BCKND_DEVICE .eq. 1) then

       call conv1(temp%x, vx%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fx%x_d, coef%B_d, temp%x_d, n)
       call conv1(temp%x, vy%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fy%x_d, coef%B_d, temp%x_d, n)
       if (coef%Xh%lz .eq. 1) then
          call device_rzero (temp%x_d, n)
       else
          call conv1(temp%x, vz%x, vx%x, vy%x, vz%x, Xh, coef)
          call device_subcol3(fz%x_d, coef%B_d, temp%x_d, n)
       end if
    else
       call conv1(temp%x, vx%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (fx%x, coef%B, temp%x, n)
       call conv1(temp%x, vy%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (fy%x, coef%B, temp%x, n)
       if (coef%Xh%lz .eq. 1) then
          call rzero (temp%x, n)
       else
          call conv1(temp%x, vz%x, vx%x, vy%x, vz%x, Xh, coef)
          call subcol3(fz%x, coef%B, temp%x, n)
       end if
    end if

    call neko_scratch_registry%relinquish_field(id_temp)
  end subroutine compute_advection_no_dealias

  !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to the
  !! RHS.
  !! @param this The object.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param s The scalar.
  !! @param fs The source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step, not required for this method.
  subroutine compute_scalar_advection_no_dealias(this, vx, vy, vz, s, fs, Xh, &
       coef, n, dt)
    class(adv_no_dealias_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: s
    type(field_t), intent(inout) :: fs
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt

    type(field_t), pointer :: temp
    integer :: id_temp
    call neko_scratch_registry%request_field(temp, id_temp, .false.)

    if (NEKO_BCKND_DEVICE .eq. 1) then

       call conv1(temp%x, s%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fs%x_d, coef%B_d, temp%x_d, n)
       if (coef%Xh%lz .eq. 1) then
          call device_rzero (temp%x_d, n)
       end if
    else
       ! temp will hold vx*ds/dx + vy*ds/dy + vz*ds/sz
       call conv1(temp%x, s%x, vx%x, vy%x, vz%x, Xh, coef)

       ! fs = fs - B*temp
       call subcol3 (fs%x, coef%B, temp%x, n)
       if (coef%Xh%lz .eq. 1) then
          call rzero (temp%x, n)
       end if
    end if

    call neko_scratch_registry%relinquish_field(id_temp)
  end subroutine compute_scalar_advection_no_dealias

  subroutine recompute_metrics_no_dealias(this, coef, moving_boundary)
    class(adv_no_dealias_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    logical, intent(in) :: moving_boundary
    ! no-op
  end subroutine recompute_metrics_no_dealias

  subroutine compute_ale_advection_no_dealias(this, vx, vy, vz, &
       wm_x, wm_y, wm_z, fx, fy, fz, Xh, coef, n, dt)
    class(adv_no_dealias_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: wm_x, wm_y, wm_z
    type(field_t), intent(inout) :: fx, fy, fz
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt
    type(field_t), pointer :: temp, work_x, work_y, work_z
    integer :: id(4)

    call neko_scratch_registry%request_field(temp, id(1), .false.)
    call neko_scratch_registry%request_field(work_x, id(2), .false.)
    call neko_scratch_registry%request_field(work_y, id(3), .false.)
    call neko_scratch_registry%request_field(work_z, id(4), .false.)

    ! u.wm_*
    call field_col3(work_x, vx, wm_x)
    call field_col3(work_y, vx, wm_y)
    call field_col3(work_z, vx, wm_z)
    call div(temp%x_d, work_x%x_d, work_y%x_d, work_z%x_d, coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_addcol3(fx%x_d, coef%B_d, temp%x_d, n)
    else
       call addcol3(fx%x, coef%B, temp%x, n)
    end if

    ! v.wm_*
    call field_col3(work_x, vy, wm_x)
    call field_col3(work_y, vy, wm_y)
    call field_col3(work_z, vy, wm_z)
    call div(temp%x_d, work_x%x_d, work_y%x_d, work_z%x_d, coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_addcol3(fy%x_d, coef%B_d, temp%x_d, n)
    else
       call addcol3(fy%x, coef%B, temp%x, n)
    end if

    ! w.wm_*
    call field_col3(work_x, vz, wm_x)
    call field_col3(work_y, vz, wm_y)
    call field_col3(work_z, vz, wm_z)
    call div(temp%x_d, work_x%x_d, work_y%x_d, work_z%x_d, coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_addcol3(fz%x_d, coef%B_d, temp%x_d, n)
    else
       call addcol3(fz%x, coef%B, temp%x, n)
    end if

    call neko_scratch_registry%relinquish_field(id)

  end subroutine compute_ale_advection_no_dealias
end module adv_no_dealias
