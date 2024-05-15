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
module advection_no_dealias
  use advection, only: advection_t
  use num_types, only: rp
  use math, only: vdot3, sub2, subcol3, rzero
  use utils, only: neko_error
  use space, only: space_t, GL
  use field, only: field_t
  use coefs, only: coef_t
  use neko_config, only: NEKO_BCKND_DEVICE, NEKO_BCKND_SX, NEKO_BCKND_XSMM, &
    NEKO_BCKND_OPENCL, NEKO_BCKND_CUDA, NEKO_BCKND_HIP
  use operators, only: opgrad, conv1, cdtp
  use interpolation, only: interpolator_t
  use device_math
  use device, only: device_free, device_map, device_get_ptr, device_memcpy, &
    HOST_TO_DEVICE
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR, &
    c_associated
  implicit none
  private

  !> Type encapsulating advection routines with no dealiasing applied
  type, public, extends(advection_t) :: adv_no_dealias_t
     real(kind=rp), allocatable :: temp(:)
     type(c_ptr) :: temp_d = C_NULL_PTR
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
     !> Add the advection term for perturbations
     ! could be either LNS or adjoint
     !! the RHS
     procedure, pass(this) :: compute_vector => &
       compute_vector_advection_no_dealias
  end type adv_no_dealias_t

contains

  !> Constructor
  !! @param coef The coefficients of the (space, mesh) pair.
  subroutine init_no_dealias(this, coef)
    class(adv_no_dealias_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

    allocate(this%temp(coef%dof%size()))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%temp, this%temp_d, coef%dof%size())
    end if

  end subroutine init_no_dealias

  !> Destructor
  subroutine free_no_dealias(this)
    class(adv_no_dealias_t), intent(inout) :: this

    if (allocated(this%temp)) then
       deallocate(this%temp)
    end if
    if (c_associated(this%temp_d)) then
       call device_free(this%temp_d)
    end if
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
  subroutine compute_advection_no_dealias(this, vx, vy, vz, fx, fy, fz, Xh, coef, n)
    class(adv_no_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: fx, fy, fz
    type(c_ptr) :: fx_d, fy_d, fz_d

    if (NEKO_BCKND_DEVICE .eq. 1) then
       fx_d = device_get_ptr(fx)
       fy_d = device_get_ptr(fy)
       fz_d = device_get_ptr(fz)

       call conv1(this%temp, vx%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fx_d, coef%B_d, this%temp_d, n)
       call conv1(this%temp, vy%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fy_d, coef%B_d, this%temp_d, n)
       if (coef%Xh%lz .eq. 1) then
          call device_rzero (this%temp_d, n)
       else
          call conv1(this%temp, vz%x, vx%x, vy%x, vz%x, Xh, coef)
          call device_subcol3(fz_d, coef%B_d, this%temp_d, n)
       end if
    else
       call conv1(this%temp, vx%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (fx, coef%B, this%temp, n)
       call conv1(this%temp, vy%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (fy, coef%B, this%temp, n)
       if (coef%Xh%lz .eq. 1) then
          call rzero (this%temp, n)
       else
          call conv1(this%temp, vz%x, vx%x, vy%x, vz%x, Xh, coef)
          call subcol3(fz, coef%B, this%temp, n)
       end if
    end if

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
  subroutine compute_scalar_advection_no_dealias(this, vx, vy, vz, s, fs, Xh, &
                                                 coef, n)
    class(adv_no_dealias_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: s
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: fs
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(c_ptr) :: fs_d

    if (NEKO_BCKND_DEVICE .eq. 1) then
       fs_d = device_get_ptr(fs)

       call conv1(this%temp, s%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fs_d, coef%B_d, this%temp_d, n)
       if (coef%Xh%lz .eq. 1) then
          call device_rzero (this%temp_d, n)
       end if
    else
       ! temp will hold vx*ds/dx + vy*ds/dy + vz*ds/sz
       call conv1(this%temp, s%x, vx%x, vy%x, vz%x, Xh, coef)

       ! fs = fs - B*temp
       call subcol3 (fs, coef%B, this%temp, n)
       if (coef%Xh%lz .eq. 1) then
          call rzero (this%temp, n)
       end if
    end if

  end subroutine compute_scalar_advection_no_dealias


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
  subroutine compute_vector_advection_no_dealias(this, vx, vy, vz, vxb, vyb, vzb, fx, fy, fz, Xh, coef, n)
    implicit none
    class(adv_no_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: vxb, vyb, vzb
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: fx, fy, fz

    logical, parameter :: use_adjoint = .false.

    ! Linearized advection term for the fluid
    if (use_adjoint) then
       call compute_adjoint_advection_no_dealias(this, vx, vy, vz, vxb, vyb, &
                                                 vzb, fx, fy, fz, Xh, coef, n)
    else
       call compute_LNS_advection_no_dealias(this, vx, vy, vz, vxb, vyb, vzb, &
                                             fx, fy, fz, Xh, coef, n)
    end if
  end subroutine compute_vector_advection_no_dealias




! THIS IS NOTHING
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
  subroutine compute_adjoint_advection_no_dealias(this, vx, vy, vz, vxb, vyb, vzb, fx, fy, fz, Xh, coef, n)
    implicit none
    class(adv_no_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: vxb, vyb, vzb
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: fx, fy, fz
    type(c_ptr) :: fx_d, fy_d, fz_d

    if (NEKO_BCKND_DEVICE .eq. 1) then
       fx_d = device_get_ptr(fx)
       fy_d = device_get_ptr(fy)
       fz_d = device_get_ptr(fz)

       call conv1(this%temp, vx%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fx_d, coef%B_d, this%temp_d, n)
       call conv1(this%temp, vy%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fy_d, coef%B_d, this%temp_d, n)
       if (coef%Xh%lz .eq. 1) then
          call device_rzero (this%temp_d, n)
       else
          call conv1(this%temp, vz%x, vx%x, vy%x, vz%x, Xh, coef)
          call device_subcol3(fz_d, coef%B_d, this%temp_d, n)
       end if
    else
       call conv1(this%temp, vx%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (fx, coef%B, this%temp, n)
       call conv1(this%temp, vy%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (fy, coef%B, this%temp, n)
       if (coef%Xh%lz .eq. 1) then
          call rzero (this%temp, n)
       else
          call conv1(this%temp, vz%x, vx%x, vy%x, vz%x, Xh, coef)
          call subcol3(fz, coef%B, this%temp, n)
       end if
    end if

  end subroutine compute_adjoint_advection_no_dealias

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
  subroutine compute_LNS_advection_no_dealias(this, vx, vy, vz, vxb, vyb, vzb, fx, fy, fz, Xh, coef, n)
    implicit none
    class(adv_no_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: vxb, vyb, vzb
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: fx, fy, fz
    type(c_ptr) :: fx_d, fy_d, fz_d

    call neko_error('compute_LNS_advection_no_dealias: Not implemented')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       fx_d = device_get_ptr(fx)
       fy_d = device_get_ptr(fy)
       fz_d = device_get_ptr(fz)

       call conv1(this%temp, vx%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fx_d, coef%B_d, this%temp_d, n)
       call conv1(this%temp, vy%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fy_d, coef%B_d, this%temp_d, n)
       if (coef%Xh%lz .eq. 1) then
          call device_rzero (this%temp_d, n)
       else
          call conv1(this%temp, vz%x, vx%x, vy%x, vz%x, Xh, coef)
          call device_subcol3(fz_d, coef%B_d, this%temp_d, n)
       end if
    else
       call conv1(this%temp, vx%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (fx, coef%B, this%temp, n)
       call conv1(this%temp, vy%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (fy, coef%B, this%temp, n)
       if (coef%Xh%lz .eq. 1) then
          call rzero (this%temp, n)
       else
          call conv1(this%temp, vz%x, vx%x, vy%x, vz%x, Xh, coef)
          call subcol3(fz, coef%B, this%temp, n)
       end if
    end if

  end subroutine compute_LNS_advection_no_dealias
end module advection_no_dealias
