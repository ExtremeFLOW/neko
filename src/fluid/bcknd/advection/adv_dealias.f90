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
module adv_dealias
  use advection, only : advection_t
  use num_types, only : rp
  use math, only : vdot3, sub2
  use space, only : space_t, GL
  use field, only : field_t
  use coefs, only : coef_t
  use device_math, only : device_vdot3, device_sub2
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX, NEKO_BCKND_XSMM, &
       NEKO_BCKND_OPENCL, NEKO_BCKND_CUDA, NEKO_BCKND_HIP
  use utils, only : neko_error
  use operators, only : opgrad
  use interpolation, only : interpolator_t
  use device, only : device_map, device_get_ptr, device_free
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated
  implicit none
  private

  !> Type encapsulating advection routines with dealiasing
  type, public, extends(advection_t) :: adv_dealias_t
     !> Coeffs of the higher-order space
     type(coef_t) :: coef_GL
     !> Coeffs of the original space in the simulation
     type(coef_t), pointer :: coef_GLL
     !> Interpolator between the original and higher-order spaces
     type(interpolator_t) :: GLL_to_GL
     !> The additional higher-order space used in dealiasing
     type(space_t) :: Xh_GL
     !> The original space used in the simulation
     type(space_t), pointer :: Xh_GLL
     real(kind=rp), allocatable :: temp(:), tbf(:)
     !> Temporary arrays
     real(kind=rp), allocatable :: tx(:), ty(:), tz(:)
     real(kind=rp), allocatable :: vr(:), vs(:), vt(:)
     !> Device pointer for `temp`
     type(c_ptr) :: temp_d = C_NULL_PTR
     !> Device pointer for `tbf`
     type(c_ptr) :: tbf_d = C_NULL_PTR
     !> Device pointer for `tx`
     type(c_ptr) :: tx_d = C_NULL_PTR
     !> Device pointer for `ty`
     type(c_ptr) :: ty_d = C_NULL_PTR
     !> Device pointer for `tz`
     type(c_ptr) :: tz_d = C_NULL_PTR
     !> Device pointer for `vr`
     type(c_ptr) :: vr_d = C_NULL_PTR
     !> Device pointer for `vs`
     type(c_ptr) :: vs_d = C_NULL_PTR
     !> Device pointer for `vt`
     type(c_ptr) :: vt_d = C_NULL_PTR

   contains
     !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
     !! the RHS.
     procedure, pass(this) :: compute => compute_advection_dealias
     !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to
     !! the RHS.
     procedure, pass(this) :: compute_scalar => compute_scalar_advection_dealias
     !> Add the advection term in ALE framework.
     procedure, pass(this) :: compute_ale => compute_ale_advection_dealias
     !> Update any metrics needed for the advection computation in ALE.
     procedure, pass(this) :: recompute_metrics => recompute_metrics_dealias
     !> Constructor
     procedure, pass(this) :: init => init_dealias
     !> Destructor
     procedure, pass(this) :: free => free_dealias
  end type adv_dealias_t

contains

  !> Constructor
  !! @param lxd The polynomial order of the space used in the dealiasing.
  !! @param coef The coefficients of the (space, mesh) pair.
  subroutine init_dealias(this, lxd, coef)
    class(adv_dealias_t), target, intent(inout) :: this
    integer, intent(in) :: lxd
    type(coef_t), intent(inout), target :: coef
    integer :: nel, n_GL, n

    call this%Xh_GL%init(GL, lxd, lxd, lxd)
    this%Xh_GLL => coef%Xh
    this%coef_GLL => coef
    call this%GLL_to_GL%init(this%Xh_GL, this%Xh_GLL)

    call this%coef_GL%init(this%Xh_GL, coef%msh)

    nel = coef%msh%nelv
    n_GL = nel*this%Xh_GL%lxyz
    n = nel*coef%Xh%lxyz
    call this%GLL_to_GL%map(this%coef_GL%drdx, coef%drdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdx, coef%dsdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdx, coef%dtdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%drdy, coef%drdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdy, coef%dsdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdy, coef%dtdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%drdz, coef%drdz, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdz, coef%dsdz, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdz, coef%dtdz, nel, this%Xh_GL)
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1) .or. (NEKO_BCKND_SX .eq. 1) .or. &
         (NEKO_BCKND_XSMM .eq. 1)) then
       allocate(this%temp(n_GL))
       allocate(this%tbf(n_GL))
       allocate(this%tx(n_GL))
       allocate(this%ty(n_GL))
       allocate(this%tz(n_GL))
       allocate(this%vr(n_GL))
       allocate(this%vs(n_GL))
       allocate(this%vt(n_GL))
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%temp, this%temp_d, n_GL)
       call device_map(this%tbf, this%tbf_d, n_GL)
       call device_map(this%tx, this%tx_d, n_GL)
       call device_map(this%ty, this%ty_d, n_GL)
       call device_map(this%tz, this%tz_d, n_GL)
       call device_map(this%vr, this%vr_d, n_GL)
       call device_map(this%vs, this%vs_d, n_GL)
       call device_map(this%vt, this%vt_d, n_GL)
    end if

  end subroutine init_dealias

  !> Destructor
  subroutine free_dealias(this)
    class(adv_dealias_t), intent(inout) :: this

    if (allocated(this%temp)) then
       deallocate(this%temp)
    end if

    if (allocated(this%tbf)) then
       deallocate(this%tbf)
    end if
    if (allocated(this%tx)) then
       deallocate(this%tx)
    end if
    if (allocated(this%ty)) then
       deallocate(this%ty)
    end if
    if (allocated(this%tz)) then
       deallocate(this%tz)
    end if
    if (allocated(this%vr)) then
       deallocate(this%vr)
    end if
    if (allocated(this%vs)) then
       deallocate(this%vs)
    end if
    if (allocated(this%vt)) then
       deallocate(this%vt)
    end if

    if (c_associated(this%temp_d)) then
       call device_free(this%temp_d)
    end if

    if (c_associated(this%tbf_d)) then
       call device_free(this%tbf_d)
    end if

    if (c_associated(this%tx_d)) then
       call device_free(this%tx_d)
    end if

    if (c_associated(this%ty_d)) then
       call device_free(this%ty_d)
    end if

    if (c_associated(this%tz_d)) then
       call device_free(this%tz_d)
    end if

    if (c_associated(this%vr_d)) then
       call device_free(this%vr_d)
    end if

    if (c_associated(this%vs_d)) then
       call device_free(this%vs_d)
    end if

    if (c_associated(this%vt_d)) then
       call device_free(this%vt_d)
    end if

    call this%coef_GL%free()
    call this%GLL_to_GL%free()
    call this%Xh_GL%free()

    nullify(this%Xh_GLL)
    nullify(this%coef_GLL)

  end subroutine free_dealias


  !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
  !! the RHS.
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
  subroutine compute_advection_dealias(this, vx, vy, vz, fx, fy, fz, Xh, &
       coef, n, dt)
    class(adv_dealias_t), intent(inout) :: this
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: fx, fy, fz
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt

    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tx, ty, tz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tfx, tfy, tfz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: vr, vs, vt
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: tempx, tempy, tempz
    integer :: e, i, idx, nel, n_GL

    nel = coef%msh%nelv
    n_GL = nel * this%Xh_GL%lxyz

    !This is extremely primitive and unoptimized  on the device //Karp
    associate(c_GL => this%coef_GL)
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call this%GLL_to_GL%map(this%tx, vx%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%ty, vy%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tz, vz%x, nel, this%Xh_GL)

         call opgrad(this%vr, this%vs, this%vt, this%tx, c_GL)
         call device_vdot3(this%tbf_d, this%vr_d, this%vs_d, this%vt_d, &
              this%tx_d, this%ty_d, this%tz_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call device_sub2(fx%x_d, this%temp_d, n)


         call opgrad(this%vr, this%vs, this%vt, this%ty, c_GL)
         call device_vdot3(this%tbf_d, this%vr_d, this%vs_d, this%vt_d, &
              this%tx_d, this%ty_d, this%tz_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call device_sub2(fy%x_d, this%temp_d, n)

         call opgrad(this%vr, this%vs, this%vt, this%tz, c_GL)
         call device_vdot3(this%tbf_d, this%vr_d, this%vs_d, this%vt_d, &
              this%tx_d, this%ty_d, this%tz_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call device_sub2(fz%x_d, this%temp_d, n)

      else if ((NEKO_BCKND_SX .eq. 1) .or. (NEKO_BCKND_XSMM .eq. 1)) then

         call this%GLL_to_GL%map(this%tx, vx%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%ty, vy%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tz, vz%x, nel, this%Xh_GL)

         call opgrad(this%vr, this%vs, this%vt, this%tx, c_GL)
         call vdot3(this%tbf, this%vr, this%vs, this%vt, &
              this%tx, this%ty, this%tz, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call sub2(fx%x, this%temp, n)


         call opgrad(this%vr, this%vs, this%vt, this%ty, c_GL)
         call vdot3(this%tbf, this%vr, this%vs, this%vt, &
              this%tx, this%ty, this%tz, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call sub2(fy%x, this%temp, n)

         call opgrad(this%vr, this%vs, this%vt, this%tz, c_GL)
         call vdot3(this%tbf, this%vr, this%vs, this%vt, &
              this%tx, this%ty, this%tz, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call sub2(fz%x, this%temp, n)

      else
         !$omp parallel do private(e, i, idx)
         do e = 1, coef%msh%nelv
            call this%GLL_to_GL%map(tx, vx%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(ty, vy%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(tz, vz%x(1,1,1,e), 1, this%Xh_GL)

            call opgrad(vr, vs, vt, tx, c_GL, e, e)
            do i = 1, this%Xh_GL%lxyz
               tfx(i) = tx(i)*vr(i) + ty(i)*vs(i) + tz(i)*vt(i)
            end do

            call opgrad(vr, vs, vt, ty, c_GL, e, e)
            do i = 1, this%Xh_GL%lxyz
               tfy(i) = tx(i)*vr(i) + ty(i)*vs(i) + tz(i)*vt(i)
            end do

            call opgrad(vr, vs, vt, tz, c_GL, e, e)
            do i = 1, this%Xh_GL%lxyz
               tfz(i) = tx(i)*vr(i) + ty(i)*vs(i) + tz(i)*vt(i)
            end do

            call this%GLL_to_GL%map(tempx, tfx, 1, this%Xh_GLL)
            call this%GLL_to_GL%map(tempy, tfy, 1, this%Xh_GLL)
            call this%GLL_to_GL%map(tempz, tfz, 1, this%Xh_GLL)

            idx = (e-1)*this%Xh_GLL%lxyz+1
            do concurrent (i = 0:this%Xh_GLL%lxyz-1)
               fx%x(i+idx,1,1,1) = fx%x(i+idx,1,1,1) - tempx(i+1)
               fy%x(i+idx,1,1,1) = fy%x(i+idx,1,1,1) - tempy(i+1)
               fz%x(i+idx,1,1,1) = fz%x(i+idx,1,1,1) - tempz(i+1)
            end do
         end do
         !$omp end parallel do
      end if
    end associate

  end subroutine compute_advection_dealias

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
  subroutine compute_scalar_advection_dealias(this, vx, vy, vz, s, fs, Xh, &
       coef, n, dt)
    class(adv_dealias_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: s
    type(field_t), intent(inout) :: fs
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt

    real(kind=rp), dimension(this%Xh_GL%lxyz) :: vx_GL, vy_GL, vz_GL, s_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: dsdx, dsdy, dsdz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: f_GL
    integer :: e, i, idx, nel, n_GL
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: temp

    nel = coef%msh%nelv
    n_GL = nel * this%Xh_GL%lxyz

    associate(c_GL => this%coef_GL)
      if (NEKO_BCKND_DEVICE .eq. 1) then

         ! Map advecting velocity onto the higher-order space
         call this%GLL_to_GL%map(this%tx, vx%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%ty, vy%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tz, vz%x, nel, this%Xh_GL)

         ! Map the scalar onto the high-order space
         call this%GLL_to_GL%map(this%temp, s%x, nel, this%Xh_GL)

         ! Compute the scalar gradient in the high-order space
         call opgrad(this%vr, this%vs, this%vt, this%temp, c_GL)

         ! Compute the convective term, i.e dot the velocity with the scalar grad
         call device_vdot3(this%tbf_d, this%vr_d, this%vs_d, this%vt_d, &
              this%tx_d, this%ty_d, this%tz_d, n_GL)

         ! Map back to the original space (we reuse this%temp)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)

         ! Update the source term
         call device_sub2(fs%x_d, this%temp_d, n)

      else if ((NEKO_BCKND_SX .eq. 1) .or. (NEKO_BCKND_XSMM .eq. 1)) then

         ! Map advecting velocity onto the higher-order space
         call this%GLL_to_GL%map(this%tx, vx%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%ty, vy%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tz, vz%x, nel, this%Xh_GL)

         ! Map the scalar onto the high-order space
         call this%GLL_to_GL%map(this%temp, s%x, nel, this%Xh_GL)

         ! Compute the scalar gradient in the high-order space
         call opgrad(this%vr, this%vs, this%vt, this%temp, c_GL)

         ! Compute the convective term, i.e dot the velocity with the scalar grad
         call vdot3(this%tbf, this%vr, this%vs, this%vt, &
              this%tx, this%ty, this%tz, n_GL)

         ! Map back to the original space (we reuse this%temp)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)

         ! Update the source term
         call sub2(fs%x, this%temp, n)

      else
         !$omp parallel do private (e, i, idx)
         do e = 1, coef%msh%nelv
            ! Map advecting velocity onto the higher-order space
            call this%GLL_to_GL%map(vx_GL, vx%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(vy_GL, vy%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(vz_GL, vz%x(1,1,1,e), 1, this%Xh_GL)

            ! Map scalar onto the higher-order space
            call this%GLL_to_GL%map(s_GL, s%x(1,1,1,e), 1, this%Xh_GL)

            ! Gradient of s in the higher-order space
            call opgrad(dsdx, dsdy, dsdz, s_GL, c_GL, e, e)

            ! vx * ds/dx + vy * ds/dy + vz * ds/dz for each point in the element
            do i = 1, this%Xh_GL%lxyz
               f_GL(i) = vx_GL(i)*dsdx(i) + vy_GL(i)*dsdy(i) + vz_GL(i)*dsdz(i)
            end do

            ! Map back the contructed operator to the original space
            call this%GLL_to_GL%map(temp, f_GL, 1, this%Xh_GLL)

            idx = (e-1)*this%Xh_GLL%lxyz + 1

            call sub2(fs%x(idx, 1, 1, 1), temp, this%Xh_GLL%lxyz)
         end do
         !$omp end parallel do
      end if
    end associate

  end subroutine compute_scalar_advection_dealias


  !!> Add the advection term in ALE framework.
  !! @param this The object.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param wm_x The x component of mesh velocity.
  !! @param wm_y The y component of mesh velocity.
  !! @param wm_z The z component of mesh velocity.
  !! @param fx The x component of source term.
  !! @param fy The y component of source term.
  !! @param fz The z component of source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step, not required for this method.
  !! Here, we compute: - div ( u_i * wm ).
  !! Based on Ho, L.W. A Legendre spectral element method for simulation of incompressible
  !! unsteady viscous free-surface flows.
  !! Ph.D. thesis, Massachusetts Institute of Technology, 1989.
  !! Note: In Nek5000, dealiasing is not done for this term.
  subroutine compute_ale_advection_dealias(this, vx, vy, vz, wm_x, wm_y, wm_z, &
       fx, fy, fz, Xh, coef, n, dt)
    class(adv_dealias_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: wm_x, wm_y, wm_z
    type(field_t), intent(inout) :: fx, fy, fz
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: vx_GL, vy_GL, vz_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: wm_x_GL, wm_y_GL, wm_z_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: flux_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: grad_x, grad_y, grad_z
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: total_div_GL
    integer :: e, i, idx, nel, n_GL
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: temp_x, temp_y, temp_z
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt

    nel = coef%msh%nelv
    n_GL = nel * this%Xh_GL%lxyz

    associate(c_GL => this%coef_GL)
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call neko_error("ALE advection with dealiasing not " // &
              "implemented yet for device")
      else if ((NEKO_BCKND_SX .eq. 1) .or. (NEKO_BCKND_XSMM .eq. 1)) then
         call neko_error("ALE advection with dealiasing not " // &
              "implemented yet for device")
      else
         !$omp parallel do private(e, i, flux_GL, total_div_GL, idx)
         do e = 1, coef%msh%nelv
            ! Map advecting velocity and mesh velocity onto the higher-order space
            call this%GLL_to_GL%map(vx_GL, vx%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(vy_GL, vy%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(vz_GL, vz%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(wm_x_GL, wm_x%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(wm_y_GL, wm_y%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(wm_z_GL, wm_z%x(1,1,1,e), 1, this%Xh_GL)

            ! x-momentum
            total_div_GL = 0.0_rp
            ! I think below can be written more efficiently. Will fix it later.
            ! This works for now.

            ! div(u * wm_*) = d/dx (u * wm_x) + d/dy (u * wm_y) + d/dz (u * wm_z)

            flux_GL = vx_GL * wm_x_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_x
            flux_GL = vx_GL * wm_y_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_y
            flux_GL = vx_GL * wm_z_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_z

            ! Map back the contructed operator to the original space
            call this%GLL_to_GL%map(temp_x, total_div_GL, 1, this%Xh_GLL)

            ! y-momentum
            total_div_GL = 0.0_rp
            ! div(v * wm_*) = d/dx (v * wm_x) + d/dy (v * wm_y) + d/dz (v * wm_z)

            flux_GL = vy_GL * wm_x_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_x
            flux_GL = vy_GL * wm_y_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_y
            flux_GL = vy_GL * wm_z_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_z

            ! Map back the contructed operator to the original space
            call this%GLL_to_GL%map(temp_y, total_div_GL, 1, this%Xh_GLL)

            ! z-momentum
            total_div_GL = 0.0_rp
            ! div(w * wm_*) = d/dx (w * wm_x) + d/dy (w * wm_y) + d/dz (w * wm_z)

            flux_GL = vz_GL * wm_x_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_x
            flux_GL = vz_GL * wm_y_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_y
            flux_GL = vz_GL * wm_z_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_z

            ! Map back the contructed operator to the original space
            call this%GLL_to_GL%map(temp_z, total_div_GL, 1, this%Xh_GLL)

            ! Note we add (+) here since the ALE advection term is
            ! - div(u * wm) on the LHS. So on the RHS it will be + div(u * wm)

            idx = (e-1)*this%Xh_GLL%lxyz+1
            do concurrent (i = 0:this%Xh_GLL%lxyz-1)
               fx%x(i+idx,1,1,1) = fx%x(i+idx,1,1,1) + temp_x(i+1)
               fy%x(i+idx,1,1,1) = fy%x(i+idx,1,1,1) + temp_y(i+1)
               fz%x(i+idx,1,1,1) = fz%x(i+idx,1,1,1) + temp_z(i+1)
            end do
         end do
        !$omp end parallel do
      end if
    end associate

  end subroutine compute_ale_advection_dealias

  subroutine recompute_metrics_dealias(this, coef, moving_boundary)
    class(adv_dealias_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    logical, intent(in) :: moving_boundary
    integer :: nel

    if (.not. moving_boundary) return

    nel = coef%msh%nelv
    call this%GLL_to_GL%map(this%coef_GL%drdx, coef%drdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdx, coef%dsdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdx, coef%dtdx, nel, this%Xh_GL)

    call this%GLL_to_GL%map(this%coef_GL%drdy, coef%drdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdy, coef%dsdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdy, coef%dtdy, nel, this%Xh_GL)

    call this%GLL_to_GL%map(this%coef_GL%drdz, coef%drdz, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdz, coef%dsdz, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdz, coef%dtdz, nel, this%Xh_GL)

  end subroutine recompute_metrics_dealias


end module adv_dealias
