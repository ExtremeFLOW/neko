! Copyright (c) 2021-2026, The Neko Authors
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
  use device_math, only : device_vdot3, device_sub2, device_col3, device_add2
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX, NEKO_BCKND_XSMM, &
       NEKO_BCKND_OPENCL, NEKO_BCKND_CUDA, NEKO_BCKND_HIP, NEKO_BCKND_METAL
  use utils, only : neko_error
  use operators, only : opgrad
  use interpolation, only : interpolator_t
  use device, only : device_map, device_get_ptr, device_unmap
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  implicit none
  private

  !> Shared, reference-counted state for every dealiased advection instance
  !! built on the same (lxd, coef) pair.
  !!
  !! The fluid and every scalar construct their own advection object, but the
  !! GL space and the geometry mapped onto it are identical between them, and
  !! the work arrays are only live inside a single compute() call. One node
  !! therefore serves all of them: at lx = 8 that is 30 field units of
  !! duplicated geometry and 27 of duplicated scratch reclaimed per scalar.
  !!
  !! Nodes are individually allocated and held in a linked list so that the
  !! pointers handed out to instances stay valid as the list grows.
  type :: adv_dealias_scratch_t
     !> Order of the dealiasing space, part of the lookup key.
     integer :: lxd = -1
     !> The fine-level coefficients this node was built from, the other half
     !! of the lookup key. Also the only coef whose geometry it holds.
     type(coef_t), pointer :: key_coef => null()
     !> Number of adv_dealias_t instances pointing at this node.
     integer :: refcount = 0
     !> Set while a compute() is using the work arrays. Sharing them is only
     !! sound because the schemes are advanced one at a time; this turns a
     !! violation of that into an error instead of a silent data race.
     logical :: busy = .false.
     !> `key_coef%metrics_version` at the last geometry refresh, so that a
     !! moving mesh re-maps the derivative arrays once rather than once per
     !! scheme sharing this node.
     integer :: metrics_version = -1
     !> The additional higher-order space used in dealiasing
     type(space_t) :: Xh_GL
     !> Coeffs of the higher-order space
     type(coef_t) :: coef_GL
     !> Interpolator between the original and higher-order spaces
     type(interpolator_t) :: GLL_to_GL
     real(kind=rp), allocatable :: temp(:), tbf(:)
     !> Temporary arrays
     real(kind=rp), allocatable :: tx(:), ty(:), tz(:)
     real(kind=rp), allocatable :: vr(:), vs(:), vt(:)
     type(c_ptr) :: temp_d = C_NULL_PTR
     type(c_ptr) :: tbf_d = C_NULL_PTR
     type(c_ptr) :: tx_d = C_NULL_PTR
     type(c_ptr) :: ty_d = C_NULL_PTR
     type(c_ptr) :: tz_d = C_NULL_PTR
     type(c_ptr) :: vr_d = C_NULL_PTR
     type(c_ptr) :: vs_d = C_NULL_PTR
     type(c_ptr) :: vt_d = C_NULL_PTR
     !> Next node in the pool.
     type(adv_dealias_scratch_t), pointer :: next => null()
  end type adv_dealias_scratch_t

  !> Head of the shared scratch pool.
  type(adv_dealias_scratch_t), pointer :: adv_dealias_pool => null()

  !> Type encapsulating advection routines with dealiasing
  type, public, extends(advection_t) :: adv_dealias_t
     !> The shared node backing every pointer below.
     type(adv_dealias_scratch_t), pointer :: scratch => null()
     !> Coeffs of the higher-order space
     type(coef_t), pointer :: coef_GL => null()
     !> Coeffs of the original space in the simulation
     type(coef_t), pointer :: coef_GLL => null()
     !> Interpolator between the original and higher-order spaces
     type(interpolator_t), pointer :: GLL_to_GL => null()
     !> The additional higher-order space used in dealiasing
     type(space_t), pointer :: Xh_GL => null()
     !> The original space used in the simulation
     type(space_t), pointer :: Xh_GLL => null()
     real(kind=rp), pointer, contiguous :: temp(:) => null(), tbf(:) => null()
     !> Temporary arrays
     real(kind=rp), pointer, contiguous :: tx(:) => null(), ty(:) => null(), &
          tz(:) => null()
     real(kind=rp), pointer, contiguous :: vr(:) => null(), vs(:) => null(), &
          vt(:) => null()
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

    call this%free()

    ! The GL space, its geometry and the work arrays are identical for every
    ! scheme built on this (lxd, coef), so take a reference to a shared node
    ! rather than building a private copy.
    this%scratch => adv_dealias_scratch_acquire(lxd, coef)

    this%Xh_GL => this%scratch%Xh_GL
    this%coef_GL => this%scratch%coef_GL
    this%GLL_to_GL => this%scratch%GLL_to_GL
    this%Xh_GLL => coef%Xh
    this%coef_GLL => coef

    ! The plain CPU path builds its temporaries per element on the stack and
    ! the node allocates nothing, so there is nothing to point at there.
    if (allocated(this%scratch%temp)) then
       this%temp => this%scratch%temp
       this%tbf => this%scratch%tbf
       this%tx => this%scratch%tx
       this%ty => this%scratch%ty
       this%tz => this%scratch%tz
       this%vr => this%scratch%vr
       this%vs => this%scratch%vs
       this%vt => this%scratch%vt

       this%temp_d = this%scratch%temp_d
       this%tbf_d = this%scratch%tbf_d
       this%tx_d = this%scratch%tx_d
       this%ty_d = this%scratch%ty_d
       this%tz_d = this%scratch%tz_d
       this%vr_d = this%scratch%vr_d
       this%vs_d = this%scratch%vs_d
       this%vt_d = this%scratch%vt_d
    end if

    ! No-op for every instance after the first on this node.
    call adv_dealias_scratch_refresh(this%scratch, coef)

  end subroutine init_dealias

  !> Destructor
  subroutine free_dealias(this)
    class(adv_dealias_t), intent(inout) :: this

    nullify(this%temp)
    nullify(this%tbf)
    nullify(this%tx)
    nullify(this%ty)
    nullify(this%tz)
    nullify(this%vr)
    nullify(this%vs)
    nullify(this%vt)

    this%temp_d = C_NULL_PTR
    this%tbf_d = C_NULL_PTR
    this%tx_d = C_NULL_PTR
    this%ty_d = C_NULL_PTR
    this%tz_d = C_NULL_PTR
    this%vr_d = C_NULL_PTR
    this%vs_d = C_NULL_PTR
    this%vt_d = C_NULL_PTR

    nullify(this%coef_GL)
    nullify(this%GLL_to_GL)
    nullify(this%Xh_GL)
    nullify(this%Xh_GLL)
    nullify(this%coef_GLL)

    ! Releases the shared node once the last instance lets go of it.
    call adv_dealias_scratch_release(this%scratch)

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
    call adv_dealias_scratch_claim(this%scratch)

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
         !$omp parallel do private(e, i, idx, tempx, tempy, tempz), &
         !$omp& private(tx, ty, tz, vr, vs, vt, tfx, tfy, tfz)
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
    call adv_dealias_scratch_yield(this%scratch)

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

    call adv_dealias_scratch_claim(this%scratch)

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

         ! Compute the convective term, i.e dot the velocity with the
         ! scalar grad
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

         ! Compute the convective term, i.e dot the velocity with the
         ! scalar grad
         call vdot3(this%tbf, this%vr, this%vs, this%vt, &
              this%tx, this%ty, this%tz, n_GL)

         ! Map back to the original space (we reuse this%temp)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)

         ! Update the source term
         call sub2(fs%x, this%temp, n)

      else
         !$omp parallel do private (e, i, idx, vx_GL, vy_GL, vz_GL, s_GL), &
         !$omp& private (f_GL, temp, dsdx, dsdy, dsdz)
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
    call adv_dealias_scratch_yield(this%scratch)

  end subroutine compute_scalar_advection_dealias


  !!> Add the advection term in ALE framework using dealiasing.
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
  !! Based on Ho, L.W. A Legendre spectral element method for
  !! simulation of incompressible unsteady viscous free-surface flows.
  !! Ph.D. thesis, Massachusetts Institute of Technology, 1989.
  !! Note: In Nek5000, dealiasing is not done for this term.
  subroutine compute_ale_advection_dealias(this, vx, vy, vz, &
       wm_x, wm_y, wm_z, fx, fy, fz, Xh, coef, n, dt)
    class(adv_dealias_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: wm_x, wm_y, wm_z
    type(field_t), intent(inout) :: fx, fy, fz
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: vx_GL, vy_GL, vz_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: wm_x_GL, wm_y_GL, wm_z_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: flux_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: grad_x, grad_y, grad_z
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: total_div_GL
    integer :: e, i, idx, nel, n_GL
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: temp_x, temp_y, temp_z

    nel = coef%msh%nelv
    n_GL = nel * this%Xh_GL%lxyz

    call adv_dealias_scratch_claim(this%scratch)

    associate(c_GL => this%coef_GL)
      if (NEKO_BCKND_DEVICE .eq. 1) then

         ! Map mesh velocity (wm) to the GL space
         call this%GLL_to_GL%map(this%vr, wm_x%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%vs, wm_y%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%vt, wm_z%x, nel, this%Xh_GL)

         ! --------------------- X-Momentum
         ! vx * wm_x
         call this%GLL_to_GL%map(this%temp, vx%x, nel, this%Xh_GL)
         call device_col3(this%temp_d, this%temp_d, this%vr_d, n_GL)
         call opgrad(this%tz, this%tx, this%ty, this%temp, c_GL)

         ! vx * wm_y
         call this%GLL_to_GL%map(this%temp, vx%x, nel, this%Xh_GL)
         call device_col3(this%temp_d, this%temp_d, this%vs_d, n_GL)
         call opgrad(this%tx, this%tbf, this%ty, this%temp, c_GL)
         call device_add2(this%tz_d, this%tbf_d, n_GL)

         ! vx * wm_z
         call this%GLL_to_GL%map(this%temp, vx%x, nel, this%Xh_GL)
         call device_col3(this%temp_d, this%temp_d, this%vt_d, n_GL)
         call opgrad(this%tx, this%ty, this%tbf, this%temp, c_GL)
         call device_add2(this%tz_d, this%tbf_d, n_GL)

         ! Map divergence back to GLL space and add to RHS
         call this%GLL_to_GL%map(this%temp, this%tz, nel, this%Xh_GLL)
         call device_add2(fx%x_d, this%temp_d, n)

         ! --------------------- Y-Momentum
         ! vy * wm_x
         call this%GLL_to_GL%map(this%temp, vy%x, nel, this%Xh_GL)
         call device_col3(this%temp_d, this%temp_d, this%vr_d, n_GL)
         call opgrad(this%tz, this%tx, this%ty, this%temp, c_GL)

         ! vy * wm_y
         call this%GLL_to_GL%map(this%temp, vy%x, nel, this%Xh_GL)
         call device_col3(this%temp_d, this%temp_d, this%vs_d, n_GL)
         call opgrad(this%tx, this%tbf, this%ty, this%temp, c_GL)
         call device_add2(this%tz_d, this%tbf_d, n_GL)

         ! vy * wm_z
         call this%GLL_to_GL%map(this%temp, vy%x, nel, this%Xh_GL)
         call device_col3(this%temp_d, this%temp_d, this%vt_d, n_GL)
         call opgrad(this%tx, this%ty, this%tbf, this%temp, c_GL)
         call device_add2(this%tz_d, this%tbf_d, n_GL)

         call this%GLL_to_GL%map(this%temp, this%tz, nel, this%Xh_GLL)
         call device_add2(fy%x_d, this%temp_d, n)

         ! --------------------- Z-Momentum
         ! wz * wm_x
         call this%GLL_to_GL%map(this%temp, vz%x, nel, this%Xh_GL)
         call device_col3(this%temp_d, this%temp_d, this%vr_d, n_GL)
         call opgrad(this%tz, this%tx, this%ty, this%temp, c_GL)

         ! wz * wm_y
         call this%GLL_to_GL%map(this%temp, vz%x, nel, this%Xh_GL)
         call device_col3(this%temp_d, this%temp_d, this%vs_d, n_GL)
         call opgrad(this%tx, this%tbf, this%ty, this%temp, c_GL)
         call device_add2(this%tz_d, this%tbf_d, n_GL)

         ! wz * wm_z
         call this%GLL_to_GL%map(this%temp, vz%x, nel, this%Xh_GL)
         call device_col3(this%temp_d, this%temp_d, this%vt_d, n_GL)
         call opgrad(this%tx, this%ty, this%tbf, this%temp, c_GL)
         call device_add2(this%tz_d, this%tbf_d, n_GL)

         call this%GLL_to_GL%map(this%temp, this%tz, nel, this%Xh_GLL)
         call device_add2(fz%x_d, this%temp_d, n)

      else
         !$omp parallel do private (e, i, idx, vx_GL, vy_GL, vz_GL), &
         !$omp& private (wm_x_GL, wm_y_GL, wm_z_GL, flux_GL, total_div_GL), &
         !$omp& private (grad_x, grad_y, grad_z, temp_x, temp_y, temp_z)
         do e = 1, coef%msh%nelv
            ! Map advecting velocity and mesh velocity onto the
            ! higher-order space
            call this%GLL_to_GL%map(vx_GL, vx%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(vy_GL, vy%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(vz_GL, vz%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(wm_x_GL, wm_x%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(wm_y_GL, wm_y%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(wm_z_GL, wm_z%x(1,1,1,e), 1, this%Xh_GL)

            ! I think below can be written more efficiently. Will fix it later.
            ! This works for now.
            ! --------------------- X-Momentum
            total_div_GL = 0.0_rp
            ! div(u * wm_*) = d/dx (u * wm_x) + d/dy (u * wm_y) +
            ! d/dz (u * wm_z)

            flux_GL = vx_GL * wm_x_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_x
            flux_GL = vx_GL * wm_y_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_y
            flux_GL = vx_GL * wm_z_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_z

            ! Map back the constructed operator to the original space
            call this%GLL_to_GL%map(temp_x, total_div_GL, 1, this%Xh_GLL)

            ! --------------------- Y-Momentum
            total_div_GL = 0.0_rp
            ! div(v * wm_*) = d/dx (v * wm_x) + d/dy (v * wm_y) +
            ! d/dz (v * wm_z)

            flux_GL = vy_GL * wm_x_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_x
            flux_GL = vy_GL * wm_y_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_y
            flux_GL = vy_GL * wm_z_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_z

            ! Map back the constructed operator to the original space
            call this%GLL_to_GL%map(temp_y, total_div_GL, 1, this%Xh_GLL)

            ! --------------------- Z-Momentum
            total_div_GL = 0.0_rp
            ! div(w * wm_*) = d/dx (w * wm_x) + d/dy (w * wm_y) +
            ! d/dz (w * wm_z)

            flux_GL = vz_GL * wm_x_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_x
            flux_GL = vz_GL * wm_y_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_y
            flux_GL = vz_GL * wm_z_GL
            call opgrad(grad_x, grad_y, grad_z, flux_GL, c_GL, e, e)
            total_div_GL = total_div_GL + grad_z

            ! Map back the constructed operator to the original space
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
    call adv_dealias_scratch_yield(this%scratch)
  end subroutine compute_ale_advection_dealias

  subroutine recompute_metrics_dealias(this, coef, moving_boundary)
    class(adv_dealias_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    logical, intent(in) :: moving_boundary

    if (.not. moving_boundary) return

    ! The node is keyed on this coef, so every scheme sharing it asks for the
    ! same refresh. The version guard inside does it once.
    call adv_dealias_scratch_refresh(this%scratch, coef)

  end subroutine recompute_metrics_dealias

  !> Find, or build, the shared node for a (lxd, coef) pair and take a
  !! reference to it.
  !! @param lxd The polynomial order of the dealiasing space.
  !! @param coef The fine-level coefficients the space is built over.
  function adv_dealias_scratch_acquire(lxd, coef) result(s)
    integer, intent(in) :: lxd
    type(coef_t), intent(inout), target :: coef
    type(adv_dealias_scratch_t), pointer :: s
    integer :: nel, n_GL

    s => adv_dealias_pool
    do while (associated(s))
       if (s%lxd .eq. lxd .and. associated(s%key_coef, coef)) then
          s%refcount = s%refcount + 1
          return
       end if
       s => s%next
    end do

    allocate(s)
    s%lxd = lxd
    s%key_coef => coef
    s%refcount = 1

    call s%Xh_GL%init(GL, lxd, lxd, lxd)
    call s%GLL_to_GL%init(s%Xh_GL, coef%Xh)
    call s%coef_GL%init(s%Xh_GL, coef%msh)

    nel = coef%msh%nelv
    n_GL = nel * s%Xh_GL%lxyz

    ! Only the backends that work on whole-field buffers need these; the plain
    ! CPU path builds its temporaries per element on the stack.
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1) .or. (NEKO_BCKND_METAL .eq. 1) .or. &
         (NEKO_BCKND_SX .eq. 1) .or. (NEKO_BCKND_XSMM .eq. 1)) then
       allocate(s%temp(n_GL))
       allocate(s%tbf(n_GL))
       allocate(s%tx(n_GL))
       allocate(s%ty(n_GL))
       allocate(s%tz(n_GL))
       allocate(s%vr(n_GL))
       allocate(s%vs(n_GL))
       allocate(s%vt(n_GL))

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(s%temp, s%temp_d, n_GL)
          call device_map(s%tbf, s%tbf_d, n_GL)
          call device_map(s%tx, s%tx_d, n_GL)
          call device_map(s%ty, s%ty_d, n_GL)
          call device_map(s%tz, s%tz_d, n_GL)
          call device_map(s%vr, s%vr_d, n_GL)
          call device_map(s%vs, s%vs_d, n_GL)
          call device_map(s%vt, s%vt_d, n_GL)
       end if
    end if

    s%next => adv_dealias_pool
    adv_dealias_pool => s

  end function adv_dealias_scratch_acquire

  !> Drop a reference to a shared node, freeing it once nothing points at it.
  !! @param s The node; nullified on return.
  subroutine adv_dealias_scratch_release(s)
    type(adv_dealias_scratch_t), pointer, intent(inout) :: s
    type(adv_dealias_scratch_t), pointer :: p

    if (.not. associated(s)) return

    s%refcount = s%refcount - 1
    if (s%refcount .gt. 0) then
       nullify(s)
       return
    end if

    if (associated(adv_dealias_pool, s)) then
       adv_dealias_pool => s%next
    else
       p => adv_dealias_pool
       do while (associated(p))
          if (associated(p%next, s)) then
             p%next => s%next
             exit
          end if
          p => p%next
       end do
    end if

    if (allocated(s%temp)) then
       if (NEKO_BCKND_DEVICE .eq. 1) call device_unmap(s%temp, s%temp_d)
       deallocate(s%temp)
    end if
    if (allocated(s%tbf)) then
       if (NEKO_BCKND_DEVICE .eq. 1) call device_unmap(s%tbf, s%tbf_d)
       deallocate(s%tbf)
    end if
    if (allocated(s%tx)) then
       if (NEKO_BCKND_DEVICE .eq. 1) call device_unmap(s%tx, s%tx_d)
       deallocate(s%tx)
    end if
    if (allocated(s%ty)) then
       if (NEKO_BCKND_DEVICE .eq. 1) call device_unmap(s%ty, s%ty_d)
       deallocate(s%ty)
    end if
    if (allocated(s%tz)) then
       if (NEKO_BCKND_DEVICE .eq. 1) call device_unmap(s%tz, s%tz_d)
       deallocate(s%tz)
    end if
    if (allocated(s%vr)) then
       if (NEKO_BCKND_DEVICE .eq. 1) call device_unmap(s%vr, s%vr_d)
       deallocate(s%vr)
    end if
    if (allocated(s%vs)) then
       if (NEKO_BCKND_DEVICE .eq. 1) call device_unmap(s%vs, s%vs_d)
       deallocate(s%vs)
    end if
    if (allocated(s%vt)) then
       if (NEKO_BCKND_DEVICE .eq. 1) call device_unmap(s%vt, s%vt_d)
       deallocate(s%vt)
    end if

    call s%coef_GL%free()
    call s%GLL_to_GL%free()
    call s%Xh_GL%free()
    nullify(s%key_coef)
    nullify(s%next)

    deallocate(s)
    nullify(s)

  end subroutine adv_dealias_scratch_release

  !> Map the fine-level derivative arrays onto the node's GL space.
  !!
  !! Runs once per `metrics_version` of the keying coefficients, so on a moving
  !! mesh the nine interpolations happen once however many schemes share the
  !! node.
  !! @param s The shared node.
  !! @param coef The fine-level coefficients the node was keyed on.
  subroutine adv_dealias_scratch_refresh(s, coef)
    type(adv_dealias_scratch_t), pointer, intent(inout) :: s
    type(coef_t), intent(in) :: coef
    integer :: nel

    if (.not. associated(s)) return
    if (s%metrics_version .eq. coef%metrics_version) return

    nel = coef%msh%nelv
    call s%GLL_to_GL%map(s%coef_GL%drdx, coef%drdx, nel, s%Xh_GL)
    call s%GLL_to_GL%map(s%coef_GL%dsdx, coef%dsdx, nel, s%Xh_GL)
    call s%GLL_to_GL%map(s%coef_GL%dtdx, coef%dtdx, nel, s%Xh_GL)

    call s%GLL_to_GL%map(s%coef_GL%drdy, coef%drdy, nel, s%Xh_GL)
    call s%GLL_to_GL%map(s%coef_GL%dsdy, coef%dsdy, nel, s%Xh_GL)
    call s%GLL_to_GL%map(s%coef_GL%dtdy, coef%dtdy, nel, s%Xh_GL)

    call s%GLL_to_GL%map(s%coef_GL%drdz, coef%drdz, nel, s%Xh_GL)
    call s%GLL_to_GL%map(s%coef_GL%dsdz, coef%dsdz, nel, s%Xh_GL)
    call s%GLL_to_GL%map(s%coef_GL%dtdz, coef%dtdz, nel, s%Xh_GL)

    s%metrics_version = coef%metrics_version

  end subroutine adv_dealias_scratch_refresh

  !> Take the shared work arrays for the duration of one compute() call.
  !!
  !! The arrays are shared across every scheme on the node, which is sound only
  !! because the schemes are advanced one at a time. Claiming makes a violation
  !! of that an error rather than a silent data race.
  !! @param s The shared node.
  subroutine adv_dealias_scratch_claim(s)
    type(adv_dealias_scratch_t), pointer, intent(inout) :: s

    if (.not. associated(s)) return
    if (.not. allocated(s%temp)) return

    if (s%busy) then
       call neko_error('adv_dealias: the shared dealiasing work arrays are ' // &
            'already in use. Two advection operators cannot compute at the ' // &
            'same time.')
    end if
    s%busy = .true.

  end subroutine adv_dealias_scratch_claim

  !> Give the shared work arrays back.
  !! @param s The shared node.
  subroutine adv_dealias_scratch_yield(s)
    type(adv_dealias_scratch_t), pointer, intent(inout) :: s

    if (.not. associated(s)) return
    s%busy = .false.

  end subroutine adv_dealias_scratch_yield

end module adv_dealias
