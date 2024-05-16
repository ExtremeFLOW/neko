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
module adv_lin_dealias
  use advection, only: advection_lin_t
  use num_types, only : rp
  use math, only : vdot3, sub2, subcol3, rzero
  use utils, only : neko_error
  use space, only : space_t, GL
  use field, only : field_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX, NEKO_BCKND_XSMM, &
    NEKO_BCKND_OPENCL, NEKO_BCKND_CUDA, NEKO_BCKND_HIP
  use operators, only : opgrad, conv1, cdtp
  use interpolation, only : interpolator_t
  use device_math
  use device, only : device_free, device_map, device_get_ptr, device_memcpy, &
    HOST_TO_DEVICE
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, &
    c_associated
  implicit none
  private

  !> Type encapsulating advection routines with dealiasing
  type, public, extends(advection_lin_t) :: adv_lin_dealias_t
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

     ! Tim,
     ! -----------------------------------------------
     ! I wanted one more set of temp arrays on GLL
     ! maybe there's a smarter way to construct the operator
     ! and reuse the previous temp arrays, but for now
     ! this will work
     ! Harry
     real(kind=rp), allocatable :: txb(:), tyb(:), tzb(:)
     !> Device pointer for `tx`
     type(c_ptr) :: txb_d = C_NULL_PTR
     !> Device pointer for `ty`
     type(c_ptr) :: tyb_d = C_NULL_PTR
     !> Device pointer for `tz`
     type(c_ptr) :: tzb_d = C_NULL_PTR

     ! revised
     ! THIS IS THE WORST WAY OF DOING IT
     ! but I want to check the operator
     ! I'll use less memory in the future
     real(kind=rp), allocatable :: duxb(:), duyb(:), duzb(:)
     real(kind=rp), allocatable :: dvxb(:), dvyb(:), dvzb(:)
     real(kind=rp), allocatable :: dwxb(:), dwyb(:), dwzb(:)
     type(c_ptr) :: duxb_d = C_NULL_PTR
     type(c_ptr) :: duyb_d = C_NULL_PTR
     type(c_ptr) :: duzb_d = C_NULL_PTR
     type(c_ptr) :: dvxb_d = C_NULL_PTR
     type(c_ptr) :: dvyb_d = C_NULL_PTR
     type(c_ptr) :: dvzb_d = C_NULL_PTR
     type(c_ptr) :: dwxb_d = C_NULL_PTR
     type(c_ptr) :: dwyb_d = C_NULL_PTR
     type(c_ptr) :: dwzb_d = C_NULL_PTR

   contains
     !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
     !! the RHS.
     procedure, pass(this) :: compute_linear => compute_linear_advection_dealias
     !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
     !! the RHS.
     procedure, pass(this) :: compute_adjoint => compute_adjoint_advection_dealias
     !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to
     !! the RHS.
     procedure, pass(this) :: compute_scalar => compute_scalar_advection_dealias
     !> Constructor
     procedure, pass(this) :: init => init_dealias
     !> Destructor
     procedure, pass(this) :: free => free_dealias
  end type adv_lin_dealias_t

contains

  !> Constructor
  !! @param lxd The polynomial order of the space used in the dealiasing.
  !! @param coef The coefficients of the (space, mesh) pair.
  subroutine init_dealias(this, lxd, coef)
    class(adv_lin_dealias_t), target, intent(inout) :: this
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
    ! HARRY
    ! I also want the det(jacobian) and the mass matrix
    ! but this shouldn't be needed if we modify cdtp to use w3 for the weights
    call this%GLL_to_GL%map(this%coef_GL%jac, coef%jac, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%B, coef%B, nel, this%Xh_GL)

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
       ! HARRY DO SOMETHING SMARTER!!!!
       allocate(this%duxb(n_GL))
       allocate(this%duyb(n_GL))
       allocate(this%duzb(n_GL))
       allocate(this%dvxb(n_GL))
       allocate(this%dvyb(n_GL))
       allocate(this%dvzb(n_GL))
       allocate(this%dwxb(n_GL))
       allocate(this%dwyb(n_GL))
       allocate(this%dwzb(n_GL))
       allocate(this%txb(n_GL))
       allocate(this%tyb(n_GL))
       allocate(this%tzb(n_GL))
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
       ! HARRY DO SOMETHING SMARTER!!!!
       call device_map(this%duxb, this%duxb_d, n_GL)
       call device_map(this%duyb, this%duyb_d, n_GL)
       call device_map(this%duzb, this%duzb_d, n_GL)
       call device_map(this%dvxb, this%dvxb_d, n_GL)
       call device_map(this%dvyb, this%dvyb_d, n_GL)
       call device_map(this%dvzb, this%dvzb_d, n_GL)
       call device_map(this%dwxb, this%dwxb_d, n_GL)
       call device_map(this%dwyb, this%dwyb_d, n_GL)
       call device_map(this%dwzb, this%dwzb_d, n_GL)
       call device_map(this%txb, this%txb_d, n_GL)
       call device_map(this%tyb, this%tyb_d, n_GL)
       call device_map(this%tzb, this%tzb_d, n_GL)
    end if

  end subroutine init_dealias

  !> Destructor
  subroutine free_dealias(this)
    class(adv_lin_dealias_t), intent(inout) :: this
  end subroutine free_dealias


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
  subroutine compute_scalar_advection_dealias(this, vx, vy, vz, s, fs, Xh, &
                                              coef, n)
    class(adv_lin_dealias_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: s
    type(field_t), intent(inout) :: fs
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: n

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
      end if
    end associate

  end subroutine compute_scalar_advection_dealias

  !> Add the linearized advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
  !! the RHS.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param vxb The x component of baseflow velocity.
  !! @param vyb The y component of baseflow velocity.
  !! @param vzb The z component of baseflow velocity.
  !! @param fx The x component of source term.
  !! @param fy The y component of source term.
  !! @param fz The z component of source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  subroutine compute_adjoint_advection_dealias(this, vx, vy, vz, vxb, vyb, vzb, fx, fy, fz, Xh, coef, n)
    !! HARRY added vxb etc for baseflow
    implicit none
    class(adv_lin_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    ! HARRY ---
    ! Baseflow
    type(field_t), intent(inout) :: vxb, vyb, vzb
    ! --------
    type(field_t), intent(inout) :: fx, fy, fz
    integer, intent(in) :: n
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tfx, tfy, tfz

!    ! these are u and U_b on dealiased mesh (one element)
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tx, ty, tz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: txb, tyb, tzb

    ! these are gradients of U_b on dealiased mesh (one element)
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: duxb, dvxb, dwxb
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: duyb, dvyb, dwyb
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: duzb, dvzb, dwzb

    real(kind=rp), dimension(this%Xh_GL%lxyz) :: vr, vs, vt
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: tempx, tempy, tempz
    integer :: e, i, idx, nel, n_GL
    real :: fac
    nel = coef%msh%nelv
    n_GL = nel * this%Xh_GL%lxyz
    associate(c_GL => this%coef_GL)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         ! HARRY
         ! These are the baseflow
         call this%GLL_to_GL%map(this%txb, vxb%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tyb, vyb%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tzb, vzb%x, nel, this%Xh_GL)

         ! These are velocity
         call this%GLL_to_GL%map(this%tx, vx%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%ty, vy%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tz, vz%x, nel, this%Xh_GL)


         ! u . grad U_b^T
         !-----------------------------
         ! take all the gradients
         call opgrad(this%duxb, this%duyb, this%duzb, this%txb, c_GL)
         call opgrad(this%dvxb, this%dvyb, this%dvzb, this%tyb, c_GL)
         call opgrad(this%dwxb, this%dwyb, this%dwzb, this%tzb, c_GL)

         ! traspose and multiply
         call device_vdot3(this%vr_d, this%tx_d, this%ty_d, this%tz_d, &
                           this%duxb_d, this%dvxb_d, this%dwxb_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%vr, nel, this%Xh_GLL)
         call device_sub2(fx%x_d, this%temp_d, n)

         call device_vdot3(this%vr_d, this%tx_d, this%ty_d, this%tz_d, &
                           this%duyb_d, this%dvyb_d, this%dwyb_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%vr, nel, this%Xh_GLL)
         call device_sub2(fy%x_d, this%temp_d, n)

         call device_vdot3(this%vr_d, this%tx_d, this%ty_d, this%tz_d, &
                           this%duzb_d, this%dvzb_d, this%dwzb_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%vr, nel, this%Xh_GLL)
         call device_sub2(fz%x_d, this%temp_d, n)

         ! \int \grad v . U_b ^ u    with ^ an outer product

         ! (x)
         ! duxb,duyb,duzb are temporary arrays
         call device_col3(this%duxb_d, this%tx_d, this%txb_d, n_GL)
         call device_col3(this%duyb_d, this%tx_d, this%tyb_d, n_GL)
         call device_col3(this%duzb_d, this%tx_d, this%tzb_d, n_GL)

         ! D^T
         ! vr,vs,vt are temporary arrays
         call cdtp(this%vr, this%duxb, c_GL%drdx, c_GL%dsdx, c_GL%dtdx, c_GL)
         call cdtp(this%vs, this%duyb, c_GL%drdy, c_GL%dsdy, c_GL%dtdy, c_GL)
         call cdtp(this%vt, this%duzb, c_GL%drdz, c_GL%dsdz, c_GL%dtdz, c_GL)

         ! reuse duxb as a temp
         call device_add4(this%duxb_d, this%vr_d, this%vs_d, this%vt_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%duxb, nel, this%Xh_GLL)
         call device_sub2(fx%x_d, this%temp_d, n)


         ! (y)
         ! duxb,duyb,duzb are temporary arrays
         call device_col3(this%duxb_d, this%ty_d, this%txb_d, n_GL)
         call device_col3(this%duyb_d, this%ty_d, this%tyb_d, n_GL)
         call device_col3(this%duzb_d, this%ty_d, this%tzb_d, n_GL)

         ! D^T
         ! vr,vs,vt are temporary arrays
         call cdtp(this%vr, this%duxb, c_GL%drdx, c_GL%dsdx, c_GL%dtdx, c_GL)
         call cdtp(this%vs, this%duyb, c_GL%drdy, c_GL%dsdy, c_GL%dtdy, c_GL)
         call cdtp(this%vt, this%duzb, c_GL%drdz, c_GL%dsdz, c_GL%dtdz, c_GL)

         ! reuse duxb as a temp
         call device_add4(this%duxb_d, this%vr_d, this%vs_d, this%vt_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%duxb, nel, this%Xh_GLL)
         call device_sub2(fy%x_d, this%temp_d, n)

         ! (z)
         ! duxb,duyb,duzb are temporary arrays
         call device_col3(this%duxb_d, this%tz_d, this%txb_d, n_GL)
         call device_col3(this%duyb_d, this%tz_d, this%tyb_d, n_GL)
         call device_col3(this%duzb_d, this%tz_d, this%tzb_d, n_GL)

         ! D^T
         ! vr,vs,vt are temporary arrays
         call cdtp(this%vr, this%duxb, c_GL%drdx, c_GL%dsdx, c_GL%dtdx, c_GL)
         call cdtp(this%vs, this%duyb, c_GL%drdy, c_GL%dsdy, c_GL%dtdy, c_GL)
         call cdtp(this%vt, this%duzb, c_GL%drdz, c_GL%dsdz, c_GL%dtdz, c_GL)

         ! reuse duxb as a temp
         call device_add4(this%duxb_d, this%vr_d, this%vs_d, this%vt_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%duxb, nel, this%Xh_GLL)
         call device_sub2(fz%x_d, this%temp_d, n)
      else if ((NEKO_BCKND_SX .eq. 1) .or. (NEKO_BCKND_XSMM .eq. 1)) then
         !TODO
      else



         do e = 1, coef%msh%nelv
            ! These are the baseflow
            call this%GLL_to_GL%map(txb, vxb%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(tyb, vyb%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(tzb, vzb%x(1,1,1,e), 1, this%Xh_GL)

            ! These are velocity
            call this%GLL_to_GL%map(tx, vx%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(ty, vy%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(tz, vz%x(1,1,1,e), 1, this%Xh_GL)


            ! u . grad U_b^T
            !-----------------------------
            call opgrad(duxb, duyb, duzb, txb, c_GL, e, e)
            call opgrad(dvxb, dvyb, dvzb, tyb, c_GL, e, e)
            call opgrad(dwxb, dwyb, dwzb, tzb, c_GL, e, e)

            ! traspose and multiply
            do i = 1, this%Xh_GL%lxyz
               tfx(i) = tx(i)*duxb(i) + ty(i)*dvxb(i) + tz(i)*dwxb(i)
               tfy(i) = tx(i)*duyb(i) + ty(i)*dvyb(i) + tz(i)*dwyb(i)
               tfz(i) = tx(i)*duzb(i) + ty(i)*dvzb(i) + tz(i)*dwzb(i)
            end do

            ! map back to GLL
            call this%GLL_to_GL%map(tempx, tfx, 1, this%Xh_GLL)
            call this%GLL_to_GL%map(tempy, tfy, 1, this%Xh_GLL)
            call this%GLL_to_GL%map(tempz, tfz, 1, this%Xh_GLL)

            ! accumulate
            idx = (e-1)*this%Xh_GLL%lxyz+1
            call sub2(fx%x(idx, 1, 1, 1), tempx, this%Xh_GLL%lxyz)
            call sub2(fy%x(idx, 1, 1, 1), tempy, this%Xh_GLL%lxyz)
            call sub2(fz%x(idx, 1, 1, 1), tempz, this%Xh_GLL%lxyz)

            ! In most literature this term is -(grad . U_b) u + BDRY
            !
            ! We prefer the weak form
            ! \int \grad v . U_b ^ u    with ^ an outer product
            ! avoiding:
            ! - an extra boundary term
            ! - the assumption of \grad . U_b pointwise
            !
            ! we're going to reuse uxb, uyb and uzb for the outer product
            !-----------------------------
            ! here we have index i, with free index j------------------------------x
            do i = 1, this%Xh_GL%lxyz
               !fac = this%Xh_GL%w3(i,1,1)*c_GL%jac(i,1,1,e)/c_GL%B(i,1,1,e)
               ! now this is correct
               fac = 1.0
               duxb(i) = tx(i)*txb(i)*fac
               duyb(i) = tx(i)*tyb(i)*fac
               duzb(i) = tx(i)*tzb(i)*fac
            end do

            ! D^T
            call cdtp(tfx, duxb, c_GL%drdx, c_GL%dsdx, c_GL%dtdx, c_GL, e, e)
            call cdtp(tfy, duyb, c_GL%drdy, c_GL%dsdy, c_GL%dtdy, c_GL, e, e)
            call cdtp(tfz, duzb, c_GL%drdz, c_GL%dsdz, c_GL%dtdz, c_GL, e, e)

            ! sum them
            do i = 1, this%Xh_GL%lxyz
               tfx(i) = tfx(i) + tfy(i) + tfz(i)
            end do

            ! map back to GLL
            call this%GLL_to_GL%map(tempx, tfx, 1, this%Xh_GLL)
            call sub2(fx%x(idx, 1, 1, 1), tempx, this%Xh_GLL%lxyz)

            ! here we have index i, with free index j------------------------------y
            do i = 1, this%Xh_GL%lxyz
               !fac = this%Xh_GL%w3(i,1,1)*c_GL%jac(i,1,1,e)/c_GL%B(i,1,1,e)
               ! now this is correct
               fac = 1.0
               duxb(i) = ty(i)*txb(i)*fac
               duyb(i) = ty(i)*tyb(i)*fac
               duzb(i) = ty(i)*tzb(i)*fac
            end do

            ! D^T
            call cdtp(tfx, duxb, c_GL%drdx, c_GL%dsdx, c_GL%dtdx, c_GL, e, e)
            call cdtp(tfy, duyb, c_GL%drdy, c_GL%dsdy, c_GL%dtdy, c_GL, e, e)
            call cdtp(tfz, duzb, c_GL%drdz, c_GL%dsdz, c_GL%dtdz, c_GL, e, e)

            ! sum them
            do i = 1, this%Xh_GL%lxyz
               tfx(i) = tfx(i) + tfy(i) + tfz(i)
            end do

            ! map back to GLL
            call this%GLL_to_GL%map(tempx, tfx, 1, this%Xh_GLL)
            call sub2(fy%x(idx, 1, 1, 1), tempx, this%Xh_GLL%lxyz)

            ! here we have index i, with free index j------------------------------z
            do i = 1, this%Xh_GL%lxyz
               !fac = this%Xh_GL%w3(i,1,1)*c_GL%jac(i,1,1,e)/c_GL%B(i,1,1,e)
               ! now this is correct
               fac = 1.0
               duxb(i) = tz(i)*txb(i)*fac
               duyb(i) = tz(i)*tyb(i)*fac
               duzb(i) = tz(i)*tzb(i)*fac
            end do
            ! D^T
            call cdtp(tfx, duxb, c_GL%drdx, c_GL%dsdx, c_GL%dtdx, c_GL, e, e)
            call cdtp(tfy, duyb, c_GL%drdy, c_GL%dsdy, c_GL%dtdy, c_GL, e, e)
            call cdtp(tfz, duzb, c_GL%drdz, c_GL%dsdz, c_GL%dtdz, c_GL, e, e)

            ! sum them
            do i = 1, this%Xh_GL%lxyz
               tfx(i) = tfx(i) + tfy(i) + tfz(i)
            end do

            ! map back to GLL
            call this%GLL_to_GL%map(tempx, tfx, 1, this%Xh_GLL)
            call sub2(fz%x(idx, 1, 1, 1), tempx, this%Xh_GLL%lxyz)

         enddo


      end if
    end associate

  end subroutine compute_adjoint_advection_dealias

  !> Add the linearized advection term for the fluid, i.e.
  !! \f$u \cdot \nabla u \f$, to the RHS.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param vxb The x component of baseflow velocity.
  !! @param vyb The y component of baseflow velocity.
  !! @param vzb The z component of baseflow velocity.
  !! @param fx The x component of source term.
  !! @param fy The y component of source term.
  !! @param fz The z component of source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  subroutine compute_linear_advection_dealias(this, vx, vy, vz, vxb, vyb, vzb, &
                                              fx, fy, fz, Xh, coef, n)
    !! HARRY added vxb etc for baseflow
    implicit none
    class(adv_lin_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    ! HARRY ---
    ! Baseflow
    type(field_t), intent(inout) :: vxb, vyb, vzb
    ! --------
    integer, intent(in) :: n
    type(field_t), intent(inout) :: fx, fy, fz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tx, ty, tz
    ! HARRY ---
    ! Baseflow temporary arrays
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: txb, tyb, tzb
    ! --------
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tfx, tfy, tfz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: vr, vs, vt
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: tempx, tempy, tempz

    integer :: e, i, idx, nel, n_GL
    nel = coef%msh%nelv
    n_GL = nel * this%Xh_GL%lxyz

    !This is extremely primitive and unoptimized  on the device //Karp
    associate(c_GL => this%coef_GL)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         ! HARRY
         ! These are the baseflow
         call this%GLL_to_GL%map(this%txb, vxb%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tyb, vyb%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tzb, vzb%x, nel, this%Xh_GL)

         ! These are velocity
         call this%GLL_to_GL%map(this%tx, vx%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%ty, vy%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tz, vz%x, nel, this%Xh_GL)

         ! OK this will be du.grad Ub
         call opgrad(this%vr, this%vs, this%vt, this%txb, c_GL)
         call device_vdot3(this%tbf_d, this%vr_d, this%vs_d, this%vt_d, &
                           this%tx_d, this%ty_d, this%tz_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call device_sub2(fx%x_d, this%temp_d, n)


         call opgrad(this%vr, this%vs, this%vt, this%tyb, c_GL)
         call device_vdot3(this%tbf_d, this%vr_d, this%vs_d, this%vt_d, &
                           this%tx_d, this%ty_d, this%tz_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call device_sub2(fy%x_d, this%temp_d, n)

         call opgrad(this%vr, this%vs, this%vt, this%tzb, c_GL)
         call device_vdot3(this%tbf_d, this%vr_d, this%vs_d, this%vt_d, &
                           this%tx_d, this%ty_d, this%tz_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call device_sub2(fz%x_d, this%temp_d, n)

         ! HARRY
         ! this will be U.grad du
         call opgrad(this%vr, this%vs, this%vt, this%tx, c_GL)
         call device_vdot3(this%tbf_d, this%vr_d, this%vs_d, this%vt_d, &
                           this%txb_d, this%tyb_d, this%tzb_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call device_sub2(fx%x_d, this%temp_d, n)


         call opgrad(this%vr, this%vs, this%vt, this%ty, c_GL)
         call device_vdot3(this%tbf_d, this%vr_d, this%vs_d, this%vt_d, &
                           this%txb_d, this%tyb_d, this%tzb_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call device_sub2(fy%x_d, this%temp_d, n)

         call opgrad(this%vr, this%vs, this%vt, this%tz, c_GL)
         call device_vdot3(this%tbf_d, this%vr_d, this%vs_d, this%vt_d, &
                           this%txb_d, this%tyb_d, this%tzb_d, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call device_sub2(fz%x_d, this%temp_d, n)

      else if ((NEKO_BCKND_SX .eq. 1) .or. (NEKO_BCKND_XSMM .eq. 1)) then
         ! HARRY
         ! These are the baseflow
         call this%GLL_to_GL%map(this%txb, vxb%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tyb, vyb%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tzb, vzb%x, nel, this%Xh_GL)

         ! These are velocity
         call this%GLL_to_GL%map(this%tx, vx%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%ty, vy%x, nel, this%Xh_GL)
         call this%GLL_to_GL%map(this%tz, vz%x, nel, this%Xh_GL)

         ! OK this will be du.grad Ub
         call opgrad(this%vr, this%vs, this%vt, this%txb, c_GL)
         call vdot3(this%tbf, this%vr, this%vs, this%vt, &
                    this%tx, this%ty, this%tz, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call sub2(fx%x, this%temp, n)


         call opgrad(this%vr, this%vs, this%vt, this%tyb, c_GL)
         call vdot3(this%tbf, this%vr, this%vs, this%vt, &
                    this%tx, this%ty, this%tz, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call sub2(fy%x, this%temp, n)

         call opgrad(this%vr, this%vs, this%vt, this%tzb, c_GL)
         call vdot3(this%tbf, this%vr, this%vs, this%vt, &
                    this%tx, this%ty, this%tz, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call sub2(fz%x, this%temp, n)

         ! OK this will be Ub.grad du
         call opgrad(this%vr, this%vs, this%vt, this%tx, c_GL)
         call vdot3(this%tbf, this%vr, this%vs, this%vt, &
                    this%txb, this%tyb, this%tzb, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call sub2(fx%x, this%temp, n)


         call opgrad(this%vr, this%vs, this%vt, this%ty, c_GL)
         call vdot3(this%tbf, this%vr, this%vs, this%vt, &
                    this%txb, this%tyb, this%tzb, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call sub2(fy%x, this%temp, n)

         call opgrad(this%vr, this%vs, this%vt, this%tz, c_GL)
         call vdot3(this%tbf, this%vr, this%vs, this%vt, &
                    this%txb, this%tyb, this%tzb, n_GL)
         call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
         call sub2(fz%x, this%temp, n)
      else

         do e = 1, coef%msh%nelv
            ! HARRY
            ! These are the baseflow
            call this%GLL_to_GL%map(txb, vxb%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(tyb, vyb%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(tzb, vzb%x(1,1,1,e), 1, this%Xh_GL)
            ! These are velocity
            call this%GLL_to_GL%map(tx, vx%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(ty, vy%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(tz, vz%x(1,1,1,e), 1, this%Xh_GL)

            ! OK this will be du.grad Ub
            call opgrad(vr, vs, vt, txb, c_GL, e, e)
            do i = 1, this%Xh_GL%lxyz
               tfx(i) = tx(i)*vr(i) + ty(i)*vs(i) + tz(i)*vt(i)
            end do

            call opgrad(vr, vs, vt, tyb, c_GL, e, e)
            do i = 1, this%Xh_GL%lxyz
               tfy(i) = tx(i)*vr(i) + ty(i)*vs(i) + tz(i)*vt(i)
            end do

            call opgrad(vr, vs, vt, tzb, c_GL, e, e)
            do i = 1, this%Xh_GL%lxyz
               tfz(i) = tx(i)*vr(i) + ty(i)*vs(i) + tz(i)*vt(i)
            end do

            call this%GLL_to_GL%map(tempx, tfx, 1, this%Xh_GLL)
            call this%GLL_to_GL%map(tempy, tfy, 1, this%Xh_GLL)
            call this%GLL_to_GL%map(tempz, tfz, 1, this%Xh_GLL)

            idx = (e-1)*this%Xh_GLL%lxyz+1
            call sub2(fx%x(idx, 1, 1, 1), tempx, this%Xh_GLL%lxyz)
            call sub2(fy%x(idx, 1, 1, 1), tempy, this%Xh_GLL%lxyz)
            call sub2(fz%x(idx, 1, 1, 1), tempz, this%Xh_GLL%lxyz)

            ! OK this will be Ub.grad du
            call opgrad(vr, vs, vt, tx, c_GL, e, e)
            do i = 1, this%Xh_GL%lxyz
               tfx(i) = txb(i)*vr(i) + tyb(i)*vs(i) + tzb(i)*vt(i)
            end do

            call opgrad(vr, vs, vt, ty, c_GL, e, e)
            do i = 1, this%Xh_GL%lxyz
               tfy(i) = txb(i)*vr(i) + tyb(i)*vs(i) + tzb(i)*vt(i)
            end do

            call opgrad(vr, vs, vt, tz, c_GL, e, e)
            do i = 1, this%Xh_GL%lxyz
               tfz(i) = txb(i)*vr(i) + tyb(i)*vs(i) + tzb(i)*vt(i)
            end do

            call this%GLL_to_GL%map(tempx, tfx, 1, this%Xh_GLL)
            call this%GLL_to_GL%map(tempy, tfy, 1, this%Xh_GLL)
            call this%GLL_to_GL%map(tempz, tfz, 1, this%Xh_GLL)

            idx = (e-1)*this%Xh_GLL%lxyz+1
            call sub2(fx%x(idx, 1, 1, 1), tempx, this%Xh_GLL%lxyz)
            call sub2(fy%x(idx, 1, 1, 1), tempy, this%Xh_GLL%lxyz)
            call sub2(fz%x(idx, 1, 1, 1), tempz, this%Xh_GLL%lxyz)
         end do
      end if
    end associate

  end subroutine compute_linear_advection_dealias

end module adv_lin_dealias
