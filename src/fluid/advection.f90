! Copyright (c) 2021-2022, The Neko Authors
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
!> Subroutines to apply advection to RHS
module advection
  use num_types
  use math
  use utils
  use space
  use field
  use coefs
  use neko_config
  use operators
  use interpolation
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  type, public, abstract :: advection_t
   contains
     procedure(apply_adv), pass(this), deferred :: apply
  end type advection_t

  type, public, extends(advection_t) :: adv_no_dealias_t
     real(kind=rp), allocatable :: temp(:)
     type(c_ptr) :: temp_d = C_NULL_PTR
   contains
     procedure, pass(this) :: apply => advab
  end type adv_no_dealias_t

  type, public, extends(advection_t) :: adv_dealias_t
     type(coef_t) :: coef_GL
     type(coef_t), pointer :: coef_GLL
     type(interpolator_t) :: GLL_to_GL
     type(space_t) :: Xh_GL
     type(space_t), pointer :: Xh_GLL
     real(kind=rp), allocatable :: temp(:), tbf(:)
     type(c_ptr) :: temp_d = C_NULL_PTR
     type(c_ptr) :: tbf_d = C_NULL_PTR
     real(kind=rp), allocatable :: tx(:), ty(:), tz(:)
     real(kind=rp), allocatable :: vr(:), vs(:), vt(:)
     type(c_ptr) :: tx_d = C_NULL_PTR
     type(c_ptr) :: ty_d = C_NULL_PTR
     type(c_ptr) :: tz_d = C_NULL_PTR
     type(c_ptr) :: vr_d = C_NULL_PTR
     type(c_ptr) :: vs_d = C_NULL_PTR
     type(c_ptr) :: vt_d = C_NULL_PTR


   contains
     procedure, pass(this) :: apply => apply_adv_dealias
     procedure, pass(this) :: init => init_dealias
  end type adv_dealias_t

  abstract interface
     subroutine apply_adv(this, vx, vy, vz, bfx, bfy, bfz, Xh, coef, n)
       import :: advection_t
       import :: coef_t
       import :: space_t
       import :: field_t
       import :: rp
       class(advection_t), intent(inout) :: this
       type(space_t), intent(inout) :: Xh
       type(coef_t), intent(inout) :: coef
       type(field_t), intent(inout) :: vx, vy, vz
       integer, intent(inout) :: n
       real(kind=rp), intent(inout), dimension(n) :: bfx, bfy, bfz
     end subroutine apply_adv
  end interface

  public :: advection_factory

contains
  
  subroutine advection_factory(this, coef, dealias, lxd)
    class(advection_t), allocatable, intent(inout) :: this
    type(coef_t), target :: coef
    logical, intent(in) :: dealias
    integer, intent(in) :: lxd

    if (allocated(this)) then
       select type(adv => this)
       type is(adv_no_dealias_t)
          if (allocated(adv%temp)) then
             deallocate(adv%temp)
          end if
          if (c_associated(adv%temp_d)) then
             call device_free(adv%temp_d)
          end if
       end select
       deallocate(this)
    end if

    if (dealias) then
       allocate(adv_dealias_t::this)
    else
       allocate(adv_no_dealias_t::this)
    end if
    
    select type(adv => this)
    type is(adv_dealias_t)
       if (lxd .gt. 0) then
          call init_dealias(adv, lxd, coef) 
       else
          call init_dealias(adv, coef%Xh%lx * 3/2,  coef)
       end if
    type is(adv_no_dealias_t)
       call init_no_dealias(adv, coef)
    end select

  end subroutine advection_factory

  subroutine init_no_dealias(this, coef)
    class(adv_no_dealias_t) :: this
    type(coef_t) :: coef

    allocate(this%temp(coef%dof%n_dofs))

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_map(this%temp, this%temp_d, coef%dof%n_dofs)
    end if

  end subroutine init_no_dealias

  subroutine init_dealias(this, lxd, coef)
    class(adv_dealias_t), target, intent(inout) :: this
    integer, intent(in) :: lxd
    type(coef_t), intent(inout), target :: coef
    integer :: nel, n_GL, n

    call space_init(this%Xh_GL, GL, lxd, lxd, lxd)
    this%Xh_GLL => coef%Xh
    this%coef_GLL => coef
    call this%GLL_to_GL%init(this%Xh_GL, this%Xh_GLL)

    call coef_init(this%coef_GL, this%Xh_GL, coef%msh)

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
    
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then       
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
  
  !> Eulerian scheme, add convection term to forcing function
  !! at current time step.
  subroutine apply_adv_dealias(this, vx, vy, vz, bfx, bfy, bfz, Xh, coef, n)
    class(adv_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    integer, intent(inout) :: n
    real(kind=rp), intent(inout), dimension(n) :: bfx, bfy, bfz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tx, ty, tz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tbfx, tbfy, tbfz 
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: vr, vs, vt
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: tempx, tempy, tempz
    type(c_ptr) :: bfx_d, bfy_d, bfz_d
    integer :: e, i, idx, nel, n_GL
    nel = coef%msh%nelv
    n_GL = nel * this%Xh_GL%lxyz
    !This is extremely primitive and unoptimized  on the device //Karp
    associate(c_GL => this%coef_GL)
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       bfx_d = device_get_ptr(bfx)
       bfy_d = device_get_ptr(bfy)
       bfz_d = device_get_ptr(bfz)
       call this%GLL_to_GL%map(this%tx, vx%x, nel, this%Xh_GL)
       call this%GLL_to_GL%map(this%ty, vy%x, nel, this%Xh_GL)
       call this%GLL_to_GL%map(this%tz, vz%x, nel, this%Xh_GL)

       call opgrad(this%vr, this%vs, this%vt, this%tx, c_GL)
       call device_col3(this%tbf_d, this%vr_d, this%tx_d, n_GL)
       call device_addcol3(this%tbf_d, this%vs_d, this%ty_d, n_GL)
       call device_addcol3(this%tbf_d, this%vt_d, this%tz_d, n_GL)
       call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
       call device_sub2(bfx_d, this%temp_d, n)


       call opgrad(this%vr, this%vs, this%vt, this%ty, c_GL)
       call device_col3(this%tbf_d, this%vr_d, this%tx_d, n_GL)
       call device_addcol3(this%tbf_d, this%vs_d, this%ty_d, n_GL)
       call device_addcol3(this%tbf_d, this%vt_d, this%tz_d, n_GL)
       call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
       call device_sub2(bfy_d, this%temp_d, n)

       call opgrad(this%vr, this%vs, this%vt, this%tz, c_GL)
       call device_col3(this%tbf_d, this%vr_d, this%tx_d, n_GL)
       call device_addcol3(this%tbf_d, this%vs_d, this%ty_d, n_GL)
       call device_addcol3(this%tbf_d, this%vt_d, this%tz_d, n_GL)
       call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
       call device_sub2(bfz_d, this%temp_d, n)

    else if ((NEKO_BCKND_SX .eq. 1) .or. (NEKO_BCKND_XSMM .eq. 1)) then

       call this%GLL_to_GL%map(this%tx, vx%x, nel, this%Xh_GL)
       call this%GLL_to_GL%map(this%ty, vy%x, nel, this%Xh_GL)
       call this%GLL_to_GL%map(this%tz, vz%x, nel, this%Xh_GL)

       call opgrad(this%vr, this%vs, this%vt, this%tx, c_GL)
       call col3(this%tbf, this%vr, this%tx, n_GL)
       call addcol3(this%tbf, this%vs, this%ty, n_GL)
       call addcol3(this%tbf, this%vt, this%tz, n_GL)
       call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
       call sub2(bfx, this%temp, n)


       call opgrad(this%vr, this%vs, this%vt, this%ty, c_GL)
       call col3(this%tbf, this%vr, this%tx, n_GL)
       call addcol3(this%tbf, this%vs, this%ty, n_GL)
       call addcol3(this%tbf, this%vt, this%tz, n_GL)
       call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
       call sub2(bfy, this%temp, n)

       call opgrad(this%vr, this%vs, this%vt, this%tz, c_GL)
       call col3(this%tbf, this%vr, this%tx, n_GL)
       call addcol3(this%tbf, this%vs, this%ty, n_GL)
       call addcol3(this%tbf, this%vt, this%tz, n_GL)
       call this%GLL_to_GL%map(this%temp, this%tbf, nel, this%Xh_GLL)
       call sub2(bfz, this%temp, n)
       
    else

       do e = 1, coef%msh%nelv
          call this%GLL_to_GL%map(tx, vx%x(1,1,1,e), 1, this%Xh_GL)
          call this%GLL_to_GL%map(ty, vy%x(1,1,1,e), 1, this%Xh_GL)
          call this%GLL_to_GL%map(tz, vz%x(1,1,1,e), 1, this%Xh_GL)

          call opgrad(vr, vs, vt, tx, c_GL, e, e)
          do i = 1, this%Xh_GL%lxyz
             tbfx(i) = tx(i)*vr(i) + ty(i)*vs(i) + tz(i)*vt(i)
          end do

          call opgrad(vr, vs, vt, ty, c_GL, e, e)
          do i = 1, this%Xh_GL%lxyz
             tbfy(i) = tx(i)*vr(i) + ty(i)*vs(i) + tz(i)*vt(i)
          end do

          call opgrad(vr, vs, vt, tz, c_GL, e, e)
          do i = 1, this%Xh_GL%lxyz
             tbfz(i) = tx(i)*vr(i) + ty(i)*vs(i) + tz(i)*vt(i)
          end do

          call this%GLL_to_GL%map(tempx, tbfx, 1, this%Xh_GLL)
          call this%GLL_to_GL%map(tempy, tbfy, 1, this%Xh_GLL)
          call this%GLL_to_GL%map(tempz, tbfz, 1, this%Xh_GLL)

          idx = (e-1)*this%Xh_GLL%lxyz+1
          call sub2(bfx(idx), tempx, this%Xh_GLL%lxyz)
          call sub2(bfy(idx), tempy, this%Xh_GLL%lxyz)
          call sub2(bfz(idx), tempz, this%Xh_GLL%lxyz)
       end do
    end if
    end associate

  end subroutine apply_adv_dealias



  !> Eulerian scheme, add convection term to forcing function
  !! at current time step.
  subroutine advab(this, vx, vy, vz, bfx, bfy, bfz, Xh, coef, n)
    class(adv_no_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    integer, intent(inout) :: n
    real(kind=rp), intent(inout), dimension(n) :: bfx, bfy, bfz
    type(c_ptr) :: bfx_d, bfy_d, bfz_d

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       bfx_d = device_get_ptr(bfx)
       bfy_d = device_get_ptr(bfy)
       bfz_d = device_get_ptr(bfz)
       
       call conv1(this%temp, vx%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (bfx_d, coef%B_d, this%temp_d, n)
       call conv1(this%temp, vy%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (bfy_d, coef%B_d, this%temp_d, n)
       if (coef%Xh%lz .eq. 1) then
          call device_rzero (this%temp_d, n)
       else
          call conv1(this%temp, vz%x, vx%x, vy%x, vz%x, Xh, coef)
          call device_subcol3(bfz_d, coef%B_d, this%temp_d, n)
       end if
    else
       call conv1(this%temp, vx%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (bfx, coef%B, this%temp, n)
       call conv1(this%temp, vy%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (bfy, coef%B, this%temp, n)
       if (coef%Xh%lz .eq. 1) then
          call rzero (this%temp, n)
       else
          call conv1(this%temp, vz%x, vx%x, vy%x, vz%x, Xh, coef)
          call subcol3(bfz, coef%B, this%temp, n)
       end if
    end if

  end subroutine advab

end module advection
