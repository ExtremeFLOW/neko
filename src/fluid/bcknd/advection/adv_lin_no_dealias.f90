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
module adv_lin_no_dealias
  use advection, only: advection_lin_t
  use num_types, only: rp
  use math, only: vdot3, sub2, subcol3, rzero
  use utils, only: neko_error
  use space, only: space_t, GL
  use field, only: field_t
  use coefs, only: coef_t
  use scratch_registry, only : neko_scratch_registry
  use neko_config, only: NEKO_BCKND_DEVICE, NEKO_BCKND_SX, NEKO_BCKND_XSMM, &
    NEKO_BCKND_OPENCL, NEKO_BCKND_CUDA, NEKO_BCKND_HIP
  use operators, only: opgrad, conv1, cdtp
  use interpolation, only: interpolator_t
  use device_math
  use device, only: device_free, device_map, device_memcpy, device_get_ptr, &
    HOST_TO_DEVICE
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR, &
    c_associated
  implicit none
  private

  !> Type encapsulating advection routines with no dealiasing applied
  type, public, extends(advection_lin_t) :: adv_lin_no_dealias_t
     real(kind=rp), allocatable :: temp(:)
     type(c_ptr) :: temp_d = C_NULL_PTR
   contains
     !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
     !! the RHS.
     procedure, pass(this) :: compute_linear => linear_advection_no_dealias
     !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
     !! the RHS.
     procedure, pass(this) :: compute_adjoint => adjoint_advection_no_dealias
     !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to
     !! the RHS.
     procedure, pass(this) :: compute_scalar => scalar_advection_no_dealias
     !> Constructor
     procedure, pass(this) :: init => init_no_dealias
     !> Destructor
     procedure, pass(this) :: free => free_no_dealias
  end type adv_lin_no_dealias_t

contains

  !> Constructor
  !! @param coef The coefficients of the (space, mesh) pair.
  subroutine init_no_dealias(this, coef)
    class(adv_lin_no_dealias_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

    allocate(this%temp(coef%dof%size()))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%temp, this%temp_d, coef%dof%size())
    end if

  end subroutine init_no_dealias

  !> Destructor
  subroutine free_no_dealias(this)
    class(adv_lin_no_dealias_t), intent(inout) :: this

    if (allocated(this%temp)) then
       deallocate(this%temp)
    end if
    if (c_associated(this%temp_d)) then
       call device_free(this%temp_d)
    end if
  end subroutine free_no_dealias

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
  subroutine scalar_advection_no_dealias(this, vx, vy, vz, s, fs, Xh, &
                                         coef, n)
    class(adv_lin_no_dealias_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: s
    type(field_t), intent(inout) :: fs
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: n
    type(c_ptr) :: fs_d

    if (NEKO_BCKND_DEVICE .eq. 1) then
       fs_d = fs%x_d

       call conv1(this%temp, s%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fs_d, coef%B_d, this%temp_d, n)
       if (coef%Xh%lz .eq. 1) then
          call device_rzero (this%temp_d, n)
       end if
    else
       ! temp will hold vx*ds/dx + vy*ds/dy + vz*ds/sz
       call conv1(this%temp, s%x, vx%x, vy%x, vz%x, Xh, coef)

       ! fs = fs - B*temp
       call subcol3 (fs%x, coef%B, this%temp, n)
       if (coef%Xh%lz .eq. 1) then
          call rzero (this%temp, n)
       end if
    end if

  end subroutine scalar_advection_no_dealias

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
  subroutine adjoint_advection_no_dealias(this, vx, vy, vz, vxb, vyb, vzb, fx, fy, fz, Xh, coef, n)
    implicit none
    class(adv_lin_no_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: vxb, vyb, vzb
    integer, intent(in) :: n
    type(field_t), intent(inout) :: fx, fy, fz
    type(c_ptr) :: fx_d, fy_d, fz_d
    type(c_ptr) :: vx_d, vy_d, vz_d

    ! these are gradients of U_b (one element)
    real(kind=rp), dimension(Xh%lxyz) :: duxb, dvxb, dwxb
    real(kind=rp), dimension(Xh%lxyz) :: duyb, dvyb, dwyb
    real(kind=rp), dimension(Xh%lxyz) :: duzb, dvzb, dwzb
    ! temporary arrays
    real(kind=rp), dimension(Xh%lxyz) :: tfx, tfy, tfz
    integer :: e, i, idx, idxx


    ! Tim F
    ! I saw this in pnpn_prs_res and it looks like we can get temp arrays this way
    ! I'm not sure if we can do the same for dealiased, but this looks clean for
    ! the no dealiased.
    !
    ! Harry
    type(field_t), pointer :: tduxb, tdvxb, tdwxb, tduyb, tdvyb, tdwyb, tduzb, tdvzb, tdwzb
    integer :: temp_indices(9)


    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_scratch_registry%request_field(tduxb, temp_indices(1))
       call neko_scratch_registry%request_field(tdvxb, temp_indices(2))
       call neko_scratch_registry%request_field(tdwxb, temp_indices(3))
       call neko_scratch_registry%request_field(tduyb, temp_indices(4))
       call neko_scratch_registry%request_field(tdvyb, temp_indices(5))
       call neko_scratch_registry%request_field(tdwyb, temp_indices(6))
       call neko_scratch_registry%request_field(tduzb, temp_indices(7))
       call neko_scratch_registry%request_field(tdvzb, temp_indices(8))
       call neko_scratch_registry%request_field(tdwzb, temp_indices(9))
       fx_d = fx%x_d
       fy_d = fy%x_d
       fz_d = fz%x_d
       vx_d = vx%x_d
       vy_d = vy%x_d
       vz_d = vz%x_d

       ! u . grad U_b^T
       !-----------------------------
       ! take all the gradients
       call opgrad(tduxb%x, tduyb%x, tduzb%x, vxb%x, coef)
       call opgrad(tdvxb%x, tdvyb%x, tdvzb%x, vyb%x, coef)
       call opgrad(tdwxb%x, tdwyb%x, tdwzb%x, vzb%x, coef)


       ! traspose and multiply
       call device_vdot3(this%temp_d, vx_d, vy_d, vz_d, &
                         tduxb%x_d, tdvxb%x_d, tdwxb%x_d, n)
       call device_sub2(fx_d, this%temp_d, n)

       call device_vdot3(this%temp_d, vx_d, vy_d, vz_d, &
                         tduyb%x_d, tdvyb%x_d, tdwyb%x_d, n)
       call device_sub2(fy_d, this%temp_d, n)

       call device_vdot3(this%temp_d, vx_d, vy_d, vz_d, &
                         tduzb%x_d, tdvzb%x_d, tdwzb%x_d, n)
       call device_sub2(fz_d, this%temp_d, n)

       ! \int \grad v . U_b ^ u
       ! with '^' an outer product
       call adjoint_weak_no_dealias_device(fx_d, vx_d, &
                                           vxb%x, vyb%x, vzb%x, &
                                           coef, Xh, n, &
                                           tduxb, tdvxb, tdwxb, &
                                           tduyb, tdvyb, tdwyb)

       call adjoint_weak_no_dealias_device(fy_d, vy_d, &
                                           vxb%x, vyb%x, vzb%x, &
                                           coef, Xh, n, &
                                           tduxb, tdvxb, tdwxb, &
                                           tduyb, tdvyb, tdwyb)

       call adjoint_weak_no_dealias_device(fz_d, vz_d, &
                                           vxb%x, vyb%x, vzb%x, &
                                           coef, Xh, n, &
                                           tduxb, tdvxb, tdwxb, &
                                           tduyb, tdvyb, tdwyb)

       call neko_scratch_registry%relinquish_field(temp_indices)
    else
       do e = 1, coef%msh%nelv
          ! u . grad U_b^T
          !-----------------------------
          call opgrad(duxb, duyb, duzb, vxb%x, coef, e, e)
          call opgrad(dvxb, dvyb, dvzb, vyb%x, coef, e, e)
          call opgrad(dwxb, dwyb, dwzb, vzb%x, coef, e, e)

          ! traspose and multiply
          idx = (e - 1)*Xh%lxyz + 1
          do i = 1, Xh%lxyz
             idxx = idx + i
             fx%x(idxx, 1, 1, 1) = fx%x(idxx, 1, 1, 1) - ( &
               & vx%x(i,1,1,e)*duxb(i) + &
               & vy%x(i,1,1,e)*dvxb(i) + &
               & vz%x(i,1,1,e)*dwxb(i) )

             fy%x(idxx, 1, 1, 1) = fy%x(idxx, 1, 1, 1) - ( &
               & vx%x(i,1,1,e)*duyb(i) + &
               & vy%x(i,1,1,e)*dvyb(i) + &
               & vz%x(i,1,1,e)*dwyb(i))

             fz%x(idxx, 1, 1, 1) = fz%x(idxx, 1, 1, 1) - ( &
               & vx%x(i,1,1,e)*duzb(i) + &
               & vy%x(i,1,1,e)*dvzb(i) + &
               & vz%x(i,1,1,e)*dwzb(i))
          end do

          ! \int \grad v . U_b ^ u
          ! with ^ an outer product
          call adjoint_weak_no_dealias_cpu( &
            & fx%x(:,:,:,e), vx%x(1,1,1,e), &
            & vxb%x(1,1,1,e), vyb%x(1,1,1,e), vzb%x(1,1,1,e), &
            & e, coef, Xh, Xh%lxyz, &
            & duxb, dvxb, dwxb, duyb, dvyb, dwyb)

          call adjoint_weak_no_dealias_cpu( &
            & fy%x(:,:,:,e), vy%x(1,1,1,e), &
            & vxb%x(1,1,1,e), vyb%x(1,1,1,e), vzb%x(1,1,1,e), &
            & e, coef, Xh, Xh%lxyz, &
            & duxb, dvxb, dwxb, duyb, dvyb, dwyb)

          call adjoint_weak_no_dealias_cpu( &
            & fz%x(:,:,:,e), vz%x(1,1,1,e), &
            & vxb%x(1,1,1,e), vyb%x(1,1,1,e), vzb%x(1,1,1,e), &
            & e, coef, Xh, Xh%lxyz, &
            & duxb, dvxb, dwxb, duyb, dvyb, dwyb)
       enddo

    end if

  end subroutine adjoint_advection_no_dealias

  subroutine adjoint_weak_no_dealias_device(f_d, u_i_d, ub, vb, wb, coef, Xh, n, &
                                            work1, work2, work3, w1, w2, w3)
    implicit none
    type(c_ptr), intent(inout) :: f_d
    type(c_ptr), intent(in) :: u_i_d
    real(kind=rp), intent(inout), dimension(n) :: ub, vb, wb
    type(field_t), intent(inout) :: w1, w2, w3
    type(field_t), intent(inout) :: work1, work2, work3
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: n
    type(c_ptr) :: ub_d, vb_d, wb_d
    type(c_ptr) :: work1_d, work2_d, work3_d, w1_d, w2_d, w3_d
    integer :: i

    work1_d = work1%x_d
    work2_d = work2%x_d
    work3_d = work3%x_d
    w1_d = w1%x_d
    w2_d = w2%x_d
    w3_d = w3%x_d
    ub_d = device_get_ptr(ub)
    vb_d = device_get_ptr(vb)
    wb_d = device_get_ptr(wb)

    ! outer product
    call device_col3(work1_d, u_i_d, ub_d, n)
    call device_col3(work2_d, u_i_d, vb_d, n)
    call device_col3(work3_d, u_i_d, wb_d, n)
    ! D^T
    call cdtp(w1%x, work1%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call cdtp(w2%x, work2%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call cdtp(w3%x, work3%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    ! sum them
    call device_add4(work1_d, w1_d, w2_d, w3_d , n)
    call device_sub2(f_d, work1_d, n)
  end subroutine adjoint_weak_no_dealias_device

  subroutine adjoint_weak_no_dealias_cpu(f, u_i, ub, vb, wb, e, coef, Xh, n, work1, work2, work3, w1, w2, w3)
    implicit none
    integer, intent(in) :: e, n
    integer :: i
    real(kind=rp), intent(inout), dimension(n) :: f
    real(kind=rp), intent(inout), dimension(n) :: u_i
    real(kind=rp), intent(inout), dimension(n) :: ub, vb, wb
    real(kind=rp), intent(inout), dimension(n) :: w1, w2, w3
    real(kind=rp), intent(inout), dimension(n) :: work1, work2, work3
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef

    ! outer product
    do i = 1, Xh%lxyz
       work1(i) = u_i(i)*ub(i)
       work2(i) = u_i(i)*vb(i)
       work3(i) = u_i(i)*wb(i)
    end do

    ! D^T
    call cdtp(w1, work1, coef%drdx, coef%dsdx, coef%dtdx, coef, e ,e)
    call cdtp(w2, work2, coef%drdy, coef%dsdy, coef%dtdy, coef, e ,e)
    call cdtp(w3, work3, coef%drdz, coef%dsdz, coef%dtdz, coef, e, e)

    ! sum them
    do i = 1, Xh%lxyz
       f(i) = f(i) - w1(i) + w2(i) + w3(i)
    end do
  end subroutine adjoint_weak_no_dealias_cpu



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
  subroutine linear_advection_no_dealias(this, vx, vy, vz, vxb, vyb, vzb, fx, fy, fz, Xh, coef, n)
    implicit none
    class(adv_lin_no_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: vxb, vyb, vzb
    integer, intent(in) :: n
    type(field_t), intent(inout) :: fx, fy, fz
    type(c_ptr) :: fx_d, fy_d, fz_d

    if (NEKO_BCKND_DEVICE .eq. 1) then
       fx_d = fx%x_d
       fy_d = fy%x_d
       fz_d = fz%x_d

       ! (x)
       call conv1(this%temp, vx%x, vxb%x, vyb%x, vzb%x, Xh, coef)
       call device_subcol3 (fx_d, coef%B_d, this%temp_d, n)
       call conv1(this%temp, vxb%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fx_d, coef%B_d, this%temp_d, n)

       ! (y)
       call conv1(this%temp, vy%x, vxb%x, vyb%x, vzb%x, Xh, coef)
       call device_subcol3 (fy_d, coef%B_d, this%temp_d, n)
       call conv1(this%temp, vyb%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fy_d, coef%B_d, this%temp_d, n)

       ! (z)
       if (coef%Xh%lz .eq. 1) then
          call device_rzero (this%temp_d, n)
       else
          call conv1(this%temp, vz%x, vxb%x, vyb%x, vzb%x, Xh, coef)
          call device_subcol3(fz_d, coef%B_d, this%temp_d, n)
          call conv1(this%temp, vzb%x, vx%x, vy%x, vz%x, Xh, coef)
          call device_subcol3(fz_d, coef%B_d, this%temp_d, n)
       end if
    else
       ! (x)
       call conv1(this%temp, vx%x, vxb%x, vyb%x, vzb%x, Xh, coef)
       call subcol3 (fx%x, coef%B, this%temp, n)
       call conv1(this%temp, vxb%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (fx%x, coef%B, this%temp, n)

       ! (y)
       call conv1(this%temp, vy%x, vxb%x, vyb%x, vzb%x, Xh, coef)
       call subcol3 (fy%x, coef%B, this%temp, n)
       call conv1(this%temp, vyb%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (fy%x, coef%B, this%temp, n)

       ! (z)
       if (coef%Xh%lz .eq. 1) then
          call rzero (this%temp, n)
       else
          call conv1(this%temp, vz%x, vxb%x, vyb%x, vzb%x, Xh, coef)
          call subcol3(fz%x, coef%B, this%temp, n)
          call conv1(this%temp, vzb%x, vx%x, vy%x, vz%x, Xh, coef)
          call subcol3(fz%x, coef%B, this%temp, n)
       end if
    end if

  end subroutine linear_advection_no_dealias
end module adv_lin_no_dealias
