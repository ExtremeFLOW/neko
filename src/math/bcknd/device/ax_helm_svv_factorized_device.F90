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
!> Device implementation of the factorized SVV Helmholtz operator.
module ax_helm_svv_factorized_device
  use ax_helm_svv, only : ax_helm_svv_t
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use device_math, only : device_addcol4
  use device, only : device_get_ptr
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  implicit none
  private

  !> Device matrix-vector product for a factorized SVV Helmholtz problem.
  type, public, extends(ax_helm_svv_t) :: ax_helm_svv_factorized_device_t
   contains
     !> Compute the product.
     procedure, pass(this) :: compute => ax_helm_svv_factorized_device_compute
  end type ax_helm_svv_factorized_device_t

#ifdef HAVE_HIP
  interface
     subroutine hip_ax_helm_svv_factorized(w_d, u_d, dx_d, dy_d, dz_d, h1_d, &
          Br_d, Bs_d, Bt_d, svv_h1_d, G11_d, G22_d, G33_d, &
          G12_d, G13_d, G23_d, nelv, lx) &
          bind(c, name='hip_ax_helm_svv_factorized')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       type(c_ptr), value :: w_d, u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, h1_d
       type(c_ptr), value :: Br_d, Bs_d, Bt_d, svv_h1_d
       type(c_ptr), value :: G11_d, G22_d, G33_d
       type(c_ptr), value :: G12_d, G13_d, G23_d
       integer(c_int) :: nelv, lx
     end subroutine hip_ax_helm_svv_factorized
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_ax_helm_svv_factorized(w_d, u_d, dx_d, dy_d, dz_d, h1_d, &
          Br_d, Bs_d, Bt_d, svv_h1_d, G11_d, G22_d, G33_d, &
          G12_d, G13_d, G23_d, nelv, lx) &
          bind(c, name='cuda_ax_helm_svv_factorized')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       type(c_ptr), value :: w_d, u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, h1_d
       type(c_ptr), value :: Br_d, Bs_d, Bt_d, svv_h1_d
       type(c_ptr), value :: G11_d, G22_d, G33_d
       type(c_ptr), value :: G12_d, G13_d, G23_d
       integer(c_int) :: nelv, lx
     end subroutine cuda_ax_helm_svv_factorized
  end interface
#endif

contains

  !> Compute the factorized SVV Helmholtz product.
  !! @param this Device factorized SVV Helmholtz operator.
  !! @param w Result.
  !! @param u Input field.
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space.
  subroutine ax_helm_svv_factorized_device_compute(this, w, u, coef, msh, Xh)
    class(ax_helm_svv_factorized_device_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    type(c_ptr) :: u_d, w_d

    u_d = device_get_ptr(u)
    w_d = device_get_ptr(w)

#ifdef HAVE_HIP
    call hip_ax_helm_svv_factorized(w_d, u_d, Xh%dx_d, Xh%dy_d, Xh%dz_d, &
         coef%h1_d, this%svv%Br_d, this%svv%Bs_d, this%svv%Bt_d, &
         this%svv%h1_d, coef%G11_d, coef%G22_d, coef%G33_d, &
         coef%G12_d, coef%G13_d, coef%G23_d, &
         msh%nelv, Xh%lx)
#elif HAVE_CUDA
    call cuda_ax_helm_svv_factorized(w_d, u_d, Xh%dx_d, Xh%dy_d, Xh%dz_d, &
         coef%h1_d, this%svv%Br_d, this%svv%Bs_d, this%svv%Bt_d, &
         this%svv%h1_d, coef%G11_d, coef%G22_d, coef%G33_d, &
         coef%G12_d, coef%G13_d, coef%G23_d, &
         msh%nelv, Xh%lx)
#elif HAVE_OPENCL
    call neko_error('OPENCL is not implemented for factorized SVV')
#endif

    if (coef%ifh2) then
       call device_addcol4(w_d, coef%h2_d, coef%B_d, u_d, coef%dof%size())
    end if

  end subroutine ax_helm_svv_factorized_device_compute

end module ax_helm_svv_factorized_device
