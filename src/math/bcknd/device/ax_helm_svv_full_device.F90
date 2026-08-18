! Copyright (c) 2025-2026, The Neko Authors
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
module ax_helm_svv_full_device
  use ax_helm_svv_full, only : ax_helm_svv_full_t
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

  !> Device asymmetric full-stress SVV Helmholtz operator.
  type, public, extends(ax_helm_svv_full_t) :: ax_helm_svv_full_device_t
   contains
     !> Compute the coupled vector product.
     procedure, pass(this) :: compute_vector => &
          ax_helm_svv_full_device_compute_vector
  end type ax_helm_svv_full_device_t

#ifdef HAVE_HIP
  interface
     subroutine hip_ax_helm_svv_full(au_d, av_d, aw_d, u_d, v_d, w_d, &
          dx_d, dy_d, dz_d, h1_d, &
          drdx_d, drdy_d, drdz_d, dsdx_d, dsdy_d, dsdz_d, &
          dtdx_d, dtdy_d, dtdz_d, jacinv_d, w3_d, svv_h1_d, &
          filter_r_d, filter_s_d, filter_t_d, nelv, lx) &
          bind(c, name='hip_ax_helm_svv_full')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       type(c_ptr), value :: au_d, av_d, aw_d
       type(c_ptr), value :: u_d, v_d, w_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, h1_d
       type(c_ptr), value :: drdx_d, drdy_d, drdz_d
       type(c_ptr), value :: dsdx_d, dsdy_d, dsdz_d
       type(c_ptr), value :: dtdx_d, dtdy_d, dtdz_d
       type(c_ptr), value :: jacinv_d, w3_d, svv_h1_d
       type(c_ptr), value :: filter_r_d, filter_s_d, filter_t_d
       integer(c_int) :: nelv, lx
     end subroutine hip_ax_helm_svv_full
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_ax_helm_svv_full(au_d, av_d, aw_d, u_d, v_d, w_d, &
          dx_d, dy_d, dz_d, h1_d, &
          drdx_d, drdy_d, drdz_d, dsdx_d, dsdy_d, dsdz_d, &
          dtdx_d, dtdy_d, dtdz_d, jacinv_d, w3_d, svv_h1_d, &
          filter_r_d, filter_s_d, filter_t_d, nelv, lx) &
          bind(c, name='cuda_ax_helm_svv_full')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       type(c_ptr), value :: au_d, av_d, aw_d
       type(c_ptr), value :: u_d, v_d, w_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, h1_d
       type(c_ptr), value :: drdx_d, drdy_d, drdz_d
       type(c_ptr), value :: dsdx_d, dsdy_d, dsdz_d
       type(c_ptr), value :: dtdx_d, dtdy_d, dtdz_d
       type(c_ptr), value :: jacinv_d, w3_d, svv_h1_d
       type(c_ptr), value :: filter_r_d, filter_s_d, filter_t_d
       integer(c_int) :: nelv, lx
     end subroutine cuda_ax_helm_svv_full
  end interface
#endif

contains

  !> Compute the coupled asymmetric full-stress SVV Helmholtz product.
  !! @param this Device asymmetric full-stress SVV operator.
  !! @param au Result for the first component.
  !! @param av Result for the second component.
  !! @param aw Result for the third component.
  !! @param u First input component.
  !! @param v Second input component.
  !! @param w Third input component.
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space.
  subroutine ax_helm_svv_full_device_compute_vector(this, au, av, aw, &
       u, v, w, coef, msh, Xh)
    class(ax_helm_svv_full_device_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    type(c_ptr) :: au_d, av_d, aw_d, u_d, v_d, w_d
    type(c_ptr) :: filter_r_d, filter_s_d, filter_t_d

    au_d = device_get_ptr(au)
    av_d = device_get_ptr(av)
    aw_d = device_get_ptr(aw)
    u_d = device_get_ptr(u)
    v_d = device_get_ptr(v)
    w_d = device_get_ptr(w)

    if (index(this%svv%direction, "r") > 0) then
       filter_r_d = this%svv%filter%fh_d
    else
       filter_r_d = this%svv%ident_d
    end if

    if (index(this%svv%direction, "s") > 0) then
       filter_s_d = this%svv%filter%fht_d
    else
       filter_s_d = this%svv%ident_d
    end if

    if (index(this%svv%direction, "t") > 0) then
       filter_t_d = this%svv%filter%fht_d
    else
       filter_t_d = this%svv%ident_d
    end if

#ifdef HAVE_HIP
    call hip_ax_helm_svv_full(au_d, av_d, aw_d, u_d, v_d, w_d, &
         Xh%dx_d, Xh%dy_d, Xh%dz_d, coef%h1_d, &
         coef%drdx_d, coef%drdy_d, coef%drdz_d, &
         coef%dsdx_d, coef%dsdy_d, coef%dsdz_d, &
         coef%dtdx_d, coef%dtdy_d, coef%dtdz_d, coef%jacinv_d, &
         Xh%w3_d, this%svv%h1_d, filter_r_d, filter_s_d, filter_t_d, &
         msh%nelv, Xh%lx)
#elif HAVE_CUDA
    call cuda_ax_helm_svv_full(au_d, av_d, aw_d, u_d, v_d, w_d, &
         Xh%dx_d, Xh%dy_d, Xh%dz_d, coef%h1_d, &
         coef%drdx_d, coef%drdy_d, coef%drdz_d, &
         coef%dsdx_d, coef%dsdy_d, coef%dsdz_d, &
         coef%dtdx_d, coef%dtdy_d, coef%dtdz_d, coef%jacinv_d, &
         Xh%w3_d, this%svv%h1_d, filter_r_d, filter_s_d, filter_t_d, &
         msh%nelv, Xh%lx)
#elif HAVE_OPENCL
    call neko_error('OPENCL is not implemented for full SVV')
#endif

    if (coef%ifh2) then
       call device_addcol4(au_d, coef%h2_d, coef%B_d, u_d, coef%dof%size())
       call device_addcol4(av_d, coef%h2_d, coef%B_d, v_d, coef%dof%size())
       call device_addcol4(aw_d, coef%h2_d, coef%B_d, w_d, coef%dof%size())
    end if

  end subroutine ax_helm_svv_full_device_compute_vector

end module ax_helm_svv_full_device
