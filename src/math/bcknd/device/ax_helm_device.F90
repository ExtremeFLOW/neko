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
module ax_helm_device
  use ax_product, only : ax_t
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use field, only : field_t
  use mesh, only : mesh_t
  use device_math, only : device_addcol4
  use device, only : device_get_ptr
  use num_types, only : rp
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  implicit none
  private
  
  type, public, extends(ax_t) :: ax_helm_device_t
   contains
     procedure, nopass :: compute => ax_helm_device_compute
  end type ax_helm_device_t

#ifdef HAVE_HIP
    interface
     subroutine hip_ax_helm(w_d, u_d, &
          dx_d, dy_d, dz_d, dxt_d, dyt_d, dzt_d, &
          h1_d, g11_d, g22_d, g33_d, g12_d, g13_d, g23_d, nelv, lx) &
          bind(c, name='hip_ax_helm')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: w_d, u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: dxt_d, dyt_d, dzt_d
       type(c_ptr), value :: h1_d, g11_d, g22_d, g33_d, g12_d, g13_d, g23_d
       integer(c_int) :: nel, lx
     end subroutine hip_ax_helm
  end interface
#elif HAVE_CUDA
    interface
     subroutine cuda_ax_helm(w_d, u_d, &
          dx_d, dy_d, dz_d, dxt_d, dyt_d, dzt_d,&
          h1_d, g11_d, g22_d, g33_d, g12_d, g13_d, g23_d, nelv, lx) &
          bind(c, name='cuda_ax_helm')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: w_d, u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: dxt_d, dyt_d, dzt_d
       type(c_ptr), value :: h1_d, g11_d, g22_d, g33_d, g12_d, g13_d, g23_d
       integer(c_int) :: nel, lx
     end subroutine cuda_ax_helm
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_ax_helm(w_d, u_d, &
          dx_d, dy_d, dz_d, dxt_d, dyt_d, dzt_d, &
          h1_d, g11_d, g22_d, g33_d, g12_d, g13_d, g23_d, nelv, lx) &
          bind(c, name='opencl_ax_helm')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: w_d, u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: dxt_d, dyt_d, dzt_d
       type(c_ptr), value :: h1_d, g11_d, g22_d, g33_d, g12_d, g13_d, g23_d
       integer(c_int) :: nel, lx
     end subroutine opencl_ax_helm
  end interface
#endif

contains

  subroutine ax_helm_device_compute(w, u, coef, msh, Xh)
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    type(c_ptr) :: u_d, w_d

    u_d = device_get_ptr(u)
    w_d = device_get_ptr(w)

#ifdef HAVE_HIP
    call hip_ax_helm(w_d, u_d, Xh%dx_d, Xh%dy_d, Xh%dz_d, &
         Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, coef%h1_d, &
         coef%G11_d, coef%G22_d, coef%G33_d, &
         coef%G12_d, coef%G13_d, coef%G23_d, &
         msh%nelv, Xh%lx)
#elif HAVE_CUDA
    call cuda_ax_helm(w_d, u_d, Xh%dx_d, Xh%dy_d, Xh%dz_d, &
         Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, coef%h1_d, &
         coef%G11_d, coef%G22_d, coef%G33_d, &
         coef%G12_d, coef%G13_d, coef%G23_d, &
         msh%nelv, Xh%lx)
#elif HAVE_OPENCL
    call opencl_ax_helm(w_d, u_d, Xh%dx_d, Xh%dy_d, Xh%dz_d, &
         Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, coef%h1_d, &
         coef%G11_d, coef%G22_d, coef%G33_d, &
         coef%G12_d, coef%G13_d, coef%G23_d, &
         msh%nelv, Xh%lx)
#endif

    if (coef%ifh2) then
       call device_addcol4(w_d ,coef%h2_d, coef%B_d, u_d, coef%dof%size())
    end if
    
  end subroutine ax_helm_device_compute
  
end module ax_helm_device


