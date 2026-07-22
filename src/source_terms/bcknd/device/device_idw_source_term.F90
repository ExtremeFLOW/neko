! Copyright (c) 2024-2026, The Neko Authors
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
!> Device backend dispatch for the IDW immersed-boundary source term
module device_idw_source_term
  use num_types, only : rp, c_rp
  use utils, only : neko_error
  use device, only : glb_cmd_queue
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  implicit none
  private

#ifdef HAVE_OPENCL
  interface
     subroutine opencl_idw_gather_one_sided(fu, fv, fw, &
          fu_ib, fv_ib, fw_ib, fum_ib, fvm_ib, fwm_ib, &
          x, y, z, ds, pmsk, w, wm, lpx, lpy, lpz, &
          active_el, el_off, el_lag, n_active, lx3, &
          dt, rmax, pwr, eps, wtol, cmd_queue) &
          bind(c, name = 'opencl_idw_gather_one_sided')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: fu, fv, fw
       type(c_ptr), value :: fu_ib, fv_ib, fw_ib
       type(c_ptr), value :: fum_ib, fvm_ib, fwm_ib
       type(c_ptr), value :: x, y, z, ds, pmsk, w, wm
       type(c_ptr), value :: lpx, lpy, lpz
       type(c_ptr), value :: active_el, el_off, el_lag
       integer(c_int) :: n_active, lx3
       real(c_rp) :: dt, rmax, pwr, eps, wtol
       type(c_ptr), value :: cmd_queue
     end subroutine opencl_idw_gather_one_sided
  end interface
#elif HAVE_METAL
  interface
     subroutine metal_idw_gather_one_sided(fu, fv, fw, &
          fu_ib, fv_ib, fw_ib, fum_ib, fvm_ib, fwm_ib, &
          x, y, z, ds, pmsk, w, wm, lpx, lpy, lpz, &
          active_el, el_off, el_lag, n_active, lx3, &
          dt, rmax, pwr, eps, wtol) &
          bind(c, name = 'metal_idw_gather_one_sided')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: fu, fv, fw
       type(c_ptr), value :: fu_ib, fv_ib, fw_ib
       type(c_ptr), value :: fum_ib, fvm_ib, fwm_ib
       type(c_ptr), value :: x, y, z, ds, pmsk, w, wm
       type(c_ptr), value :: lpx, lpy, lpz
       type(c_ptr), value :: active_el, el_off, el_lag
       integer(c_int) :: n_active, lx3
       real(c_rp) :: dt, rmax, pwr, eps, wtol
     end subroutine metal_idw_gather_one_sided
  end interface
#endif

  public :: device_idw_gather

contains

  !> Assemble the IDW forcing on the device (atomic-free gather kernel).
  !! Argument order mirrors the OpenCL kernel; see idw_kernel.cl.
  subroutine device_idw_gather(fu, fv, fw, fu_ib, fv_ib, fw_ib, &
       fum_ib, fvm_ib, fwm_ib, x, y, z, ds, pmsk, w, wm, &
       lpx, lpy, lpz, active_el, el_off, el_lag, n_active, lx3, &
       dt, rmax, pwr, eps, wtol)
    type(c_ptr) :: fu, fv, fw
    type(c_ptr) :: fu_ib, fv_ib, fw_ib
    type(c_ptr) :: fum_ib, fvm_ib, fwm_ib
    type(c_ptr) :: x, y, z, ds, pmsk, w, wm
    type(c_ptr) :: lpx, lpy, lpz
    type(c_ptr) :: active_el, el_off, el_lag
    integer, intent(in) :: n_active, lx3
    real(kind=rp), intent(in) :: dt, rmax, pwr, eps, wtol

#ifdef HAVE_HIP
    call neko_error('IDW source term: HIP kernel not implemented')
#elif HAVE_CUDA
    call neko_error('IDW source term: CUDA kernel not implemented')
#elif HAVE_OPENCL
    call opencl_idw_gather_one_sided(fu, fv, fw, &
         fu_ib, fv_ib, fw_ib, fum_ib, fvm_ib, fwm_ib, &
         x, y, z, ds, pmsk, w, wm, lpx, lpy, lpz, &
         active_el, el_off, el_lag, n_active, lx3, &
         dt, rmax, pwr, eps, wtol, glb_cmd_queue)
#elif HAVE_METAL
    call metal_idw_gather_one_sided(fu, fv, fw, &
         fu_ib, fv_ib, fw_ib, fum_ib, fvm_ib, fwm_ib, &
         x, y, z, ds, pmsk, w, wm, lpx, lpy, lpz, &
         active_el, el_off, el_lag, n_active, lx3, &
         dt, rmax, pwr, eps, wtol)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_idw_gather

end module device_idw_source_term
