! Copyright (c) 2022, The Neko Authors
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
module scalar_residual_device
  use scalar_residual
  use gather_scatter
  use operators
  use device_math
  use device_mathops
  use, intrinsic :: iso_c_binding
  implicit none
  private
 
  type, public, extends(scalar_residual_t) :: scalar_residual_device_t
   contains
     procedure, nopass :: compute => scalar_residual_device_compute
  end type scalar_residual_device_t

#ifdef HAVE_HIP
  interface
     subroutine scalar_residual_update_hip(s_res_d,f_s_d, n) &
          bind(c, name='scalar_residual_update_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: s_res_d
       type(c_ptr), value :: f_s_d
       integer(c_int) :: n
     end subroutine scalar_residual_update_hip
  end interface
#elif HAVE_CUDA
  
  interface
     subroutine scalar_residual_update_cuda(s_res_d,f_s_d, n) &
          bind(c, name='scalar_residual_update_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: s_res_d
       type(c_ptr), value :: f_s_d
       integer(c_int) :: n
     end subroutine scalar_residual_update_cuda
  end interface
#elif HAVE_OPENCL
  
  interface
     subroutine scalar_residual_update_opencl(s_res_d,f_s_d, n) &
          bind(c, name='scalar_residual_update_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: s_res_d
       type(c_ptr), value :: f_s_d
       integer(c_int) :: n
     end subroutine scalar_residual_update_opencl
  end interface
#endif

  
contains


  subroutine scalar_residual_device_compute(Ax, s, s_res, f_Xh, c_Xh, msh, Xh, &
             lambda, rhocp, bd, dt, n)
    class(ax_t), intent(in) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh    
    type(field_t), intent(inout) :: s
    type(field_t), intent(inout) :: s_res
    type(source_scalar_t), intent(inout) :: f_Xh
    type(coef_t), intent(inout) :: c_Xh
    real(kind=rp), intent(in) :: lambda
    real(kind=rp), intent(in) :: rhocp
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n
    
    call device_cfill(c_Xh%h1_d, lambda, n)
    call device_cfill(c_Xh%h2_d, rhocp * (bd / dt), n)
    c_Xh%ifh2 = .true.
    
    call Ax%compute(s_res%x, s%x, c_Xh, msh, Xh)

#ifdef HAVE_HIP
    call scalar_residual_update_hip(s_res%x_d, f_Xh%s_d, n)
#elif HAVE_CUDA
    call scalar_residual_update_cuda(s_res%x_d, f_Xh%s_d, n)
#elif HAVE_OPENCL
    call scalar_residual_update_opencl(s_res%x_d, f_Xh%s_d, n)
#endif
    
  end subroutine scalar_residual_device_compute
  
end module scalar_residual_device
