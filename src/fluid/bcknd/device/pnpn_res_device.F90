! Copyright (c) 2022-2023, The Neko Authors
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
module pnpn_res_device
  use gather_scatter, only : gs_t, GS_OP_ADD
  use operators
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use facet_normal, only : facet_normal_t
  use mesh, only : mesh_t
  use num_types, only : rp, c_rp
  use space, only : space_t
  use device_math
  use device_mathops
  use pnpn_residual, only : pnpn_prs_res_t, pnpn_vel_res_t
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  use scratch_registry, only : neko_scratch_registry
  implicit none
  private
 
  type, public, extends(pnpn_prs_res_t) :: pnpn_prs_res_device_t
   contains
     procedure, nopass :: compute => pnpn_prs_res_device_compute
  end type pnpn_prs_res_device_t

  type, public, extends(pnpn_vel_res_t) :: pnpn_vel_res_device_t
   contains
     procedure, nopass :: compute => pnpn_vel_res_device_compute
  end type pnpn_vel_res_device_t

#ifdef HAVE_HIP
    interface
     subroutine pnpn_prs_res_part1_hip(ta1_d, ta2_d, ta3_d, &
          wa1_d, wa2_d, wa3_d, f_u_d, f_v_d, f_w_d, &
          B_d, h1_d, mu, rho, n) &
          bind(c, name='pnpn_prs_res_part1_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: wa1_d, wa2_d, wa3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       type(c_ptr), value :: B_d, h1_d
       real(c_rp) :: mu, rho
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part1_hip
  end interface

  interface
     subroutine pnpn_prs_res_part2_hip(p_res_d, wa1_d, wa2_d, wa3_d, n) &
          bind(c, name='pnpn_prs_res_part2_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: p_res_d, wa1_d, wa2_d, wa3_d
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part2_hip
  end interface

  interface
     subroutine pnpn_prs_res_part3_hip(p_res_d, ta1_d, ta2_d, ta3_d, dtbd, n) &
          bind(c, name='pnpn_prs_res_part3_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_res_d, ta1_d, ta2_d, ta3_d
       real(c_rp) :: dtbd
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part3_hip
  end interface
  
  interface
     subroutine pnpn_vel_res_update_hip(u_res_d, v_res_d, w_res_d, &
          ta1_d, ta2_d, ta3_d, f_u_d, f_v_d, f_w_d, n) &
          bind(c, name='pnpn_vel_res_update_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: u_res_d, v_res_d, w_res_d
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       integer(c_int) :: n
     end subroutine pnpn_vel_res_update_hip
  end interface
#elif HAVE_CUDA
  interface
     subroutine pnpn_prs_res_part1_cuda(ta1_d, ta2_d, ta3_d, &
          wa1_d, wa2_d, wa3_d, f_u_d, f_v_d, f_w_d, &
          B_d, h1_d, mu, rho, n) &
          bind(c, name='pnpn_prs_res_part1_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: wa1_d, wa2_d, wa3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       type(c_ptr), value :: B_d, h1_d
       real(c_rp) :: mu, rho
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part1_cuda
  end interface

  interface
     subroutine pnpn_prs_res_part2_cuda(p_res_d, wa1_d, wa2_d, wa3_d, n) &
          bind(c, name='pnpn_prs_res_part2_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: p_res_d, wa1_d, wa2_d, wa3_d
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part2_cuda
  end interface

  interface
     subroutine pnpn_prs_res_part3_cuda(p_res_d, ta1_d, ta2_d, ta3_d, dtbd, n) &
          bind(c, name='pnpn_prs_res_part3_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_res_d, ta1_d, ta2_d, ta3_d
       real(c_rp) :: dtbd
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part3_cuda
  end interface
  
  interface
     subroutine pnpn_vel_res_update_cuda(u_res_d, v_res_d, w_res_d, &
          ta1_d, ta2_d, ta3_d, f_u_d, f_v_d, f_w_d, n) &
          bind(c, name='pnpn_vel_res_update_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: u_res_d, v_res_d, w_res_d
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       integer(c_int) :: n
     end subroutine pnpn_vel_res_update_cuda
  end interface
#elif HAVE_OPENCL
  interface
     subroutine pnpn_prs_res_part1_opencl(ta1_d, ta2_d, ta3_d, &
          wa1_d, wa2_d, wa3_d, f_u_d, f_v_d, f_w_d, &
          B_d, h1_d, mu, rho, n) &
          bind(c, name='pnpn_prs_res_part1_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: wa1_d, wa2_d, wa3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       type(c_ptr), value :: B_d, h1_d
       real(c_rp) :: mu, rho
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part1_opencl
  end interface

  interface
     subroutine pnpn_prs_res_part2_opencl(p_res_d, wa1_d, wa2_d, wa3_d, n) &
          bind(c, name='pnpn_prs_res_part2_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: p_res_d, wa1_d, wa2_d, wa3_d
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part2_opencl
  end interface

  interface
     subroutine pnpn_prs_res_part3_opencl(p_res_d, ta1_d, ta2_d, ta3_d, dtbd, n) &
          bind(c, name='pnpn_prs_res_part3_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_res_d, ta1_d, ta2_d, ta3_d
       real(c_rp) :: dtbd
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part3_opencl
  end interface
  
  interface
     subroutine pnpn_vel_res_update_opencl(u_res_d, v_res_d, w_res_d, &
          ta1_d, ta2_d, ta3_d, f_u_d, f_v_d, f_w_d, n) &
          bind(c, name='pnpn_vel_res_update_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: u_res_d, v_res_d, w_res_d
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       integer(c_int) :: n
     end subroutine pnpn_vel_res_update_opencl
  end interface
#endif

  
contains

  subroutine pnpn_prs_res_device_compute(p, p_res, u, v, w, u_e, v_e, w_e, &
       f_x, f_y, f_z, c_Xh, gs_Xh, bc_prs_surface, bc_sym_surface, Ax, bd, dt,&
       mu, rho)
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(inout) :: u_e, v_e, w_e
    type(field_t), intent(inout) :: p_res
    type(field_t), intent(inout) :: f_x, f_y, f_z
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(facet_normal_t), intent(inout) :: bc_prs_surface
    type(facet_normal_t), intent(inout) :: bc_sym_surface
    class(Ax_t), intent(inout) :: Ax
    real(kind=rp), intent(inout) :: bd
    real(kind=rp), intent(in) :: dt
    real(kind=rp), intent(in) :: mu
    real(kind=rp), intent(in) :: rho
    real(kind=rp) :: dtbd
    integer :: n, gdim
    type(field_t), pointer :: ta1, ta2, ta3, wa1, wa2, wa3, work1, work2
    integer :: temp_indices(8)

    call neko_scratch_registry%request_field(ta1, temp_indices(1))
    call neko_scratch_registry%request_field(ta2, temp_indices(2))
    call neko_scratch_registry%request_field(ta3, temp_indices(3))
    call neko_scratch_registry%request_field(wa1, temp_indices(4))
    call neko_scratch_registry%request_field(wa2, temp_indices(5))
    call neko_scratch_registry%request_field(wa3, temp_indices(6))
    call neko_scratch_registry%request_field(work1, temp_indices(7))
    call neko_scratch_registry%request_field(work2, temp_indices(8))

    n = u%dof%size()
    gdim = c_Xh%msh%gdim
    
    call curl(ta1, ta2, ta3, u_e, v_e, w_e, work1, work2, c_Xh)
    call curl(wa1, wa2, wa3, ta1, ta2, ta3, work1, work2, c_Xh)

#ifdef HAVE_HIP
    call pnpn_prs_res_part1_hip(ta1%x_d, ta2%x_d, ta3%x_d, &
         wa1%x_d, wa2%x_d, wa3%x_d, f_x%x_d, f_y%x_d, f_z%x_d, &
         c_Xh%B_d, c_Xh%h1_d, mu, rho, n) 

#elif HAVE_CUDA
    call pnpn_prs_res_part1_cuda(ta1%x_d, ta2%x_d, ta3%x_d, &
         wa1%x_d, wa2%x_d, wa3%x_d, f_x%x_d, f_y%x_d, f_z%x_d, &
         c_Xh%B_d, c_Xh%h1_d, mu, rho, n) 
#elif HAVE_OPENCL
    call pnpn_prs_res_part1_opencl(ta1%x_d, ta2%x_d, ta3%x_d, &
         wa1%x_d, wa2%x_d, wa3%x_d, f_x%x_d, f_z%x_d, f_z%x_d, &
         c_Xh%B_d, c_Xh%h1_d, mu, rho, n) 
#endif
     c_Xh%ifh2 = .false.
         
    call gs_Xh%op(ta1, GS_OP_ADD) 
    call gs_Xh%op(ta2, GS_OP_ADD) 
    call gs_Xh%op(ta3, GS_OP_ADD) 

    call device_opcolv(ta1%x_d, ta2%x_d, ta3%x_d, c_Xh%Binv_d, gdim, n)

    call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call cdtp(wa3%x, ta3%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)

    call Ax%compute(p_res%x,p%x,c_Xh,p%msh,p%Xh)

#ifdef HAVE_HIP
    call pnpn_prs_res_part2_hip(p_res%x_d, wa1%x_d, wa2%x_d, wa3%x_d, n);
#elif HAVE_CUDA
    call pnpn_prs_res_part2_cuda(p_res%x_d, wa1%x_d, wa2%x_d, wa3%x_d, n);
#elif HAVE_OPENCL
    call pnpn_prs_res_part2_opencl(p_res%x_d, wa1%x_d, wa2%x_d, wa3%x_d, n);
#endif

    !
    ! Surface velocity terms
    call device_rzero(wa1%x_d, n)
    call device_rzero(wa2%x_d, n)
    call device_rzero(wa3%x_d, n)
    dtbd = 1.0_rp

    call bc_sym_surface%apply_surfvec_dev(wa1%x_d, wa2%x_d, wa3%x_d, ta1%x_d, ta2%x_d, ta3%x_d)

#ifdef HAVE_HIP
    call pnpn_prs_res_part3_hip(p_res%x_d, wa1%x_d, wa2%x_d, wa3%x_d, dtbd, n);
#elif HAVE_CUDA
    call pnpn_prs_res_part3_cuda(p_res%x_d, wa1%x_d, wa2%x_d, wa3%x_d, dtbd, n);
#elif HAVE_OPENCL
    call pnpn_prs_res_part3_opencl(p_res%x_d, wa1%x_d, wa2%x_d, wa3%x_d, dtbd, n);
#endif
   !
    dtbd = bd / dt

    call device_rzero(ta1%x_d, n)
    call device_rzero(ta2%x_d, n)
    call device_rzero(ta3%x_d, n)

    call bc_prs_surface%apply_surfvec_dev(ta1%x_d, ta2%x_d, ta3%x_d, &
         u%x_d, v%x_d, w%x_d)

#ifdef HAVE_HIP
    call pnpn_prs_res_part3_hip(p_res%x_d, ta1%x_d, ta2%x_d, ta3%x_d, dtbd, n);
#elif HAVE_CUDA
    call pnpn_prs_res_part3_cuda(p_res%x_d, ta1%x_d, ta2%x_d, ta3%x_d, dtbd, n);
#elif HAVE_OPENCL
    call pnpn_prs_res_part3_opencl(p_res%x_d, ta1%x_d, ta2%x_d, ta3%x_d, dtbd, n);
#endif
        
  call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine pnpn_prs_res_device_compute

  subroutine pnpn_vel_res_device_compute(Ax, u, v, w, u_res, v_res, w_res, &
       p, f_x, f_y, f_z, c_Xh, msh, Xh, mu, rho, bd, dt, n)
    class(ax_t), intent(in) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh    
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(inout) :: u_res, v_res, w_res
    type(field_t), intent(inout) :: f_x, f_y, f_z
    type(coef_t), intent(inout) :: c_Xh
    real(kind=rp), intent(in) :: mu 
    real(kind=rp), intent(in) :: rho
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n
    integer :: temp_indices(3)
    type(field_t), pointer :: ta1, ta2, ta3

    call device_cfill(c_Xh%h1_d, mu, n)
    call device_cfill(c_Xh%h2_d, rho * (bd / dt), n)
    c_Xh%ifh2 = .true.
    
    call Ax%compute(u_res%x, u%x, c_Xh, msh, Xh)
    call Ax%compute(v_res%x, v%x, c_Xh, msh, Xh)
    call Ax%compute(w_res%x, w%x, c_Xh, msh, Xh)

    call neko_scratch_registry%request_field(ta1, temp_indices(1))
    call neko_scratch_registry%request_field(ta2, temp_indices(2))
    call neko_scratch_registry%request_field(ta3, temp_indices(3))

    call opgrad(ta1%x, ta2%x, ta3%x, p%x, c_Xh)

#ifdef HAVE_HIP
    call pnpn_vel_res_update_hip(u_res%x_d, v_res%x_d, w_res%x_d, &
         ta1%x_d, ta2%x_d, ta3%x_d, f_x%x_d, f_y%x_d, f_z%x_d, n)
#elif HAVE_CUDA
    call pnpn_vel_res_update_cuda(u_res%x_d, v_res%x_d, w_res%x_d, &
         ta1%x_d, ta2%x_d, ta3%x_d, f_x%x_d, f_y%x_d, f_z%x_d, n)
#elif HAVE_OPENCL
    call pnpn_vel_res_update_opencl(u_res%x_d, v_res%x_d, w_res%x_d, &
         ta1%x_d, ta2%x_d, ta3%x_d, f_x%x_d, f_y%x_d, f_z%x_d, n)
#endif
    
    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine pnpn_vel_res_device_compute

end module pnpn_res_device
