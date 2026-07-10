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
!> Device dispatch for LPT elastic wall-collision handling.
module lpt_wall_collision_device
  use vector, only : vector_t
  use utils, only : neko_error
  use mesh, only : mesh_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use device, only : glb_cmd_queue
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr
  implicit none
  private

#ifdef HAVE_HIP
  interface
     !> HIP kernel entry point for elastic wall-collision reflection.
     subroutine hip_lpt_handle_elastic_wall_collisions(wall_facet_mask, &
          el_list, x_old, y_old, z_old, x, y, z, d, u, v, w, u_lag, &
          v_lag, w_lag, u_laglag, v_laglag, w_laglag, acc_xlag, &
          acc_ylag, acc_zlag, acc_xlaglag, acc_ylaglag, acc_zlaglag, &
          u_old, v_old, w_old, acc_x, acc_y, acc_z, dm_x, dm_y, dm_z, &
          nx, ny, nz, n, gdim, nelv, lx, ly, lz, lag_len, strm) &
          bind(c, name = 'hip_lpt_handle_elastic_wall_collisions')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: wall_facet_mask, el_list
       type(c_ptr), value :: x_old, y_old, z_old
       type(c_ptr), value :: x, y, z, d, u, v, w
       type(c_ptr), value :: u_lag, v_lag, w_lag
       type(c_ptr), value :: u_laglag, v_laglag, w_laglag
       type(c_ptr), value :: acc_xlag, acc_ylag, acc_zlag
       type(c_ptr), value :: acc_xlaglag, acc_ylaglag, acc_zlaglag
       type(c_ptr), value :: u_old, v_old, w_old
       type(c_ptr), value :: acc_x, acc_y, acc_z
       type(c_ptr), value :: dm_x, dm_y, dm_z
       type(c_ptr), value :: nx, ny, nz, strm
       integer(c_int) :: n, gdim, nelv, lx, ly, lz, lag_len
     end subroutine hip_lpt_handle_elastic_wall_collisions
  end interface
#elif HAVE_CUDA
  interface
     !> CUDA kernel entry point for elastic wall-collision reflection.
     subroutine cuda_lpt_handle_elastic_wall_collisions(wall_facet_mask, &
          el_list, x_old, y_old, z_old, x, y, z, d, u, v, w, u_lag, &
          v_lag, w_lag, u_laglag, v_laglag, w_laglag, acc_xlag, &
          acc_ylag, acc_zlag, acc_xlaglag, acc_ylaglag, acc_zlaglag, &
          u_old, v_old, w_old, acc_x, acc_y, acc_z, dm_x, dm_y, dm_z, &
          nx, ny, nz, n, gdim, nelv, lx, ly, lz, lag_len, strm) &
          bind(c, name = 'cuda_lpt_handle_elastic_wall_collisions')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: wall_facet_mask, el_list
       type(c_ptr), value :: x_old, y_old, z_old
       type(c_ptr), value :: x, y, z, d, u, v, w
       type(c_ptr), value :: u_lag, v_lag, w_lag
       type(c_ptr), value :: u_laglag, v_laglag, w_laglag
       type(c_ptr), value :: acc_xlag, acc_ylag, acc_zlag
       type(c_ptr), value :: acc_xlaglag, acc_ylaglag, acc_zlaglag
       type(c_ptr), value :: u_old, v_old, w_old
       type(c_ptr), value :: acc_x, acc_y, acc_z
       type(c_ptr), value :: dm_x, dm_y, dm_z
       type(c_ptr), value :: nx, ny, nz, strm
       integer(c_int) :: n, gdim, nelv, lx, ly, lz, lag_len
     end subroutine cuda_lpt_handle_elastic_wall_collisions
  end interface
#endif

  public :: lpt_handle_elastic_wall_collisions_device

contains

  !> Launch device kernel for elastic wall-collision reflection.
  !! @param msh Mesh containing candidate wall facets.
  !! @param dm_Xh Coordinate dofmap used for wall geometry.
  !! @param coef Coefficients used for facet normals.
  !! @param wall_facet_mask_d Device mask of elastic wall facets.
  !! @param el_list_d Device list of local owner elements.
  !! @param n Number of local particles.
  subroutine lpt_handle_elastic_wall_collisions_device(msh, dm_Xh, coef, &
       wall_facet_mask_d, el_list_d, x_old, y_old, z_old, x, y, z, d, &
       u, v, w, u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, &
       acc_xlag, acc_ylag, acc_zlag, acc_xlaglag, acc_ylaglag, &
       acc_zlaglag, u_old, v_old, w_old, acc_x, acc_y, acc_z, lag_len, &
       n, strm)
    type(mesh_t), intent(in) :: msh
    type(dofmap_t), intent(in) :: dm_Xh
    type(coef_t), intent(in) :: coef
    type(c_ptr), intent(in) :: wall_facet_mask_d
    type(c_ptr), intent(in) :: el_list_d
    type(vector_t), intent(in) :: x_old, y_old, z_old
    type(vector_t), intent(inout) :: x, y, z
    type(vector_t), intent(in) :: d
    type(vector_t), intent(inout) :: u, v, w
    type(vector_t), intent(inout) :: u_lag, v_lag, w_lag
    type(vector_t), intent(inout) :: u_laglag, v_laglag, w_laglag
    type(vector_t), intent(inout) :: acc_xlag, acc_ylag, acc_zlag
    type(vector_t), intent(inout) :: acc_xlaglag, acc_ylaglag
    type(vector_t), intent(inout) :: acc_zlaglag
    type(vector_t), intent(inout) :: u_old, v_old, w_old
    type(vector_t), intent(inout) :: acc_x, acc_y, acc_z
    integer, intent(in) :: lag_len
    integer, intent(in) :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_
    integer(c_int) :: n_
    integer(c_int) :: gdim_
    integer(c_int) :: nelv_
    integer(c_int) :: lx_
    integer(c_int) :: ly_
    integer(c_int) :: lz_
    integer(c_int) :: lag_len_

    if (n .lt. 1) return

    strm_ = glb_cmd_queue
    if (present(strm)) strm_ = strm

    n_ = int(n, c_int)
    gdim_ = int(msh%gdim, c_int)
    nelv_ = int(msh%nelv, c_int)
    lx_ = int(dm_Xh%Xh%lx, c_int)
    ly_ = int(dm_Xh%Xh%ly, c_int)
    lz_ = int(dm_Xh%Xh%lz, c_int)
    lag_len_ = int(lag_len, c_int)

#ifdef HAVE_HIP
    call hip_lpt_handle_elastic_wall_collisions(wall_facet_mask_d, &
         el_list_d, x_old%x_d, y_old%x_d, z_old%x_d, x%x_d, y%x_d, z%x_d, &
         d%x_d, u%x_d, v%x_d, w%x_d, u_lag%x_d, v_lag%x_d, w_lag%x_d, &
         u_laglag%x_d, v_laglag%x_d, w_laglag%x_d, acc_xlag%x_d, &
         acc_ylag%x_d, acc_zlag%x_d, acc_xlaglag%x_d, acc_ylaglag%x_d, &
         acc_zlaglag%x_d, u_old%x_d, v_old%x_d, w_old%x_d, acc_x%x_d, &
         acc_y%x_d, acc_z%x_d, dm_Xh%x_d, dm_Xh%y_d, dm_Xh%z_d, &
         coef%nx_d, coef%ny_d, coef%nz_d, n_, gdim_, nelv_, lx_, ly_, &
         lz_, lag_len_, strm_)
#elif HAVE_CUDA
    call cuda_lpt_handle_elastic_wall_collisions(wall_facet_mask_d, &
         el_list_d, x_old%x_d, y_old%x_d, z_old%x_d, x%x_d, y%x_d, z%x_d, &
         d%x_d, u%x_d, v%x_d, w%x_d, u_lag%x_d, v_lag%x_d, w_lag%x_d, &
         u_laglag%x_d, v_laglag%x_d, w_laglag%x_d, acc_xlag%x_d, &
         acc_ylag%x_d, acc_zlag%x_d, acc_xlaglag%x_d, acc_ylaglag%x_d, &
         acc_zlaglag%x_d, u_old%x_d, v_old%x_d, w_old%x_d, acc_x%x_d, &
         acc_y%x_d, acc_z%x_d, dm_Xh%x_d, dm_Xh%y_d, dm_Xh%z_d, &
         coef%nx_d, coef%ny_d, coef%nz_d, n_, gdim_, nelv_, lx_, ly_, &
         lz_, lag_len_, strm_)
#else
    call neko_error('LPT wall collision device handling requires CUDA or HIP')
#endif
  end subroutine lpt_handle_elastic_wall_collisions_device

end module lpt_wall_collision_device
