! Copyright (c) 2020-2023, The Neko Authors
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

module legendre_rst_finder
  use num_types, only: rp, dp, xp, i8
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : HOST_TO_DEVICE, DEVICE_TO_HOST
  use space, only: space_t
  use utils, only: neko_error, neko_warning
  use vector, only: vector_t
  use matrix, only: matrix_t
  use tensor, only: tnsr3d
  use math, only: NEKO_EPS, matinv39
  use tensor_cpu, only: tnsr3d_cpu, tnsr3d_el_cpu
  use device_local_interpolation, only: device_find_rst_legendre
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, &
       c_sizeof, c_loc, c_bool
  use device, only: device_alloc, device_free, device_memcpy, &
       device_get_ptr, glb_cmd_queue
  use device_math, only: device_vlsc3
  implicit none
  private

  !> Type to compute local element (rst) coordinates
  !! for a gives set points in physical (xyz) space on a SEM grid.
  type, public :: legendre_rst_finder_t
     type(vector_t) :: x_hat, y_hat, z_hat
     type(space_t), pointer :: Xh => Null()
     integer :: nelv = 0
     real(kind=rp) :: tol = NEKO_EPS
     integer :: max_iter = 10
   contains
     procedure, pass(this) :: init => legendre_rst_finder_init
     procedure, pass(this) :: free => legendre_rst_finder_free
     procedure, pass(this) :: find => legendre_rst_finder_find
  end type legendre_rst_finder_t

contains

  subroutine legendre_rst_finder_init(this, x, y, z, nelv, Xh, tol, max_iter)
    class(legendre_rst_finder_t), intent(inout) :: this
    type(space_t), target, intent(in) :: Xh
    integer, intent(in) :: nelv
    real(kind=rp), intent(in), dimension(nelv*Xh%lxyz) :: x, y, z
    real(kind=rp), intent(in), optional :: tol
    integer, intent(in), optional :: max_iter

    call this%free()

    if (present(tol)) then
       if (tol < NEKO_EPS) then
          call neko_error('Tolerance for the Legendre finder is too small')
       end if
       this%tol = tol

    else
       this%tol = NEKO_EPS
    end if
    if (present(max_iter)) then
       if (max_iter < 1) then
          call neko_error('Max iterations for the Legendre finder is too small')
       end if
       this%max_iter = max_iter
    else
       this%max_iter = 10
    end if

    this%Xh => Xh
    this%nelv = nelv

    call this%x_hat%init(nelv*Xh%lxyz)
    call this%y_hat%init(nelv*Xh%lxyz)
    call this%z_hat%init(nelv*Xh%lxyz)

    call tnsr3d_cpu(this%x_hat%x, Xh%lx, x, &
         Xh%lx, Xh%vinv, &
         Xh%vinvt, Xh%vinvt, nelv)
    call tnsr3d_cpu(this%y_hat%x, Xh%lx, y, &
         Xh%lx, Xh%vinv, &
         Xh%vinvt, Xh%vinvt, nelv)
    call tnsr3d_cpu(this%z_hat%x, Xh%lx, z, &
         Xh%lx, Xh%vinv, &
         Xh%vinvt, Xh%vinvt, nelv)

    !> Copy the data to the device (if device exists)
    call this%x_hat%copy_from(HOST_TO_DEVICE,.false.)
    call this%y_hat%copy_from(HOST_TO_DEVICE,.false.)
    call this%z_hat%copy_from(HOST_TO_DEVICE,.false.)

  end subroutine legendre_rst_finder_init

  subroutine legendre_rst_finder_free(this)
    class(legendre_rst_finder_t), intent(inout) :: this

    call this%x_hat%free()
    call this%y_hat%free()
    call this%z_hat%free()

    this%Xh => Null()

  end subroutine legendre_rst_finder_free

!> Given a set of element candidates containing
!! the given points and computes the local RST coordinates for those points.
!! This subroutine supports both CPU and device (GPU) execution, depending
!! on the configugered backend of `NEKO_BCKND_DEVICE`.
!! @param this An instance of the `legendre_rst_finder_t` class.
!! @param rst_local_cand A matrix to store the local RST coordinates
!! of the candidate points.
!! @param x_t A vector containing the x-coordinates of the target points.
!! @param y_t A vector containing the y-coordinates of the target points.
!! @param z_t A vector containing the z-coordinates of the target points.
!! @param el_cands An array with the indices of the candidate elements.
!! @param n_point_cand The number of candidate points to process.
!! @param resx A vector to store the residuals in the x-direction.
!! @param resy A vector to store the residuals in the y-direction.
!! @param resz A vector to store the residuals in the z-direction.
  subroutine legendre_rst_finder_find(this, rst_local_cand, x_t, y_t, z_t, &
       el_cands, n_point_cand, resx, resy, resz)
    class(legendre_rst_finder_t), intent(inout) :: this
    type(matrix_t), intent(inout) :: rst_local_cand
    type(vector_t), intent(in) :: x_t, y_t, z_t
    integer, intent(in) :: n_point_cand
    integer, intent(inout) :: el_cands(:)
    type(vector_t), intent(inout) :: resx, resy, resz
    type(c_ptr) :: el_cands_d = c_null_ptr

    !> Find the elements that contain the points
    if (n_point_cand .lt. 1) return

    rst_local_cand = 0.0_rp

    if (NEKO_BCKND_DEVICE .eq. 1) then
       el_cands_d = device_get_ptr(el_cands)
       call find_rst_legendre_device(this, rst_local_cand%x_d, &
            x_t%x_d, y_t%x_d, z_t%x_d, &
            el_cands_d, n_point_cand, &
            resx%x_d, resy%x_d, resz%x_d)

    else
       call find_rst_legendre_cpu(this, rst_local_cand%x, x_t%x, y_t%x, z_t%x, &
            el_cands, n_point_cand, &
            resx%x, resy%x, resz%x)

    end if

  end subroutine legendre_rst_finder_find

  !> Using the Legendre polynomials to find the rst coordinates on GPU.
  !! @param rst holds the computed rst values
  !! @param pt_{x,y,z} are the xyz coords of the points
  !! @param el_list are the elements in which we look for rst
  !! @param n_pts the number of points
  !! @param res{x,y,z} are the difference between pt_{xyz} and xyz
  subroutine find_rst_legendre_device(this, rst, pt_x, pt_y, pt_z, &
       el_list, n_pts, resx, resy, resz)
    type(legendre_rst_finder_t), intent(inout) :: this
    type(c_ptr), intent(inout) :: rst
    type(c_ptr), intent(in) :: pt_x, pt_y, pt_z
    type(c_ptr), intent(in) :: el_list
    type(c_ptr), intent(inout) :: resx, resy, resz
    integer, intent(in) :: n_pts
    type(vector_t) :: conv_pts
    integer :: i, iter
    logical :: converged
    real(kind=rp) :: conv_sum

    if (n_pts .eq. 0) return

    call conv_pts%init(n_pts)

    conv_pts = 1.0_rp

    iter = 0
    converged = .false.
    !Iterate until found, not heavily optimized
    do while (.not. converged)
       call device_find_rst_legendre(rst, pt_x, pt_y, pt_z, &
            this%x_hat%x_d, this%y_hat%x_d, this%z_hat%x_d, &
            resx, resy, resz, &
            this%Xh%lx,el_list, n_pts, this%tol, &
            conv_pts%x_d)
       !This can be made more approriate... avoid memcpy at least
       conv_sum = device_vlsc3(conv_pts%x_d,conv_pts%x_d,conv_pts%x_d,n_pts)
       converged = conv_sum .lt. 0.5
       print *, conv_sum
       if( iter .ge. this%max_iter) converged = .true.
    end do

    call conv_pts%free()
  end subroutine find_rst_legendre_device

  !> Using the Legendre polynomials to find the rst coordinates.
  !! @param rst holds the computed rst values
  !! @param pt_{x,y,z} are the xyz coords of the points
  !! @param el_list are the elements in which we look for rst
  !! @param n_pts the number of points
  !! @param res{x,y,z} are the difference between pt_{xyz} and xyz
  subroutine find_rst_legendre_cpu(this, rst, pt_x, pt_y, pt_z, &
       el_list, n_pts, resx, resy, resz)
    type(legendre_rst_finder_t), intent(inout) :: this
    integer, intent(in) :: n_pts
    real(kind=rp), intent(inout) :: rst(3, n_pts)
    real(kind=rp), intent(in) :: pt_x(n_pts)
    real(kind=rp), intent(in) :: pt_y(n_pts)
    real(kind=rp), intent(in) :: pt_z(n_pts)
    real(kind=rp), intent(inout) :: resx(n_pts)
    real(kind=rp), intent(inout) :: resy(n_pts)
    real(kind=rp), intent(inout) :: resz(n_pts)
    integer, intent(in) :: el_list(n_pts)
    real(kind=rp) :: r_legendre(this%Xh%lx)
    real(kind=rp) :: s_legendre(this%Xh%lx)
    real(kind=rp) :: t_legendre(this%Xh%lx)
    real(kind=rp) :: dr_legendre(this%Xh%lx)
    real(kind=rp) :: ds_legendre(this%Xh%lx)
    real(kind=rp) :: dt_legendre(this%Xh%lx)
    real(kind=rp) :: jac(3,3)
    real(kind=xp) :: tmp(this%Xh%lx), tmp2(this%Xh%lx)
    real(kind=xp) :: rst_d(3), jacinv(3,3)
    integer :: conv_pts
    logical :: converged
    integer :: i, j, e, iter, lx, lx2



    lx = this%Xh%lx
    if (n_pts .lt. 1) return

    rst = 0.0_rp
    ! If performance critical we should do multiple points at the time
    ! Currently we do one point at the time
    do i = 1, n_pts
       iter = 0
       converged = .false.
       do while (.not. converged)
          iter = iter + 1
          ! Compute legendre polynomials in this rst coordinate
          r_legendre(1) = 1.0
          r_legendre(2) = rst(1,i)
          s_legendre(1) = 1.0
          s_legendre(2) = rst(2,i)
          t_legendre(1) = 1.0
          t_legendre(2) = rst(3,i)
          dr_legendre(1) = 0.0
          dr_legendre(2) = 1.0
          ds_legendre(1) = 0.0
          ds_legendre(2) = 1.0
          dt_legendre(1) = 0.0
          dt_legendre(2) = 1.0
          do j = 2, lx-1
             r_legendre(j+1) = ((2.0_xp*(j-1.0_xp)+1.0_xp) * rst(1,i) &
                  * r_legendre(j) - (j-1.0_xp) &
                  * r_legendre(j-1)) / (real(j,xp))
             s_legendre(j+1) = ((2.0_xp*(j-1.0_xp)+1.0_xp) * rst(2,i) &
                  * s_legendre(j) - (j-1.0_xp) &
                  * s_legendre(j-1))/(real(j,xp))
             t_legendre(j+1) = ((2.0_xp*(j-1.0_xp)+1.0_xp) * rst(3,i) &
                  * t_legendre(j) - (j-1.0_xp) &
                  * t_legendre(j-1))/(real(j,xp))
             dr_legendre(j+1) = ((j-1.0_xp)+1.0_xp) * r_legendre(j) &
                  + rst(1,i)*dr_legendre(j)
             ds_legendre(j+1) = ((j-1.0_xp)+1.0_xp) * s_legendre(j) &
                  + rst(2,i)*ds_legendre(j)
             dt_legendre(j+1) = ((j-1.0_xp)+1.0_xp) * t_legendre(j) &
                  + rst(3,i)*dt_legendre(j)
          end do
          e = (el_list(i))*this%Xh%lxyz + 1
          ! Compute the current xyz value
          call tnsr3d_el_cpu(resx(i), 1, this%x_hat%x(e), lx, &
               r_legendre, s_legendre, t_legendre)
          call tnsr3d_el_cpu(resy(i), 1, this%y_hat%x(e), lx, &
               r_legendre, s_legendre, t_legendre)
          call tnsr3d_el_cpu(resz(i), 1, this%z_hat%x(e), lx, &
               r_legendre, s_legendre, t_legendre)
          ! This should in principle be merged into some larger kernel
          ! Compute the jacobian
          call tnsr3d_el_cpu(jac(1,1), 1, this%x_hat%x(e), lx, &
               dr_legendre, s_legendre(1), t_legendre)
          call tnsr3d_el_cpu(jac(1,2), 1, this%y_hat%x(e), lx, &
               dr_legendre, s_legendre(1), t_legendre)
          call tnsr3d_el_cpu(jac(1,3), 1, this%z_hat%x(e), lx, &
               dr_legendre, s_legendre, t_legendre)
          call tnsr3d_el_cpu(jac(2,1), 1, this%x_hat%x(e), lx, &
               r_legendre, ds_legendre, t_legendre)
          call tnsr3d_el_cpu(jac(2,2), 1, this%y_hat%x(e), lx, &
               r_legendre, ds_legendre, t_legendre)
          call tnsr3d_el_cpu(jac(2,3), 1, this%z_hat%x(e), lx, &
               r_legendre, ds_legendre, t_legendre)
          call tnsr3d_el_cpu(jac(3,1), 1, this%x_hat%x(e), lx, &
               r_legendre, s_legendre, dt_legendre)
          call tnsr3d_el_cpu(jac(3,2), 1, this%y_hat%x(e), lx, &
               r_legendre, s_legendre, dt_legendre)
          call tnsr3d_el_cpu(jac(3,3), 1, this%z_hat%x(e), lx, &
               r_legendre, s_legendre, dt_legendre)
          resx(i) = pt_x(i) - resx(i)
          resy(i) = pt_y(i) - resy(i)
          resz(i) = pt_z(i) - resz(i)
          ! Jacobian inverse
          jacinv = matinv39(jac(1,1), jac(1,2), jac(1,3),&
               jac(2,1), jac(2,2), jac(2,3),&
               jac(3,1), jac(3,2), jac(3,3))
          ! Update direction
          rst_d(1) = (resx(i)*jacinv(1,1) &
               + jacinv(2,1)*resy(i) &
               + jacinv(3,1)*resz(i))
          rst_d(2) = (resx(i)*jacinv(1,2) &
               + jacinv(2,2)*resy(i) &
               + jacinv(3,2)*resz(i))
          rst_d(3) = (resx(i)*jacinv(1,3) &
               + jacinv(2,3)*resy(i) &
               + jacinv(3,3)*resz(i))

          conv_pts = 0
          if (norm2(real(rst_d,xp)) .le. this%tol) then
             conv_pts = 1
          end if
          if (norm2(real(rst_d,xp)) .gt. 4.0) then
             conv_pts = 1
          end if
          ! Update rst coordinates
          rst(1,i) = rst(1,i) + rst_d(1)
          rst(2,i) = rst(2,i) + rst_d(2)
          rst(3,i) = rst(3,i) + rst_d(3)

          converged = conv_pts .eq. 1
          if (iter .ge. this%max_iter) converged = .true.
       end do
    end do
  end subroutine find_rst_legendre_cpu
end module legendre_rst_finder
