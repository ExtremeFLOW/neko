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
  use math, only: NEKO_EPS
  use tensor_cpu, only: tnsr3d_cpu, tnsr3d_el_cpu
  use device_local_interpolation, only: device_find_rst_legendre
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, &
      c_sizeof, c_loc, c_bool
  use device, only: device_alloc, device_free, device_memcpy, &
      device_get_ptr, glb_cmd_queue
  implicit none
  private


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

  type, public :: rst_finder_monitor_t


  end type rst_finder_monitor_t


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
    call this%x_hat%copyto(HOST_TO_DEVICE,.false.)
    call this%y_hat%copyto(HOST_TO_DEVICE,.false.)
    call this%z_hat%copyto(HOST_TO_DEVICE,.false.)

  end subroutine legendre_rst_finder_init

  subroutine legendre_rst_finder_free(this)
    class(legendre_rst_finder_t), intent(inout) :: this

    call this%x_hat%free()
    call this%y_hat%free()
    call this%z_hat%free()

    this%Xh => Null()

  end subroutine legendre_rst_finder_free

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
       call find_rst_legendre_device(this, rst_local_cand%x_d, x_t%x_d, y_t%x_d, &
            z_t%x_d, el_cands_d, n_point_cand, &
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
    logical(kind=c_bool) , allocatable, target :: conv_pts(:)
    type(c_ptr) :: conv_pts_d = c_null_ptr
    integer :: i, iter
    logical :: converged
    integer(kind=i8) :: bytes

    if (allocated(conv_pts)) then
       deallocate(conv_pts)
    end if

    allocate(conv_pts(n_pts))

    conv_pts = .false.
    bytes = n_pts*c_sizeof(conv_pts(1))
    call device_alloc(conv_pts_d,bytes)
    call device_memcpy(c_loc(conv_pts), conv_pts_d, bytes, HOST_TO_DEVICE, &
    .false., glb_cmd_queue)
    iter = 0
    converged = .false.
    !Iterate until found, not heavily optimized
    do while (.not. converged)
       call device_find_rst_legendre(rst, pt_x, pt_y, pt_z, &
            this%x_hat%x_d, this%y_hat%x_d, this%z_hat%x_d, &
            resx, resy, resz, &
            this%Xh%lx,el_list, n_pts, this%tol, &
            conv_pts_d)
       call device_memcpy(c_loc(conv_pts), conv_pts_d, bytes, DEVICE_TO_HOST, &
             .true., glb_cmd_queue)
       converged = .true.
       iter = iter + 1
       do i = 1, n_pts
          converged = converged .and. conv_pts(i)
       end do
       if( iter .ge. this%max_iter) converged = .true.
    end do

    call device_free(conv_pts_d)
    deallocate(conv_pts)

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
    real(kind=rp), intent(in)  :: pt_x(n_pts)
    real(kind=rp), intent(in)  :: pt_y(n_pts)
    real(kind=rp), intent(in)  :: pt_z(n_pts)
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
    do i = 1, n_pts
       iter = 0
       converged = .false.
       do while (.not. converged)
          iter  = iter + 1
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
             r_legendre(j+1) = ((2.0_xp*(j-1.0_xp)+1.0_xp)*rst(1,i)*r_legendre(j) - (j-1.0_xp)*r_legendre(j-1))/(real(j,xp))
             s_legendre(j+1) = ((2.0_xp*(j-1.0_xp)+1.0_xp)*rst(2,i)*s_legendre(j) - (j-1.0_xp)*s_legendre(j-1))/(real(j,xp))
             t_legendre(j+1) = ((2.0_xp*(j-1.0_xp)+1.0_xp)*rst(3,i)*t_legendre(j) - (j-1.0_xp)*t_legendre(j-1))/(real(j,xp))
             dr_legendre(j+1) = ((j-1.0_xp)+1.0_xp) * r_legendre(j) + rst(1,i)*dr_legendre(j)
             ds_legendre(j+1) = ((j-1.0_xp)+1.0_xp) * s_legendre(j) + rst(2,i)*ds_legendre(j)
             dt_legendre(j+1) = ((j-1.0_xp)+1.0_xp) * t_legendre(j) + rst(3,i)*dt_legendre(j)
          end do
          e = (el_list(i))*this%Xh%lxyz + 1
          call tnsr3d_el_cpu(resx(i), 1, this%x_hat%x(e), lx, &
               r_legendre, s_legendre, t_legendre)
          call tnsr3d_el_cpu(resy(i), 1, this%y_hat%x(e), lx, &
               r_legendre, s_legendre, t_legendre)
          call tnsr3d_el_cpu(resz(i), 1, this%z_hat%x(e), lx, &
               r_legendre, s_legendre, t_legendre)

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
          jacinv = matinv39(jac(1,1),jac(1,2),jac(1,3),&
               jac(2,1),jac(2,2),jac(2,3),&
               jac(3,1),jac(3,2),jac(3,3))
          rst_d(1) = (resx(i)*jacinv(1,1)+jacinv(2,1)*resy(i)+jacinv(3,1)*resz(i))
          rst_d(2) = (resx(i)*jacinv(1,2)+jacinv(2,2)*resy(i)+jacinv(3,2)*resz(i))
          rst_d(3) = (resx(i)*jacinv(1,3)+jacinv(2,3)*resy(i)+jacinv(3,3)*resz(i))

          conv_pts = 0
          if (norm2(real(rst_d,xp)) .le. this%tol)then
             conv_pts = 1
          end if
          if (norm2(real(rst_d,xp)) .gt. 4.0)then
             conv_pts = 1
          end if
          rst(1,i) = rst(1,i) + rst_d(1)
          rst(2,i) = rst(2,i) + rst_d(2)
          rst(3,i) = rst(3,i) + rst_d(3)

          converged = conv_pts .eq. 1
          if (iter .ge. this%max_iter) converged = .true.
       end do
    end do
  end subroutine find_rst_legendre_cpu

!> Compares two sets of rst coordinates and checks whether rst2 is better than rst1 given a tolerance
!> res1 and res2 are the distances to the interpolated xyz coordinate and true xyz coord
  function rst_cmp(rst1, rst2,res1, res2, tol) result(rst2_better)
    real(kind=rp) :: rst1(3), res1(3)
    real(kind=rp) :: rst2(3), res2(3)
    real(kind=rp) :: tol
    logical :: rst2_better
!If rst1 is invalid and rst2 is valid, take rst2
! If both invalidl, take smallest residual
    if (abs(rst1(1)) .gt. 1.0_xp+tol .or. &
         abs(rst1(2)) .gt. 1.0_xp+tol .or. &
         abs(rst1(3)) .gt. 1.0_xp+tol) then
       if (abs(rst2(1)) .le. 1.0_xp+tol .and. &
            abs(rst2(2)) .le. 1.0_xp+tol .and. &
            abs(rst2(3)) .le. 1.0_xp+tol) then
          rst2_better = .true.
       else
          rst2_better = (norm2(real(res2,xp)) .lt. norm2(real(res1,xp)))
       end if
    else
!> Else we check rst2 is inside and has a smaller distance
       rst2_better = (norm2(real(res2,xp)) .lt. norm2(real(res1,xp)) .and.&
            abs(rst2(1)) .le. 1.0_xp+tol .and. &
            abs(rst2(2)) .le. 1.0_xp+tol .and. &
            abs(rst2(3)) .le. 1.0_xp+tol)
    end if
  end function rst_cmp

  subroutine jacobian_inverse(jacinv, rst, x, y, z, n_pts, Xh)
    integer :: n_pts
    real(kind=rp), intent(inout) :: rst(3, n_pts)
    type(space_t), intent(inout) :: Xh
    real(kind=rp), intent(inout) :: x(Xh%lx, Xh%ly, Xh%lz, n_pts)
    real(kind=rp), intent(inout) :: y(Xh%lx, Xh%ly, Xh%lz, n_pts)
    real(kind=rp), intent(inout) :: z(Xh%lx, Xh%ly, Xh%lz, n_pts)
    real(kind=rp), intent(out) :: jacinv(3,3, n_pts)
    real(kind=rp) :: tmp(3,3)
    integer :: i

    call jacobian(jacinv, rst, x, y, z, n_pts, Xh)

    do i = 1, n_pts
       tmp = matinv3(real(jacinv(:,:,3),xp))
       jacinv(:,:,i) = tmp
    end do

  end subroutine jacobian_inverse
  ! M33INV and M44INV by David G. Simpson pure function version from
  ! https://fortranwiki.org/fortran/show/Matrix+inversion
  ! Invert 3x3 matrix
  function matinv39(a11, a12, a13, a21, a22, a23, a31, a32, a33) &
      result(B)
    real(kind=rp), intent(in) :: a11, a12, a13, a21, a22, a23, a31, a32, a33
    real(xp) :: A(3,3)   !! Matrix
    real(rp) :: B(3,3)   !! Inverse matrix
    A(1,1) = a11
    A(1,2) = a12
    A(1,3) = a13
    A(2,1) = a21
    A(2,2) = a22
    A(2,3) = a23
    A(3,1) = a31
    A(3,2) = a32
    A(3,3) = a33
    B = matinv3(A)
  end function matinv39
  !! Performs a direct calculation of the inverse of a 3×3 matrix.
  !M33INV and M44INV by David G. Simpson pure function version from https://fortranwiki.org/fortran/show/Matrix+inversion
  ! Invert 3x3 matrix
  function matinv3(A) result(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(xp), intent(in) :: A(3,3)   !! Matrix
    real(xp) :: B(3,3)   !! Inverse matrix
    real(xp) :: detinv

    ! Calculate the inverse determinant of the matrix
    !first index x,y,z, second r, s, t
    detinv = 1.0_xp/real(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1),xp)
    !print *, "detinv",detinv
    ! Calculate the inverse of the matrix
    ! first index r, s, t, second x, y, z
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  end function matinv3


end module legendre_rst_finder
