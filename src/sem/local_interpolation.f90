! Copyright (c) 2021-2023, The Neko Authors
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
!> Routines to obtain interpolated values on a set of points with known
!! rst coordinates in elements local to this process.
module local_interpolation
  use tensor, only: triple_tensor_product, tnsr3d_el_list, tnsr3d
  use space, only: space_t, GL, GLL
  use num_types, only: rp, xp
  use point, only: point_t
  use math, only: abscmp, NEKO_EPS
  use speclib
  use fast3d, only: fd_weights_full, setup_intp
  use utils, only: neko_error
  use field, only: field_t
  use field_list, only: field_list_t
  use device
  use tensor_cpu
  use device_math, only: device_rzero
  use neko_config, only: NEKO_BCKND_DEVICE
  implicit none

  !> Interpolation on a set of points with known rst coordinates in elements local
  !! to this process.
  !! Similar to point_interpolator, but prioritizes performance
  !! Only works with arrays of coordinates
  !! Performs interpolation with the configured NEKO_BCKND
  type :: local_interpolator_t
     !> First space.
     type(space_t), pointer :: Xh => null()
     !> Number of points to interpolate on
     integer :: n_points
     !> Weights for local interpolation
     real(kind=rp), allocatable :: weights_r(:,:)
     real(kind=rp), allocatable :: weights_s(:,:)
     real(kind=rp), allocatable :: weights_t(:,:)
     type(c_ptr) :: weights_r_d = c_null_ptr
     type(c_ptr) :: weights_s_d = c_null_ptr
     type(c_ptr) :: weights_t_d = c_null_ptr
   contains
     !> Constructor.
     procedure, pass(this) :: init => local_interpolator_init
     !> Destructor.
     procedure, pass(this) :: free => local_interpolator_free
     !> Interpolates the scalar field \f$ X \f$ on the specified coordinates
     procedure, pass(this) :: evaluate => local_interpolator_evaluate
     !> COmputes weights based on rst coordinates
     procedure, pass(this) :: compute_weights => local_interpolator_compute_weights

  end type local_interpolator_t
   
  public :: find_rst_legendre

contains

  !> Initialization of point interpolation.
  !! @param xh Function space.
  subroutine local_interpolator_init(this, Xh, r, s, t, n_points)
    class(local_interpolator_t), intent(inout), target :: this
    type(space_t), intent(in), target :: Xh
    integer, intent(in) :: n_points
    real(kind=rp) :: r(n_points), s(n_points), t(n_points)
    integer :: size_weights
    call this%free()
    if ((Xh%t .eq. GL) .or. (Xh%t .eq. GLL)) then
    else
       call neko_error('Unsupported interpolation')
    end if

    this%Xh => Xh
    this%n_points = n_points
    allocate(this%weights_r(Xh%lx,n_points))
    allocate(this%weights_s(Xh%ly,n_points))
    allocate(this%weights_t(Xh%lz,n_points))
    call this%compute_weights(r, s, t)
    size_weights = Xh%lx * n_points

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%weights_r, this%weights_r_d, size_weights)
       call device_map(this%weights_s, this%weights_s_d, size_weights)
       call device_map(this%weights_t, this%weights_t_d, size_weights)
       call device_memcpy(this%weights_r, this%weights_r_d,&
            size_weights, HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%weights_s, this%weights_s_d,&
            size_weights, HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%weights_t, this%weights_t_d,&
            size_weights, HOST_TO_DEVICE, sync = .true.)
    end if

  end subroutine local_interpolator_init

  !> Free pointers
  subroutine local_interpolator_free(this)
    class(local_interpolator_t), intent(inout) :: this

    if (associated(this%Xh)) this%Xh => null()

    if(allocated(this%weights_r)) deallocate(this%weights_r)
    if(allocated(this%weights_s)) deallocate(this%weights_s)
    if(allocated(this%weights_t)) deallocate(this%weights_t)
    if (c_associated(this%weights_r_d)) then
       call device_free(this%weights_r_d)
    end if
    if (c_associated(this%weights_s_d)) then
       call device_free(this%weights_s_d)
    end if
    if (c_associated(this%weights_t_d)) then
       call device_free(this%weights_t_d)
    end if
  end subroutine local_interpolator_free

  !> Computes interpolation weights \f$ w_r, w_s, w_t \f$ for a
  !! list of points.
  !! @param r local r-coordinates.
  !! @param s local s-coordinates.
  !! @param t local t-coordinates.
  !! @param wr Weights in the r-direction.
  !! @param ws Weights in the s-direction.
  !! @param wt Weights in the t-direction.
  !! @note `wr`, `ws` and `wt` must be arrays of dimensions `(lx, N)` where `N`
  !! is the number of points (size of the `r`,`s`,`t` arrays).
  subroutine local_interpolator_compute_weights(this, r, s, t)
    class(local_interpolator_t), intent(inout) :: this
    real(kind=rp), intent(in) :: r(:), s(:), t(:)

    integer :: N, i, lx
    lx = this%Xh%lx
    N = size(r)

    do i = 1, N
       call fd_weights_full(r(i), this%Xh%zg(:,1), lx-1, 0, this%weights_r(:,i))
       call fd_weights_full(s(i), this%Xh%zg(:,2), lx-1, 0, this%weights_s(:,i))
       call fd_weights_full(t(i), this%Xh%zg(:,3), lx-1, 0, this%weights_t(:,i))
    end do

  end subroutine local_interpolator_compute_weights

  !> Interpolates a list of fields based on a set of element ids.
  !! @param rst r,s,t coordinates.
  !! @param el_owners Array of element ids that "own" a given point `i`.
  !! @param sampled_fields_list A list of fields to interpolate.
  !! @param wr Weights in the r-direction of shape `(lx, N)` where `N` is the
  !! number of points to interpolate.
  !! @param ws Weights in the s-direction of shape `(lx, N)` where `N` is the
  !! number of points to interpolate.
  !! @param wt Weights in the t-direction of shape `(lx, N)` where `N` is the
  !! number of points to interpolate.
  !! @note The weights can be generated with the subroutine `compute_weights`.
  !! Assumes weights have been computed for these points.
  subroutine local_interpolator_evaluate(this, interp_values, el_list, field,nel, on_host)
    class(local_interpolator_t), intent(inout) :: this
    integer, intent(in) :: el_list(this%n_points)
    integer, intent(in) :: nel
    real(kind=rp), intent(inout) :: interp_values(this%n_points)
    real(kind=rp), intent(inout) :: field(this%Xh%lxyz, nel)
    logical, intent(in) :: on_host

    call tnsr3d_el_list(interp_values, 1, field, this%Xh%lx, &
         this%weights_r, this%weights_s, this%weights_t, el_list, this%n_points, on_host)

  end subroutine local_interpolator_evaluate

  !> Constructs the Jacobian, returns a 3-by-3 times number of points where
  !! \f$ [J(\mathbf{r}]_{ij} = \frac{d\mathbf{x}_i}{d\mathbf{r}_j}\f$.
  !! @param rst r,s,t coordinates.
  !! @param X Values of the field \f$ X \f$ at points.
  !! @param Y Values of the field \f$ Y \f$ at points.
  !! @param Z Values of the field \f$ Z \f$ at points.
  subroutine jacobian(jac, rst, x, y, z, n_pts, Xh)
    integer :: n_pts
    real(kind=rp), intent(inout) :: rst(3, n_pts)
    type(space_t), intent(inout) :: Xh
    real(kind=rp), intent(inout) :: x(Xh%lx, Xh%ly, Xh%lz, n_pts)
    real(kind=rp), intent(inout) :: y(Xh%lx, Xh%ly, Xh%lz, n_pts)
    real(kind=rp), intent(inout) :: z(Xh%lx, Xh%ly, Xh%lz, n_pts)
    real(kind=rp), intent(out) :: jac(3,3, n_pts)

    real(kind=rp) :: tmp(3)

    real(kind=rp) :: hr(Xh%lx, 2), hs(Xh%ly, 2), ht(Xh%lz, 2)
    integer :: lx, ly, lz, i
    lx = Xh%lx
    ly = Xh%ly
    lz = Xh%lz
    
    do i = 1, n_pts 
       ! Weights
       call fd_weights_full(rst(1,i), Xh%zg(:,1), lx-1, 1, hr)
       call fd_weights_full(rst(2,i), Xh%zg(:,2), ly-1, 1, hs)
       call fd_weights_full(rst(3,i), Xh%zg(:,3), lz-1, 1, ht)

       ! d(x,y,z)/dr
       call triple_tensor_product(tmp(1), X(:,:,:,i), lx, hr(:,2), hs(:,1), ht(:,1))
       jac(1,1,i) = tmp(1)
       call triple_tensor_product(tmp(1), Y(:,:,:,i), lx, hr(:,2), hs(:,1), ht(:,1))
       jac(1,2,i) = tmp(1)
       call triple_tensor_product(tmp(1), Z(:,:,:,i), lx, hr(:,2), hs(:,1), ht(:,1))
       jac(1,3,i) = tmp(1)

       ! d(x,y,z)/ds
       call triple_tensor_product(tmp, X(:,:,:,i), Y(:,:,:,i), Z(:,:,:,i), lx, hr(:,1), hs(:,2), ht(:,1))
       jac(2,:,i) = tmp

       ! d(x,y,z)/dt
       call triple_tensor_product(tmp, X(:,:,:,i), Y(:,:,:,i), Z(:,:,:,i), lx, hr(:,1), hs(:,1), ht(:,2))
       jac(3,:,i) = tmp
    end do

  end subroutine jacobian

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
       tmp = matinv3(jacinv(:,:,3))
       jacinv(:,:,i) = tmp
    end do

  end subroutine jacobian_inverse
 !M33INV and M44INV by David G. Simpson pure function version from https://fortranwiki.org/fortran/show/Matrix+inversion
 ! Invert 3x3 matrix
  pure function matinv39(a11,a12,a13,a21,a22,a23,a31,a32,a33) result(B)
    real(rp), intent(in) :: a11, a12, a13,a21, a22, a23,a31, a32, a33
    real(rp) :: A(3,3)   !! Matrix
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
  pure function matinv3(A) result(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(rp), intent(in) :: A(3,3)   !! Matrix
    real(rp) :: B(3,3)   !! Inverse matrix
    real(xp) :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.0_xp/real(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1),xp)

    ! Calculate the inverse of the matrix
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

!  !> Using the Legendre-Perez method to find the rst coordinates.
!  subroutine find_rst_legendre_bad(rst, pt_x, pt_y, pt_z, Xh, x, y, z, &
!                                   el_list, n_pts, nelv, &
!                                   resx, resy, resz, tol)
!    integer, intent(in) :: n_pts, nelv
!    real(kind=rp), intent(inout) :: rst(3, n_pts)
!    type(space_t), intent(inout) :: Xh
!    real(kind=rp), intent(inout) :: x(Xh%lx, Xh%ly, Xh%lz, nelv)
!    real(kind=rp), intent(inout) :: y(Xh%lx, Xh%ly, Xh%lz, nelv)
!    real(kind=rp), intent(inout) :: z(Xh%lx, Xh%ly, Xh%lz, nelv)
!    real(kind=rp), intent(inout) :: resx(n_pts)
!    real(kind=rp), intent(inout) :: resy(n_pts)
!    real(kind=rp), intent(inout) :: resz(n_pts)
!    integer, intent(in) :: el_list(n_pts)
!    real(kind=rp), intent(in) :: tol 
!    real(kind=rp), allocatable :: x_hat(:, :, :, :)
!    real(kind=rp), allocatable :: y_hat(:, :, :, :)
!    real(kind=rp), allocatable :: z_hat(:, :, :, :)
!    real(kind=rp), allocatable :: pt_x(:)
!    real(kind=rp), allocatable :: pt_y(:)
!    real(kind=rp), allocatable :: pt_z(:)
!    real(kind=rp), allocatable :: xt(:)
!    real(kind=rp), allocatable :: yt(:)
!    real(kind=rp), allocatable :: zt(:)
!    real(kind=rp), allocatable :: r_legendre(:,:) 
!    real(kind=rp), allocatable :: s_legendre(:,:) 
!    real(kind=rp), allocatable :: t_legendre(:,:) 
!    real(kind=rp), allocatable :: dr_legendre(:,:) 
!    real(kind=rp), allocatable :: ds_legendre(:,:) 
!    real(kind=rp), allocatable :: dt_legendre(:,:) 
!    real(kind=rp), allocatable :: jac11(:), jac12(:), jac13(:)
!    real(kind=rp), allocatable :: jac21(:), jac22(:), jac23(:)
!    real(kind=rp), allocatable :: jac31(:), jac32(:), jac33(:)
!    real(kind=rp) :: jacinv(3,3), tmp(Xh%lx), tmp2(Xh%lx), rst_d(3), avgx, avgy, avgz
!    integer, allocatable :: conv_pts(:)
!    logical :: converged
!    integer :: i, j, k, iter, lx, n_pts1
!
!    allocate(x_hat(xh%lx, xh%ly, xh%lz, nelv))
!    allocate(y_hat(xh%lx, xh%ly, xh%lz, nelv))
!    allocate(z_hat(xh%lx, xh%ly, xh%lz, nelv))
!    allocate(r_legendre(xh%lx, n_pts))
!    allocate(s_legendre(xh%lx, n_pts))
!    allocate(t_legendre(xh%lx, n_pts))
!    allocate(dr_legendre(xh%lx, n_pts))
!    allocate(ds_legendre(xh%lx, n_pts))
!    allocate(dt_legendre(xh%lx, n_pts))
!    allocate(xt(n_pts))
!    allocate(yt(n_pts))
!    allocate(zt(n_pts))
!    allocate(conv_pts(n_pts))
!    allocate(jac11(n_pts), jac12(n_pts), jac13(n_pts))
!    allocate(jac21(n_pts), jac22(n_pts), jac23(n_pts))
!    allocate(jac31(n_pts), jac32(n_pts), jac33(n_pts))
!    
!    lx = Xh%lx
!    if (n_pts .lt. 1) return
!    
!    !> Transform into legendre space
!    !> hats are in Legendre space
!    call tnsr3d(x_hat, Xh%lx, x, &
!                Xh%lx, Xh%vinv, &
!                Xh%vinvt, Xh%vinvt, nelv)
!
!    call tnsr3d(y_hat, Xh%lx, x_hat, &
!                Xh%lx, Xh%v, &
!                Xh%vt, Xh%vt, nelv)
!    call tnsr3d(y_hat, Xh%lx, y, &
!                Xh%lx, Xh%vinv, &
!                Xh%vinvt, Xh%vinvt, nelv)
!    call tnsr3d(z_hat, Xh%lx, z, &
!                Xh%lx, Xh%vinv, &
!                Xh%vinvt, Xh%vinvt, nelv)
!
!    do i = 1, n_pts
!    end do
!    n_pts1 = n_pts
! 
!    do k = 1, n_pts1
!       i = k 
!       avgx = 0.0
!       avgy = 0.0
!       avgz = 0.0
!       do j = 1, lx*lx*lx
!          avgx = avgx + x(j,1,1,el_list(i)+1)
!          avgy = avgy + y(j,1,1,el_list(i)+1)
!          avgz = avgz + z(j,1,1,el_list(i)+1)
!       end do
!       rst(1,i) = 1.7/(lx*lx*lx)*(pt_x(i)*(lx*lx*lx) - avgx)/&
!  (maxval(x(:,:,:,el_list(i)+1))-minval(x(:,:,:,el_list(i)+1)))
!       rst(2,i) = 1.7/(lx*lx*lx)*(pt_y(i)*(lx*lx*lx) - avgy)/&
!       (maxval(y(:,:,:,el_list(i)+1))-minval(y(:,:,:,el_list(i)+1)))
!       rst(3,i) = 1.7/(lx*lx*lx)*(pt_z(i)*(lx*lx*lx) - avgz)/&
!(maxval(z(:,:,:,el_list(i)+1))-minval(z(:,:,:,el_list(i)+1)))
!       
!       call legendre_poly(r_legendre(1,1),rst(1,i),lx-1)     
!       call legendre_poly(s_legendre(1,1),rst(2,i),lx-1)     
!       call legendre_poly(t_legendre(1,1),rst(3,i),lx-1)     
!       do j = 0, lx-1
!          dr_legendre(j+1,1) = PNDLEG(rst(1,i),j)
!          ds_legendre(j+1,1) = PNDLEG(rst(2,i),j)
!          dt_legendre(j+1,1) = PNDLEG(rst(3,i),j)
!       end do
!    iter = 0 
!    converged = .false.
!    do while (.not. converged)
!       iter  = iter + 1
!
!       !xyz = local_interpolator_evaluate(rst)
!       !Compute xyz based on coord
!       call tnsr3d_el_list(xt, 1, x_hat, Xh%lx, &
!            r_legendre, s_legendre, t_legendre, el_list(k), 1)
!       call tnsr3d_el_list(yt, 1, y_hat, Xh%lx, &
!            r_legendre, s_legendre, t_legendre, el_list(k), 1)
!       call tnsr3d_el_list(zt, 1, z_hat, Xh%lx, &
!            r_legendre, s_legendre, t_legendre, el_list(k), 1)
!       !check where to go
!       xt(1) = pt_x(k) - xt(1)
!       yt(1) = pt_y(k) - yt(1)
!       zt(1) = pt_z(k) - zt(1)
!       
!       ! compute jacobian
!       call tnsr3d_el_list(jac11, 1, x_hat, Xh%lx, &
!            dr_legendre, s_legendre, t_legendre, el_list(k),1)
!       call tnsr3d_el_list(jac12, 1, y_hat, Xh%lx, &
!            dr_legendre, s_legendre, t_legendre, el_list(k),1)
!       call tnsr3d_el_list(jac13, 1, z_hat, Xh%lx, &
!            dr_legendre, s_legendre, t_legendre, el_list(k),1)
!       call tnsr3d_el_list(jac21, 1, x_hat, Xh%lx, &
!            r_legendre, ds_legendre, t_legendre, el_list(k),1)
!       call tnsr3d_el_list(jac22, 1, y_hat, Xh%lx, &
!            r_legendre, ds_legendre, t_legendre, el_list(k),1)
!       call tnsr3d_el_list(jac23, 1, z_hat, Xh%lx, &
!            r_legendre, ds_legendre, t_legendre, el_list(k),1)
!       call tnsr3d_el_list(jac31, 1, x_hat, Xh%lx, &
!            r_legendre, s_legendre, dt_legendre, el_list(k),1)
!       call tnsr3d_el_list(jac32, 1, y_hat, Xh%lx, &
!            r_legendre, s_legendre, dt_legendre, el_list(k),1)
!       call tnsr3d_el_list(jac33, 1, z_hat, Xh%lx, &
!            r_legendre, s_legendre, dt_legendre, el_list(k),1)
!       !do jac inverse
!       conv_pts = 1
!          jacinv = matinv39(jac11(1),jac12(1),jac13(1),&
!                            jac21(1),jac22(1),jac23(1),&
!                            jac31(1),jac32(1),jac33(1))
!          rst_d(1) = (xt(1)*jacinv(1,1)+jacinv(2,1)*yt(1)+jacinv(3,1)*zt(1))
!          rst_d(2) = (xt(1)*jacinv(1,2)+jacinv(2,2)*yt(1)+jacinv(3,2)*zt(1))
!          rst_d(3) = (xt(1)*jacinv(1,3)+jacinv(2,3)*yt(1)+jacinv(3,3)*zt(1))
!
!          rst(1,k) = rst(1,k) + rst_d(1)
!          rst(2,k) = rst(2,k) + rst_d(2)
!          rst(3,k) = rst(3,k) + rst_d(3)
!          if ( (rst_d(1) > rst(1,k)*tol .and. &
!               rst_d(2) >  rst(2,k)*tol .and. &
!               rst_d(3) >  rst(3,k)*tol ) .and. &
!               norm2((/xt(1),yt(1),zt(1)/)) > tol) then
!               if (norm2((/xt(1),yt(1),zt(1)/)) < 10) then
!                  if (norm2((/xt(1),yt(1),zt(1)/)) < 10) then
!                     if( abs(rst(1,k)) .lt. 1.0+tol .and. &
!                         abs(rst(2,k)) .lt. 1.0+tol .and. &
!                         abs(rst(3,k)) .lt. 1.0+tol) then
!                         conv_pts(1) = 0
!                    end if
!                 end if   
!             end if
!          end if
!           
!       !call jacobian_inverse(jacinv,xyz)
!       !rst = rst -(jacinv*dir)
!       !converged = (maxval(abs(xt)) < NEKO_EPS*10) .and. &
!       !            (maxval(abs(yt)) < NEKO_EPS*10) .and. &
!       !            (maxval(abs(zt)) < NEKO_EPS*10) 
!       resx(k) = xt(1)
!       resy(k) = yt(1)
!       resz(k) = zt(1)
!       converged = conv_pts(1) .eq. 1
!       !if (converged) print *, 'are we true?',converged
!       !if (converged) print *,'maxabsnorm', sum(conv_pts), n_pts,'iter', iter
!       ! Need some kind of reasonable convergence criterion  
!
!       if (iter .gt. 50) exit
!       i = 1
!          call legendre_poly(r_legendre(1,i),rst(1,k),lx-1)     
!          call legendre_poly(s_legendre(1,i),rst(2,k),lx-1)     
!          call legendre_poly(t_legendre(1,i),rst(3,k),lx-1)     
!          do j = 0, lx-1
!             dr_legendre(j+1,i) = PNDLEG(rst(1,k),j)
!             ds_legendre(j+1,i) = PNDLEG(rst(2,k),j)
!             dt_legendre(j+1,i) = PNDLEG(rst(3,k),j)
!          end do
!    end do
!    end do
!
!  end subroutine find_rst_legendre_bad
 
  !> Using the Legendre-Perez method to find the rst coordinates.
  subroutine find_rst_legendre(rst, pt_x, pt_y, pt_z, Xh, x, y, z, &
                               el_list, n_pts, nelv, resx, resy, resz, &
                               tol)
    integer, intent(in) :: n_pts, nelv
    real(kind=rp), intent(inout) :: rst(3, n_pts)
    type(space_t), intent(inout) :: Xh
    real(kind=rp), allocatable :: pt_x(:)
    real(kind=rp), allocatable :: pt_y(:)
    real(kind=rp), allocatable :: pt_z(:)
    real(kind=rp), intent(inout) :: x(Xh%lx, Xh%ly, Xh%lz, nelv)
    real(kind=rp), intent(inout) :: y(Xh%lx, Xh%ly, Xh%lz, nelv)
    real(kind=rp), intent(inout) :: z(Xh%lx, Xh%ly, Xh%lz, nelv)
    real(kind=rp), intent(inout) :: resx(n_pts)
    real(kind=rp), intent(inout) :: resy(n_pts)
    real(kind=rp), intent(inout) :: resz(n_pts)
    integer, intent(in) :: el_list(n_pts)
    real(kind=rp), intent(in) :: tol 
    real(kind=rp), allocatable :: x_hat(:, :, :, :)
    real(kind=rp), allocatable :: y_hat(:, :, :, :)
    real(kind=rp), allocatable :: z_hat(:, :, :, :)
    real(kind=rp) :: r_legendre(Xh%lx) 
    real(kind=rp) :: s_legendre(Xh%lx) 
    real(kind=rp) :: t_legendre(Xh%lx) 
    real(kind=rp) :: dr_legendre(Xh%lx) 
    real(kind=rp) :: ds_legendre(Xh%lx) 
    real(kind=rp) :: dt_legendre(Xh%lx) 
    real(kind=rp) :: jac(3,3), jacinv(3,3)
    real(kind=rp) :: tmp(Xh%lx), tmp2(Xh%lx), rst_d(3), avgx, avgy, avgz
    integer :: conv_pts
    logical :: converged
    integer :: i, j, e, iter, lx, lx2

    allocate(x_hat(xh%lx, xh%ly, xh%lz, nelv))
    allocate(y_hat(xh%lx, xh%ly, xh%lz, nelv))
    allocate(z_hat(xh%lx, xh%ly, xh%lz, nelv))
    
    lx = Xh%lx
    if (n_pts .lt. 1) return
    
    !> Transform into legendre space
    !> hats are in Legendre space
    call tnsr3d_cpu(x_hat, Xh%lx, x, &
                Xh%lx, Xh%vinv, &
                Xh%vinvt, Xh%vinvt, nelv)

    call tnsr3d_cpu(y_hat, Xh%lx, x_hat, &
                Xh%lx, Xh%v, &
                Xh%vt, Xh%vt, nelv)
    call tnsr3d_cpu(y_hat, Xh%lx, y, &
                Xh%lx, Xh%vinv, &
                Xh%vinvt, Xh%vinvt, nelv)
    call tnsr3d_cpu(z_hat, Xh%lx, z, &
                Xh%lx, Xh%vinv, &
                Xh%vinvt, Xh%vinvt, nelv)
    rst = 0.0 
    do i = 1, n_pts
       iter = 0 
       converged = .false.
       do while (.not. converged)
          iter  = iter + 1
          call legendre_poly(r_legendre,rst(1,i),lx-1)     
          call legendre_poly(s_legendre,rst(2,i),lx-1)     
          call legendre_poly(t_legendre,rst(3,i),lx-1)     
          do j = 0, lx-1
             dr_legendre(j+1) = PNDLEG(rst(1,i),j)
             ds_legendre(j+1) = PNDLEG(rst(2,i),j)
             dt_legendre(j+1) = PNDLEG(rst(3,i),j)
          end do
          e = el_list(i)+1
          call tnsr3d_el_cpu(resx(i), 1, x_hat(1,1,1,e), Xh%lx, &
               r_legendre, s_legendre, t_legendre)
          call tnsr3d_el_cpu(resy(i), 1, y_hat(1,1,1,e), Xh%lx, &
               r_legendre, s_legendre, t_legendre)
          call tnsr3d_el_cpu(resz(i), 1, z_hat(1,1,1,e), Xh%lx, &
               r_legendre, s_legendre, t_legendre)
          call tnsr3d_el_cpu(jac(1,1), 1, x_hat(1,1,1,e), Xh%lx, &
               dr_legendre, s_legendre(1), t_legendre)
          call tnsr3d_el_cpu(jac(1,2), 1, y_hat(1,1,1,e), Xh%lx, &
               dr_legendre, s_legendre(1), t_legendre)
          call tnsr3d_el_cpu(jac(1,3), 1, z_hat(1,1,1,e), Xh%lx, &
               dr_legendre, s_legendre, t_legendre)
          call tnsr3d_el_cpu(jac(2,1), 1, x_hat(1,1,1,e), Xh%lx, &
               r_legendre, ds_legendre, t_legendre)
          call tnsr3d_el_cpu(jac(2,2), 1, y_hat(1,1,1,e), Xh%lx, &
               r_legendre, ds_legendre, t_legendre)
          call tnsr3d_el_cpu(jac(2,3), 1, z_hat(1,1,1,e), Xh%lx, &
               r_legendre, ds_legendre, t_legendre)
          call tnsr3d_el_cpu(jac(3,1), 1, x_hat(1,1,1,e), Xh%lx, &
               r_legendre, s_legendre, dt_legendre)
          call tnsr3d_el_cpu(jac(3,2), 1, y_hat(1,1,1,e), Xh%lx, &
               r_legendre, s_legendre, dt_legendre)
          call tnsr3d_el_cpu(jac(3,3), 1, z_hat(1,1,1,e), Xh%lx, &
               r_legendre, s_legendre, dt_legendre)
          !do jac inverse
          resx(i) = pt_x(i) - resx(i)
          resy(i) = pt_y(i) - resy(i)
          resz(i) = pt_z(i) - resz(i)
          jacinv = matinv39(jac(1,1),jac(1,2),jac(1,3),&
                            jac(2,1),jac(2,2),jac(2,3),&
                            jac(3,1),jac(3,2),jac(3,3))
          rst_d(1) = (resx(i)*jacinv(1,1)+jacinv(2,1)*resy(i)+jacinv(3,1)*resz(i))
          rst_d(2) = (resx(i)*jacinv(1,2)+jacinv(2,2)*resy(i)+jacinv(3,2)*resz(i))
          rst_d(3) = (resx(i)*jacinv(1,3)+jacinv(2,3)*resy(i)+jacinv(3,3)*resz(i))

          rst(1,i) = rst(1,i) + rst_d(1)
          rst(2,i) = rst(2,i) + rst_d(2)
          rst(3,i) = rst(3,i) + rst_d(3)
          conv_pts = 0
          if (norm2(rst_d) < tol) conv_pts = 1
          
          if (norm2(rst_d) > 4) then
             conv_pts = 1 
          end if
          converged = conv_pts .eq. 1
          if (iter .ge. 50) converged = .true.
       end do
    end do
  end subroutine find_rst_legendre
  !> Compares two sets of rst coordinates and checks whether rst2 is better than rst1  
  !> res1 and res2 are the distances to the interpolated xyz coordinate and true xyz coord
  function rst_cmp(rst1, rst2,res1, res2, tol) result(rst2_better)
    real(kind=rp) :: rst1(3), res1(3)
    real(kind=rp) :: rst2(3), res2(3)
    real(kind=rp) :: tol
    logical :: rst2_better
    rst2_better = (norm2(res2) .lt. norm2(res1) .and. &
                  abs(rst2(1)) .lt. 1.0+tol .and. &
                  abs(rst2(2)) .lt. 1.0+tol .and. &
                  abs(rst2(3)) .lt. 1.0+tol)

  end function rst_cmp
 

end module local_interpolation
