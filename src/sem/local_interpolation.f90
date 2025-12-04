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
  use math, only : matinv3, matinv39
  use tensor_cpu
  use device_math, only: device_rzero
  use neko_config, only: NEKO_BCKND_DEVICE
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Interpolation on a set of points with known rst coordinates in elements local
  !! to this process.
  !! Similar to point_interpolator, but prioritizes performance
  !! Only works with arrays of coordinates
  !! Performs interpolation with the configured NEKO_BCKND
  type, public :: local_interpolator_t
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
     procedure, pass(this) :: init_3arrays => local_interpolator_init_3arrays
     procedure, pass(this) :: init_1array => local_interpolator_init_1array
     !> Destructor.
     procedure, pass(this) :: free => local_interpolator_free
     !> Interpolates the scalar field \f$ X \f$ on the specified coordinates
     procedure, pass(this) :: evaluate => local_interpolator_evaluate
     !> COmputes weights based on rst coordinates
     procedure, pass(this) :: compute_weights => local_interpolator_compute_weights
     generic :: init => init_3arrays, init_1array

  end type local_interpolator_t

contains

  !> Initialization of point interpolation.
  !! @param xh Function space.
  subroutine local_interpolator_init_3arrays(this, Xh, r, s, t, n_points)
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

  end subroutine local_interpolator_init_3arrays

  !> Initialization of point interpolation.
  !! @param xh Function space.
  subroutine local_interpolator_init_1array(this, Xh, rst, n_points)
    class(local_interpolator_t), intent(inout), target :: this
    type(space_t), intent(in), target :: Xh
    integer, intent(in) :: n_points
    real(kind=rp), intent(in) :: rst(3,n_points)
    real(kind=rp), allocatable :: r(:), s(:), t(:)
    integer :: i

    if (allocated(r)) deallocate(r)
    allocate(r(n_points))
    if (allocated(s)) deallocate(s)
    allocate(s(n_points))
    if (allocated(t)) deallocate(t)
    allocate(t(n_points))

    do i = 1, n_points
       r(i) = rst(1,i)
       s(i) = rst(2,i)
       t(i) = rst(3,i)
    end do

    call this%init_3arrays(Xh, r, s, t, n_points)

    deallocate(r,s,t)

  end subroutine local_interpolator_init_1array


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
       if ((r(i) <= 1.1_rp .and. r(i) >= -1.1_rp) .and. &
            (s(i) <= 1.1_rp .and. s(i) >= -1.1_rp) .and. &
            (t(i) <= 1.1_rp .and. t(i) >= -1.1_rp)) then
          call fd_weights_full(r(i), this%Xh%zg(:,1), lx-1, 0, this%weights_r(:,i))
          call fd_weights_full(s(i), this%Xh%zg(:,2), lx-1, 0, this%weights_s(:,i))
          call fd_weights_full(t(i), this%Xh%zg(:,3), lx-1, 0, this%weights_t(:,i))
       else
          this%weights_r(:,i) = 0.0_rp
          this%weights_s(:,i) = 0.0_rp
          this%weights_t(:,i) = 0.0_rp
       end if

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
  subroutine local_interpolator_evaluate(this, interp_values, el_list, field, &
       nel, on_host)
    class(local_interpolator_t), intent(inout) :: this
    integer, intent(in) :: el_list(this%n_points)
    integer, intent(in) :: nel
    real(kind=rp), intent(inout) :: interp_values(this%n_points)
    real(kind=rp), intent(inout) :: field(this%Xh%lxyz, nel)
    logical, intent(in) :: on_host

    call tnsr3d_el_list(interp_values, 1, field, this%Xh%lx, &
         this%weights_r, this%weights_s, this%weights_t, el_list, &
         this%n_points, on_host)

  end subroutine local_interpolator_evaluate

  !> Constructs the Jacobian, returns a 3-by-3 times number of points where
  !! \f$ [J(\mathbf{r}]_{ij} = \frac{d\mathbf{x}_i}{d\mathbf{r}_j}\f$.
  !! @param rst r,s,t coordinates.
  !! @param X Values of the field \f$ X \f$ at points.
  !! @param Y Values of the field \f$ Y \f$ at points.
  !! @param Z Values of the field \f$ Z \f$ at points.
  subroutine jacobian(jac, rst, x, y, z, n_pts, Xh)
    integer, intent(in) :: n_pts
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
       call triple_tensor_product(tmp(1), X(:,:,:,i), lx, hr(:,2), hs(:,1), &
            ht(:,1))
       jac(1,1,i) = tmp(1)
       call triple_tensor_product(tmp(1), Y(:,:,:,i), lx, hr(:,2), hs(:,1), &
            ht(:,1))
       jac(1,2,i) = tmp(1)
       call triple_tensor_product(tmp(1), Z(:,:,:,i), lx, hr(:,2), hs(:,1), &
            ht(:,1))
       jac(1,3,i) = tmp(1)

       ! d(x,y,z)/ds
       call triple_tensor_product(tmp, X(:,:,:,i), Y(:,:,:,i), Z(:,:,:,i), lx, &
            hr(:,1), hs(:,2), ht(:,1))
       jac(2,:,i) = tmp

       ! d(x,y,z)/dt
       call triple_tensor_product(tmp, X(:,:,:,i), Y(:,:,:,i), Z(:,:,:,i), lx, &
            hr(:,1), hs(:,1), ht(:,2))
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
       tmp = matinv3(real(jacinv(:,:,3),xp))
       jacinv(:,:,i) = tmp
    end do

  end subroutine jacobian_inverse

end module local_interpolation
