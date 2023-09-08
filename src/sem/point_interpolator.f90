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
!> Routines to interpolate scalar/vector fields
module point_interpolator
  use tensor, only: triple_tensor_product, tnsr3d_el_list
  use space, only: space_t, GL, GLL
  use num_types, only: rp
  use point, only: point_t
  use math, only: abscmp
  use fast3d, only: fd_weights_full, setup_intp
  use utils, only: neko_error
  use field, only: field_t
  use field_list, only: field_list_t
  use device
  use device_math, only: device_rzero
  use neko_config, only: NEKO_BCKND_DEVICE
  use, intrinsic :: iso_c_binding
  implicit none

  !> Field interpolator to artbitrary points within an element.
  type :: point_interpolator_t
     !> First space.
     type(space_t), pointer :: Xh => null()
     !> Weights in the r direction on the host
     real(kind=rp), allocatable :: weights_r(:,:)
     !> Weights in the s direction on the host
     real(kind=rp), allocatable :: weights_s(:,:)
     !> Weights in the s direction on the host
     real(kind=rp), allocatable :: weights_t(:,:)
     !> Weights in the r direction on the device
     type(c_ptr) :: weights_r_d = C_NULL_PTR
     !> Weights in the s direction on the device
     type(c_ptr) :: weights_s_d = C_NULL_PTR
     !> Weights in the t direction on the device
     type(c_ptr) :: weights_t_d = C_NULL_PTR
     !> To know if we need to recompute the weights
     logical :: update_weights = .true.
   contains
     !> Constructor.
     procedure, pass(this) :: init => point_interpolator_init
     !> Destructor.
     procedure, pass(this) :: free => point_interpolator_free
     !> Interpolates a scalar field \f$ X \f$ on a set of points.
     procedure, pass(this) :: point_interpolator_interpolate_scalar
     !> Interpolates a vector field \f$ \vec f = (X,Y,Z) \f$ on a set of points.
     procedure, pass(this) :: point_interpolator_interpolate_vector
     !> Constructs the Jacobian at a single point.
     procedure, pass(this) :: point_interpolator_interpolate_jacobian
     !> Interpolates a vector field and builds the Jacobian at a single point.
     procedure, pass(this) :: point_interpolator_interpolate_vector_jacobian
     !> Interpolates a scalar field \f$ X \f$ from a `field_t` object and element id
     procedure, pass(this) :: point_interpolator_interpolate_field_el
     procedure, pass(this) :: point_interpolator_interpolate_fields
     !> Interpolates a scalar or vector field on a set of points.
     generic :: interpolate => point_interpolator_interpolate_scalar, &
          point_interpolator_interpolate_vector, &
          point_interpolator_interpolate_field_el, &
          point_interpolator_interpolate_fields
     !> Constructs the Jacobian for a point \f$ (r,s,t) \f$.
     generic :: jacobian => point_interpolator_interpolate_jacobian, &
          point_interpolator_interpolate_vector_jacobian

  end type point_interpolator_t

contains

  !> Initialization of point interpolation.
  !! @param xh Function space.
  subroutine point_interpolator_init(this, Xh)
    class(point_interpolator_t), intent(inout), target :: this
    type(space_t), intent(in), target :: Xh

    if ((Xh%t .eq. GL) .or. (Xh%t .eq. GLL)) then
    else
       call neko_error('Unsupported interpolation')
    end if

    this%Xh => Xh

  end subroutine point_interpolator_init

  !> Free pointers
  subroutine point_interpolator_free(this)
    class(point_interpolator_t), intent(inout) :: this

    if (associated(this%Xh)) this%Xh => null()
    if (allocated(this%weights_r)) deallocate(this%weights_r)
    if (allocated(this%weights_s)) deallocate(this%weights_s)
    if (allocated(this%weights_t)) deallocate(this%weights_t)

  end subroutine point_interpolator_free

  !> Interpolates a scalar field \f$ X \f$ on a set of \f$ N \f$ points
  !! \f$ \mathbf{r}_i , i\in[1,N]\f$. Returns a vector of N coordinates
  !! \f$ [x_i(\mathbf{r}_i)], i\in[1,N]\f$.
  !! @param rst r,s,t coordinates.
  !! @param X Values of the field \f$ X \f$ at GLL points.
  function point_interpolator_interpolate_scalar(this, rst, X) result(res)
    class(point_interpolator_t), intent(in) :: this
    type(point_t), intent(in) :: rst(:)
    real(kind=rp), intent(inout) :: X(this%Xh%lx, this%Xh%ly, this%Xh%lz)
    real(kind=rp), allocatable :: res(:)

    real(kind=rp) :: hr(this%Xh%lx), hs(this%Xh%ly), ht(this%Xh%lz)
    integer :: lx,ly,lz, i
    integer :: N
    lx = this%Xh%lx
    ly = this%Xh%ly
    lz = this%Xh%lz

    N = size(rst)
    allocate(res(N))

    !
    ! Compute weights and perform interpolation for the first point
    !
    call fd_weights_full(real(rst(1)%x(1), rp), this%Xh%zg(:,1), lx-1, 0, hr)
    call fd_weights_full(real(rst(1)%x(2), rp), this%Xh%zg(:,2), ly-1, 0, hs)
    call fd_weights_full(real(rst(1)%x(3), rp), this%Xh%zg(:,3), lz-1, 0, ht)

    ! And interpolate!
    call triple_tensor_product(res(1),X,lx,hr,hs,ht)

    if (N .eq. 1) return

    !
    ! Loop through the rest of the points
    !
    do i = 2, N

       ! If the coordinates are different, then recompute weights
       if ( .not. abscmp(rst(i)%x(1), rst(i-1)%x(1)) ) then
          call fd_weights_full(real(rst(i)%x(1), rp), &
                               this%Xh%zg(:,1), lx-1, 0, hr)
       end if
       if ( .not. abscmp(rst(i)%x(2), rst(i-1)%x(2)) ) then
          call fd_weights_full(real(rst(i)%x(2), rp), &
                               this%Xh%zg(:,2), ly-1, 0, hs)
       end if
       if ( .not. abscmp(rst(i)%x(3), rst(i-1)%x(3)) ) then
          call fd_weights_full(real(rst(i)%x(3), rp), &
                               this%Xh%zg(:,3), lz-1, 0, ht)
       end if

       ! And interpolate!
       call triple_tensor_product(res(i), X, lx, hr, hs, ht)

    end do

  end function point_interpolator_interpolate_scalar

  !> Interpolates a vector field \f$ \vec f = (X,Y,Z) \f$ on a set of \f$ N \f$ points
  !! \f$ \mathbf{r}_i \f$. Returns an array of N points
  !! \f$ [x(\mathbf{r}_i), y(\mathbf{r}_i), z(\mathbf{r}_i)], i\in[1,N]\f$.
  !! @param N number of points (use `1` to interpolate a scalar).
  !! @param X Values of the field \f$ X \f$ at GLL points.
  !! @param Y Values of the field \f$ Y \f$ at GLL points.
  !! @param Z Values of the field \f$ Z \f$ at GLL points.
  function point_interpolator_interpolate_vector(this, rst, X, Y, Z) result(res)
    class(point_interpolator_t), intent(in) :: this
    type(point_t), intent(in) :: rst(:)
    real(kind=rp), intent(inout) :: X(this%Xh%lx, this%Xh%ly, this%Xh%lz)
    real(kind=rp), intent(inout) :: Y(this%Xh%lx, this%Xh%ly, this%Xh%lz)
    real(kind=rp), intent(inout) :: Z(this%Xh%lx, this%Xh%ly, this%Xh%lz)

    type(point_t), allocatable :: res(:)
    real(kind=rp), allocatable :: tmp(:,:)
    real(kind=rp) :: hr(this%Xh%lx), hs(this%Xh%ly), ht(this%Xh%lz)
    integer :: lx,ly,lz, i
    integer :: N
    lx = this%Xh%lx
    ly = this%Xh%ly
    lz = this%Xh%lz
    
    N = size(rst)
    allocate(res(N))
    allocate(tmp(3, N))

    !
    ! Compute weights and perform interpolation for the first point
    !
    call fd_weights_full(real(rst(1)%x(1), rp), this%Xh%zg(:,1), lx-1, 0, hr)
    call fd_weights_full(real(rst(1)%x(2), rp), this%Xh%zg(:,2), ly-1, 0, hs)
    call fd_weights_full(real(rst(1)%x(3), rp), this%Xh%zg(:,3), lz-1, 0, ht)

    ! And interpolate!
    call triple_tensor_product(tmp(:,1), X, Y, Z, lx, hr, hs, ht)

    if (N .eq. 1) then
       res(1)%x = tmp(:, 1)
       return
    end if
    

    !
    ! Loop through the rest of the points
    !
    do i = 2, N

       ! If the coordinates are different, then recompute weights
       if ( .not. abscmp(rst(i)%x(1), rst(i-1)%x(1)) ) then
          call fd_weights_full(real(rst(i)%x(1), rp), &
                               this%Xh%zg(:,1), lx-1, 0, hr)
       end if
       if ( .not. abscmp(rst(i)%x(2), rst(i-1)%x(2)) ) then
          call fd_weights_full(real(rst(i)%x(2), rp), &
                               this%Xh%zg(:,2), ly-1, 0, hs)
       end if
       if ( .not. abscmp(rst(i)%x(3), rst(i-1)%x(3)) ) then
          call fd_weights_full(real(rst(i)%x(3), rp), &
                               this%Xh%zg(:,3), lz-1, 0, ht)
       end if

       ! And interpolate!
       call triple_tensor_product(tmp(:,i), X, Y, Z, lx, hr, hs, ht)
    end do

    ! Cast result to point_t dp
    do i = 1, N
       res(i)%x = tmp(:,i)
    end do

  end function point_interpolator_interpolate_vector

  !> Interpolates a vector field \f$ \vec f = (X,Y,Z) \f$ and constructs the Jacobian
  !! at a point \f$ (r,s,t) \f$. Returns a vector
  !! \f$ [x(\mathbf{r}_i), y(\mathbf{r}_i), z(\mathbf{r}_i)], i\in[1,N]\f$.
  !! @param jac Jacobian.
  !! @param rst r,s,t coordinates;
  !! @param X Values of the field \f$ X \f$ at GLL points.
  !! @param Y Values of the field \f$ Y \f$ at GLL points.
  !! @param Z Values of the field \f$ Z \f$ at GLL points.
  function point_interpolator_interpolate_vector_jacobian(this, jac, rst, X, Y, Z) result(res)
    class(point_interpolator_t), intent(in) :: this
    real(kind=rp), intent(inout) :: jac(3,3)
    type(point_t), intent(in) :: rst
    real(kind=rp), intent(inout) :: X(this%Xh%lx, this%Xh%ly, this%Xh%lz)
    real(kind=rp), intent(inout) :: Y(this%Xh%lx, this%Xh%ly, this%Xh%lz)
    real(kind=rp), intent(inout) :: Z(this%Xh%lx, this%Xh%ly, this%Xh%lz)

    real(kind=rp) :: hr(this%Xh%lx, 2), hs(this%Xh%ly, 2), ht(this%Xh%lz, 2)
    type(point_t) :: res
    real(kind=rp) :: tmp(3)

    integer :: lx,ly,lz, i
    lx = this%Xh%lx
    ly = this%Xh%ly
    lz = this%Xh%lz

    !
    ! Compute weights
    !
    call fd_weights_full(real(rst%x(1), rp), this%Xh%zg(:,1), lx-1, 1, hr)
    call fd_weights_full(real(rst%x(2), rp), this%Xh%zg(:,2), ly-1, 1, hs)
    call fd_weights_full(real(rst%x(3), rp),this%Xh%zg(:,3), lz-1, 1, ht)

    !
    ! Interpolate
    !
    call triple_tensor_product(tmp, X, Y, Z, lx, hr(:,1), hs(:,1), ht(:,1))
    res%x = dble(tmp)! Cast from rp -> point_t dp 

    !
    ! Build jacobian
    !

    ! d(x,y,z)/dr
    call triple_tensor_product(jac(1,:), X,Y,Z, lx, hr(:,2), hs(:,1), ht(:,1))

    ! d(x,y,z)/ds
    call triple_tensor_product(jac(2,:), X,Y,Z, lx, hr(:,1), hs(:,2), ht(:,1))

    ! d(x,y,z)/dt
    call triple_tensor_product(jac(3,:), X,Y,Z, lx, hr(:,1), hs(:,1), ht(:,2))


  end function point_interpolator_interpolate_vector_jacobian

  !> Constructs the Jacobian, returns a 3-by-3 array where
  !! \f$ [J(\mathbf{r}]_{ij} = \frac{d\mathbf{x}_i}{d\mathbf{r}_j}\f$.
  !! @param rst r,s,t coordinates.
  !! @param X Values of the field \f$ X \f$ at GLL points.
  !! @param Y Values of the field \f$ Y \f$ at GLL points.
  !! @param Z Values of the field \f$ Z \f$ at GLL points.
  function point_interpolator_interpolate_jacobian(this, rst, X,Y,Z) result(jac)
    class(point_interpolator_t), intent(in) :: this
    type(point_t), intent(in) :: rst
    real(kind=rp), intent(inout) :: X(this%Xh%lx, this%Xh%ly, this%Xh%lz)
    real(kind=rp), intent(inout) :: Y(this%Xh%lx, this%Xh%ly, this%Xh%lz)
    real(kind=rp), intent(inout) :: Z(this%Xh%lx, this%Xh%ly, this%Xh%lz)

    real(kind=rp) :: jac(3,3)

    real(kind=rp) :: hr(this%Xh%lx, 2), hs(this%Xh%ly, 2), ht(this%Xh%lz, 2)
    integer :: lx, ly, lz
    lx = this%Xh%lx
    ly = this%Xh%ly
    lz = this%Xh%lz

    ! Weights
    call fd_weights_full(real(rst%x(1), rp), this%Xh%zg(:,1), lx-1, 1, hr)
    call fd_weights_full(real(rst%x(2), rp), this%Xh%zg(:,2), ly-1, 1, hs)
    call fd_weights_full(real(rst%x(3), rp), this%Xh%zg(:,3), lz-1, 1, ht)

    ! d(x,y,z)/dr
    call triple_tensor_product(jac(1,:), X, Y, Z, lx, hr(:,2), hs(:,1), ht(:,1))

    ! d(x,y,z)/ds
    call triple_tensor_product(jac(2,:), X, Y, Z, lx, hr(:,1), hs(:,2), ht(:,1))

    ! d(x,y,z)/dt
    call triple_tensor_product(jac(3,:), X, Y, Z, lx, hr(:,1), hs(:,1), ht(:,2))
  end function point_interpolator_interpolate_jacobian

  !> Interpolates a `field_t` \f$ X \f$ for a given element on a set of
  !! \f$ N \f$ points \f$ \mathbf{r}_i , i\in[1,N]\f$. Returns a vector
  !! of N coordinates \f$ [x_i(\mathbf{r}_i)], i\in[1,N]\f$.
  !! @param rst r,s,t coordinates.
  !! @param X Values of the field \f$ X \f$ at GLL points.
  function point_interpolator_interpolate_field_el(this, rst, f, el) result(res)
    class(point_interpolator_t), intent(in) :: this
    type(point_t), intent(in) :: rst(:)
    type(field_t), intent(inout) :: f
    integer, intent(in) :: el
    real(kind=rp), allocatable :: res(:)

    res = point_interpolator_interpolate_scalar(this, rst, f%x(:,:,:,el))

  end function point_interpolator_interpolate_field_el

  function point_interpolator_interpolate_fields(this, rst, el_owners, sampled_fields_list) result(res)
    class(point_interpolator_t), intent(inout) :: this
    type(point_t), intent(inout), allocatable :: rst(:)
    integer, intent(in), allocatable :: el_owners(:)
    type(field_list_t), intent(in) :: sampled_fields_list
    real(kind=rp), allocatable :: res(:,:)

    integer :: n_points, n_fields, lx, n, i,j
    type(c_ptr) :: tmp_d = C_NULL_PTR
    real(kind=rp), allocatable :: tmp(:)

    lx = this%Xh%lx
    n_points = size(rst)
    n_fields = size(sampled_fields_list%fields)

    !
    ! Only update weights if necessary
    !
    if (this%update_weights) then

       if (.not. allocated(this%weights_r)) allocate(this%weights_r(lx,n_points))
       if (.not. allocated(this%weights_s)) allocate(this%weights_s(lx,n_points))
       if (.not. allocated(this%weights_t)) allocate(this%weights_t(lx,n_points))

       !
       ! Build weights
       !
       do i = 1, n_points
          call fd_weights_full(rst(i)%x(1), this%Xh%zg(:,1), &
               lx-1, 0, this%weights_r(:,i))
          call fd_weights_full(rst(i)%x(2), this%Xh%zg(:,2), &
               lx-1, 0, this%weights_s(:,i))
          call fd_weights_full(rst(i)%x(3), this%Xh%zg(:,3), &
               lx-1, 0, this%weights_t(:,i))
       end do

       this%update_weights = .false.

    end if

    allocate(res(n_points, n_fields))
    allocate(tmp(n_points))
    n = n_points*lx

    !
    ! Associate device pointers
    !
    if (NEKO_BCKND_DEVICE .eq. 1) then

       call device_map(this%weights_r, this%weights_r_d, n)
       call device_map(this%weights_s, this%weights_s_d, n)
       call device_map(this%weights_t, this%weights_t_d, n)
       call device_memcpy(this%weights_r, this%weights_r_d, n, HOST_TO_DEVICE)
       call device_memcpy(this%weights_s, this%weights_s_d, n, HOST_TO_DEVICE)
       call device_memcpy(this%weights_t, this%weights_t_d, n, HOST_TO_DEVICE)

       call device_map(tmp, tmp_d, n_points)
       call device_rzero(tmp_d, n_points)

    end if

    !
    ! Interpolate each field
    !
    do i = 1, n_fields

       call tnsr3d_el_list(tmp, 1, sampled_fields_list%fields(i)%f%x, lx, &
            this%weights_r, this%weights_s, this%weights_t, &
            el_owners, n_points)
       res(:,i) = tmp
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) call device_free(tmp_d)

  end function point_interpolator_interpolate_fields

end module point_interpolator
