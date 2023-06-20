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
!> Routines to interpolate between different spaces
module interpolation
  use speclib
  use device
  use utils
  use math
  use fast3d
  use tensor
  use tensor_cpu
  use space
  use device
  use mxm_wrapper, only: mxm
  use coefs, only: coef_t
  use point, only: point_t
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Interpolation between two \ref space::space_t.
  !! @details
  !! This type implements functionality to interpolate between a pair of spaces.
  !! Simply put, given some data of form (lx1, lx1, lx1, nelem) we can map it to
  !! (lx2, lx2, lx2, nelem), corresponding to a different polynomial order in
  !! each element.
  type, public :: interpolator_t
     !> First space
     type(space_t), pointer :: Xh
     !> Second space
     type(space_t), pointer :: Yh
     !> Interpolation weights from Xh to Yh
     real(kind=rp), allocatable :: Xh_to_Yh(:,:), Xh_to_YhT(:,:)
     !> Interpolation weights from Yh to Xh
     real(kind=rp), allocatable :: Yh_to_Xh(:,:), Yh_to_XhT(:,:)
     !> Device pointer for Xh_to_Yh
     type(c_ptr) :: Xh_Yh_d = C_NULL_PTR
     !> Device pointer for Xh_to_YhT
     type(c_ptr) :: Xh_YhT_d = C_NULL_PTR
     !> Device pointer for Yh_to_Xh
     type(c_ptr) :: Yh_Xh_d = C_NULL_PTR
     !> Device pointer for Yh_to_XhT
     type(c_ptr) :: Yh_XhT_d = C_NULL_PTR
     !> Pointer to coefficients (for interpolation)
     type(coef_t), pointer :: coef
   contains

     !> Constructor.
     procedure, pass(this) :: init => interpolator_init
     !> Destructor.
     procedure, pass(this) :: free => interpolator_free
     !> Interpolate an array to one of Xh or Yh.
     procedure, pass(this) :: map => interpolator_map_twospaces
     !> Interpolate an array to one of Xh or Yh on the host.
     procedure, pass(this) :: map_host => interpolator_map_host
     !> Interpolates a scalar field \f$ X \f on a set of \f$ N \f$ points.
     procedure, pass(this) :: interp1d => interpolator_interp1d
     !> Interpolates a vector field \f$ \vec f = (X,Y,Z) \f$ on a set of \f$ N \f$ points.
     procedure, pass(this) :: interp3d => interpolator_interp3d
     !> Constructs the Jacobian for a point \f$ (r,s,t) \f$.
     procedure, pass(this) :: jacobian => interpolator_jacobian
     !> Interpolates a vector field and builds the Jacobian at a single point.
     procedure, pass(this) :: interp3d_jacobian => interpolator_interp3d_jacobian

  end type interpolator_t

contains


  !> Constructor for interpolation between two spaces if `Yh` is specifed,
  !! or regular spectral interpolation if `coef` is specified.
  !! @param Xh The first space.
  !! @param Xh The second space.
  !! @param coef Coefficients.
  !! @NOTE: If `coef` is present, `Xh` and `Yh` will be disregarded and the
  !! initialization will be performed with the space in `coef%xh`.
  subroutine interpolator_init(this, Xh, Yh, coef)
    class(interpolator_t), intent(inout), target :: this
    type(space_t), intent(inout), target :: Xh
    type(space_t), intent(inout), target, optional :: Yh
    type(coef_t), intent(inout), target, optional :: coef

    if (present(coef)) then
       call interpolator_init_onespace(this,coef)
    else if (present(Yh)) then
       call interpolator_init_twospaces(this,Xh,Yh)
    end if

  end subroutine interpolator_init

  !> Constructor to initialize with two different spaces.
  !> @param Xh The first space.
  !> @param Xh The second space.
  subroutine interpolator_init_twospaces(this, Xh, Yh)
    class(interpolator_t), intent(inout), target :: this
    type(space_t), intent(inout), target :: Xh
    type(space_t), intent(inout), target, optional :: Yh
    integer :: deg_derivate

    call this%free()

    allocate(this%Xh_to_Yh(Yh%lx,Xh%lx))
    allocate(this%Xh_to_YhT(Xh%lx,Yh%lx))
    allocate(this%Yh_to_Xh(Xh%lx,Yh%lx))
    allocate(this%Yh_to_XhT(Yh%lx,Xh%lx))

    if (Xh%t .eq. GLL .and. Yh%t .eq. GLL) then
    else if ((Xh%t .eq. GL .and. Yh%t .eq. GLL) .or. &
         (Yh%t .eq. GL .and. Xh%t .eq. GLL)) then
    else
       call neko_error('Unsupported interpolation')
    end if
    deg_derivate = 0
    call setup_intp(this%Xh_to_Yh, this%Xh_to_YhT, &
         Yh%zg, Xh%zg, Yh%lx, Xh%lx, deg_derivate)
    call setup_intp(this%Yh_to_Xh, this%Yh_to_XhT, &
         Xh%zg, Yh%zg, Xh%lx, Yh%lx, deg_derivate)

    this%Xh => Xh
    this%Yh => Yh
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%Xh_to_Yh, this%Xh_Yh_d, Yh%lx*Xh%lx)
       call device_map(this%Xh_to_YhT, this%Xh_YhT_d, Yh%lx*Xh%lx)
       call device_map(this%Yh_to_Xh, this%Yh_Xh_d, Yh%lx*Xh%lx)
       call device_map(this%Yh_to_XhT, this%Yh_XhT_d, Yh%lx*Xh%lx)
       call device_memcpy(this%Xh_to_Yh, this%Xh_Yh_d, Yh%lx*Xh%lx, HOST_TO_DEVICE)
       call device_memcpy(this%Xh_to_YhT, this%Xh_YhT_d, Yh%lx*Xh%lx, HOST_TO_DEVICE)
       call device_memcpy(this%Yh_to_Xh, this%Yh_Xh_d, Yh%lx*Xh%lx, HOST_TO_DEVICE)
       call device_memcpy(this%Yh_to_XhT, this%Yh_XhT_d, Yh%lx*Xh%lx, HOST_TO_DEVICE)
    end if

  end subroutine interpolator_init_twospaces

  !> Constructor with only one space, for "regular" interpolation.
  !! @param Xh Space with which to initialize.
  !! @param coef Coefficients.
  subroutine interpolator_init_onespace(this, coef)
    class(interpolator_t), intent(inout), target :: this
    type(coef_t), intent(inout), target :: coef

    if ((coef%xh%t .eq. GL) .or. (coef%xh%t .eq. GLL)) then
    else
       call neko_error('Unsupported interpolation')
    end if

    this%Xh => coef%xh
    this%coef => coef

  end subroutine interpolator_init_onespace

  subroutine interpolator_free(this)
    class(interpolator_t), intent(inout) :: this

    if (allocated(this%Xh_to_Yh)) then
       deallocate(this%Xh_to_Yh)
    end if
    if (allocated(this%Xh_to_YhT)) then
       deallocate(this%Xh_to_YhT)
    end if
    if (allocated(this%Yh_to_Xh)) then
       deallocate(this%Yh_to_Xh)
    end if
    if (allocated(this%Yh_to_XhT)) then
       deallocate(this%Yh_to_XhT)
    end if
    if (c_associated(this%Yh_Xh_d)) then
       call device_free(this%Yh_Xh_d)
    end if
    if (c_associated(this%Yh_XhT_d)) then
       call device_free(this%Yh_XhT_d)
    end if
    if (c_associated(this%Xh_Yh_d)) then
       call device_free(this%Xh_Yh_d)
    end if
    if (c_associated(this%Xh_YhT_d)) then
       call device_free(this%Xh_YhT_d)
    end if

  end subroutine interpolator_free

  !> Interpolates an array to one of Xh or Yh.
  !! @param x Original array.
  !! @param y Interpolated array.
  !! @param nel Number of elements in the mesh.
  !! @param to_space The space to interpolate to, must be either Xh or Yh.
  subroutine interpolator_map_twospaces(this, y, x, nel, to_space)
    class(interpolator_t), intent(inout) :: this
    integer :: nel
    type(space_t) :: to_space
    real(kind=rp), intent(inout) :: x(1,nel)
    real(kind=rp), intent(inout) :: y(1,nel)
    if (to_space .eq. this%Yh) then
       call tnsr3d(y, this%Yh%lx, x, &
                   this%Xh%lx,this%Yh_to_XhT, &
                   this%Yh_to_Xh, this%Yh_to_Xh, nel)
    else if (to_space .eq. this%Xh) then
       call tnsr3d(y, this%Xh%lx, x, &
                   this%Yh%lx,this%Yh_to_Xh, &
                   this%Yh_to_XhT, this%Yh_to_XhT, nel)
    else
       call neko_error('Invalid interpolation')
    end if
  end subroutine interpolator_map_twospaces

  !> Interpolates a scalar field \f$ X \f$ on a set of \f$ N \f$ points
  !! \f$ \mathbf{r}_i , i\in[1,N]\f$. Returns a vector of N coordinates
  !! \f$ [x_i(\mathbf{r}_i)], i\in[1,N]\f$.
  !! @param N number of points (use `1` to interpolate a scalar).
  !! @param rst r,s,t coordinates.
  !! @param X Values of the field \f$ X \f$ at GLL points.
  function interpolator_interp1d(this, N, rst, X) result(res)
    class(interpolator_t), intent(in) :: this
    integer, intent(in) :: N
    type(point_t), intent(in) :: rst(N)
    real(kind=rp), intent(inout) :: X(this%xh%lx*this%xh%ly*this%xh%lz)
    real(kind=rp) :: res(N)

    real(kind=rp) :: hr(this%xh%lx), hs(this%xh%ly), ht(this%xh%lz)
    integer :: lx,ly,lz, i
    lx = this%xh%lx
    ly = this%xh%ly
    lz = this%xh%lz

    !
    ! Compute weights and perform interpolation for the first point
    !

    call fd_weights_full(rst(1)%x(1), this%xh%zg(:,1), lx-1, 0, hr)
    call fd_weights_full(rst(1)%x(2), this%xh%zg(:,2), ly-1, 0, hs)
    call fd_weights_full(rst(1)%x(3), this%xh%zg(:,3), lz-1, 0, ht)

    ! And interpolate!
    call tensor_scalar1(res(1),X,lx,hr,hs,ht)

    if (N .eq. 1) return

    !
    ! Loop through the rest of the points
    !
    do i = 2, N

       ! If the coordinates are different, then recompute weights
       if ( .not. abscmp(rst(i)%x(1), rst(i-1)%x(1)) ) then
          call fd_weights_full(rst(i)%x(1), this%xh%zg(:,1), lx-1, 0, hr)
       end if
       if ( .not. abscmp(rst(i)%x(2), rst(i-1)%x(2)) ) then
          call fd_weights_full(rst(i)%x(2), this%xh%zg(:,2), ly-1, 0, hs)
       end if
       if ( .not. abscmp(rst(i)%x(3), rst(i-1)%x(3)) ) then
          call fd_weights_full(rst(i)%x(3), this%xh%zg(:,3), lz-1, 0, ht)
       end if

       ! And interpolate!
       call tensor_scalar1(res(i),X,lx,hr,hs,ht)

    end do

  end function interpolator_interp1d

  !> Interpolates a vector field \f$ \vec f = (X,Y,Z) \f$ on a set of \f$ N \f$ points
  !! \f$ \mathbf{r}_i \f$. Returns an array of N points
  !! \f$ [x(\mathbf{r}_i), y(\mathbf{r}_i), z(\mathbf{r}_i)], i\in[1,N]\f$.
  !! @param N number of points (use `1` to interpolate a scalar).
  !! @param rst Array of `N` (r,s,t) coordinates.
  !! @param X Values of the field \f$ X \f$ at GLL points.
  !! @param Y Values of the field \f$ Y \f$ at GLL points.
  !! @param Z Values of the field \f$ Z \f$ at GLL points.
  function interpolator_interp3d(this, N, rst, X, Y, Z) result(res)
    class(interpolator_t), intent(in) :: this
    integer, intent(in) :: N
    type(point_t), intent(in) :: rst(N)
    real(kind=rp), intent(inout) :: X(this%xh%lx*this%xh%ly*this%xh%lz)
    real(kind=rp), intent(inout) :: Y(this%xh%lx*this%xh%ly*this%xh%lz)
    real(kind=rp), intent(inout) :: Z(this%xh%lx*this%xh%ly*this%xh%lz)

    type(point_t):: res(N)
    real(kind=rp) :: hr(this%xh%lx), hs(this%xh%ly), ht(this%xh%lz)
    integer :: lx,ly,lz, i
    lx = this%xh%lx
    ly = this%xh%ly
    lz = this%xh%lz

    !
    ! Compute weights and perform interpolation for the first point
    !
    call fd_weights_full(rst(1)%x(1), this%xh%zg(:,1), lx-1, 0, hr)
    call fd_weights_full(rst(1)%x(2), this%xh%zg(:,2), ly-1, 0, hs)
    call fd_weights_full(rst(1)%x(3), this%xh%zg(:,3), lz-1, 0, ht)

    ! And interpolate!
    call tensor_scalar1(res(1)%x(1),X,lx,hr,hs,ht)
    call tensor_scalar1(res(1)%x(2),Y,ly,hr,hs,ht)
    call tensor_scalar1(res(1)%x(3),Z,lz,hr,hs,ht)

    if (N .eq. 1) return

    !
    ! Loop through the rest of the points
    !
    do i = 2, N

       ! If the coordinates are different, then recompute weights
       if ( .not. abscmp(rst(i)%x(1), rst(i-1)%x(1)) ) then
          call fd_weights_full(rst(i)%x(1), this%xh%zg(:,1), lx-1, 0, hr)
       end if
       if ( .not. abscmp(rst(i)%x(2), rst(i-1)%x(2)) ) then
          call fd_weights_full(rst(i)%x(2), this%xh%zg(:,2), ly-1, 0, hs)
       end if
       if ( .not. abscmp(rst(i)%x(3), rst(i-1)%x(3)) ) then
          call fd_weights_full(rst(i)%x(3), this%xh%zg(:,3), lz-1, 0, ht)
       end if

       ! And interpolate!
       call tensor_scalar3(res(i)%x, X, Y, Z, lx, hr, hs, ht)
    end do

  end function interpolator_interp3d

  !> Interpolates a vector field \f$ \vec f = (X,Y,Z) \f$ and constructs the Jacobian
  !! at a point \f$ (r,s,t) \f$. Returns a vector
  !! \f$ [x(\mathbf{r}_i), y(\mathbf{r}_i), z(\mathbf{r}_i)], i\in[1,N]\f$.
  !! @param jac Jacobian.
  !! @param rst r,s,t coordinates;
  !! @param X Values of the field \f$ X \f$ at GLL points.
  !! @param Y Values of the field \f$ Y \f$ at GLL points.
  !! @param Z Values of the field \f$ Z \f$ at GLL points.
  function interpolator_interp3d_jacobian(this, jac, rst, X, Y, Z) result(res)
    class(interpolator_t), intent(in) :: this
    real(kind=rp), intent(inout) :: jac(3,3)
    type(point_t), intent(in) :: rst
    real(kind=rp), intent(inout) :: X(this%xh%lx*this%xh%ly*this%xh%lz)
    real(kind=rp), intent(inout) :: Y(this%xh%lx*this%xh%ly*this%xh%lz)
    real(kind=rp), intent(inout) :: Z(this%xh%lx*this%xh%ly*this%xh%lz)

    real(kind=rp) :: hr(this%xh%lx, 2), hs(this%xh%ly, 2), ht(this%xh%lz, 2)
    type(point_t) :: res

    integer :: lx,ly,lz, i
    lx = this%xh%lx
    ly = this%xh%ly
    lz = this%xh%lz

    !
    ! Compute weights
    !
    call fd_weights_full(rst%x(1), this%xh%zg(:,1), lx-1, 1, hr)
    call fd_weights_full(rst%x(2), this%xh%zg(:,2), ly-1, 1, hs)
    call fd_weights_full(rst%x(3), this%xh%zg(:,3), lz-1, 1, ht)

    !
    ! Interpolate
    !
    call tensor_scalar1(res%x(1), X, lx, hr(:,1), hs(:,1), ht(:,1))
    call tensor_scalar1(res%x(2), Y, ly, hr(:,1), hs(:,1), ht(:,1))
    call tensor_scalar1(res%x(3), Z, lz, hr(:,1), hs(:,1), ht(:,1))

    !
    ! Build jacobian
    !

    ! d(x,y,z)/dr
    call tensor_scalar3(jac(1,:), X,Y,Z, lx, hr(:,2), hs(:,1), ht(:,1))

    ! d(x,y,z)/ds
    call tensor_scalar3(jac(2,:), X,Y,Z, lx, hr(:,1), hs(:,2), ht(:,1))

    ! d(x,y,z)/dt
    call tensor_scalar3(jac(3,:), X,Y,Z, lx, hr(:,1), hs(:,1), ht(:,2))


  end function interpolator_interp3d_jacobian

  !> Constructs the Jacobian, returns a 3-by-3 array where
  !! \f$ [J(\mathbf{r}]_{ij} = \frac{d\mathbf{x}_i}{d\mathbf{r}_j}\f$.
  !! @param rst r,s,t coordinates.
  !! @param X Values of the field \f$ X \f$ at GLL points.
  !! @param Y Values of the field \f$ Y \f$ at GLL points.
  !! @param Z Values of the field \f$ Z \f$ at GLL points.
  function interpolator_jacobian(this, rst, X,Y,Z) result(jac)
    class(interpolator_t), intent(in) :: this
    type(point_t), intent(in) :: rst
    real(kind=rp), intent(inout) :: X(this%xh%lx*this%xh%ly*this%xh%lz)
    real(kind=rp), intent(inout) :: Y(this%xh%lx*this%xh%ly*this%xh%lz)
    real(kind=rp), intent(inout) :: Z(this%xh%lx*this%xh%ly*this%xh%lz)

    real(kind=rp) :: jac(3,3)

    real(kind=rp) :: hr(this%xh%lx, 2), hs(this%xh%ly, 2), ht(this%xh%lz, 2)
    integer :: lx, ly, lz
    lx = this%xh%lx
    ly = this%xh%ly
    lz = this%xh%lz

    ! Weights
    call fd_weights_full(rst%x(1), this%xh%zg(:,1), lx-1, 1, hr)
    call fd_weights_full(rst%x(2), this%xh%zg(:,2), ly-1, 1, hs)
    call fd_weights_full(rst%x(3), this%xh%zg(:,3), lz-1, 1, ht)

    ! d(x,y,z)/dr
    call tensor_scalar3(jac(1,:), X, Y, Z, lx, hr(:,2), hs(:,1), ht(:,1))

    ! d(x,y,z)/ds
    call tensor_scalar3(jac(2,:), X, Y, Z, lx, hr(:,1), hs(:,2), ht(:,1))

    ! d(x,y,z)/dt
    call tensor_scalar3(jac(3,:), X, Y, Z, lx, hr(:,1), hs(:,1), ht(:,2))
  end function interpolator_jacobian


  !> Interpolates an array to one of Xh or Yh on host.
  !! @param x Original array.
  !! @param y Interpolated array.
  !! @param nel Number of elements in the mesh.
  !! @param to_space The space to interpolate to, must be either Xh or Yh.
  !! @note Not optimized for performance, should only be used during init.
  subroutine interpolator_map_host(this, y, x, nel, to_space)
    class(interpolator_t), intent(inout) :: this
    integer :: nel
    type(space_t) :: to_space
    real(kind=rp), intent(inout) :: x(1,nel)
    real(kind=rp), intent(inout) :: y(1,nel)
    if (to_space .eq. this%Yh) then
       call tnsr3d_cpu(y, this%Yh%lx, x, &
                   this%Xh%lx,this%Yh_to_XhT, &
                   this%Yh_to_Xh, this%Yh_to_Xh, nel)
    else if (to_space .eq. this%Xh) then
       call tnsr3d_cpu(y, this%Xh%lx, x, &
                   this%Yh%lx,this%Yh_to_Xh, &
                   this%Yh_to_XhT, this%Yh_to_XhT, nel)
    else
       call neko_error('Invalid interpolation')
    end if
  end subroutine interpolator_map_host


end module interpolation
