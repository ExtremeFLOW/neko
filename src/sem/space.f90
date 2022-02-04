! Copyright (c) 2019-2022, The Neko Authors
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
!> Defines a function space
module space
  use neko_config
  use num_types
  use speclib
  use device
  use utils
  use fast3d
  use math
  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter :: GL = 0, GLL = 1, GJ = 2

  type space_t
     integer :: t               !< Space type (GL, GLL, GJ, ...)
     integer :: lx              !< Polynomial dimension in x-direction
     integer :: ly              !< Polynomial dimension in y-direction
     integer :: lz              !< Polynomial dimension in z-direction
     integer :: lxy             !< Number of points in xy-plane
     integer :: lyz             !< Number of points in yz-plane
     integer :: lxz             !< Number of points in xz-plane
     integer :: lxyz            !< Number of points in xyz-block
     
     real(kind=rp), allocatable :: zg(:,:) !< Quadrature points
     
     real(kind=rp), allocatable :: dr_inv(:) !< 1/dist quadrature points
     real(kind=rp), allocatable :: ds_inv(:) !< 1/dist quadrature points
     real(kind=rp), allocatable :: dt_inv(:) !< 1/dist quadrature points

     real(kind=rp), allocatable :: wx(:)   !< Quadrature weights
     real(kind=rp), allocatable :: wy(:)   !< Quadrature weights
     real(kind=rp), allocatable :: wz(:)   !< Quadrature weights

     real(kind=rp), allocatable :: w3(:,:,:)

     !> Derivative operator \f$ D_1 \f$
     real(kind=rp), allocatable :: dx(:,:)
     !> Derivative operator \f$ D_2 \f$
     real(kind=rp), allocatable :: dy(:,:)
     !> Derivative operator \f$ D_3 \f$
     real(kind=rp), allocatable :: dz(:,:)

     !> Transposed derivative operator \f$ D_1^T \f$
     real(kind=rp), allocatable :: dxt(:,:)
     !> Transposed derivative operator \f$ D_2^T \f$
     real(kind=rp), allocatable :: dyt(:,:)
     !> Transposed derivative operator \f$ D_3^T \f$
     real(kind=rp), allocatable :: dzt(:,:)

     !
     ! Device pointers (if present)
     !
     type(c_ptr) :: dr_inv_d = C_NULL_PTR
     type(c_ptr) :: ds_inv_d = C_NULL_PTR
     type(c_ptr) :: dt_inv_d = C_NULL_PTR
     type(c_ptr) :: dxt_d = C_NULL_PTR
     type(c_ptr) :: dyt_d = C_NULL_PTR
     type(c_ptr) :: dzt_d = C_NULL_PTR
     type(c_ptr) :: dx_d = C_NULL_PTR
     type(c_ptr) :: dy_d = C_NULL_PTR
     type(c_ptr) :: dz_d = C_NULL_PTR
     type(c_ptr) :: wx_d = C_NULL_PTR
     type(c_ptr) :: wy_d = C_NULL_PTR
     type(c_ptr) :: wz_d = C_NULL_PTR
     type(c_ptr) :: zg_d = C_NULL_PTR
     type(c_ptr) :: w3_d = C_NULL_PTR

  end type space_t

  interface operator(.eq.)
     module procedure space_eq
  end interface operator(.eq.)

  interface operator(.ne.)
     module procedure space_ne
  end interface operator(.ne.)
  
contains

  !> Initialize a function space @a s with given polynomial dimensions
  subroutine space_init(s, t, lx, ly, lz)
    type(space_t), intent(inout) :: s
    integer, intent(in) :: t            !< Quadrature type
    integer, intent(in) :: lx           !< Polynomial dimension in x-direction
    integer, intent(in) :: ly           !< Polynomial dimension in y-direction
    integer, optional, intent(in) :: lz !< Polynomial dimension in z-direction
    integer :: ix, iy, iz
 
    call space_free(s)

    s%lx = lx
    s%ly = ly
    s%t = t
    if (present(lz)) then
       if (lz .ne. 1) then
          s%lz = lz
          if (lx .ne. ly .or. lx .ne. lz) then
             call neko_error("Unsupported polynomial dimension")
          end if
       end if
    else
       if (lx .ne. ly) then
          call neko_error("Unsupported polynomial dimension")
       end if
       s%lz = 1
    end if
    s%lxy = s%ly*s%lx
    s%lyz = s%ly*s%lz
    s%lxz = s%lx*s%lz
    s%lxyz = s%lx*s%ly*s%lz

    allocate(s%zg(lx, 3))

    allocate(s%wx(s%lx))
    allocate(s%wy(s%ly))
    allocate(s%wz(s%lz))
    
    allocate(s%dr_inv(s%lx))
    allocate(s%ds_inv(s%ly))
    allocate(s%dt_inv(s%lz))

    allocate(s%w3(s%lx, s%ly, s%lz))

    allocate(s%dx(s%lx, s%lx))
    allocate(s%dy(s%ly, s%ly))
    allocate(s%dz(s%lz, s%lz))

    allocate(s%dxt(s%lx, s%lx))
    allocate(s%dyt(s%ly, s%ly))
    allocate(s%dzt(s%lz, s%lz))
    
    if (t .eq. GLL) then
       call zwgll(s%zg(1,1), s%wx, s%lx)
       call zwgll(s%zg(1,2), s%wy, s%ly)
       if (s%lz .gt. 1) then
          call zwgll(s%zg(1,3), s%wz, s%lz)
       else
          s%zg(:,3) = 0d0
          s%wz = 1d0
       end if
    else if (t .eq. GL) then
       call zwgl(s%zg(1,1), s%wx, s%lx)
       call zwgl(s%zg(1,2), s%wy, s%ly)
       if (s%lz .gt. 1) then
          call zwgl(s%zg(1,3), s%wz, s%lz)
       else
          s%zg(:,3) = 0d0
          s%wz = 1d0
       end if
    else
       call neko_error("Invalid quadrature rule")
    end if

    do iz = 1, s%lz
       do iy = 1, s%ly
          do ix = 1, s%lx
             s%w3(ix, iy, iz) = s%wx(ix) * s%wy(iy) * s%wz(iz)
          end do
       end do
    end do
    !> Setup derivative matrices
    if (t .eq. GLL) then
        call dgll(s%dx, s%dxt, s%zg(1,1), s%lx, s%lx)
        call dgll(s%dy, s%dyt, s%zg(1,2), s%ly, s%ly)
        if (s%lz .gt. 1) then
           call dgll(s%dz, s%dzt, s%zg(1,3), s%lz, s%lz)
        else
           s%dz = 0d0
           s%dzt = 0d0
        end if
    else if (t .eq. GL) then
       call setup_intp(s%dx,s%dxt,s%zg(1,1),s%zg(1,1),s%lx,s%lx,1)
       call setup_intp(s%dy,s%dyt,s%zg(1,2),s%zg(1,2),s%ly,s%ly,1)
        if (s%lz .gt. 1) then
           call setup_intp(s%dz,s%dzt,s%zg(1,3),s%zg(1,3),s%lz,s%lz,1)
        else
           s%dz = 0d0
           s%dzt = 0d0
        end if
     else
        call neko_error("Invalid quadrature rule")
     end if
    
    call space_compute_dist(s%dr_inv, s%zg(1,1), s%lx)
    call space_compute_dist(s%ds_inv, s%zg(1,2), s%ly)
    if (s%lz .gt. 1) then
       call space_compute_dist(s%dt_inv, s%zg(1,3), s%lz)
    else
       s%dt_inv = 0d0
    end if

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
        (NEKO_BCKND_OPENCL .eq. 1)) then 
       call device_map(s%dr_inv, s%dr_inv_d, s%lx)
       call device_map(s%ds_inv, s%ds_inv_d, s%lx)
       call device_map(s%dt_inv, s%dt_inv_d, s%lx)
       call device_map(s%wx, s%wx_d, s%lx)
       call device_map(s%wy, s%wy_d, s%lx)
       call device_map(s%wz, s%wz_d, s%lx)
       call device_map(s%dx, s%dx_d, s%lxy)
       call device_map(s%dy, s%dy_d, s%lxy)
       call device_map(s%dz, s%dz_d, s%lxy)
       call device_map(s%dxt, s%dxt_d, s%lxy)
       call device_map(s%dyt, s%dyt_d, s%lxy)
       call device_map(s%dzt, s%dzt_d, s%lxy)
       call device_map(s%w3, s%w3_d, s%lxyz)

       call device_memcpy(s%dr_inv, s%dr_inv_d, s%lx, HOST_TO_DEVICE)
       call device_memcpy(s%ds_inv, s%ds_inv_d, s%lx, HOST_TO_DEVICE)
       call device_memcpy(s%dt_inv, s%dt_inv_d, s%lx, HOST_TO_DEVICE)
       call device_memcpy(s%wx, s%wx_d, s%lx, HOST_TO_DEVICE)
       call device_memcpy(s%wy, s%wy_d, s%lx, HOST_TO_DEVICE)
       call device_memcpy(s%wz, s%wz_d, s%lx, HOST_TO_DEVICE)
       call device_memcpy(s%dx, s%dx_d, s%lxy, HOST_TO_DEVICE)
       call device_memcpy(s%dy, s%dy_d, s%lxy, HOST_TO_DEVICE)
       call device_memcpy(s%dz, s%dz_d, s%lxy, HOST_TO_DEVICE)
       call device_memcpy(s%dxt, s%dxt_d, s%lxy, HOST_TO_DEVICE)
       call device_memcpy(s%dyt, s%dyt_d, s%lxy, HOST_TO_DEVICE)
       call device_memcpy(s%dzt, s%dzt_d, s%lxy, HOST_TO_DEVICE)
       call device_memcpy(s%w3, s%w3_d, s%lxyz, HOST_TO_DEVICE)

       ix = s%lx * 3
       call device_map(s%zg, s%zg_d, ix)
       call device_memcpy(s%zg, s%zg_d, ix, HOST_TO_DEVICE)
    end if

  end subroutine space_init
   
  !> Deallocate a space @a s
  subroutine space_free(s)
    type(space_t), intent(inout) :: s

    if (allocated(s%zg)) then
       deallocate(s%zg)
    end if

    if (allocated(s%wx)) then
       deallocate(s%wx)
    end if

    if (allocated(s%wy)) then
       deallocate(s%wy)
    end if

    if (allocated(s%wz)) then
       deallocate(s%wz)
    end if

    if (allocated(s%w3)) then
       deallocate(s%w3)
    end if

    if (allocated(s%dx)) then
       deallocate(s%dx)
    end if

    if (allocated(s%dy)) then
       deallocate(s%dy)
    end if

    if (allocated(s%dz)) then
       deallocate(s%dz)
    end if

    if (allocated(s%dxt)) then
       deallocate(s%dxt)
    end if

    if (allocated(s%dyt)) then
       deallocate(s%dyt)
    end if

    if (allocated(s%dzt)) then
       deallocate(s%dzt)
    end if
    
    if (allocated(s%dr_inv)) then
       deallocate(s%dr_inv)
    end if
    
    if (allocated(s%ds_inv)) then
       deallocate(s%ds_inv)
    end if
    
    if (allocated(s%dt_inv)) then
       deallocate(s%dt_inv)
    end if

    !
    ! Cleanup the device (if present)
    !
    
    if (c_associated(s%dr_inv_d)) then
       call device_free(s%dr_inv_d)
    end if

    if (c_associated(s%ds_inv_d)) then
       call device_free(s%ds_inv_d)
    end if

    if (c_associated(s%dt_inv_d)) then
       call device_free(s%dt_inv_d)
    end if

    if (c_associated(s%dxt_d)) then
       call device_free(s%dxt_d)
    end if

    if (c_associated(s%dyt_d)) then
       call device_free(s%dyt_d)
    end if

    if (c_associated(s%dzt_d)) then
       call device_free(s%dzt_d)
    end if
    
    if (c_associated(s%dx_d)) then
       call device_free(s%dx_d)
    end if

    if (c_associated(s%dy_d)) then
       call device_free(s%dy_d)
    end if

    if (c_associated(s%dz_d)) then
       call device_free(s%dz_d)
    end if

    if (c_associated(s%wx_d)) then
       call device_free(s%wx_d)
    end if

    if (c_associated(s%wy_d)) then
       call device_free(s%wy_d)
    end if

    if (c_associated(s%wz_d)) then
       call device_free(s%wz_d)
    end if

    if (c_associated(s%w3_d)) then
       call device_free(s%w3_d)
    end if

    if (c_associated(s%zg_d)) then
       call device_free(s%zg_d)
    end if

  end subroutine space_free

  !> Check if \f$ X_h = Y_H \f$
  !! @note this only checks the polynomial dimensions
  pure function space_eq(Xh, Yh) result(res)
    type(space_t), intent(in) :: Xh
    type(space_t), intent(in) :: Yh
    logical :: res

    if ( (Xh%lx .eq. Yh%lx) .and. &
         (Xh%ly .eq. Yh%ly) .and. &
         (Xh%lz .eq. Yh%lz) ) then
       res = .true.
    else
       res = .false.
    end if
    
  end function space_eq

  !> Check if \f$ X_h \ne Y_H \f$
  !! @note this only checks the polynomial dimensions
  pure function space_ne(Xh, Yh) result(res)
    type(space_t), intent(in) :: Xh
    type(space_t), intent(in) :: Yh
    logical :: res

    if ( (Xh%lx .eq. Yh%lx) .and. &
         (Xh%ly .eq. Yh%ly) .and. &
         (Xh%lz .eq. Yh%lz) ) then
       res = .false.
    else
       res = .true.
    end if
    
  end function space_ne
  
  subroutine space_compute_dist(dx, x, lx)
    integer, intent(in) :: lx
    real(kind=rp), intent(inout) :: dx(lx), x(lx)
    integer :: i
    dx(1) = x(2) - x(1)
    do i = 2, lx - 1
       dx(i) = 0.5*(x(i+1) - x(i-1))
    enddo
    dx(lx) = x(lx) - x(lx-1)
    call invcol1(dx, lx)
  end subroutine space_compute_dist

end module space
