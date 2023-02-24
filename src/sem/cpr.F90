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
!> Compression
module cpr
  use gather_scatter
  use neko_config
  use num_types
  use field  
  use space  
  use math
  use mesh
  use coefs
  use tensor
  use mxm_wrapper
  use speclib
  use device
  use utils
  use device_cpr

  use, intrinsic :: iso_c_binding
  implicit none

  !> include information needed for compressing fields
  type :: cpr_t
     real(kind=rp), allocatable :: v(:,:) !< Transformation matrix
     real(kind=rp), allocatable :: vt(:,:) !< Transformation matrix transposed
     real(kind=rp), allocatable :: vinv(:,:) !< Transformation matrix inversed 
     real(kind=rp), allocatable :: vinvt(:,:) !< Transformation matrix
     !! inversed and transposed 
     real(kind=rp), allocatable :: w(:,:) !< Diagonal matrix with weights
     real(kind=rp), allocatable :: specmat(:,:) !< Transformation matrix
     real(kind=rp), allocatable :: specmatt(:,:) !< Transformation matrix
     type(field_t), pointer :: fld  => null()
     type(field_t) :: fldhat
     type(field_t) :: wk
     type(space_t), pointer :: Xh => null()
     type(mesh_t), pointer :: msh => null()
     type(dofmap_t), pointer :: dof => null()

     !
     ! Device pointers (if present)
     !
     type(c_ptr) :: v_d = C_NULL_PTR
     type(c_ptr) :: vt_d = C_NULL_PTR
     type(c_ptr) :: vinv_d = C_NULL_PTR
     type(c_ptr) :: vinvt_d = C_NULL_PTR
     type(c_ptr) :: w_d = C_NULL_PTR
     type(c_ptr) :: specmat_d = C_NULL_PTR
     type(c_ptr) :: specmatt_d = C_NULL_PTR

  end type cpr_t

  interface cpr_init
     module procedure cpr_init_all
  end interface cpr_init

  public :: cpr_init, cpr_free

contains

  !> Initialize cpr
  subroutine cpr_init_all(cpr, u)
    type(cpr_t), intent(inout) :: cpr
    type(field_t), intent(in), target :: u

    call cpr_free(cpr)

    cpr%fld => u
    cpr%fldhat = u
    cpr%wk = u
    cpr%msh => u%msh
    cpr%Xh => u%Xh
    cpr%dof => u%dof


    ! Allocate arrays
    allocate(cpr%v(cpr%Xh%lx, cpr%Xh%lx))
    allocate(cpr%vt(cpr%Xh%lx, cpr%Xh%lx))
    allocate(cpr%vinv(cpr%Xh%lx, cpr%Xh%lx))
    allocate(cpr%vinvt(cpr%Xh%lx, cpr%Xh%lx))
    allocate(cpr%w(cpr%Xh%lx, cpr%Xh%lx))
    allocate(cpr%specmat(cpr%Xh%lx, cpr%Xh%lx))
    allocate(cpr%specmatt(cpr%Xh%lx, cpr%Xh%lx))

    ! Initialize all the matrices
    call cpr_generate_specmat(cpr)

    ! Generate the uhat field (legendre coeff)

    call cpr_goto_space(cpr,'spec') !< 'spec' / 'phys'

  end subroutine cpr_init_all


  !> Deallocate coefficients
  subroutine cpr_free(cpr)
    type(cpr_t), intent(inout) :: cpr

    if(allocated(cpr%v)) then
       deallocate(cpr%v)
    end if

    if(allocated(cpr%vt)) then
       deallocate(cpr%vt)
    end if

    if(allocated(cpr%vinv)) then
       deallocate(cpr%vinv)
    end if

    if(allocated(cpr%vinvt)) then
       deallocate(cpr%vinvt)
    end if

    if(allocated(cpr%w)) then
       deallocate(cpr%w)
    end if
    
    if(allocated(cpr%specmat)) then
       deallocate(cpr%specmat)
    end if
    
    if(allocated(cpr%specmatt)) then
       deallocate(cpr%specmatt)
    end if

    call field_free(cpr%fldhat)

    call field_free(cpr%wk)
    
    nullify(cpr%fld)
    nullify(cpr%msh)
    nullify(cpr%Xh)
    nullify(cpr%dof)


    !
    ! Cleanup the device (if present)
    !
    
    if (c_associated(cpr%v_d)) then
       call device_free(cpr%v_d)
    end if

    if (c_associated(cpr%vt_d)) then
       call device_free(cpr%vt_d)
    end if

    if (c_associated(cpr%vinv_d)) then
       call device_free(cpr%vinv_d)
    end if

    if (c_associated(cpr%vinvt_d)) then
       call device_free(cpr%vinvt_d)
    end if

    if (c_associated(cpr%w_d)) then
       call device_free(cpr%w_d)
    end if
    
    if (c_associated(cpr%specmat_d)) then
       call device_free(cpr%specmat_d)
    end if

    if (c_associated(cpr%specmatt_d)) then
       call device_free(cpr%specmatt_d)
    end if

  end subroutine cpr_free


  !> Generate spectral tranform matrices
  subroutine cpr_generate_specmat(cpr)
    type(cpr_t), intent(inout) :: cpr
    real(kind=rp) :: L(0:cpr%Xh%lx-1)
    real(kind=rp) :: delta(cpr%Xh%lx)
    integer :: i, kj, j, j2, kk
    character(len=LOG_SIZE) :: log_buf 

    associate(Xh => cpr%Xh, v=> cpr%v, vt => cpr%vt, &
         vinv => cpr%vinv, vinvt => cpr%vinvt, w => cpr%w)
      ! Get the Legendre polynomials for each point
      ! Then proceed to compose the transform matrix
      kj = 0
      do j = 1, Xh%lx
         L(0) = 1.
         L(1) = Xh%zg(j,1)
         do j2 = 2, Xh%lx-1
            L(j2) = ( (2*j2-1) * Xh%zg(j,1) * L(j2-1) &
                 - (j2-1) * L(j2-2) ) / j2 
         end do
         do kk = 1, Xh%lx
            kj = kj+1
            v(kj,1) = L(KK-1)
         end do
      end do
      
      ! transpose the matrix
      call trsp1(v, Xh%lx) !< non orthogonal wrt weights

      ! Calculate the nominal scaling factors
      do i = 1, Xh%lx
         delta(i) = 2.0_rp / (2*(i-1)+1)
      end do
      ! modify last entry  
      delta(Xh%lx) = 2.0_rp / (Xh%lx-1)
      
      ! calculate the inverse to multiply the matrix
      do i = 1, Xh%lx
         delta(i) = sqrt(1.0_rp / delta(i))
      end do
      ! scale the matrix      
      do i = 1, Xh%lx
         do j = 1, Xh%lx
            v(i,j) = v(i,j) * delta(j) ! orthogonal wrt weights
         end do
      end do
    
      ! get the trasposed
      call copy(vt, v, Xh%lx * Xh%lx)
      call trsp1(vt, Xh%lx)
      
      !populate the mass matrix
      kk = 1
      do i = 1, Xh%lx
         do j = 1, Xh%lx
            if (i .eq. j) then
               w(i,j) = Xh%wx(kk)
               kk = kk+1
            else
               cpr%w(i,j) = 0
            end if
         end do
      end do
      
      !Get the inverse of the transform matrix
      call mxm(vt, Xh%lx, w, Xh%lx, vinv, Xh%lx)
      
      !get the transposed of the inverse
      call copy(vinvt, vinv, Xh%lx * Xh%lx)
      call trsp1(vinvt, Xh%lx)
    end associate

    ! Copy the data to the GPU
    ! Move all this to space.f90 to for next version 
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
       (NEKO_BCKND_OPENCL .eq. 1)) then 
       call device_map(cpr%v,     cpr%v_d,     cpr%Xh%lxy)
       call device_map(cpr%vt,    cpr%vt_d,    cpr%Xh%lxy)
       call device_map(cpr%vinv,  cpr%vinv_d,  cpr%Xh%lxy)
       call device_map(cpr%vinvt, cpr%vinvt_d, cpr%Xh%lxy)
       call device_map(cpr%w,     cpr%w_d,     cpr%Xh%lxy)
       !Map the following pointers but do not copy data for them
       call device_map(cpr%specmat,  cpr%specmat_d,  cpr%Xh%lxy)
       call device_map(cpr%specmatt, cpr%specmatt_d, cpr%Xh%lxy)


       call device_memcpy(cpr%v,     cpr%v_d,     cpr%Xh%lxy, &
                          HOST_TO_DEVICE)
       call device_memcpy(cpr%vt,    cpr%vt_d,    cpr%Xh%lxy, &
                          HOST_TO_DEVICE)
       call device_memcpy(cpr%vinv,  cpr%vinv_d,  cpr%Xh%lxy, &
                          HOST_TO_DEVICE)
       call device_memcpy(cpr%vinvt, cpr%vinvt_d, cpr%Xh%lxy, &
                          HOST_TO_DEVICE)
       call device_memcpy(cpr%w,     cpr%w_d,     cpr%Xh%lxy, &
                          HOST_TO_DEVICE)

    end if

  end subroutine cpr_generate_specmat


  !> Tranform to spectral space (using tensor product)
  !the result of the transform is given in fldhat
  subroutine cpr_goto_space(cpr, space)
    type(cpr_t), intent(inout) :: cpr
    integer :: i, j, k, e, nxyz, nelv, n
    character(len=LOG_SIZE) :: log_buf 
    character(len=4) :: space 

    ! define some constants
    nxyz = cpr%Xh%lx*cpr%Xh%lx*cpr%Xh%lx
    nelv = cpr%msh%nelv
    n    = nxyz*nelv

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
       (NEKO_BCKND_OPENCL .eq. 1)) then 

       write(*,*) 'Transform in the GPU'

       ! Define the matrix according to which transform to do 
       if (space .eq. 'spec') then
          call device_copy(cpr%specmat_d,  cpr%vinv_d,  cpr%Xh%lxy)
          call device_copy(cpr%specmatt_d, cpr%vinvt_d, cpr%Xh%lxy)
          call device_copy(cpr%wk%x_d, cpr%fld%x_d, n)
       endif
       if (space .eq. 'phys') then
          call device_copy(cpr%specmat_d,  cpr%v_d,  cpr%Xh%lxy)
          call device_copy(cpr%specmatt_d, cpr%vt_d, cpr%Xh%lxy)
          call device_copy(cpr%wk%x_d, cpr%fldhat%x_d, n)
       endif

    else
       
       write(*,*) 'Transform in the CPU'

       ! Define the matrix according to which transform to do 
       if (space .eq. 'spec') then
          call copy(cpr%specmat, cpr%vinv, cpr%Xh%lx*cpr%Xh%lx)
          call copy(cpr%specmatt, cpr%vinvt, cpr%Xh%lx*cpr%Xh%lx)
          call copy(cpr%wk%x,cpr%fld%x,n)
       endif
       if (space .eq. 'phys') then
          call copy(cpr%specmat, cpr%v, cpr%Xh%lx*cpr%Xh%lx)
          call copy(cpr%specmatt, cpr%vt, cpr%Xh%lx*cpr%Xh%lx)
          call copy(cpr%wk%x,cpr%fldhat%x,n)
       endif

    end if

    call tnsr3d(cpr%fldhat%x, cpr%Xh%lx, cpr%wk%x, &
                cpr%Xh%lx,cpr%specmat, &
                cpr%specmatt, cpr%specmatt, nelv)

  end subroutine cpr_goto_space

  !> Truncate the frequencies
  subroutine cpr_truncate_wn(cpr, coef)
    type(cpr_t), intent(inout) :: cpr
    type(coef_t), intent(inout) :: coef
    real(kind=rp) :: w2(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: w1(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: vsort(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: vtrunc(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: vtemp(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: errvec(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx) 
    real(kind=rp) :: fx(cpr%Xh%lx, cpr%Xh%lx) 
    real(kind=rp) :: fy(cpr%Xh%lx, cpr%Xh%lx) 
    real(kind=rp) :: fz(cpr%Xh%lx, cpr%Xh%lx) 
    real(kind=rp) :: ident(cpr%Xh%lx, cpr%Xh%lx) 
    real(kind=rp) :: l2norm, oldl2norm, targeterr
    integer :: i, j, k, e, nxyz, nelv, n
    integer :: targetkut, comp_counter
    real(kind=rp) :: restemp,res, restemp2, res2, vol2
    integer :: kut, kutx, kuty, kutz, nx
    character(len=LOG_SIZE) :: log_buf 

    ! For sorting
    integer :: isort(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: fldhat_sort(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx, cpr%msh%nelv)
    integer :: key(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx, cpr%msh%nelv)
    integer :: dummykey(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx, cpr%msh%nelv)
    integer :: key_sort(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx, cpr%msh%nelv)



    ! Definitions for the l2norm in device
    real(kind=rp) :: ev(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx,cpr%msh%nelv) 
    real(kind=rp) :: tmp(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx,cpr%msh%nelv) 
    real(kind=rp) :: res_e(cpr%msh%nelv) 
    real(kind=rp) :: vol_e(cpr%msh%nelv) 
    real(kind=rp) :: compress_e(cpr%msh%nelv) 
    
     !
     ! Device pointers (if present)
     !
     type(c_ptr) :: fx_d = C_NULL_PTR
     type(c_ptr) :: fy_d = C_NULL_PTR
     type(c_ptr) :: fz_d = C_NULL_PTR
     type(c_ptr) :: ident_d = C_NULL_PTR
     type(c_ptr) :: ev_d = C_NULL_PTR
     type(c_ptr) :: tmp_d = C_NULL_PTR
     type(c_ptr) :: res_e_d = C_NULL_PTR
     type(c_ptr) :: vol_e_d = C_NULL_PTR
     type(c_ptr) :: compress_e_d = C_NULL_PTR
      
     type(c_ptr) :: fldhat_sort_d = C_NULL_PTR
     type(c_ptr) :: key_d = C_NULL_PTR
     type(c_ptr) :: dummykey_d = C_NULL_PTR
     type(c_ptr) :: key_sort_d = C_NULL_PTR


    ! define some constants
    nx = cpr%Xh%lx
    nxyz = cpr%Xh%lx*cpr%Xh%lx*cpr%Xh%lx
    nelv = cpr%msh%nelv
    n    = nxyz*nelv
    targeterr = 1e-3_rp
    targetkut = 4

    ! Initialize a vector that tells wich element to truncate
    ! Initially truncate all the elemtns so all entries set to 1
    do i= 1, nelv
       compress_e(i) = 1
    end do

    ! fill up key sort
    do e= 1, nelv
       do i = 1, nxyz
         key(i,1,1,e)= i
         key_sort(i,1,1,e)= 0
       end do
    end do


    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
       (NEKO_BCKND_OPENCL .eq. 1)) then 
       

       write(*,*) 'GPU truncation'

       ! Map the pointers
       call device_map(fx,  fx_d,  cpr%Xh%lxy)
       call device_map(fy,  fy_d,  cpr%Xh%lxy)
       call device_map(fz,  fz_d,  cpr%Xh%lxy)
       call device_map(ident,  ident_d,  cpr%Xh%lxy)
       call device_map(ev,  ev_d,  n)
       call device_map(tmp,  tmp_d,  n)
       call device_map(res_e,  res_e_d,  nelv)
       call device_map(vol_e,  vol_e_d,  nelv)
       call device_map(compress_e,  compress_e_d,  nelv)
       
       
       call device_map(fldhat_sort,  fldhat_sort_d,  n)
       call device_map(key,  key_d,  n)
       call device_map(dummykey, dummykey_d,  n)
       call device_map(key_sort,  key_sort_d,  n)


       ! Transfer the key to the device
       call device_memcpy(key,     key_d,  n, &
                          HOST_TO_DEVICE)

       comp_counter=2 
       kut = 0

       !Copy the untruncated coefficients to a temporal array
       call device_copy(tmp_d, cpr%fldhat%x_d, n)
       
       ! Sort coefficients before truncating them
       call device_lcsort_abs(fldhat_sort_d,key_sort_d,tmp_d, &
                              key_d, cpr%Xh%lx,nelv)
       
       do while (comp_counter .ge. 1 .and. kut .le. (cpr%Xh%lx-1)*3)
             
          ! advance kut
          kut = kut+1
      
          ! create identity matrices in cpu (bypass for truncation)
          call build_filter_tf(fx, fy, fz, 0, cpr%Xh%lx)
          call copy(ident, fx, cpr%Xh%lxy)

          ! create filters on cpu
          call build_filter_tf(fx, fy, fz, kut, cpr%Xh%lx)

          ! Transfer the filter to the GPU
          call device_memcpy(fx,     fx_d,     cpr%Xh%lxy, &
                             HOST_TO_DEVICE)
          call device_memcpy(fy,     fy_d,     cpr%Xh%lxy, &
                             HOST_TO_DEVICE)
          call device_memcpy(fz,     fz_d,     cpr%Xh%lxy, &
                             HOST_TO_DEVICE)
          call device_memcpy(ident,     ident_d,     cpr%Xh%lxy, &
                             HOST_TO_DEVICE)

          ! Transfer compression keys to GPU (keys says which element is to be compressed)
          call device_memcpy(compress_e,     compress_e_d, nelv, &
                             HOST_TO_DEVICE)

          ! Copy the spectral coefficients to the working array in GPU
          call device_copy(cpr%wk%x_d, fldhat_sort_d, n)

          ! Perform the truncation in the GPU
          !call tnsr3d(cpr%fldhat%x, cpr%Xh%lx, cpr%wk%x, &
          !            cpr%Xh%lx,fx, &
          !            fy, fz, nelv)

          ! Perform the truncation in the GPU
          call device_lctnsr3d(fldhat_sort_d, cpr%Xh%lx, cpr%wk%x_d, &
                               cpr%Xh%lx,fx_d, &
                               fy_d, fz_d, ident_d, compress_e_d, nelv)

          ! Unsort the coefficients before evaluating norms       
          call device_lcsort_bykey(cpr%fldhat%x_d,dummykey_d,&
                                   fldhat_sort_d, &
                                   key_sort_d, cpr%Xh%lx,nelv)

          ! Get error vector 
          call device_sub3(ev_d,cpr%fldhat%x_d,tmp_d,n)

          ! Get local inner products 
          call device_lcsc3(res_e_d,ev_d,coef%jac_d, &
                            ev_d, cpr%Xh%lx,nelv) 

          ! Get local volumes
          call device_lcsum(vol_e_d,coef%B_d, &
                            cpr%Xh%lx,nelv) 

          ! Send the local values to the host
          call device_memcpy(res_e, res_e_d,     nelv, &
                             DEVICE_TO_HOST)

          call device_memcpy(vol_e, vol_e_d,     nelv, &
                             DEVICE_TO_HOST)
              
          !evaluate the element error and decide if shoudl keep compressing
          res2 = 0
          comp_counter=0       
          do e= 1, nelv
             res2 = sqrt(res_e(e))/sqrt(vol_e(e))
             if (res2.ge.targeterr) then
                compress_e(e)=0
             else
                compress_e(e)=1
             end if
             comp_counter=comp_counter+compress_e(e)   
          end do

       
       end do ! This is the end do of the while loop

       !================debugging info
       !call device_memcpy(cpr%fldhat%x, cpr%fldhat%x_d,     n, &
       !                   DEVICE_TO_HOST)
       !call device_memcpy(key_sort, key_sort_d,     n, &
       !                   DEVICE_TO_HOST)
       !write(*,*) 'Checking sorted back coefficients'
       !do i = 1, 10
       !  write(log_buf, '(A,E15.7,A,I5)') &
       !        'u hat value:', cpr%fldhat%x(i,1,1,10), &
       !        ' index:', &
       !        key_sort(i,1,1,10)
       !  call neko_log%message(log_buf)
       !enddo
       !do i = 500, nxyz
       !  write(log_buf, '(A,E15.7,A,I5)') &
       !        'u hat value:', cpr%fldhat%x(i,1,1,10), &
       !        ' index:', &
       !        key_sort(i,1,1,10)
       !  call neko_log%message(log_buf)
       !===============


       ! Print the final l2 norm of the error
       ! Get global inner product 
       restemp = device_glsc3(ev_d,coef%jac_d,ev_d,n) 
       ! Get the general l2norm 
       res= sqrt(restemp)/sqrt(coef%volume)
       write(*,*) 'l2 norm of the error is: ', res


       ! Free memory after performing actions
       if (c_associated(fx_d)) then
          call device_free(fx_d)
       end if
       if (c_associated(fy_d)) then
          call device_free(fy_d)
       end if
       if (c_associated(fz_d)) then
          call device_free(fz_d)
       end if
       if (c_associated(ev_d)) then
          call device_free(ev_d)
       end if
       if (c_associated(tmp_d)) then
          call device_free(tmp_d)
       end if
       if (c_associated(res_e_d)) then
          call device_free(res_e_d)
       end if
       if (c_associated(vol_e_d)) then
          call device_free(vol_e_d)
       end if
       if (c_associated(ident_d)) then
          call device_free(ident_d)
       end if
       if (c_associated(compress_e_d)) then
          call device_free(compress_e_d)
       end if
       
       if (c_associated(fldhat_sort_d)) then
          call device_free(fldhat_sort_d)
       end if
       if (c_associated(key_d)) then
          call device_free(key_d)
       end if
       if (c_associated(dummykey_d)) then
          call device_free(dummykey_d)
       end if
       if (c_associated(key_sort_d)) then
          call device_free(key_sort_d)
       end if
 
       !if(allocated(temp)) then
       !   deallocate(temp)
       !end if
       !if(allocated(ev)) then
       !   deallocate(ev)
       !end if

    else

       !truncate for every element
       do e = 1, nelv
          ! create temp vector where the trunc coeff will be
          call copy(vtemp, cpr%fldhat%x(1,1,1,e), nxyz)
          call copy(vtrunc, cpr%fldhat%x(1,1,1,e), nxyz)
          ! sort the coefficients by absolute value
          call sortcoeff(vsort, cpr%fldhat%x(1,1,1,e), isort, nxyz) 
          ! initialize values for iterative procedure
          l2norm = 0.0_rp
          kut = 0

          do while (l2norm .le. targeterr .and. kut .le. (cpr%Xh%lx-1)*3)
             ! save value from prev it since it met requirements to enter
             call copy(vtrunc, vtemp, nxyz)
             oldl2norm = l2norm

             ! advance kut
             kut = kut+1

             ! create filters
             call build_filter_tf(fx, fy, fz, kut, cpr%Xh%lx)

             ! truncate, remember that transposed of diag mat is the same
             call copy(w2, vsort, nxyz)
             ! apply the operator to the x direction, result in w1
             call mxm(fx, nx, w2, nx, w1, nx*nx)
             ! apply matrix to y direction, result in w2
             do k=1,nx
                call mxm(w1(1,1,k), nx, fy, nx, w2(1,1,k), nx)
             enddo
             ! apply matrix to z direction, result in vtemp
             call mxm (w2, nx*nx, fz, nx, vtemp, nx)

             ! vtemp is sorted, so bring it back to real order
             ! by inverting the swap operation in sortcoeff
             call reord(vtemp, isort, nxyz)

             ! calculate the error vector
             call sub3(errvec, cpr%fldhat%x(1,1,1,e), vtemp, nxyz)

             ! get the norm of the error
             l2norm = get_elem_l2norm(errvec, coef, 'spec' ,e)

          end do

          ! copy the truncated field back to the main object
          call copy(cpr%fldhat%x(1,1,1,e), vtrunc, nxyz)


          !================debugging info

          !just to check that all is good, resort vsort
          !if (e .eq. 50) then
          !   write(log_buf, '(A)') &
          !        'Debugging info for e=50. Do not forget to delete'
          !   call neko_log%message(log_buf)

          !   call reord(vsort,isort,nxyz)

          !   write(log_buf, '(A)') &
          !     'Norm'
          !   call neko_log%message(log_buf)
          !   !chech that the copy is fine in one entry
          !   write(log_buf, '(A,E15.7)') &
          !        'l2norm:', oldl2norm
          !   call neko_log%message(log_buf)
          !   write(log_buf, '(A)') &
          !        'cut off level'
          !   call neko_log%message(log_buf)
          !   !chech that the copy is fine in one entry
          !   write(log_buf, '(A,I5)') &
          !        'kut:', kut
          !   call neko_log%message(log_buf)
          !   write(log_buf, '(A)') &
          !        'spectral coefficients'
          !   call neko_log%message(log_buf)
          !   do i = 1, 10
          !      write(log_buf, '(A,E15.7,A,E15.7)') &
          !           'org coeff:', vsort(i,1,1), ' truncated coeff', &
          !           cpr%fldhat%x(i,1,1,e)
          !      call neko_log%message(log_buf)
          !   end do
          !end if
       end do

    end if

  end subroutine cpr_truncate_wn

  !> Sort the spectral coefficient in descending order 
  !! array vsort. The original indices are stored in the isort vector.
  subroutine sortcoeff(vsort, v, isort, nxyz) 
    integer, intent(in) :: nxyz      
    real(kind=rp), intent(inout) :: vsort(nxyz)
    real(kind=rp), intent(inout) :: v(nxyz)
    integer, intent(inout) :: isort(nxyz)
    real(kind=rp) :: wrksort(nxyz)
    integer :: wrkisort(nxyz)
    integer :: i, j, k

    ! copy absolute values to sort by magnitude
    do i =1, nxyz
       vsort(i) = abs(v(i))
       isort(i) = i
    end do

    ! sort the absolute values of the vectors, here the index is the
    ! important part (isort)
    call sort(vsort, isort, nxyz)

    ! Flip the indices so they are in a descending order
    call flipv(vsort, isort, nxyz)

    ! now re order the fields (non abs value) based on the indices
    call copy(vsort, v, nxyz)
    call swap(vsort, isort, nxyz)

  end subroutine sortcoeff

  !> Flip vector b and ind 
  subroutine flipv(b, ind, n)
    integer, intent(in) :: n      
    real(kind=rp), intent(inout) :: b(n)
    integer, intent(inout) :: ind(n)
    real(kind=rp) :: temp(n)
    integer :: tempind(n)
    integer :: i, jj

    do i = 1, n
       jj = n+1-i
       temp(jj) = b(i)
       tempind(jj) = ind(i)
    end do
    do i = 1,n
       b(i) = temp(i)
       ind(i) = tempind(i)
    end do

  end subroutine flipv

  !> sort the array acording to ind vector
  subroutine swap(b, ind, n)
    integer, intent(in) :: n      
    real(kind=rp), intent(inout) :: b(n)
    integer, intent(inout) :: ind(n)
    real(kind=rp) :: temp(n)
    integer :: i, jj

    do i = 1, n
       temp(i)=b(i)
    end do
    do i = 1, n
       jj=ind(i)
       b(i)=temp(jj)
    end do

  end subroutine swap

  !> reorder the array - inverse of swap
  subroutine reord(b, ind, n)
    integer, intent(in) :: n      
    real(kind=rp), intent(inout) :: b(n)
    integer, intent(inout) :: ind(n)
    real(kind=rp) :: temp(n)
    integer :: i, jj

    do i = 1, n
       temp(i) = b(i)
    end do
    do i = 1, n
       jj = ind(i)
       b(jj) = temp(i)
    end do

  end subroutine reord

  !> create filter transfer function 
  subroutine build_filter_tf(fx, fy, fz, kut, lx) 
    integer, intent(in) :: lx      
    integer, intent(in) :: kut      
    real(kind=rp), intent(inout) :: fx(lx,lx)
    real(kind=rp), intent(inout) :: fy(lx,lx)
    real(kind=rp), intent(inout) :: fz(lx,lx)
    integer :: i, j, k, k0, kk, kutx, kuty, kutz

    ! set the values acording to kut
    if (kut .eq. 0) then
       kutz = 0
       kuty = 0
       kutx = 0
    end if
    if (kut .gt. 0 .and. kut .le. (lx-1)) then
       kutz = kut
       kuty = 0
       kutx = 0
    end if
    if (kut .gt. (lx-1) .and. kut .le. (lx-1)*2) then
       kutz = lx-1
       kuty = kut-(lx-1)
       kutx = 0
    end if
    if (kut .gt. (lx-1)*2 .and. kut .le. (lx-1)*3) then
       kutz = lx-1
       kuty = lx-1
       kutx = kut-(lx-1)*2
    end if

    ! create identity matrices
    do i = 1, lx
       do j = 1, lx
          if (i .eq. j) then
             fx(i,j) = 1
             fy(i,j) = 1
             fz(i,j) = 1
          else
             fx(i,j) = 0
             fy(i,j) = 0
             fz(i,j) = 0
          end if
       end do
    end do

    ! truncate the matrices acording to input kut
    k0 = lx-kutx
    do k = k0+1, lx
       kk = k+lx*(k-1)
       fx(kk,1) = 0
    end do
    k0 = lx-kuty
    do k = k0+1, lx
       kk = k+lx*(k-1)
       fy(kk,1) = 0
    end do
    k0 = lx-kutz
    do k = k0+1, lx
       kk = k+lx*(k-1)
       fz(kk,1) = 0
    end do

  end subroutine build_filter_tf

  !< Get l2 norm of the element
  function get_elem_l2norm(elemdata, coef, space, e) result(l2e)
    type(coef_t), intent(in) :: coef
    real(kind=rp) :: elemdata(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx)
    real(kind=rp) :: vole, suma, l2e
    integer i, e, eg, nxyz
    character(len=4) :: space 

    ! Get the volume of the element
    nxyz = coef%Xh%lx*coef%Xh%lx*coef%Xh%lx
    vole=0
    do i = 1, nxyz
       vole = vole+coef%B(i,1,1,e)
    end do

    ! Get the weighted l2 norm of the element
    suma = 0
    if (space .eq. 'spec') then
       do i = 1, nxyz
          suma = suma+elemdata(i,1,1)*elemdata(i,1,1)*coef%jac(i,1,1,e)
       end do
    end if
    if (space .eq. 'phys') then
       do i = 1, nxyz
          suma = suma+elemdata(i,1,1)*elemdata(i,1,1)*coef%B(i,1,1,e)
       end do
    end if

    l2e = sqrt(suma)/sqrt(vole)

  end function get_elem_l2norm


end module cpr
