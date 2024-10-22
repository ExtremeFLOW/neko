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
  use field, only : field_t
  use space, only : space_t
  use math
  use mesh, only : mesh_t
  use coefs, only : coef_t
  use tensor
  use mxm_wrapper
  use logger
  use dofmap, only : dofmap_t
  implicit none
  private

  !> include information needed for compressing fields
  type, public  :: cpr_t
     real(kind=rp), allocatable :: v(:,:) !< Transformation matrix

     real(kind=rp), allocatable :: vt(:,:) !< Transformation matrix transposed
     real(kind=rp), allocatable :: vinv(:,:) !< Transformation matrix inversed
     real(kind=rp), allocatable :: vinvt(:,:) !< Transformation matrix
     !! inversed and transposed
     real(kind=rp), allocatable :: w(:,:) !< Diagonal matrix with weights

     real(kind=rp), allocatable :: fldhat(:,:,:,:) !< transformed Field data

     type(field_t), pointer :: fld  => null()
     type(space_t), pointer :: Xh => null()
     type(mesh_t), pointer :: msh => null()
     type(dofmap_t), pointer :: dof => null()

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
    cpr%msh => u%msh
    cpr%Xh => u%Xh
    cpr%dof => u%dof


    ! Allocate arrays
    allocate(cpr%v(cpr%Xh%lx, cpr%Xh%lx))
    allocate(cpr%vt(cpr%Xh%lx, cpr%Xh%lx))
    allocate(cpr%vinv(cpr%Xh%lx, cpr%Xh%lx))
    allocate(cpr%vinvt(cpr%Xh%lx, cpr%Xh%lx))
    allocate(cpr%w(cpr%Xh%lx, cpr%Xh%lx))
    allocate(cpr%fldhat(cpr%Xh%lx, cpr%Xh%ly, cpr%Xh%lz, cpr%msh%nelv))

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

    if(allocated(cpr%fldhat)) then
       deallocate(cpr%fldhat)
    end if

    nullify(cpr%fld)
    nullify(cpr%msh)
    nullify(cpr%Xh)
    nullify(cpr%dof)


  end subroutine cpr_free


  !> Generate spectral tranform matrices
  subroutine cpr_generate_specmat(cpr)
    type(cpr_t), intent(inout) :: cpr
    real(kind=rp) :: L(0:cpr%Xh%lx-1)
    real(kind=rp) :: delta(cpr%Xh%lx)
    integer :: i, kj, j, j2, kk

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

  end subroutine cpr_generate_specmat


  !> Tranform to spectral space (using tensor product)
  !the result of the transform is given in fldhat
  subroutine cpr_goto_space(cpr, space)
    type(cpr_t), intent(inout) :: cpr
    real(kind=rp) :: w2(cpr%Xh%lx,cpr%Xh%lx,cpr%Xh%lx)
    real(kind=rp) :: w1(cpr%Xh%lx,cpr%Xh%lx,cpr%Xh%lx)
    real(kind=rp) :: specmat(cpr%Xh%lx,cpr%Xh%lx)
    real(kind=rp) :: specmatt(cpr%Xh%lx,cpr%Xh%lx)
    integer :: k, e, nxyz, nelv
    character(len=4) :: space

    ! define some constants
    nxyz = cpr%Xh%lx*cpr%Xh%lx*cpr%Xh%lx
    nelv = cpr%msh%nelv

    ! Define the matrix according to which transform to do
    if (space .eq. 'spec') then
       call copy(specmat, cpr%vinv, cpr%Xh%lx*cpr%Xh%lx)
       call copy(specmatt, cpr%vinvt, cpr%Xh%lx*cpr%Xh%lx)
    endif
    if (space .eq. 'phys') then
       call copy(specmat, cpr%v, cpr%Xh%lx*cpr%Xh%lx)
       call copy(specmatt, cpr%vt, cpr%Xh%lx*cpr%Xh%lx)
    endif

    ! Apply the operator (transform to given space)
    do e=1,nelv
       if (space .eq. 'spec') then
          call copy(w2, cpr%fld%x(1,1,1,e), nxyz) ! start from phys field
       else
          call copy(w2, cpr%fldhat(1,1,1,e), nxyz) ! start from spec coeff
       endif
       ! apply the operator to the x direction, result in w1
       call mxm(specmat, cpr%Xh%lx, w2, cpr%Xh%lx, w1, cpr%Xh%lx*cpr%Xh%lx)
       ! apply matrix to y direction, result in w2
       do k=1,cpr%Xh%lx
          call mxm(w1(1,1,k),cpr%Xh%lx,specmatt,cpr%Xh%lx,w2(1,1,k),cpr%Xh%lx)
       enddo
       ! apply matrix to z direction, result always in fldhat
       call mxm (w2, cpr%Xh%lx*cpr%Xh%lx, specmatt, cpr%Xh%lx,&
            cpr%fldhat(1,1,1,e), cpr%Xh%lx)
    enddo

  end subroutine cpr_goto_space

  !> Truncate the frequencies
  subroutine cpr_truncate_wn(cpr, coef)
    type(cpr_t), intent(inout) :: cpr
    type(coef_t), intent(inout) :: coef
    real(kind=rp) :: w2(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: w1(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: vsort(cpr%Xh%lx * cpr%Xh%lx * cpr%Xh%lx)
    real(kind=rp) :: vtrunc(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: vtemp(cpr%Xh%lx * cpr%Xh%lx * cpr%Xh%lx)
    real(kind=rp) :: errvec(cpr%Xh%lx, cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: fx(cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: fy(cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: fz(cpr%Xh%lx, cpr%Xh%lx)
    real(kind=rp) :: l2norm, oldl2norm, targeterr
    integer :: isort(cpr%Xh%lx * cpr%Xh%lx * cpr%Xh%lx)
    integer :: i, j, k, e, nxyz, nelv
    integer :: kut, kutx, kuty, kutz, nx
    character(len=LOG_SIZE) :: log_buf

    ! define some constants
    nx = cpr%Xh%lx
    nxyz = cpr%Xh%lx*cpr%Xh%lx*cpr%Xh%lx
    nelv = cpr%msh%nelv
    targeterr = 1e-3_rp

    !truncate for every element
    do e = 1, nelv
       ! create temp vector where the trunc coeff will be
       call copy(vtemp, cpr%fldhat(1,1,1,e), nxyz)
       call copy(vtrunc, cpr%fldhat(1,1,1,e), nxyz)
       ! sort the coefficients by absolute value
       call sortcoeff(vsort, cpr%fldhat(1,1,1,e), isort, nxyz)
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
          call sub3(errvec, cpr%fldhat(1,1,1,e), vtemp, nxyz)

          ! get the norm of the error
          l2norm = get_elem_l2norm(errvec, coef, 'spec' ,e)

       end do

       ! copy the truncated field back to the main object
       call copy(cpr%fldhat(1,1,1,e), vtrunc, nxyz)


       !================debugging info

       !just to check that all is good, resort vsort
       if (e .eq. 50) then
          write(log_buf, '(A)') &
               'Debugging info for e=50. Do not forget to delete'
          call neko_log%message(log_buf)

          call reord(vsort,isort,nxyz)

          write(log_buf, '(A)') &
               'Norm'
          call neko_log%message(log_buf)
          !chech that the copy is fine in one entry
          write(log_buf, '(A,E15.7)') &
               'l2norm:', oldl2norm
          call neko_log%message(log_buf)
          write(log_buf, '(A)') &
               'cut off level'
          call neko_log%message(log_buf)
          !chech that the copy is fine in one entry
          write(log_buf, '(A,I5)') &
               'kut:', kut
          call neko_log%message(log_buf)
          write(log_buf, '(A)') &
               'spectral coefficients'
          call neko_log%message(log_buf)
          do i = 1, 10
             write(log_buf, '(A,E15.7,A,E15.7)') &
                  'org coeff:', vsort(i), ' truncated coeff', &
                  cpr%fldhat(i,1,1,e)
             call neko_log%message(log_buf)
          end do
       end if
    end do

  end subroutine cpr_truncate_wn

  !> Sort the spectral coefficient in descending order
  !! array vsort. The original indices are stored in the isort vector.
  subroutine sortcoeff(vsort, v, isort, nxyz)
    integer, intent(in) :: nxyz
    real(kind=rp), intent(inout) :: vsort(nxyz)
    real(kind=rp), intent(inout) :: v(nxyz)
    integer, intent(inout) :: isort(nxyz)
    integer :: i

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
    integer i, e, nxyz
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
