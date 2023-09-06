! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC. 
!
! The UChicago Argonne, LLC as Operator of Argonne National
! Laboratory holds copyright in the Software. The copyright holder
! reserves all rights except those expressly granted to licensees,
! and U.S. Government license rights.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the disclaimer below.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the disclaimer (as noted below)
! in the documentation and/or other materials provided with the
! distribution.
!
! 3. Neither the name of ANL nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF 
! ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Additional BSD Notice
! ---------------------
! 1. This notice is required to be provided under our contract with
! the U.S. Department of Energy (DOE). This work was produced at
! Argonne National Laboratory under Contract 
! No. DE-AC02-06CH11357 with the DOE.
!
! 2. Neither the United States Government nor UCHICAGO ARGONNE, 
! LLC nor any of their employees, makes any warranty, 
! express or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.
!
! 3. Also, reference herein to any specific commercial products, process, 
! or services by trade name, trademark, manufacturer or otherwise does 
! not necessarily constitute or imply its endorsement, recommendation, 
! or favoring by the United States Government or UCHICAGO ARGONNE LLC. 
! The views and opinions of authors expressed 
! herein do not necessarily state or reflect those of the United States 
! Government or UCHICAGO ARGONNE, LLC, and shall 
! not be used for advertising or product endorsement purposes.
!
!> Tensor operations.
module tensor
  use tensor_xsmm
  use tensor_cpu
  use tensor_sx
  use tensor_device
  use num_types
  use mxm_wrapper
  use neko_config
  use device
  implicit none
  private

  interface transpose
     module procedure trsp, trsp1
  end interface transpose

  interface triple_tensor_product
     module procedure triple_tensor_product_scalar, triple_tensor_product_vector
  end interface triple_tensor_product

public tensr3, transpose, trsp, trsp1, &
     tnsr2d_el, tnsr3d_el, tnsr3d, tnsr1_3d, addtnsr, &
     triple_tensor_product, tnsr3d_el_list


contains

  !> Tensor product \f$ v =(C \otimes B \otimes A) u \f$.
  subroutine tensr3(v, nv, u, nu, A, Bt, Ct, w)
    integer :: nv
    integer :: nu
    real(kind=rp), intent(inout) :: v(nv, nv, nv)
    real(kind=rp), intent(inout) :: u(nu, nu, nu)
    real(kind=rp), intent(inout) :: w(nu*nu*nv)
    real(kind=rp), intent(inout) :: A(nv, nu)
    real(kind=rp), intent(inout) :: Bt(nu, nv)
    real(kind=rp), intent(inout) :: Ct(nu, nv)
    integer :: j, k, l, nunu, nvnv, nunv
    
    nunu = nu**2
    nvnv = nv**2
    nunv = nu*nv

    !>@todo Add 2d case

    call mxm(A, nv, u, nu, v, nunu)
    k = 1
    l = 1
    do j = 1, nu
       call mxm(v(k, 1, 1), nv, Bt, nu, w(l), nv)
       k = k + nunv
       l = l + nvnv
    end do
    call mxm(w, nvnv, Ct, nu, v, nv)
    
  end subroutine tensr3

  !> Transpose of a rectangular tensor \f$ A = B^T \f$.
  subroutine trsp(a, lda, b, ldb)
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    real(kind=rp), intent(inout) :: a(lda, ldb)
    real(kind=rp), intent(in) :: b(ldb, lda)
    integer :: i, j

    do j = 1, ldb
       do i = 1, lda
          a(i, j) = b(j, i)
       end do
    end do
    
  end subroutine trsp

  !> In-place transpose of a square tensor.
  subroutine trsp1(a, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: a(n, n)
    real(kind=rp) :: tmp
    integer :: i, j

    do j = 1, n
       do i = j + 1, n
          tmp = a(i, j)
          a(i, j) = a(j, i)
          a(j, i) = tmp
       end do
    end do
    
  end subroutine trsp1

  !> Computes \f$ v = A u B^T \f$.
  subroutine tnsr2d_el(v, nv, u, nu, A, Bt)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv), u(nu*nu)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu,nv)

    if (NEKO_BCKND_SX .eq. 1) then
       call tnsr2d_el_sx(v, nv, u, nu, A, Bt)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call tnsr2d_el_xsmm(v, nv, u, nu, A, Bt)
    else
       call tnsr2d_el_cpu(v, nv, u, nu, A, Bt)
    end if
    
  end subroutine tnsr2d_el

  !> Tensor product \f$ v =(C \otimes B \otimes A) u \f$
  !! performed on a single element.
  subroutine tnsr3d_el(v, nv, u, nu, A, Bt, Ct)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv*nv), u(nu*nu*nu)
    real(kind=rp), intent(inout) :: A(nv,nu),Bt(nu, nv),Ct(nu,nv)

    if (NEKO_BCKND_SX .eq. 1) then
       call tnsr3d_el_sx(v, nv, u, nu, A, Bt, Ct)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call tnsr3d_el_xsmm(v, nv, u, nu, A, Bt, Ct)
    else
       call tnsr3d_el_cpu(v, nv, u, nu, A, Bt, Ct)
    end if
    
  end subroutine tnsr3d_el

  !> Tensor product \f$ v =(C \otimes B \otimes A) u \f$
  !! performed on a subset of the  elements.
  subroutine tnsr3d_el_list(v, nv, u, nu, A, Bt, Ct, el_list, n_pt)
    integer, intent(in) :: nv, nu, n_pt, el_list(n_pt)
    real(kind=rp), intent(inout) :: v(nv*nv*nv, n_pt), u(nu*nu*nu,1)
    real(kind=rp), intent(inout) :: A(nv,nu,n_pt),Bt(nu, nv,n_pt),Ct(nu,nv,n_pt)
    integer :: i

    if (NEKO_BCKND_SX .eq. 1) then
       do i = 1, n_pt
          call tnsr3d_el_sx(v(1,i), nv, u(1,el_list(i)), nu, A(1,1,i), Bt(1,1,i), Ct(1,1,i))
       end do
    else if (NEKO_BCKND_XSMM .eq. 1) then
       do i = 1, n_pt
          call tnsr3d_el_xsmm(v(1,i), nv, u(1,el_list(i)), nu, A(1,1,i), Bt(1,1,i), Ct(1,1,i))
       end do
    else if (NEKO_BCKND_DEVICE .eq. 1) then
    !   call tnsr3d_el_list_device(1,i), nv, u(1,i), nu, A(1,1,i), Bt(1,1,i), Ct(1,1,i))
    else
       do i = 1, n_pt
          call tnsr3d_el_cpu(v(1,i), nv, u(1,el_list(i)), nu, A(1,1,i), Bt(1,1,i), Ct(1,1,i))
       end do
    end if
    
  end subroutine tnsr3d_el_list


  !> Tensor product \f$ v =(C \otimes B \otimes A) u \f$ performed on
  !!`nelv` elements.
  subroutine tnsr3d(v, nv, u, nu, A, Bt, Ct, nelv)
    integer, intent(inout) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv,nelv), u(nu*nu*nu,nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    type(c_ptr) :: v_d, u_d, A_d, Bt_d, Ct_d
   

    if (NEKO_BCKND_SX .eq. 1) then
       call tnsr3d_sx(v, nv, u, nu, A, Bt, Ct, nelv)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call tnsr3d_xsmm(v, nv, u, nu, A, Bt, Ct, nelv)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       v_d = device_get_ptr(v)
       u_d = device_get_ptr(u)
       A_d = device_get_ptr(A)
       Bt_d = device_get_ptr(Bt)
       Ct_d = device_get_ptr(Ct)
       call tnsr3d_device(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv)
    else
       call tnsr3d_cpu(v, nv, u, nu, A, Bt, Ct, nelv)
    end if
    
  end subroutine tnsr3d

  !> In place tensor product \f$ v =(C \otimes B \otimes A) v \f$.
  subroutine tnsr1_3d(v, nv, nu, A, Bt, Ct, nelv)
    integer, intent(inout) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv*nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)

    if (NEKO_BCKND_SX .eq. 1) then
       call tnsr1_3d_sx(v, nv, nu, A, Bt, Ct, nelv)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call tnsr1_3d_xsmm(v, nv, nu, A, Bt, Ct, nelv)
    else
       call tnsr1_3d_cpu(v, nv, nu, A, Bt, Ct, nelv)
    end if

  end subroutine tnsr1_3d

  !> Maps and adds to S a tensor product form of the three functions H1,H2,H3.
  !! This is a single element routine used for deforming geometry.
  subroutine addtnsr(s, h1, h2, h3, nx, ny, nz)

    integer, intent(in) :: nx, ny, nz
    real(kind=rp), intent(in) :: h1(nx), h2(ny), h3(nz) 
    real(kind=rp), intent(inout) ::  s(nx, ny, nz)
    real(kind=rp) :: hh
    integer :: ix, iy, iz
  
    do iz = 1,nz
       do iy = 1,ny
          hh = h2(iy)*h3(iz)
          do ix = 1,nx
             s(ix,iy,iz) = s(ix,iy,iz)+hh*h1(ix)
          end do
       end do
    end do
    
  end subroutine addtnsr

  !> Computes the tensor product \f$ v =(H_t \otimes H_s \otimes H_r) u \f$.
  !! This operation is usually performed for spectral interpolation of a
  !! scalar field as defined by
  !! \f{eqnarray*}{
  !!    v(r,s,t) = \sum_{i=0}^{N}{\sum_{j=0}^{N}{
  !!      \sum_{k=0}^{N}{u_{ijk}h_i(r)h_j(s)h_k(t)}}}
  !! \f}
  !!
  !! @param v Interpolated value (scalar).
  !! @param u Field values at the GLL points (e.g. velocity in x-direction).
  !! @param nu Size of the interpolation weights (usually `lx`).
  !! @param Hr Interpolation weights in the r-direction.
  !! @param Hs Interpolation weights in the s-direction.
  !! @param Ht Interpolation weights in the t-direction.
  subroutine triple_tensor_product_scalar(v, u, nu, Hr, Hs, Ht)
    real(kind=rp), intent(inout) :: v
    integer, intent(in) :: nu
    real(kind=rp), intent(inout) :: u(nu,nu,nu)
    real(kind=rp), intent(inout) :: Hr(nu)
    real(kind=rp), intent(inout) :: Hs(nu)
    real(kind=rp), intent(inout) :: Ht(nu)

    ! Artificially reshape v into a 1-dimensional array
    ! since this is what tnsr3d_el needs as input argument
    real(kind=rp) :: vv(1)
    ! vv(1) = v

    call tnsr3d_el(vv,1,u,nu,Hr,Hs,Ht)

    v = vv(1)

  end subroutine triple_tensor_product_scalar

  !> Computes the tensor product on a vector field
  !! \f$ \mathbf{v} =(H_t \otimes H_s \otimes H_r) \mathbf{u} \f$.
  !! This operation is usually performed for spectral interpolation on
  !! a 3D vector field \f$ \mathbf{u} = (u_1,u_2,u_3) \f$ as defined by
  !! \f{eqnarray*}{
  !!    \mathbf{v}(r,s,t) = \sum_{i=0}^{N}{\sum_{j=0}^{N}{
  !!      \sum_{k=0}^{N}{\mathbf{u}_{ijk}h_i(r)h_j(s)h_k(t)}}}
  !! \f}
  !!
  !! @param v Interpolated value (scalar).
  !! @param u1 3D-array containing values at the GLL points (e.g. velocity).
  !! @param u2 3D-array containing values at the GLL points (e.g. velocity).
  !! @param u3 3D-array containing values at the GLL points (e.g. velocity).
  !! @param nu Size of the interpolation weights (usually `lx`).
  !! @param Hr Interpolation weights in the r-direction.
  !! @param Hs Interpolation weights in the s-direction.
  !! @param Ht Interpolation weights in the t-direction.
  subroutine triple_tensor_product_vector(v, u1, u2, u3, nu, Hr, Hs, Ht)
    real(kind=rp), intent(inout) :: v(3)
    integer, intent(in) :: nu
    real(kind=rp), intent(inout) :: u1(nu,nu,nu)
    real(kind=rp), intent(inout) :: u2(nu,nu,nu)
    real(kind=rp), intent(inout) :: u3(nu,nu,nu)
    real(kind=rp), intent(inout) :: Hr(nu)
    real(kind=rp), intent(inout) :: Hs(nu)
    real(kind=rp), intent(inout) :: Ht(nu)

    call triple_tensor_product_scalar(v(1), u1, nu, Hr, Hs, Ht)
    call triple_tensor_product_scalar(v(2), u2, nu, Hr, Hs, Ht)
    call triple_tensor_product_scalar(v(3), u3, nu, Hr, Hs, Ht)

  end subroutine triple_tensor_product_vector

end module tensor
