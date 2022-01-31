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
!> Tensor operations libxsmm backend
module tensor_xsmm
  use num_types
  use mxm_wrapper
  implicit none
  private

  public :: tnsr2d_el_xsmm, tnsr3d_el_xsmm, tnsr3d_xsmm, tnsr1_3d_xsmm
  
contains
  
  subroutine tnsr2d_el_xsmm(v, nv, u, nu, A, Bt)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv), u(nu*nu)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu,nv)
    real(kind=rp) :: work(0:nu**2*nv)

    call mxm(A, nv, u, nu, work, nu)
    call mxm(work, nv, Bt, nu, v, nv)
    
  end subroutine tnsr2d_el_xsmm

  subroutine tnsr3d_el_xsmm(v, nv, u, nu, A, Bt, Ct)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv*nv), u(nu*nu*nu)
    real(kind=rp), intent(inout) :: A(nv,nu),Bt(nu, nv),Ct(nu,nv)
    real(kind=rp) :: work(0:nu**2*nv), work2(0:nu*nv**2)
    integer :: i, nunu, nvnu, nvnv

    nvnu = nv * nu
    nunu = nu * nu 
    nvnv = nv * nv
    
    call mxm(A, nv, u(1), nu ,work, nunu)
    do i = 0,nu-1
       call mxm(work(nvnu*i), nv, Bt, nu, work2(nv*nv*i), nv)
    end do
    call mxm(work2, nvnv, Ct, nu, v(1), nv)
    
  end subroutine tnsr3d_el_xsmm
 
  subroutine tnsr3d_xsmm(v, nv, u, nu, A, Bt, Ct, nelv)
    integer, intent(inout) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv,nelv), u(nu*nu*nu,nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    real(kind=rp) :: work(0:nu**2*nv), work2(0:nu*nv**2)
    integer :: ie, i, nunu, nvnu, nvnv

    nvnu = nv * nu
    nunu = nu * nu 
    nvnv = nv * nv
    
    do ie = 1,nelv
       call mxm(A, nv, u(1,ie), nu, work, nunu)
       do i = 0,nu-1
          call mxm(work(nvnu*i), nv, Bt, nu, work2(nv*nv*i), nv)
       end do
       call mxm(work2, nvnv, Ct, nu, v(1,ie), nv)
    end do
    
  end subroutine tnsr3d_xsmm

  subroutine tnsr1_3d_xsmm(v, nv, nu, A, Bt, Ct, nelv) 
    integer, intent(in) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv*nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    real(kind=rp) :: work(0:nu**2*nv), work2(0:nu*nv**2)
    integer :: e, e0, ee, es, iu, iv, i, nu3, nv3

    e0 = 1
    es = 1
    ee = nelv

    if (nv.gt.nu) then
       e0 = nelv
       es = -1
       ee = 1
    endif

    nu3 = nu**3
    nv3 = nv**3

    do e = e0,ee,es
       iu = 1 + (e-1)*nu3
       iv = 1 + (e-1)*nv3
       call mxm(A, nv, v(iu), nu, work, nu*nu)
       do i = 0,nu-1
          call mxm(work(nv*nu*i), nv, Bt, nu, work2(nv*nv*i), nv)
       end do
       call mxm(work2, nv*nv, Ct, nu, v(iv), nv)
    end do
  end subroutine tnsr1_3d_xsmm

end module tensor_xsmm
