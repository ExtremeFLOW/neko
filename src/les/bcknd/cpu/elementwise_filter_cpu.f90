! Copyright (c) 2024, The Neko Authors
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
!> Implements the CPU kernel for the `elementwise_filter_t` type.
module elementwise_filter_cpu
  use num_types, only : rp
  use math, only : copy
  use speclib, only : zwgll, legendre_poly
  use matrix, only : matrix_t
  use tensor, only : trsp
  use mxm_wrapper, only : mxm
  implicit none
  private

  public :: build_1d_cpu

contains

  !> Build the 1d filter for an element on the CPU.
  !> Suppose field x is filtered into x_hat by x_hat = fh*x.
  !! @param fh The 1D filter operator.
  !! @param fht The transpose of fh.
  !! @param trnfr The transfer function containing weights for different modes.
  !! @param nx number of points, dimension of x.
  !! @param filter_type
  subroutine build_1d_cpu(fh, fht, trnsfr, nx, filter_type)
    integer, intent(in) :: nx
    real(kind=rp), intent(inout) :: fh(nx, nx), fht(nx, nx)
    real(kind=rp), intent(in) :: trnsfr(nx)
    real(kind=rp) :: diag(nx, nx), rmult(nx), Lj(nx), zpts(nx)
    type(matrix_t) :: phi, pht
    integer :: n, i, j, k
    real(kind=rp) :: z
    character(len=*), intent(in) :: filter_type

    call phi%init(nx, nx)
    call pht%init(nx, nx)

    call zwgll(zpts, rmult, nx)

    n  = nx-1
    do j = 1, nx
       z = zpts(j)
       call legendre_poly(Lj, z, n)
       select case (filter_type)
       case("Boyd")
          pht%x(1,j) = Lj(1)
          pht%x(2,j) = Lj(2)
          do k=3,nx
             pht%x(k,j) = Lj(k)-Lj(k-2)
          end do
       case("nonBoyd")
          pht%x(:,j) = Lj
       end select
    end do

    call trsp(phi%x, nx, pht%x, nx)
    pht = phi
    call pht%inverse()

    diag = 0.0_rp

    do i=1,nx
       diag(i,i) = trnsfr(i)
    end do

    call mxm  (diag, nx, pht%x, nx, fh, nx)       !          -1
    call mxm  (phi%x, nx, fh, nx, pht%x, nx)      !     V D V

    call copy      (fh, pht%x, nx*nx)
    call trsp (fht, nx, fh, nx)

  end subroutine build_1d_cpu

end module elementwise_filter_cpu

