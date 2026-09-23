! Copyright (c) 2026, The Neko Authors
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
module cpu_matrix_math
  use num_types, only : rp, xp
  use utils, only : neko_error
  implicit none
  private

  public :: cpu_matrix_inverse
contains

  !> Gauss-Jordan matrix inversion with full pivoting
  !! Num. Rec. p. 30, 2nd Ed., Fortran
  !! @param x     is an square matrix
  !! @param r     is the number of rows
  !! @param c     is the number of columns
  !! rmult is this  work array of length nrows = ncols
  subroutine cpu_matrix_inverse(x, r, c)
    integer, intent(in) :: r, c
    real(kind=rp), intent(inout) :: x(r, c)
    integer :: indr(r), indc(c), ipiv(c)
    real(kind=xp) :: rmult(r), amx, tmp, piv, eps
    integer :: i, j, k, ir, jc

    if (.not. (c .eq. r)) then
       call neko_error("Fatal error: trying to invert this matrix that is " // &
            "not square")
    end if

    eps = 1e-9_rp
    ipiv = 0

    do k = 1, r
       amx = 0.0_rp
       do i = 1, r ! Pivot search
          if (ipiv(i) .ne. 1) then
             do j = 1, r
                if (ipiv(j) .eq. 0) then
                   if (abs(x(i, j)) .ge. amx) then
                      amx = abs(x(i, j))
                      ir = i
                      jc = j
                   end if
                else if (ipiv(j) .gt. 1) then
                   return
                end if
             end do
          end if
       end do
       ipiv(jc) = ipiv(jc) + 1

       !  Swap rows
       if (ir .ne. jc) then
          do j = 1, c
             tmp = x(ir, j)
             x(ir, j) = x(jc, j)
             x(jc, j) = tmp
          end do
       end if
       indr(k) = ir
       indc(k) = jc

       if (abs(x(jc, jc)) .lt. eps) then
          call neko_error("matrix_inverse error: small Gauss Jordan Piv")
       end if
       piv = 1.0_xp/x(jc, jc)
       x(jc, jc) = 1.0_xp
       do j = 1, c
          x(jc, j) = x(jc, j)*piv
       end do

       do j = 1, c
          tmp = x(jc, j)
          x(jc, j) = x(1 , j)
          x(1 , j) = tmp
       end do
       do i = 2, r
          rmult(i) = x(i, jc)
          x(i, jc) = 0.0_rp
       end do

       do j = 1, c
          do i = 2, r
             x(i, j) = x(i, j) - rmult(i)*x(1, j)
          end do
       end do

       do j = 1, c
          tmp = x(jc, j)
          x(jc, j) = x(1 , j)
          x(1 , j) = tmp
       end do
    end do

    ! Unscramble matrix
    do j = r, 1, -1
       if (indr(j) .ne. indc(j)) then
          do i = 1, r
             tmp = x(i, indr(j))
             x(i, indr(j)) = x(i, indc(j))
             x(i, indc(j)) = tmp
          end do
       end if
    end do

  end subroutine cpu_matrix_inverse


end module cpu_matrix_math
