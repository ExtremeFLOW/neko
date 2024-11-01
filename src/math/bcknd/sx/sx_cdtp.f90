! Copyright (c) 2021, The Neko Authors
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
!> DT*X kernels for SX-Aurora
submodule (opr_sx) sx_cdtp
  use math, only : col3, invcol2
  implicit none

contains

  module subroutine opr_sx_cdtp(dtx, x, dr, ds, dt, coef)
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: dtx(coef%Xh%lxyz, coef%msh%nelv)
    real(kind=rp), intent(inout) :: x(coef%Xh%lxyz, coef%msh%nelv)
    real(kind=rp), intent(in) :: dr(coef%Xh%lxyz, coef%msh%nelv)
    real(kind=rp), intent(in) :: ds(coef%Xh%lxyz, coef%msh%nelv)
    real(kind=rp), intent(in) :: dt(coef%Xh%lxyz, coef%msh%nelv)

    associate(Xh => coef%Xh, msh => coef%msh, dof => coef%dof)
      select case (Xh%lx)
      case (14)
         call sx_cdtp_lx14(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case (13)
         call sx_cdtp_lx13(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case (12)
         call sx_cdtp_lx12(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case (11)
         call sx_cdtp_lx11(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case (10)
         call sx_cdtp_lx10(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case (9)
         call sx_cdtp_lx9(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case (8)
         call sx_cdtp_lx8(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case (7)
         call sx_cdtp_lx7(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case (6)
         call sx_cdtp_lx6(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case (5)
         call sx_cdtp_lx5(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case (4)
         call sx_cdtp_lx4(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case (3)
         call sx_cdtp_lx3(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case (2)
         call sx_cdtp_lx2(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case default
         call sx_cdtp_lx(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size(), Xh%lx)
      end select
    end associate

  end subroutine opr_sx_cdtp
  
  subroutine sx_cdtp_lx(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd, lx)
    integer, intent(in) :: nel, nd, lx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk) * ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx

  subroutine sx_cdtp_lx14(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 14
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk) * ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx14

  subroutine sx_cdtp_lx13(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 13
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk) * ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx13

  subroutine sx_cdtp_lx12(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 12
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx12

  subroutine sx_cdtp_lx11(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 11
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx11

  subroutine sx_cdtp_lx10(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 10
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk) * ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx10

  subroutine sx_cdtp_lx9(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 9
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk) * ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx9

  subroutine sx_cdtp_lx8(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 8
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk) * ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx8

  subroutine sx_cdtp_lx7(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 7
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk) * ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx7

  subroutine sx_cdtp_lx6(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 6
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk) * ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx6

  subroutine sx_cdtp_lx5(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 5
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk) * ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx5

  subroutine sx_cdtp_lx4(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 4
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk) * ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx4

  subroutine sx_cdtp_lx3(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 3
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk) * ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx3

  subroutine sx_cdtp_lx2(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 2
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt, &
         jac, B
    real(kind=rp), intent(in)  :: dxt(lx, lx), dyt(lx, lx), dzt(lx, lx)
    real(kind=rp), dimension(lx, lx, lx, nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i = 1, lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk = 1, lx
             tmp = tmp + dxt(i, kk) * ta1(kk, jj,1,1)
          end do
          dtx(i, jj,1,1) = tmp
       end do
    end do

    call col3(ta1, wx, ds, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dyt(j, kk) * ta1(i, kk,k,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

    call col3(ta1, wx, dt, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   tmp = tmp + dzt(k, kk) * ta1(i,j, kk,e)
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
    end do

  end subroutine sx_cdtp_lx2


end submodule sx_cdtp
