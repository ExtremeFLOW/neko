!> DT*X kernels for SX-Aurora
module sx_cdtp
  use num_types
  use math
  implicit none

contains

  subroutine sx_cdtp_lx12(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 12
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i=1,lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk=1,lx
             tmp = tmp + dxt(i,kk) * ta1(kk,jj,1,1)
          end do
          dtx(i,jj,1,1) = tmp
       enddo
    enddo

    call col3(ta1, wx, ds, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dyt(j, kk) * ta1(i,kk,k,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

    call col3(ta1, wx, dt, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j,kk,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

  end subroutine sx_cdtp_lx12

  subroutine sx_cdtp_lx11(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 11
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i=1,lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk=1,lx
             tmp = tmp + dxt(i,kk) * ta1(kk,jj,1,1)
          end do
          dtx(i,jj,1,1) = tmp
       enddo
    enddo

    call col3(ta1, wx, ds, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dyt(j, kk) * ta1(i,kk,k,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

    call col3(ta1, wx, dt, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j,kk,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

  end subroutine sx_cdtp_lx11

  subroutine sx_cdtp_lx10(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 10
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i=1,lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk=1,lx
             tmp = tmp + dxt(i,kk) * ta1(kk,jj,1,1)
          end do
          dtx(i,jj,1,1) = tmp
       enddo
    enddo

    call col3(ta1, wx, ds, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dyt(j, kk) * ta1(i,kk,k,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

    call col3(ta1, wx, dt, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j,kk,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

  end subroutine sx_cdtp_lx10

  subroutine sx_cdtp_lx9(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 9
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i=1,lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk=1,lx
             tmp = tmp + dxt(i,kk) * ta1(kk,jj,1,1)
          end do
          dtx(i,jj,1,1) = tmp
       enddo
    enddo

    call col3(ta1, wx, ds, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dyt(j, kk) * ta1(i,kk,k,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

    call col3(ta1, wx, dt, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j,kk,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

  end subroutine sx_cdtp_lx9

  subroutine sx_cdtp_lx8(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 8
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i=1,lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk=1,lx
             tmp = tmp + dxt(i,kk) * ta1(kk,jj,1,1)
          end do
          dtx(i,jj,1,1) = tmp
       enddo
    enddo

    call col3(ta1, wx, ds, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dyt(j, kk) * ta1(i,kk,k,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

    call col3(ta1, wx, dt, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j,kk,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

  end subroutine sx_cdtp_lx8

  subroutine sx_cdtp_lx7(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 7
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i=1,lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk=1,lx
             tmp = tmp + dxt(i,kk) * ta1(kk,jj,1,1)
          end do
          dtx(i,jj,1,1) = tmp
       enddo
    enddo

    call col3(ta1, wx, ds, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dyt(j, kk) * ta1(i,kk,k,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

    call col3(ta1, wx, dt, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j,kk,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

  end subroutine sx_cdtp_lx7

  subroutine sx_cdtp_lx6(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 6
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i=1,lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk=1,lx
             tmp = tmp + dxt(i,kk) * ta1(kk,jj,1,1)
          end do
          dtx(i,jj,1,1) = tmp
       enddo
    enddo

    call col3(ta1, wx, ds, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dyt(j, kk) * ta1(i,kk,k,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

    call col3(ta1, wx, dt, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j,kk,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

  end subroutine sx_cdtp_lx6

  subroutine sx_cdtp_lx5(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 5
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i=1,lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk=1,lx
             tmp = tmp + dxt(i,kk) * ta1(kk,jj,1,1)
          end do
          dtx(i,jj,1,1) = tmp
       enddo
    enddo

    call col3(ta1, wx, ds, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dyt(j, kk) * ta1(i,kk,k,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

    call col3(ta1, wx, dt, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j,kk,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

  end subroutine sx_cdtp_lx5

  subroutine sx_cdtp_lx4(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 4
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i=1,lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk=1,lx
             tmp = tmp + dxt(i,kk) * ta1(kk,jj,1,1)
          end do
          dtx(i,jj,1,1) = tmp
       enddo
    enddo

    call col3(ta1, wx, ds, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dyt(j, kk) * ta1(i,kk,k,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

    call col3(ta1, wx, dt, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j,kk,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

  end subroutine sx_cdtp_lx4

  subroutine sx_cdtp_lx3(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 3
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i=1,lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk=1,lx
             tmp = tmp + dxt(i,kk) * ta1(kk,jj,1,1)
          end do
          dtx(i,jj,1,1) = tmp
       enddo
    enddo

    call col3(ta1, wx, ds, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dyt(j, kk) * ta1(i,kk,k,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

    call col3(ta1, wx, dt, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j,kk,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

  end subroutine sx_cdtp_lx3

  subroutine sx_cdtp_lx2(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 2
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, kk, jj

    call col3(wx, B, x, nd)
    call invcol2(wx, jac, nd)
    call col3(ta1, wx, dr, nd)

    do i=1,lx
       do jj = 1, lx * lx * nel
          tmp = 0d0
          do kk=1,lx
             tmp = tmp + dxt(i,kk) * ta1(kk,jj,1,1)
          end do
          dtx(i,jj,1,1) = tmp
       enddo
    enddo

    call col3(ta1, wx, ds, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dyt(j, kk) * ta1(i,kk,k,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

    call col3(ta1, wx, dt, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,nel
                tmp = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   tmp = tmp + dzt(k, kk)*ta1(i,j,kk,e)
                enddo
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             enddo
          enddo
       enddo
    end do

  end subroutine sx_cdtp_lx2


end module sx_cdtp
