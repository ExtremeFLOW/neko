!> DT*X kernels for SX-Aurora
module sx_cdtp
  use num_types
  use math
  implicit none

contains

  subroutine sx_cdtp_lx10(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 10
    integer, intent(in) :: nel, nd
    real(kind=dp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=dp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=dp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=dp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=dp) :: tmp
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

  subroutine sx_cdtp_lx8(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, nd)
    integer, parameter :: lx = 8
    integer, intent(in) :: nel, nd
    real(kind=dp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=dp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=dp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=dp), dimension(lx,lx,lx,nel) :: wx, ta1
    real(kind=dp) :: tmp
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
  
end module sx_cdtp
