!> Fast Diagonalization SX-Aurora backend
module fdm_sx
  use num_types
  use tensor_sx
  implicit none

contains

  subroutine fdm_do_fast_sx(e, r, s, d, nl, ldim, nelv)
    integer, intent(in) :: nl, nelv, ldim
    real(kind=rp), intent(inout) :: e(nl**ldim, nelv)
    real(kind=rp), intent(inout) :: r(nl**ldim, nelv)
    real(kind=rp), intent(inout) :: s(nl*nl,2,ldim, nelv)
    real(kind=rp), intent(inout) :: d(nl**ldim, nelv)    
    integer ::  ie, nn, i

    nn = nl**ldim
    if(.not. ldim .eq. 3) then
       do ie = 1, nelv
          call tnsr2d_el_sx(e(1,ie), nl, r(1,ie), nl, s(1,2,1,ie), s(1,1,2,ie))
          do i = 1, nn
             r(i,ie) = d(i,ie) * e(i,ie)
          end do
          call tnsr2d_el_sx(e(1,ie), nl, r(1,ie), nl, s(1,1,1,ie), s(1,2,2,ie))
       end do
    else       
       select case (nl)
       case (14)
          call fdm_do_fast_sx_nl14(e, r, s, d, nelv)
       case (13)
          call fdm_do_fast_sx_nl13(e, r, s, d, nelv)
       case (12)
          call fdm_do_fast_sx_nl12(e, r, s, d, nelv)
       case (11)
          call fdm_do_fast_sx_nl11(e, r, s, d, nelv)
       case (10)
          call fdm_do_fast_sx_nl10(e, r, s, d, nelv)
       case (9)
          call fdm_do_fast_sx_nl9(e, r, s, d, nelv)
       case (8)
          call fdm_do_fast_sx_nl8(e, r, s, d, nelv)
       case (7)
          call fdm_do_fast_sx_nl7(e, r, s, d, nelv)
       case (6)
          call fdm_do_fast_sx_nl6(e, r, s, d, nelv)
       case (5)
          call fdm_do_fast_sx_nl5(e, r, s, d, nelv)
       case (4)
          call fdm_do_fast_sx_nl4(e, r, s, d, nelv)
       case (3)
          call fdm_do_fast_sx_nl3(e, r, s, d, nelv)
       case (2)
          call fdm_do_fast_sx_nl2(e, r, s, d, nelv)
       end select
    end if
  end subroutine fdm_do_fast_sx

  subroutine fdm_do_fast_sx_nl14(e, r, s, d, nelv)
    integer, parameter :: n = 14
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,2,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,2,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,2,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,2,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,2,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,2,1,ie) * r(8 + n * (j - 1), ie) &
                         + s(i,9,2,1,ie) * r(9 + n * (j - 1), ie) &
                         + s(i,10,2,1,ie) * r(10 + n * (j - 1), ie) &
                         + s(i,11,2,1,ie) * r(11 + n * (j - 1), ie) &
                         + s(i,12,2,1,ie) * r(12 + n * (j - 1), ie) &
                         + s(i,13,2,1,ie) * r(13 + n * (j - 1), ie) &
                         + s(i,14,2,1,ie) * r(14 + n * (j - 1), ie)    
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,1,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,1,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,1,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,1,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,1,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,1,2,ie) &
                            + wrk(l + n * (9 - 1) + nn * (i - 1), ie) &
                              * s(9,j,1,2,ie) &
                            + wrk(l + n * (10 - 1) + nn * (i - 1), ie) &
                              * s(10,j,1,2,ie) &
                            + wrk(l + n * (11 - 1) + nn * (i - 1), ie) &
                              * s(11,j,1,2,ie) &
                            + wrk(l + n * (12 - 1) + nn * (i - 1), ie) &
                              * s(12,j,1,2,ie) &
                            + wrk(l + n * (13 - 1) + nn * (i - 1), ie) &
                              * s(13,j,1,2,ie) &
                            + wrk(l + n * (14 - 1) + nn * (i - 1), ie) &
                              * s(14,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 1, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 1, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 1, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 1, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 1, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 1, 3, ie) &
                      + wrk2(i + nn * (9 - 1), ie) * s(9, j, 1, 3, ie) &
                      + wrk2(i + nn * (10 - 1), ie) * s(10, j, 1, 3, ie) &
                      + wrk2(i + nn * (11 - 1), ie) * s(11, j, 1, 3, ie) &
                      + wrk2(i + nn * (12 - 1), ie) * s(12, j, 1, 3, ie) &
                      + wrk2(i + nn * (13 - 1), ie) * s(13, j, 1, 3, ie) &
                      + wrk2(i + nn * (14 - 1), ie) * s(14, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,1,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,1,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,1,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,1,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,1,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,1,1,ie) * r(8 + n * (j - 1), ie) &
                         + s(i,9,1,1,ie) * r(9 + n * (j - 1), ie) &
                         + s(i,10,1,1,ie) * r(10 + n * (j - 1), ie) &
                         + s(i,11,1,1,ie) * r(11 + n * (j - 1), ie) &
                         + s(i,12,1,1,ie) * r(12 + n * (j - 1), ie) &
                         + s(i,13,1,1,ie) * r(13 + n * (j - 1), ie) &
                         + s(i,14,1,1,ie) * r(14 + n * (j - 1), ie)    
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,2,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,2,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,2,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,2,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,2,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,2,2,ie) &
                            + wrk(l + n * (9 - 1) + nn * (i - 1), ie) &
                              * s(9,j,2,2,ie) &
                            + wrk(l + n * (10 - 1) + nn * (i - 1), ie) &
                              * s(10,j,2,2,ie) &
                            + wrk(l + n * (11 - 1) + nn * (i - 1), ie) &
                              * s(11,j,2,2,ie) &
                            + wrk(l + n * (12 - 1) + nn * (i - 1), ie) &
                              * s(12,j,2,2,ie) &
                            + wrk(l + n * (13 - 1) + nn * (i - 1), ie) &
                              * s(13,j,2,2,ie) &
                            + wrk(l + n * (14 - 1) + nn * (i - 1), ie) &
                              * s(14,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 2, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 2, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 2, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 2, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 2, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 2, 3, ie) &
                      + wrk2(i + nn * (9 - 1), ie) * s(9, j, 2, 3, ie) &
                      + wrk2(i + nn * (10 - 1), ie) * s(10, j, 2, 3, ie) &
                      + wrk2(i + nn * (11 - 1), ie) * s(11, j, 2, 3, ie) &
                      + wrk2(i + nn * (12 - 1), ie) * s(12, j, 2, 3, ie) &
                      + wrk2(i + nn * (13 - 1), ie) * s(13, j, 2, 3, ie) &
                      + wrk2(i + nn * (14 - 1), ie) * s(14, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl14

  subroutine fdm_do_fast_sx_nl13(e, r, s, d, nelv)
    integer, parameter :: n = 13
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,2,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,2,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,2,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,2,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,2,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,2,1,ie) * r(8 + n * (j - 1), ie) &
                         + s(i,9,2,1,ie) * r(9 + n * (j - 1), ie) &
                         + s(i,10,2,1,ie) * r(10 + n * (j - 1), ie) &
                         + s(i,11,2,1,ie) * r(11 + n * (j - 1), ie) &
                         + s(i,12,2,1,ie) * r(12 + n * (j - 1), ie) &
                         + s(i,13,2,1,ie) * r(13 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,1,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,1,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,1,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,1,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,1,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,1,2,ie) &
                            + wrk(l + n * (9 - 1) + nn * (i - 1), ie) &
                              * s(9,j,1,2,ie) &
                            + wrk(l + n * (10 - 1) + nn * (i - 1), ie) &
                              * s(10,j,1,2,ie) &
                            + wrk(l + n * (11 - 1) + nn * (i - 1), ie) &
                              * s(11,j,1,2,ie) &
                            + wrk(l + n * (12 - 1) + nn * (i - 1), ie) &
                              * s(12,j,1,2,ie) &
                            + wrk(l + n * (13 - 1) + nn * (i - 1), ie) &
                              * s(13,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 1, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 1, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 1, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 1, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 1, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 1, 3, ie) &
                      + wrk2(i + nn * (9 - 1), ie) * s(9, j, 1, 3, ie) &
                      + wrk2(i + nn * (10 - 1), ie) * s(10, j, 1, 3, ie) &
                      + wrk2(i + nn * (11 - 1), ie) * s(11, j, 1, 3, ie) &
                      + wrk2(i + nn * (12 - 1), ie) * s(12, j, 1, 3, ie) &
                      + wrk2(i + nn * (13 - 1), ie) * s(13, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,1,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,1,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,1,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,1,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,1,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,1,1,ie) * r(8 + n * (j - 1), ie) &
                         + s(i,9,1,1,ie) * r(9 + n * (j - 1), ie) &
                         + s(i,10,1,1,ie) * r(10 + n * (j - 1), ie) &
                         + s(i,11,1,1,ie) * r(11 + n * (j - 1), ie) &
                         + s(i,12,1,1,ie) * r(12 + n * (j - 1), ie) &
                         + s(i,13,1,1,ie) * r(13 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,2,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,2,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,2,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,2,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,2,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,2,2,ie) &
                            + wrk(l + n * (9 - 1) + nn * (i - 1), ie) &
                              * s(9,j,2,2,ie) &
                            + wrk(l + n * (10 - 1) + nn * (i - 1), ie) &
                              * s(10,j,2,2,ie) &
                            + wrk(l + n * (11 - 1) + nn * (i - 1), ie) &
                              * s(11,j,2,2,ie) &
                            + wrk(l + n * (12 - 1) + nn * (i - 1), ie) &
                              * s(12,j,2,2,ie) &
                            + wrk(l + n * (13 - 1) + nn * (i - 1), ie) &
                              * s(13,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 2, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 2, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 2, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 2, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 2, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 2, 3, ie) &
                      + wrk2(i + nn * (9 - 1), ie) * s(9, j, 2, 3, ie) &
                      + wrk2(i + nn * (10 - 1), ie) * s(10, j, 2, 3, ie) &
                      + wrk2(i + nn * (11 - 1), ie) * s(11, j, 2, 3, ie) &
                      + wrk2(i + nn * (12 - 1), ie) * s(12, j, 2, 3, ie) &
                      + wrk2(i + nn * (13 - 1), ie) * s(13, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl13

  subroutine fdm_do_fast_sx_nl12(e, r, s, d, nelv)
    integer, parameter :: n = 12
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,2,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,2,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,2,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,2,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,2,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,2,1,ie) * r(8 + n * (j - 1), ie) &
                         + s(i,9,2,1,ie) * r(9 + n * (j - 1), ie) &
                         + s(i,10,2,1,ie) * r(10 + n * (j - 1), ie) &
                         + s(i,11,2,1,ie) * r(11 + n * (j - 1), ie) &
                         + s(i,12,2,1,ie) * r(12 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,1,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,1,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,1,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,1,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,1,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,1,2,ie) &
                            + wrk(l + n * (9 - 1) + nn * (i - 1), ie) &
                              * s(9,j,1,2,ie) &
                            + wrk(l + n * (10 - 1) + nn * (i - 1), ie) &
                              * s(10,j,1,2,ie) &
                            + wrk(l + n * (11 - 1) + nn * (i - 1), ie) &
                              * s(11,j,1,2,ie) &
                            + wrk(l + n * (12 - 1) + nn * (i - 1), ie) &
                              * s(12,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 1, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 1, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 1, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 1, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 1, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 1, 3, ie) &
                      + wrk2(i + nn * (9 - 1), ie) * s(9, j, 1, 3, ie) &
                      + wrk2(i + nn * (10 - 1), ie) * s(10, j, 1, 3, ie) &
                      + wrk2(i + nn * (11 - 1), ie) * s(11, j, 1, 3, ie) &
                      + wrk2(i + nn * (12 - 1), ie) * s(12, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,1,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,1,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,1,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,1,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,1,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,1,1,ie) * r(8 + n * (j - 1), ie) &
                         + s(i,9,1,1,ie) * r(9 + n * (j - 1), ie) &
                         + s(i,10,1,1,ie) * r(10 + n * (j - 1), ie) &
                         + s(i,11,1,1,ie) * r(11 + n * (j - 1), ie) &
                         + s(i,12,1,1,ie) * r(12 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,2,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,2,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,2,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,2,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,2,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,2,2,ie) &
                            + wrk(l + n * (9 - 1) + nn * (i - 1), ie) &
                              * s(9,j,2,2,ie) &
                            + wrk(l + n * (10 - 1) + nn * (i - 1), ie) &
                              * s(10,j,2,2,ie) &
                            + wrk(l + n * (11 - 1) + nn * (i - 1), ie) &
                              * s(11,j,2,2,ie) &
                            + wrk(l + n * (12 - 1) + nn * (i - 1), ie) &
                              * s(12,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 2, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 2, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 2, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 2, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 2, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 2, 3, ie) &
                      + wrk2(i + nn * (9 - 1), ie) * s(9, j, 2, 3, ie) &
                      + wrk2(i + nn * (10 - 1), ie) * s(10, j, 2, 3, ie) &
                      + wrk2(i + nn * (11 - 1), ie) * s(11, j, 2, 3, ie) &
                      + wrk2(i + nn * (12 - 1), ie) * s(12, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl12

  subroutine fdm_do_fast_sx_nl11(e, r, s, d, nelv)
    integer, parameter :: n = 11
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,2,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,2,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,2,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,2,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,2,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,2,1,ie) * r(8 + n * (j - 1), ie) &
                         + s(i,9,2,1,ie) * r(9 + n * (j - 1), ie) &
                         + s(i,10,2,1,ie) * r(10 + n * (j - 1), ie) &
                         + s(i,11,2,1,ie) * r(11 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,1,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,1,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,1,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,1,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,1,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,1,2,ie) &
                            + wrk(l + n * (9 - 1) + nn * (i - 1), ie) &
                              * s(9,j,1,2,ie) &
                            + wrk(l + n * (10 - 1) + nn * (i - 1), ie) &
                              * s(10,j,1,2,ie) &
                            + wrk(l + n * (11 - 1) + nn * (i - 1), ie) &
                              * s(11,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 1, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 1, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 1, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 1, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 1, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 1, 3, ie) &
                      + wrk2(i + nn * (9 - 1), ie) * s(9, j, 1, 3, ie) &
                      + wrk2(i + nn * (10 - 1), ie) * s(10, j, 1, 3, ie) &
                      + wrk2(i + nn * (11 - 1), ie) * s(11, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,1,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,1,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,1,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,1,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,1,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,1,1,ie) * r(8 + n * (j - 1), ie) &
                         + s(i,9,1,1,ie) * r(9 + n * (j - 1), ie) &
                         + s(i,10,1,1,ie) * r(10 + n * (j - 1), ie) &
                         + s(i,11,1,1,ie) * r(11 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,2,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,2,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,2,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,2,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,2,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,2,2,ie) &
                            + wrk(l + n * (9 - 1) + nn * (i - 1), ie) &
                              * s(9,j,2,2,ie) &
                            + wrk(l + n * (10 - 1) + nn * (i - 1), ie) &
                              * s(10,j,2,2,ie) &
                            + wrk(l + n * (11 - 1) + nn * (i - 1), ie) &
                              * s(11,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 2, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 2, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 2, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 2, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 2, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 2, 3, ie) &
                      + wrk2(i + nn * (9 - 1), ie) * s(9, j, 2, 3, ie) &
                      + wrk2(i + nn * (10 - 1), ie) * s(10, j, 2, 3, ie) &
                      + wrk2(i + nn * (11 - 1), ie) * s(11, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl11

  subroutine fdm_do_fast_sx_nl10(e, r, s, d, nelv)
    integer, parameter :: n = 10
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,2,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,2,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,2,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,2,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,2,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,2,1,ie) * r(8 + n * (j - 1), ie) &
                         + s(i,9,2,1,ie) * r(9 + n * (j - 1), ie) &
                         + s(i,10,2,1,ie) * r(10 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,1,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,1,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,1,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,1,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,1,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,1,2,ie) &
                            + wrk(l + n * (9 - 1) + nn * (i - 1), ie) &
                              * s(9,j,1,2,ie) &
                            + wrk(l + n * (10 - 1) + nn * (i - 1), ie) &
                              * s(10,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 1, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 1, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 1, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 1, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 1, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 1, 3, ie) &
                      + wrk2(i + nn * (9 - 1), ie) * s(9, j, 1, 3, ie) &
                      + wrk2(i + nn * (10 - 1), ie) * s(10, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,1,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,1,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,1,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,1,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,1,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,1,1,ie) * r(8 + n * (j - 1), ie) &
                         + s(i,9,1,1,ie) * r(9 + n * (j - 1), ie) &
                         + s(i,10,1,1,ie) * r(10 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,2,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,2,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,2,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,2,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,2,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,2,2,ie) &
                            + wrk(l + n * (9 - 1) + nn * (i - 1), ie) &
                              * s(9,j,2,2,ie) &
                            + wrk(l + n * (10 - 1) + nn * (i - 1), ie) &
                              * s(10,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 2, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 2, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 2, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 2, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 2, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 2, 3, ie) &
                      + wrk2(i + nn * (9 - 1), ie) * s(9, j, 2, 3, ie) &
                      + wrk2(i + nn * (10 - 1), ie) * s(10, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl10
  
  subroutine fdm_do_fast_sx_nl9(e, r, s, d, nelv)
    integer, parameter :: n = 9
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,2,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,2,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,2,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,2,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,2,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,2,1,ie) * r(8 + n * (j - 1), ie) &
                         + s(i,9,2,1,ie) * r(9 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,1,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,1,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,1,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,1,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,1,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,1,2,ie) &
                            + wrk(l + n * (9 - 1) + nn * (i - 1), ie) &
                              * s(9,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 1, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 1, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 1, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 1, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 1, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 1, 3, ie) &
                      + wrk2(i + nn * (9 - 1), ie) * s(9, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,1,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,1,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,1,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,1,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,1,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,1,1,ie) * r(8 + n * (j - 1), ie) &
                         + s(i,9,1,1,ie) * r(9 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,2,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,2,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,2,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,2,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,2,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,2,2,ie) &
                            + wrk(l + n * (9 - 1) + nn * (i - 1), ie) &
                              * s(9,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 2, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 2, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 2, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 2, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 2, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 2, 3, ie) &
                      + wrk2(i + nn * (9 - 1), ie) * s(9, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl9

  subroutine fdm_do_fast_sx_nl8(e, r, s, d, nelv)
    integer, parameter :: n = 8
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,2,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,2,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,2,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,2,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,2,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,2,1,ie) * r(8 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,1,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,1,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,1,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,1,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,1,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 1, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 1, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 1, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 1, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 1, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,1,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,1,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,1,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,1,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,1,1,ie) * r(7 + n * (j - 1), ie) &
                         + s(i,8,1,1,ie) * r(8 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,2,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,2,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,2,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,2,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,2,2,ie) &
                            + wrk(l + n * (8 - 1) + nn * (i - 1), ie) &
                              * s(8,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 2, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 2, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 2, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 2, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 2, 3, ie) &
                      + wrk2(i + nn * (8 - 1), ie) * s(8, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl8

  subroutine fdm_do_fast_sx_nl7(e, r, s, d, nelv)
    integer, parameter :: n = 7
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,2,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,2,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,2,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,2,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,2,1,ie) * r(7 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,1,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,1,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,1,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,1,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 1, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 1, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 1, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 1, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,1,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,1,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,1,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,1,1,ie) * r(6 + n * (j - 1), ie) &
                         + s(i,7,1,1,ie) * r(7 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,2,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,2,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,2,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,2,2,ie) &
                            + wrk(l + n * (7 - 1) + nn * (i - 1), ie) &
                              * s(7,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 2, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 2, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 2, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 2, 3, ie) &
                      + wrk2(i + nn * (7 - 1), ie) * s(7, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl7

  subroutine fdm_do_fast_sx_nl6(e, r, s, d, nelv)
    integer, parameter :: n = 6
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,2,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,2,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,2,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,2,1,ie) * r(6 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,1,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,1,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,1,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 1, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 1, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 1, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,1,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,1,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,1,1,ie) * r(5 + n * (j - 1), ie) &
                         + s(i,6,1,1,ie) * r(6 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,2,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,2,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,2,2,ie) &
                            + wrk(l + n * (6 - 1) + nn * (i - 1), ie) &
                              * s(6,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 2, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 2, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 2, 3, ie) &
                      + wrk2(i + nn * (6 - 1), ie) * s(6, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl6

  subroutine fdm_do_fast_sx_nl5(e, r, s, d, nelv)
    integer, parameter :: n = 5
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,2,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,2,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,2,1,ie) * r(5 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,1,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,1,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 1, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 1, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,1,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,1,1,ie) * r(4 + n * (j - 1), ie) &
                         + s(i,5,1,1,ie) * r(5 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,2,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,2,2,ie) &
                            + wrk(l + n * (5 - 1) + nn * (i - 1), ie) &
                              * s(5,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 2, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 2, 3, ie) &
                      + wrk2(i + nn * (5 - 1), ie) * s(5, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl5

  subroutine fdm_do_fast_sx_nl4(e, r, s, d, nelv)
    integer, parameter :: n = 4
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,2,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,2,1,ie) * r(4 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,1,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 1, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,1,1,ie) * r(3 + n * (j - 1), ie) &
                         + s(i,4,1,1,ie) * r(4 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,2,2,ie) &
                            + wrk(l + n * (4 - 1) + nn * (i - 1), ie) &
                              * s(4,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 2, 3, ie) &
                      + wrk2(i + nn * (4 - 1), ie) * s(4, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl4

  subroutine fdm_do_fast_sx_nl3(e, r, s, d, nelv)
    integer, parameter :: n = 3
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,2,1,ie) * r(3 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) &
                         + s(i,3,1,1,ie) * r(3 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) &
                            + wrk(l + n * (3 - 1) + nn * (i - 1), ie) &
                              * s(3,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) &
                      + wrk2(i + nn * (3 - 1), ie) * s(3, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl3

  subroutine fdm_do_fast_sx_nl2(e, r, s, d, nelv)
    integer, parameter :: n = 2
    integer, parameter :: nn = n**2
    integer, parameter :: nnn = n**3
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: e(n**3, nelv)
    real(kind=rp), intent(inout) :: r(n**3, nelv)
    real(kind=rp), intent(inout) :: s(n,n,2,3, nelv)
    real(kind=rp), intent(inout) :: d(n**3, nelv)
    real(kind=rp) :: wrk(n**3, nelv), wrk2(n**3, nelv)
    integer ::  ie, i, j, k, l, ii, jj
       
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,2,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,2,1,ie) * r(2 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,1,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,1,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 1, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 1, 3, ie) 
          end do
       end do
    end do

    do i = 1, nnn * nelv
       r(i,1) = d(i,1) * e(i,1)
    end do
   
    do j = 1, nn
       do i = 1, n
          do ie = 1, nelv
             ii = i + n * (j - 1)
             wrk(ii, ie) = s(i,1,1,1,ie) * r(1 + n * (j - 1), ie) &
                         + s(i,2,1,1,ie) * r(2 + n * (j - 1), ie) 
          end do
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             do ie = 1, nelv
                ii = l + n * (j - 1) + nn * (i - 1)
                wrk2(ii,ie) = wrk(l + n * (1 - 1) + nn * (i - 1), ie) &
                              * s(1,j,2,2,ie) &
                            + wrk(l + n * (2 - 1) + nn * (i - 1), ie) &
                              * s(2,j,2,2,ie) 
             end do
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          do ie = 1, nelv
             jj = i + nn * (j - 1)
             e(jj,ie) = wrk2(i + nn * (1 - 1), ie) * s(1, j, 2, 3, ie) &
                      + wrk2(i + nn * (2 - 1), ie) * s(2, j, 2, 3, ie) 
          end do
       end do
    end do
    

  end subroutine fdm_do_fast_sx_nl2
    
end module fdm_sx
